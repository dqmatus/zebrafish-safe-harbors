#!/bin/bash
# 07_find_consensus_peaks.sh
# Identifies consensus peaks using peak-summit approach with score normalization
# and replicate collapsing
#
# PHASE 2 REDESIGN (COMPLETE):
# 1. Extract peak signal with normalized scores (col-10 summit not retained downstream)
# 2. Collapse replicates to biosample level (median score)
# 3. Collapse biosamples to stage level (median score)
# 4. Create fixed windows (±SUMMIT_WINDOW_HALF) around each summit
# 5. Merge windows and count how many stages have a summit
# 6. Filter by consensus tiers (80%, 90%, 100% of stages)
#
# This ensures consensus reflects true biological consistency across developmental
# stages, not just technical replication depth.

set -euo pipefail

# Load configuration
if [ ! -f "../config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ../config.sh

echo "=== Finding Consensus ATAC-seq Peaks (Summit-Based, Score-Normalized) ==="
echo ""
echo "Configuration:"
echo "  Summit window: ±${SUMMIT_WINDOW_HALF}bp"
echo "  Score normalization: ${SCORE_NORMALIZATION_METHOD}"
echo "  Consensus tiers: ${CONSENSUS_MIN_PCT_80}%, ${CONSENSUS_MIN_PCT_90}%, ${CONSENSUS_MIN_PCT_100}%"
echo ""

cd ../peaks

# Check for normalized scores
if [ ! -d "normalized_scores" ]; then
    echo "ERROR: normalized_scores/ directory not found"
    echo "Run: python3 ../scripts/06b_normalize_scores.py"
    exit 1
fi

# ============================================================================
# Load accession metadata for replicate collapsing
# ============================================================================
echo "Loading DANIO-CODE accession metadata..."

METADATA_FILE="../annotation/danio_code_accession_map.txt"
if [ ! -f "$METADATA_FILE" ]; then
    echo "WARNING: $METADATA_FILE not found"
    echo "Treating each file as a unique biosample/stage"
    # Create fallback metadata: each normalized file is its own biosample/stage
    mkdir -p metadata_temp
    for f in normalized_scores/*.normalized.bed; do
        basename=$(basename "$f" .normalized.bed)
        # Extract SQ ID from filename (assume format: ATAC-seq_Lab_AS.SQ_ID...)
        sq_id=$(echo "$basename" | grep -oE 'DCD[0-9]+SQ' || echo "$basename")
        echo -e "$sq_id\t$sq_id\t$basename\tunknown_lab\tfallback" >> metadata_temp/fallback.txt
    done
    METADATA_FILE="metadata_temp/fallback.txt"
fi

# ============================================================================
# Step 1: Collapse replicates to biosample level
# ============================================================================
echo ""
echo "Step 1: Collapsing replicates to biosample level..."

mkdir -p biosamples_collapsed

# For each biosample, merge all replicates and take median score
awk '
!/^#/ && NF > 0 {
    sq_id = $1
    bs_id = $2
    stage = $3
    biosamples[bs_id] = 1
    stages[bs_id] = stage
}
END {
    for (bs in biosamples) print bs, stages[bs]
}
' "$METADATA_FILE" | sort -u > biosamples_list.txt

echo "Found $(wc -l < biosamples_list.txt) unique biosamples"

while read biosample stage; do
    # Find all SQ files for this biosample
    sqs=$(awk -v bs="$biosample" '!/^#/ && $2 == bs {print $1}' "$METADATA_FILE" | tr '\n' ' ')

    if [ -z "$sqs" ]; then
        continue
    fi

    # Cat only the normalized files belonging to this biosample's SQ accessions
    > "biosamples_collapsed/${biosample}.tmp"
    for sq in $sqs; do
        matching_file=$(ls normalized_scores/*${sq}*.normalized.bed 2>/dev/null || true)
        if [ -n "$matching_file" ]; then
            cat "$matching_file" >> "biosamples_collapsed/${biosample}.tmp"
        fi
    done

    if [ ! -s "biosamples_collapsed/${biosample}.tmp" ]; then
        rm -f "biosamples_collapsed/${biosample}.tmp"
        continue
    fi

    # Merge overlapping peaks and take median score per region
    sort -k1,1 -k2,2n "biosamples_collapsed/${biosample}.tmp" \
    | bedtools merge -i - -c 5 -o median \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $4, "."}' \
    > "biosamples_collapsed/${biosample}.bed"

    rm -f "biosamples_collapsed/${biosample}.tmp"

    num_sqs=$(echo $sqs | wc -w | tr -d ' ')
    echo "  Collapsed $biosample (stage: $stage, replicates: $num_sqs)"
done < biosamples_list.txt

echo "Biosample collapsing complete: $(ls biosamples_collapsed/ | wc -l) biosamples"

# ============================================================================
# Step 2: Collapse biosamples to stage level
# ============================================================================
echo ""
echo "Step 2: Collapsing biosamples to stage level..."

mkdir -p stages_collapsed

# Group biosamples by stage
awk '
!/^#/ && NF > 0 {
    stages[$3] = stages[$3] " " $2
}
END {
    for (stage in stages) print stage, stages[stage]
}
' "$METADATA_FILE" | sort -u > stages_list.txt

num_stages=$(wc -l < stages_list.txt)
echo "Found $num_stages unique developmental stages"

while read stage biosamples; do
    # Remove leading space
    biosamples=$(echo "$biosamples" | sed 's/^ //')

    # Merge all peaks from biosamples at this stage
    > "stages_collapsed/${stage}.tmp"
    for bs in $biosamples; do
        if [ -f "biosamples_collapsed/${bs}.bed" ] && [ -s "biosamples_collapsed/${bs}.bed" ]; then
            cat "biosamples_collapsed/${bs}.bed"
        fi
    done > "stages_collapsed/${stage}.tmp"

    if [ ! -s "stages_collapsed/${stage}.tmp" ]; then
        rm -f "stages_collapsed/${stage}.tmp"
        continue
    fi

    # Merge overlapping regions and take median score
    sort -k1,1 -k2,2n "stages_collapsed/${stage}.tmp" \
    | bedtools merge -i - -c 5 -o median \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $4, "."}' \
    > "stages_collapsed/${stage}.bed"

    rm -f "stages_collapsed/${stage}.tmp"

    echo "  Collapsed stage $stage (biosamples: $(echo $biosamples | wc -w | tr -d ' '))"
done < stages_list.txt

echo "Stage collapsing complete: $num_stages stages"

# ============================================================================
# Step 3: Create fixed windows around stage-level summits
# ============================================================================
echo ""
echo "Step 3: Creating fixed windows around stage-level summits..."

mkdir -p stage_summit_windows

while read stage_file; do
    stage=$(basename "$stage_file" .bed)
    
    # For each peak, create fixed window around summit
    awk -v window="${SUMMIT_WINDOW_HALF}" 'BEGIN {OFS="\t"} {
        chrom = $1
        start = $2
        end = $3
        score = $5
        
        # Assume summit is at peak midpoint (conservative; ideally use narrowPeak col-10)
        summit = int((start + end) / 2)
        
        # Create fixed window
        window_start = summit - window
        window_end = summit + window
        
        # Ensure window stays within peak bounds
        if (window_start < start) window_start = start
        if (window_end > end) window_end = end
        
        print chrom, window_start, window_end, $4"_window", score, "."
    }' "$stage_file" > "stage_summit_windows/${stage}.windows.bed"
    
done < <(ls stages_collapsed/*.bed)

echo "Created summit windows for all stages"

# ============================================================================
# Step 4: Merge summit windows and count stage overlap
# ============================================================================
echo ""
echo "Step 4: Merging stage summit windows and counting overlaps..."

# Combine all stage summit windows
cat stage_summit_windows/*.windows.bed | cut -f1-3 | sort -k1,1 -k2,2n > all_stage_windows_combined.bed
bedtools merge -i all_stage_windows_combined.bed > merged_stage_windows.bed

echo "Merged regions: $(wc -l < merged_stage_windows.bed)"

# Count how many stages have a window overlapping each merged region
bedtools annotate -i merged_stage_windows.bed \
    -files stage_summit_windows/*.windows.bed -counts \
    > stage_overlap_counts.txt

echo "Calculated stage overlap counts"

# ============================================================================
# Step 5: Filter by consensus tiers
# ============================================================================
echo ""
echo "Filtering by consensus tiers (% of developmental stages)..."

# Calculate minimum stage counts for each tier
min_stages_80=$(echo "scale=0; $num_stages * $CONSENSUS_MIN_PCT_80 / 100" | bc | cut -d. -f1)
min_stages_90=$(echo "scale=0; $num_stages * $CONSENSUS_MIN_PCT_90 / 100" | bc | cut -d. -f1)
min_stages_100=$num_stages

echo "Minimum stages required:"
echo "  80% consensus: $min_stages_80 / $num_stages"
echo "  90% consensus: $min_stages_90 / $num_stages"
echo "  100% consensus: $min_stages_100 / $num_stages"

# Extract consensus peaks: count stage columns (4 onwards) for non-zero entries
awk -v min="$min_stages_80" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' stage_overlap_counts.txt > consensus_peaks_80pct.bed

awk -v min="$min_stages_90" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' stage_overlap_counts.txt > consensus_peaks_90pct.bed

awk -v min="$min_stages_100" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' stage_overlap_counts.txt > consensus_peaks_100pct.bed

echo ""
echo "=== Results ==="
echo "Consensus peaks (>=80% stages, $min_stages_80/$num_stages): $(wc -l < consensus_peaks_80pct.bed)"
echo "Consensus peaks (>=90% stages, $min_stages_90/$num_stages): $(wc -l < consensus_peaks_90pct.bed)"
echo "Consensus peaks (>=100% stages, $min_stages_100/$num_stages): $(wc -l < consensus_peaks_100pct.bed)"

# ============================================================================
# Quality Control
# ============================================================================
echo ""
echo "=== Quality Control ==="

# Check region sizes
if [ -s consensus_peaks_80pct.bed ]; then
    median_size_80=$(awk '{print $3 - $2}' consensus_peaks_80pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')
    echo "Median region size (80% consensus): ${median_size_80:-0} bp"
fi

if [ -s consensus_peaks_90pct.bed ]; then
    median_size_90=$(awk '{print $3 - $2}' consensus_peaks_90pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')
    echo "Median region size (90% consensus): ${median_size_90:-0} bp"
fi

if [ -s consensus_peaks_100pct.bed ]; then
    median_size_100=$(awk '{print $3 - $2}' consensus_peaks_100pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')
    echo "Median region size (100% consensus): ${median_size_100:-0} bp"
fi

# ============================================================================
# Copy results to results folder
# ============================================================================
mkdir -p ../results
cp consensus_peaks_*.bed ../results/

echo ""
echo "=== Consensus peak finding complete ==="
echo "Output files copied to ../results/"
echo ""
echo "Summary:"
echo "  - Replicate collapsing: SQ → biosample → stage"
echo "  - Score normalization: $SCORE_NORMALIZATION_METHOD"
echo "  - Consensus method: stage-level summit windows (±${SUMMIT_WINDOW_HALF}bp)"
echo "  - Final consensus tiers: 80/90/100% of developmental stages"
