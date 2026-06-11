#!/bin/bash
# 07_find_consensus_peaks.sh
# Identifies consensus peaks using peak summit approach
# 
# PHASE 2 REDESIGN: Consensus methodology now uses peak summits (narrowPeak col-10)
# with a fixed window (parameterized ±SUMMIT_WINDOW_HALF bp). Counts biosamples/stages
# with a summit in each window, rather than merging all peaks and counting overlaps.
# This avoids inflating consensus by large merged regions that represent different
# sub-intervals from different samples.

set -euo pipefail

# Load configuration
if [ ! -f "../config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ../config.sh

echo "=== Finding Consensus ATAC-seq Peaks (Summit-Based) ==="
echo "Configuration:"
echo "  Summit window: ±${SUMMIT_WINDOW_HALF}bp"
echo "  Consensus tiers: ${CONSENSUS_MIN_PCT_80}%, ${CONSENSUS_MIN_PCT_90}%, ${CONSENSUS_MIN_PCT_100}%"
echo ""

cd peaks

# Count samples
num_samples=$(ls *.gz 2>/dev/null | wc -l | tr -d ' ')
if [ "$num_samples" -eq 0 ]; then
    echo "ERROR: No .gz peak files found in $(pwd)"
    exit 1
fi
echo "Number of samples: $num_samples"

# Decompress all files
echo "Decompressing peak files..."
for f in *.gz; do
    gunzip -k "$f" 2>/dev/null || true
done

# ============================================================================
# Extract peak summits from each sample and convert to fixed windows
# ============================================================================
echo ""
echo "Extracting peak summits and converting to fixed windows..."

mkdir -p sample_summits sample_beds

for narrowpeak in *.narrowPeak; do
    sample=$(basename "$narrowpeak" .narrowPeak)
    
    # narrowPeak format: chrom chromStart chromEnd name score strand pValue qValue fold_change peak_offset
    # Column 10 (peak_offset) is the offset from chromStart to peak summit
    # We create a fixed window around the summit: [summit - WINDOW_HALF, summit + WINDOW_HALF]
    
    awk -v window="${SUMMIT_WINDOW_HALF}" 'BEGIN{OFS="\t"} {
        chrom = $1
        start = $2
        end = $3
        name = $4
        score = $5
        strand = $6
        peak_offset = $10
        
        # Calculate summit position
        summit = start + peak_offset
        
        # Create fixed window around summit
        window_start = summit - window
        window_end = summit + window
        
        # Ensure window stays within peak bounds (sanity check)
        if (window_start < start) window_start = start
        if (window_end > end) window_end = end
        
        print chrom, window_start, window_end, name"_"sample, score, strand
    }' "$narrowpeak" | sort -k1,1 -k2,2n > "sample_summits/${sample}.summits.bed"
done

echo "Created summit windows for $num_samples samples"

# ============================================================================
# Merge summit windows across all samples to define consensus regions
# ============================================================================
echo ""
echo "Merging summit windows across all samples..."

# Combine all summit windows
cat sample_summits/*.summits.bed | cut -f1-3 | sort -k1,1 -k2,2n > all_summits_combined.bed
bedtools merge -i all_summits_combined.bed > merged_summit_regions.bed

echo "Merged regions: $(wc -l < merged_summit_regions.bed)"

# ============================================================================
# Count how many samples have a summit in each merged region
# ============================================================================
echo ""
echo "Counting sample overlaps for each consensus region..."

# For each merged region, count how many sample summit files have overlap
bedtools annotate -i merged_summit_regions.bed \
    -files sample_summits/*.summits.bed -counts \
    > peak_overlap_counts.txt

echo "Calculated overlap counts for all regions"

# ============================================================================
# Filter by consensus tiers (80%, 90%, 100%)
# ============================================================================
echo ""
echo "Filtering by consensus tiers..."

# Calculate minimum sample counts for each tier
min_samples_80=$(echo "scale=0; $num_samples * $CONSENSUS_MIN_PCT_80 / 100" | bc | cut -d. -f1)
min_samples_90=$(echo "scale=0; $num_samples * $CONSENSUS_MIN_PCT_90 / 100" | bc | cut -d. -f1)
min_samples_100=$num_samples

echo "Minimum samples required:"
echo "  80% consensus: $min_samples_80 / $num_samples"
echo "  90% consensus: $min_samples_90 / $num_samples"
echo "  100% consensus: $min_samples_100 / $num_samples"

# Extract consensus peaks (sum the count columns to get total sample count)
# Column 4 onwards contain per-sample counts; sum them to get total overlapping samples
awk -v min="$min_samples_80" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' peak_overlap_counts.txt > consensus_peaks_80pct.bed

awk -v min="$min_samples_90" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' peak_overlap_counts.txt > consensus_peaks_90pct.bed

awk -v min="$min_samples_100" 'BEGIN {OFS="\t"} {
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' peak_overlap_counts.txt > consensus_peaks_100pct.bed

echo ""
echo "=== Results ==="
echo "Consensus peaks (>=80% samples, $min_samples_80/$num_samples): $(wc -l < consensus_peaks_80pct.bed)"
echo "Consensus peaks (>=90% samples, $min_samples_90/$num_samples): $(wc -l < consensus_peaks_90pct.bed)"
echo "Consensus peaks (>=100% samples, $min_samples_100/$num_samples): $(wc -l < consensus_peaks_100pct.bed)"

# ============================================================================
# Validation & QC
# ============================================================================
echo ""
echo "=== Quality Control ==="

# Check that consensus regions have reasonable sizes
median_size_80=$(awk '{print $3 - $2}' consensus_peaks_80pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')
median_size_90=$(awk '{print $3 - $2}' consensus_peaks_90pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')
median_size_100=$(awk '{print $3 - $2}' consensus_peaks_100pct.bed | sort -n | awk '{a[NR]=$1} END {if (NR%2==0) print (a[NR/2]+a[NR/2+1])/2; else print a[(NR+1)/2]}')

echo "Median region sizes:"
echo "  80% consensus: ${median_size_80:-0} bp"
echo "  90% consensus: ${median_size_90:-0} bp"
echo "  100% consensus: ${median_size_100:-0} bp"

# ============================================================================
# Copy results to results folder
# ============================================================================
mkdir -p ../results
cp consensus_peaks_*.bed ../results/

echo ""
echo "=== Consensus peak finding complete ==="
echo "Output files copied to ../results/"
