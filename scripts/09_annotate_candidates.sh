#!/bin/bash
# 09_annotate_candidates.sh
# Phase 4: Generate rich per-candidate feature table and apply sgRNA filter
#
# Produces a comprehensive annotation table with:
# - Nearest gene + distance
# - Overlapping/nearby CREs (if ChIP data available)
# - Conservation score
# - TAD context
# - GC%
# - Repeat content
# - Stages open in
# - Best sgRNA score (if available)
# - Composite rank

set -euo pipefail

# Load configuration
if [ ! -f "./config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ./config.sh

echo "=== Phase 4: Candidate Annotation & Feature Table ==="
echo ""

cd results

# Tier parameters (defaults preserve strict pipeline):
#   HARBOR_PREFIX  — input candidate prefix (e.g. safe_harbors_relaxed)
#   FEATURE_TABLE  — output table name
HP="${HARBOR_PREFIX:-safe_harbors}"
FEATURE_TABLE="${FEATURE_TABLE:-candidate_feature_table.tsv}"

# ======================================================================
# Input: use biologically filtered candidates from Phase 3
# Fall back to sized candidates from script 08 if Phase 3 hasn't run
# ======================================================================
if [ -s "${HP}_biological_100pct_sized.bed" ]; then
    INPUT="${HP}_biological_100pct_sized.bed"
    echo "Using Phase 3 filtered candidates: $INPUT"
elif [ -s "${HP}_100pct_sized.bed" ]; then
    INPUT="${HP}_100pct_sized.bed"
    echo "WARNING: Phase 3 output not found. Using basic filtered candidates: $INPUT"
else
    echo "ERROR: No candidate files found. Run scripts 08 (and 08b) first."
    exit 1
fi

NUM_CANDIDATES=$(wc -l < "$INPUT")
echo "Annotating $NUM_CANDIDATES candidates..."
echo ""

# Sort input for bedtools
sort -k1,1 -k2,2n "$INPUT" > candidates_sorted.bed

# ======================================================================
# Step 1: Nearest gene + distance
# ======================================================================
echo "--- Step 1: Nearest gene annotation ---"

# Prepare gene file with chr prefix
awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    ../annotation/genes.bed \
| sort -k1,1 -k2,2n \
> genes_chr_sorted.bed

# Ensure consistent sort order for both files
LC_ALL=C sort -k1,1 -k2,2n candidates_sorted.bed > candidates_resort.bed
LC_ALL=C sort -k1,1 -k2,2n genes_chr_sorted.bed > genes_resort.bed
mv candidates_resort.bed candidates_sorted.bed
mv genes_resort.bed genes_chr_sorted.bed

bedtools closest \
    -a candidates_sorted.bed \
    -b genes_chr_sorted.bed \
    -d \
    -t first \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $8, $NF}' \
> nearest_gene_annotated.txt

echo "  Done. Results in nearest_gene_annotated.txt"
echo ""

# ======================================================================
# Step 2: GC content
# ======================================================================
echo "--- Step 2: GC content ---"

GENOME_FASTA="../reference/GRCz11.fa"
if [ -f "$GENOME_FASTA" ]; then
    # Strip chr prefix for genome FASTA compatibility (FASTA uses bare chrom names)
    awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); print $1, $2, $3}' candidates_sorted.bed \
    | bedtools nuc \
        -fi "$GENOME_FASTA" \
        -bed - \
    | tail -n +2 \
    | awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, $5}' \
    > gc_content.txt
    echo "  Done. GC content computed for all candidates."
else
    echo "  [SKIP] Genome FASTA not found at $GENOME_FASTA"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > gc_content.txt
fi
echo ""

# ======================================================================
# Step 3: Repeat content (fraction of candidate overlapping repeats)
# ======================================================================
echo "--- Step 3: Repeat content ---"

REPEATS="../annotation/repeats_all.bed"
if [ -f "$REPEATS" ]; then
    bedtools coverage \
        -a candidates_sorted.bed \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$REPEATS" | cut -f1-3 | sort -k1,1 -k2,2n) \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF}' \
    > repeat_content.txt
    echo "  Done. Repeat fraction computed."
else
    echo "  [SKIP] repeats_all.bed not found"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > repeat_content.txt
fi
echo ""

# ======================================================================
# Step 4: Conservation score (mean across candidate)
# ======================================================================
echo "--- Step 4: Conservation annotation ---"

if [ -f "$CONSERVATION_BED" ] && [ -s "$CONSERVATION_BED" ]; then
    bedtools coverage \
        -a candidates_sorted.bed \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$CONSERVATION_BED" | sort -k1,1 -k2,2n) \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF}' \
    > conservation_overlap.txt
    echo "  Done. Conservation overlap fraction computed."
else
    echo "  [SKIP] Conservation track not available"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > conservation_overlap.txt
fi
echo ""

# ======================================================================
# Step 5: TAD boundary distance
# ======================================================================
echo "--- Step 5: TAD boundary distance ---"

if [ -f "$TAD_BOUNDARIES_BED" ] && [ -s "$TAD_BOUNDARIES_BED" ]; then
    bedtools closest \
        -a candidates_sorted.bed \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$TAD_BOUNDARIES_BED" | sort -k1,1 -k2,2n) \
        -d \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF}' \
    > tad_distance.txt
    echo "  Done. TAD boundary distances computed."
else
    echo "  [SKIP] TAD boundary data not available"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > tad_distance.txt
fi
echo ""

# ======================================================================
# Step 6: sgRNA targetability
# ======================================================================
echo "--- Step 6: sgRNA targetability ---"

if [ -f "$SGRNA_BED" ] && [ -s "$SGRNA_BED" ]; then
    echo "  Filtering candidates by sgRNA availability (min specificity: $SGRNA_MIN_SPECIFICITY)..."

    # Find candidates with at least one good sgRNA
    bedtools intersect \
        -a candidates_sorted.bed \
        -b <(awk -v min="$SGRNA_MIN_SPECIFICITY" 'BEGIN{OFS="\t"} $5 >= min {if ($1 !~ /^chr/) $1 = "chr" $1; print}' "$SGRNA_BED" | sort -k1,1 -k2,2n) \
        -u \
    > candidates_with_sgrna.bed

    # Get best sgRNA score per candidate
    bedtools intersect \
        -a candidates_sorted.bed \
        -b <(awk -v min="$SGRNA_MIN_SPECIFICITY" 'BEGIN{OFS="\t"} $5 >= min {if ($1 !~ /^chr/) $1 = "chr" $1; print}' "$SGRNA_BED" | sort -k1,1 -k2,2n) \
        -wa -wb \
    | awk 'BEGIN{OFS="\t"} {
        key = $1 ":" $2 ":" $3
        if (!(key in best) || $NF > best[key]) best[key] = $NF
    }
    END {
        for (key in best) {
            split(key, parts, ":")
            print parts[1], parts[2], parts[3], best[key]
        }
    }' | sort -k1,1 -k2,2n > sgrna_best_score.txt

    TARGETABLE=$(wc -l < candidates_with_sgrna.bed)
    echo "  Candidates with ≥1 good sgRNA: $TARGETABLE / $NUM_CANDIDATES"
    echo "  Best sgRNA scores written to sgrna_best_score.txt"
else
    echo "  [SKIP] sgRNA data not available (run CRISPOR/FlashFry externally)"
    echo "  >> See DATA_NEEDED.md for instructions"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > sgrna_best_score.txt
    cp candidates_sorted.bed candidates_with_sgrna.bed
fi
echo ""

# ======================================================================
# Step 7: Stages open in (from consensus peak calling)
# ======================================================================
echo "--- Step 7: Stage count annotation ---"

# The consensus peaks should have stage count in column 4 from script 07
# If not, compute from the input file
awk 'BEGIN{OFS="\t"} {
    if (NF >= 4) print $1, $2, $3, $4
    else print $1, $2, $3, "NA"
}' candidates_sorted.bed > stage_counts.txt
echo "  Stage counts extracted from candidate file."
echo ""

# ======================================================================
# Step 7b: COPE/DOPE overlap flags (DANIO-CODE positive "open but inert" signal)
# A candidate overlapping a COPE is gene-distal AND independently validated by
# DANIO-CODE as constitutively open but lacking chromatin marks — a high-
# confidence positive. (Annotation only; not used in the composite rank.)
# ======================================================================
echo "--- Step 7b: COPE/DOPE overlap flags ---"

# config paths are repo-root-relative; we run from results/, so prefix ../
_cope="$COPE_BED"; case "$_cope" in /*) : ;; *) _cope="../$_cope" ;; esac
_dope="$DOPE_BED"; case "$_dope" in /*) : ;; *) _dope="../$_dope" ;; esac

flag_overlap() {  # $1=track  $2=output (chr,start,end,0/1 in candidate order)
    if [ -f "$1" ] && [ -s "$1" ]; then
        bedtools intersect -a candidates_sorted.bed \
            -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$1" | cut -f1-3 | sort -k1,1 -k2,2n) -c \
        | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ($NF>0)?1:0}' > "$2"
    else
        awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NA"}' candidates_sorted.bed > "$2"
    fi
}
flag_overlap "$_cope" cope_flag.txt
flag_overlap "$_dope" dope_flag.txt
echo "  COPE-overlapping: $(awk '$4==1' cope_flag.txt | wc -l | tr -d ' ') / $NUM_CANDIDATES; DOPE-overlapping: $(awk '$4==1' dope_flag.txt | wc -l | tr -d ' ')"
echo ""

# ======================================================================
# Step 8: Build composite feature table
# ======================================================================
echo "--- Step 8: Building composite feature table ---"
echo ""

# Header
echo -e "chrom\tstart\tend\tsize\tstage_count\tnearest_gene\tgene_distance\tgc_pct\trepeat_frac\tconservation_frac\ttad_distance\tsgrna_score\tin_cope\tin_dope\tcomposite_rank" \
    > ${FEATURE_TABLE}

# Join all annotations
paste \
    <(awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' candidates_sorted.bed) \
    <(cut -f4 stage_counts.txt) \
    <(cut -f4 nearest_gene_annotated.txt) \
    <(cut -f5 nearest_gene_annotated.txt) \
    <(cut -f4 gc_content.txt) \
    <(cut -f4 repeat_content.txt) \
    <(cut -f4 conservation_overlap.txt) \
    <(cut -f4 tad_distance.txt) \
    <(cut -f4 sgrna_best_score.txt) \
    <(cut -f4 cope_flag.txt) \
    <(cut -f4 dope_flag.txt) \
| awk 'BEGIN{OFS="\t"} {
    # Composite rank with distance penalty (in_cope/in_dope are cols 13/14,
    # annotation only — not used in the score; composite printed as last col)
    stage = ($5 == "NA" || $5 == ".") ? 0 : $5+0
    size = ($4 == "NA" || $4 == ".") ? 0 : $4+0
    gene_dist = ($7 == "NA" || $7 == "." || $7 == "" || $7+0 < 0) ? 100000 : $7+0
    cons = ($10 == "NA" || $10 == ".") ? 0 : $10+0
    rep = ($9 == "NA" || $9 == ".") ? 0 : $9+0
    sgrna = ($12 == "NA" || $12 == ".") ? 0 : $12+0

    # Distance bonus/penalty
    dist_bonus = (gene_dist > 50000) ? 50 : (gene_dist > 10000) ? 30 : (gene_dist > 0) ? 10 : -200

    # Composite: prioritize stage count, size, distance from genes
    composite = (stage * 100) + (size / 100) + (sgrna * 50) - (cons * 200) - (rep * 100) + dist_bonus
    print $0, composite
}' \
| sort -t$'\t' -k15,15 -nr \
>> ${FEATURE_TABLE}

TOTAL_ANNOTATED=$(tail -n +2 ${FEATURE_TABLE} | wc -l)
echo "Feature table complete: $TOTAL_ANNOTATED candidates annotated"
echo "Output: results/${FEATURE_TABLE}"
echo ""

# ======================================================================
# Display top candidates
# ======================================================================
echo "=== Top 20 Candidates (by composite rank) ==="
echo ""
head -21 ${FEATURE_TABLE} | column -t
echo ""

# ======================================================================
# Summary
# ======================================================================
echo "=== Annotation Summary ==="
echo "  Total candidates annotated: $TOTAL_ANNOTATED"
echo "  Feature table: ${FEATURE_TABLE}"
if [ -f "$SGRNA_BED" ] && [ -s "$SGRNA_BED" ]; then
    echo "  sgRNA-targetable candidates: $(wc -l < candidates_with_sgrna.bed)"
fi
echo ""

# Clean up temp files
rm -f candidates_sorted.bed genes_chr_sorted.bed nearest_gene_raw.txt
rm -f nearest_gene_annotated.txt gc_content.txt repeat_content.txt
rm -f conservation_overlap.txt tad_distance.txt sgrna_best_score.txt
rm -f cope_flag.txt dope_flag.txt
rm -f stage_counts.txt candidates_with_sgrna.bed

echo "=== Annotation complete ==="