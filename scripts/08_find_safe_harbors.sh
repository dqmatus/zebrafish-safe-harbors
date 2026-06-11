#!/bin/bash
# 08_find_safe_harbors.sh
# Filters consensus peaks against exclusion regions to find safe harbors

set -euo pipefail

# Load configuration
if [ ! -f "../config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ../config.sh

echo "=== Finding Safe Harbor Candidates ==="

mkdir -p results
cd results

# FIX A1: Guard against chromosome naming mismatch before filtering
# If consensus peaks and exclusions have zero chrom overlap, abort with clear error.
# This prevents silent regression where bedtools subtract silently returns empty.
PEAK_CHROMS=$(cut -f1 consensus_peaks_100pct.bed 2>/dev/null | sort -u | wc -l)
EXCL_CHROMS=$(cut -f1 ../annotation/exclusions_minimal.bed 2>/dev/null | sort -u | wc -l)
COMMON_CHROMS=$(comm -12 \
  <(cut -f1 consensus_peaks_100pct.bed 2>/dev/null | sort -u) \
  <(cut -f1 ../annotation/exclusions_minimal.bed 2>/dev/null | sort -u) | wc -l)

if [ "$COMMON_CHROMS" -eq 0 ]; then
  echo "ERROR: Chromosome naming mismatch detected!"
  echo "  Consensus peaks chroms: $(cut -f1 consensus_peaks_100pct.bed | sort -u | head -3 | tr '\n' ' ')"
  echo "  Exclusions chroms:       $(cut -f1 ../annotation/exclusions_minimal.bed | sort -u | head -3 | tr '\n' ' ')"
  echo "  (Expected: both 'chr1' or both '1')"
  exit 1
fi

# Check what we're starting with
echo "Starting consensus peaks (100%): $(wc -l < consensus_peaks_100pct.bed)"
echo "Starting consensus peaks (90%): $(wc -l < consensus_peaks_90pct.bed)"
echo "Starting consensus peaks (80%): $(wc -l < consensus_peaks_80pct.bed)"

# Filter against minimal exclusions
echo ""
echo "=== Filtering against minimal exclusions ==="
# FIX A2: Update log message to reflect gene-distal filtering
echo "(exons + essential genes + transposons + gene bodies ±10kb + telomeres)"

bedtools subtract \
    -a consensus_peaks_100pct.bed \
    -b ../annotation/exclusions_minimal.bed \
    > safe_harbors_100pct_minimal.bed

bedtools subtract \
    -a consensus_peaks_90pct.bed \
    -b ../annotation/exclusions_minimal.bed \
    > safe_harbors_90pct_minimal.bed

bedtools subtract \
    -a consensus_peaks_80pct.bed \
    -b ../annotation/exclusions_minimal.bed \
    > safe_harbors_80pct_minimal.bed

echo "After minimal filtering:"
echo "  100% consensus: $(wc -l < safe_harbors_100pct_minimal.bed)"
echo "  90% consensus: $(wc -l < safe_harbors_90pct_minimal.bed)"
echo "  80% consensus: $(wc -l < safe_harbors_80pct_minimal.bed)"

# Filter by minimum size (need room for landing pad)
echo ""
echo "=== Filtering by size (>=${MIN_HARBOR_SIZE}bp) ==="

awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_100pct_minimal.bed > safe_harbors_100pct_sized.bed
awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_90pct_minimal.bed > safe_harbors_90pct_sized.bed
awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_80pct_minimal.bed > safe_harbors_80pct_sized.bed

echo "After size filtering (>=${MIN_HARBOR_SIZE}bp):"
echo "  100% consensus: $(wc -l < safe_harbors_100pct_sized.bed)"
echo "  90% consensus: $(wc -l < safe_harbors_90pct_sized.bed)"
echo "  80% consensus: $(wc -l < safe_harbors_80pct_sized.bed)"

# Sort by size (largest first) and add size column
echo ""
echo "=== Creating final ranked lists ==="

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2, $4}' safe_harbors_100pct_sized.bed \
    | sort -k4 -nr > safe_harbors_final_100pct.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2, $4}' safe_harbors_90pct_sized.bed \
    | sort -k4 -nr > safe_harbors_final_90pct.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2, $4}' safe_harbors_80pct_sized.bed \
    | sort -k4 -nr > safe_harbors_final_80pct.bed

echo ""
echo "=== Top 20 Safe Harbor Candidates (100% consensus) ==="
echo "chrom	start	end	size	num_samples"
head -20 safe_harbors_final_100pct.bed

echo ""
echo "=== Summary ==="
echo "Final safe harbor candidates:"
echo "  100% consensus (best): $(wc -l < safe_harbors_final_100pct.bed)"
echo "  90% consensus: $(wc -l < safe_harbors_final_90pct.bed)"
echo "  80% consensus: $(wc -l < safe_harbors_final_80pct.bed)"

echo ""
echo "=== Safe harbor identification complete ==="
