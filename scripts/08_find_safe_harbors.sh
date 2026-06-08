#!/bin/bash
# 08_find_safe_harbors.sh
# Filters consensus peaks against exclusion regions to find safe harbors

set -e

echo "=== Finding Safe Harbor Candidates ==="

mkdir -p results
cd results

# Check what we're starting with
echo "Starting consensus peaks (100%): $(wc -l < consensus_peaks_100pct.bed)"
echo "Starting consensus peaks (90%): $(wc -l < consensus_peaks_90pct.bed)"
echo "Starting consensus peaks (80%): $(wc -l < consensus_peaks_80pct.bed)"

# Filter against minimal exclusions
echo ""
echo "=== Filtering against minimal exclusions ==="
echo "(exons + essential genes + transposons)"

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

# Filter by minimum size (need room for landing pad, at least 300bp)
echo ""
echo "=== Filtering by size (>=300bp) ==="

awk '($3 - $2) >= 300' safe_harbors_100pct_minimal.bed > safe_harbors_100pct_sized.bed
awk '($3 - $2) >= 300' safe_harbors_90pct_minimal.bed > safe_harbors_90pct_sized.bed
awk '($3 - $2) >= 300' safe_harbors_80pct_minimal.bed > safe_harbors_80pct_sized.bed

echo "After size filtering (>=300bp):"
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