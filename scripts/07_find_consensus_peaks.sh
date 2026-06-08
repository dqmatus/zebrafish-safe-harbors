#!/bin/bash
# 07_find_consensus_peaks.sh
# Identifies peaks present across multiple ATAC-seq samples

set -e

echo "=== Finding Consensus ATAC-seq Peaks ==="

cd peaks

# Count samples
num_samples=$(ls *.gz | wc -l | tr -d ' ')
echo "Number of samples: $num_samples"

# Decompress all files
echo "Decompressing peak files..."
for f in *.gz; do
    gunzip -k "$f" 2>/dev/null || true
done

# Combine all peaks
echo "Combining all peaks..."
cat *.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n > all_peaks_combined.bed
echo "Total peak entries: $(wc -l < all_peaks_combined.bed)"

# Merge overlapping peaks to create reference regions
echo "Merging overlapping regions..."
bedtools merge -i all_peaks_combined.bed > merged_peak_regions.bed
echo "Merged regions: $(wc -l < merged_peak_regions.bed)"

# Count how many samples overlap each merged region
echo "Counting sample overlaps (this may take a few minutes)..."

# Create a file for each sample's peaks
mkdir -p sample_beds
for f in *.narrowPeak; do
    sample=$(basename "$f" .narrowPeak)
    cut -f1-3 "$f" | sort -k1,1 -k2,2n > "sample_beds/${sample}.bed"
done

# Use bedtools to count overlaps
echo "Calculating overlap counts..."
bedtools annotate -i merged_peak_regions.bed -files sample_beds/*.bed -counts > peak_overlap_counts.txt

# Find consensus peaks (present in >80% of samples)
min_samples=$(echo "$num_samples * 0.8" | bc | cut -d. -f1)
echo "Requiring overlap in at least $min_samples samples (80%)"

# Sum the overlap counts and filter
awk -v min="$min_samples" -v ncol="$((num_samples + 3))" '
BEGIN {OFS="\t"}
{
    sum = 0
    for (i = 4; i <= NF; i++) {
        if ($i > 0) sum++
    }
    if (sum >= min) {
        print $1, $2, $3, sum
    }
}' peak_overlap_counts.txt > consensus_peaks_80pct.bed

echo ""
echo "=== Results ==="
echo "Consensus peaks (>=80% samples): $(wc -l < consensus_peaks_80pct.bed)"

# Also create 90% and 100% consensus sets
min_90=$(echo "$num_samples * 0.9" | bc | cut -d. -f1)
awk -v min="$min_90" '{
    sum = 0
    for (i = 4; i <= NF; i++) if ($i > 0) sum++
    if (sum >= min) print $1, $2, $3, sum
}' OFS='\t' peak_overlap_counts.txt > consensus_peaks_90pct.bed

awk -v min="$num_samples" '{
    sum = 0
    for (i = 4; i <= NF; i++) if ($i > 0) sum++
    if (sum >= min) print $1, $2, $3, sum
}' OFS='\t' peak_overlap_counts.txt > consensus_peaks_100pct.bed

echo "Consensus peaks (>=90% samples): $(wc -l < consensus_peaks_90pct.bed)"
echo "Consensus peaks (100% samples): $(wc -l < consensus_peaks_100pct.bed)"

# Copy results to results folder
mkdir -p ../results
cp consensus_peaks_*.bed ../results/

echo ""
echo "=== Consensus peak finding complete ==="