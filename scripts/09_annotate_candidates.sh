#!/bin/bash
# 09_annotate_candidates.sh
# Annotates safe harbor candidates with distance to nearest genes

set -euo pipefail

echo "=== Annotating Safe Harbor Candidates ==="

cd annotation

# Add chr prefix to genes for compatibility with ATAC-seq peak files
awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, $4, $5, $6}' genes.bed | sort -k1,1 -k2,2n > genes_chr_sorted.bed

cd ../results

# Get top 50 candidates
head -50 safe_harbors_final_100pct.bed > top50_candidates.bed

# Sort for bedtools
sort -k1,1 -k2,2n top50_candidates.bed > top50_sorted.bed

# Find nearest genes
echo "Finding nearest genes..."
bedtools closest \
    -a top50_sorted.bed \
    -b ../annotation/genes_chr_sorted.bed \
    -d \
    > top50_nearest_genes.txt

# Display results
echo ""
echo "=== Top Candidates by Distance to Nearest Gene ==="
echo ""
awk 'BEGIN{OFS="\t"; print "chrom","start","end","size","nearest_gene","distance"} 
     $NF>0 {print $1,$2,$3,$4,$9,$NF}' top50_nearest_genes.txt \
| sort -k6 -nr \
| head -20

# Save formatted results
awk 'BEGIN{OFS="\t"; print "chrom","start","end","size","num_samples","nearest_gene","distance"} 
     {print $1,$2,$3,$4,$5,$9,$NF}' top50_nearest_genes.txt \
> top50_candidates_annotated.txt

echo ""
echo "Full annotated list saved to: results/top50_candidates_annotated.txt"
echo ""
echo "=== Annotation complete ==="
