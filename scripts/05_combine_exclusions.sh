#!/bin/bash
# 05_combine_exclusions.sh
# Combines all exclusion regions into master file

set -e

echo "=== Combining exclusion regions ==="

cd annotation

# Add flanking to essential genes
bedtools slop \
    -i essential_genes.bed \
    -g ../reference/GRCz11.chrom.sizes \
    -b 20000 \
> essential_genes_plus_20kb.bed

# Create minimal exclusions (exons + essential genes + transposons)
cat \
    exons.bed \
    essential_genes_plus_20kb.bed \
    repeats_transposons.bed \
| cut -f1-3 \
| sort -k1,1 -k2,2n \
| bedtools merge \
> exclusions_minimal.bed

# Create master exclusions (more stringent)
cat \
    genes_plus_10kb.bed \
    repeats_transposons.bed \
    essential_genes_plus_20kb.bed \
    telomere_regions.bed \
    centromere_regions.bed \
| cut -f1-3 \
| sort -k1,1 -k2,2n \
| bedtools merge \
> master_exclusions.bed

# Statistics
echo ""
echo "=== Exclusion Statistics ==="
echo "Genome size: $(awk '{sum+=$2}END{print sum}' ../reference/GRCz11.chrom.sizes) bp"
echo ""
echo "Minimal exclusions: $(awk '{sum+=$3-$2}END{print sum}' exclusions_minimal.bed) bp"
echo "Master exclusions: $(awk '{sum+=$3-$2}END{print sum}' master_exclusions.bed) bp"

echo "=== Exclusion combination complete ==="