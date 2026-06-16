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
# FIX A1: Add chromosome name normalizer to ensure 'chr' prefix consistency
# (ATAC/UCSC repeats use chr1, but Ensembl-derived BEDs use 1)
cat \
    exons.bed \
    essential_genes_plus_20kb.bed \
    repeats_transposons.bed \
| cut -f1-3 \
| awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
| sort -k1,1 -k2,2n \
| bedtools merge \
> exclusions_minimal.bed

# Create master exclusions (more stringent)
# FIX A2: Add genes_plus_10kb and telomere_regions to minimal exclusions
# (Previously only exons + essential genes + transposons were excluded,
#  failing to enforce "intergenic / gene-distal" criterion)
# NOTE: centromere_regions not included - current zebrafish centromere model
# (chromosome midpoint) is biologically wrong and over-excludes ~4 Mb mid-chromosome
cat \
    genes_plus_10kb.bed \
    repeats_transposons.bed \
    essential_genes_plus_20kb.bed \
    telomere_regions.bed \
| cut -f1-3 \
| awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
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
