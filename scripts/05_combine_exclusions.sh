#!/bin/bash
# 05_combine_exclusions.sh
# Combines all exclusion regions into master file

set -e

echo "=== Combining exclusion regions ==="

cd annotation
# Gene flanking buffers (kb). Strict = GENE_BUFFER_KB (default 10); a relaxed
# tier (RELAXED_GENE_BUFFER_KB, default 2) is also built for an alternate,
# more-permissive candidate list. (Previously 05 hardcoded genes_plus_10kb.bed
# while config's GENE_BUFFER_KB was unused — this wires the knob in.)
[ -f ../config.sh ] && source ../config.sh
GENE_BUFFER_KB=${GENE_BUFFER_KB:-10}
RELAXED_GENE_BUFFER_KB=${RELAXED_GENE_BUFFER_KB:-2}

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

# Create master exclusions (gene flank + repeats + essential genes + telomeres).
# NOTE: centromere_regions not included - current zebrafish centromere model
# (chromosome midpoint) is biologically wrong and over-excludes ~4 Mb mid-chromosome
# Helper: build master_exclusions for a given gene buffer (kb) -> output file.
build_master() {
    local buf_kb="$1" out="$2"
    bedtools slop -i genes.bed -g ../reference/GRCz11.chrom.sizes -b "$((buf_kb * 1000))" \
        > genes_plus_${buf_kb}kb.bed
    cat \
        genes_plus_${buf_kb}kb.bed \
        repeats_transposons.bed \
        essential_genes_plus_20kb.bed \
        telomere_regions.bed \
    | cut -f1-3 \
    | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > "$out"
}

# Strict (default 10kb) — the primary candidate list
build_master "$GENE_BUFFER_KB" master_exclusions.bed
# Relaxed (default 2kb) — alternate, more-permissive candidate list
build_master "$RELAXED_GENE_BUFFER_KB" master_exclusions_relaxed.bed

# Statistics
echo ""
echo "=== Exclusion Statistics ==="
echo "Genome size: $(awk '{sum+=$2}END{print sum}' ../reference/GRCz11.chrom.sizes) bp"
echo ""
echo "Minimal exclusions: $(awk '{sum+=$3-$2}END{print sum}' exclusions_minimal.bed) bp"
echo "Master exclusions (strict, ${GENE_BUFFER_KB}kb): $(awk '{sum+=$3-$2}END{print sum}' master_exclusions.bed) bp"
echo "Master exclusions (relaxed, ${RELAXED_GENE_BUFFER_KB}kb): $(awk '{sum+=$3-$2}END{print sum}' master_exclusions_relaxed.bed) bp"

echo "=== Exclusion combination complete ==="
