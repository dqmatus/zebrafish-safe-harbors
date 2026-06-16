#!/bin/bash
# 03_process_annotations.sh
# Processes downloaded annotations into BED format

set -e

echo "=== Processing annotations ==="

cd annotation

# -----------------------------------------
# Process gene annotations
# -----------------------------------------
echo "Processing genes..."

# Extract all genes with gene names using perl
awk '$3 == "gene"' Danio_rerio.GRCz11.110.gtf | \
perl -ne '
    if (/^(\S+)\t\S+\tgene\t(\d+)\t(\d+)\t\S+\t([+-])/) {
        $chr = $1;
        $start = $2 - 1;
        $end = $3;
        $strand = $4;
        
        if (/gene_name "([^"]+)"/) {
            $name = $1;
        } elsif (/gene_id "([^"]+)"/) {
            $name = $1;
        } else {
            $name = "unknown";
        }
        
        print "$chr\t$start\t$end\t$name\t.\t$strand\n";
    }
' > genes.bed

echo "Found $(wc -l < genes.bed) genes"

# Add 10kb flanking regions
bedtools slop \
    -i genes.bed \
    -g ../reference/GRCz11.chrom.sizes \
    -b 10000 \
> genes_plus_10kb.bed

# Add 5kb flanking regions
bedtools slop \
    -i genes.bed \
    -g ../reference/GRCz11.chrom.sizes \
    -b 5000 \
> genes_plus_5kb.bed

# -----------------------------------------
# Process repeat annotations
# -----------------------------------------
echo "Processing repeats..."

# Convert RepeatMasker format to BED
# Emit repClass ($12) as column 5 in place of swScore ($2)
awk 'BEGIN{OFS="\t"} {
    print $6, $7, $8, $11, $12, $10
}' rmsk.txt \
| sort -k1,1 -k2,2n \
> repeats_all.bed

echo "Found $(wc -l < repeats_all.bed) repeat elements"

# Extract transposable elements (most problematic)
# Use repClass-based selector to match DNA|LINE|SINE|LTR|Retroposon|RC
awk 'BEGIN{OFS="\t"} $5 ~ /^(DNA|LINE|SINE|LTR|Retroposon|RC)\??$/ {print $1,$2,$3}' repeats_all.bed \
| sort -k1,1 -k2,2n \
| bedtools merge \
> repeats_transposons.bed

# -----------------------------------------
# Extract exons
# -----------------------------------------
echo "Extracting exons..."

awk '$3 == "exon"' Danio_rerio.GRCz11.110.gtf \
| grep -v "^#" \
| awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' \
| sort -k1,1 -k2,2n \
| bedtools merge \
> exons.bed

echo "=== Annotation processing complete ==="
