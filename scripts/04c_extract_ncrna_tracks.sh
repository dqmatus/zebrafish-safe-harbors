#!/bin/bash
# 04c_extract_ncrna_tracks.sh
# Extract small-ncRNA and lncRNA gene coordinates from the Ensembl GTF for the
# GSH distance criteria (08b Step 4), following the eLife GSH paper
# (Aznauryan et al. 2022): keep safe harbors away from miRNAs and lncRNAs.
#
# genes.bed already excludes ALL gene biotypes at +/-GENE_BUFFER_KB (10kb), so
# these tracks add the *extended* class-specific buffers (configured in config.sh:
# GSH_MIRNA_BUFFER_KB, GSH_LNCRNA_BUFFER_KB). Distances are scaled down from the
# paper's human values (300kb miRNA / 100kb lncRNA) for the ~1.4Gb zebrafish genome.
#
# Outputs (Ensembl chrom naming, matching genes.bed; 08b normalizes to chr):
#   annotation/mirna_ncrna_coordinates.bed  (miRNA + snoRNA + snRNA)
#   annotation/lncrna_coordinates.bed       (lincRNA + antisense + sense_intronic)

set -euo pipefail
GTF=annotation/Danio_rerio.GRCz11.110.gtf
[ -f "$GTF" ] || { echo "ERROR: $GTF not found (run 02/03 first)"; exit 1; }

echo "=== Extracting ncRNA tracks for GSH distance criteria ==="

# Small functional ncRNAs (miRNA / snoRNA / snRNA)
awk -F'\t' '$3=="gene" && $9 ~ /gene_biotype "(miRNA|snoRNA|snRNA)"/ {
    bt="ncRNA";
    if ($9 ~ /"miRNA"/)       bt="miRNA";
    else if ($9 ~ /"snoRNA"/) bt="snoRNA";
    else if ($9 ~ /"snRNA"/)  bt="snRNA";
    print $1"\t"($4-1)"\t"$5"\t"bt }' "$GTF" \
  | sort -k1,1 -k2,2n > annotation/mirna_ncrna_coordinates.bed
echo "  miRNA/snoRNA/snRNA: $(wc -l < annotation/mirna_ncrna_coordinates.bed) genes -> mirna_ncrna_coordinates.bed"

# lncRNAs (lincRNA / antisense / sense_intronic)
awk -F'\t' '$3=="gene" && $9 ~ /gene_biotype "(lincRNA|antisense|sense_intronic)"/ {
    print $1"\t"($4-1)"\t"$5"\tlncRNA" }' "$GTF" \
  | sort -k1,1 -k2,2n > annotation/lncrna_coordinates.bed
echo "  lncRNA: $(wc -l < annotation/lncrna_coordinates.bed) genes -> lncrna_coordinates.bed"

echo "=== ncRNA track extraction complete ==="
