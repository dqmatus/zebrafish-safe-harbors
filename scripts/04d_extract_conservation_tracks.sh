#!/bin/bash
# 04d_extract_conservation_tracks.sh
# Build the two conservation exclusion tracks for 08b Step 2 (eLife GSH paper +
# collaborator review #3 P1.1): safe harbors should avoid evolutionarily
# constrained sequence, especially conserved NON-coding elements (regulatory).
#
# UCSC has NO phastCons/phyloP track for danRer11, so we use:
#   CONSERVATION_BED  <- Ensembl Compara 65-fish GERP constrained elements
#                        (danio_rerio, already on GRCz11; phastCons-equivalent)
#   CNE_BED           <- Ancora HCNEs (zebrafish vs human, the canonical
#                        conserved-non-coding set), lifted danRer10 -> danRer11
#
# Requires UCSC tools on PATH: bigBedToBed, liftOver (install from
# https://hgdownload.soe.ucsc.edu/admin/exe/<platform>/).
#
# Outputs (data BEDs gitignored; THIS SCRIPT is the reproducible artifact):
#   annotation/conservation/gerp_constrained_elements.bed
#   annotation/conservation/CNEs.bed

set -euo pipefail

for t in bigBedToBed liftOver bedtools; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not on PATH"; exit 1; }
done

OUT=annotation/conservation
mkdir -p "$OUT"
norm() { awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }'; }

# Tunables
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-110}"
HCNE_SET="${HCNE_SET:-70pc_50col}"     # Ancora stringency (identity_over_window)
HCNE_VS="${HCNE_VS:-hg38}"             # comparison species (human = deeply conserved)

GERP_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/bed/ensembl-compara/65_fish.gerp_constrained_element/gerp_constrained_elements.danio_rerio.bb"
HCNE_URL="http://ancora.genereg.net/downloads/danRer10/vs_human/HCNE_danRer10_${HCNE_VS}_${HCNE_SET}.bed.gz"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/liftOver/danRer10ToDanRer11.over.chain.gz"

echo "=== 1. Ensembl GERP constrained elements (GRCz11) ==="
curl -sL --fail --max-time 120 "$GERP_URL" -o "$OUT/gerp.bb"
bigBedToBed "$OUT/gerp.bb" "$OUT/gerp_raw.bed"
norm < "$OUT/gerp_raw.bed" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$OUT/gerp_constrained_elements.bed"
echo "  $(wc -l < "$OUT/gerp_constrained_elements.bed") elements, $(awk '{s+=$3-$2}END{print s}' "$OUT/gerp_constrained_elements.bed") bp"

echo "=== 2. Ancora HCNEs ($HCNE_VS $HCNE_SET) -> liftOver danRer10->danRer11 ==="
curl -sL --fail --max-time 60 "$HCNE_URL" -o "$OUT/hcne.bed.gz"
curl -sL --fail --max-time 60 "$CHAIN_URL" -o "$OUT/d10ToD11.chain.gz"
# Strip Ancora 'track'/'browser' header lines; keep integer-coord rows only.
gzip -dc "$OUT/hcne.bed.gz" | grep -vE '^(track|browser|#)' | awk -F'\t' '$2 ~ /^[0-9]+$/' \
    | norm | cut -f1-3 > "$OUT/hcne_d10.bed"
liftOver "$OUT/hcne_d10.bed" "$OUT/d10ToD11.chain.gz" "$OUT/CNEs.raw.bed" "$OUT/cne_unmapped.txt"
sort -k1,1 -k2,2n "$OUT/CNEs.raw.bed" | bedtools merge > "$OUT/CNEs.bed"
echo "  $(wc -l < "$OUT/hcne_d10.bed") danRer10 -> $(wc -l < "$OUT/CNEs.bed") merged danRer11 CNEs, $(awk '{s+=$3-$2}END{print s}' "$OUT/CNEs.bed") bp"

# Tidy intermediates (keep the two final BEDs)
rm -f "$OUT/gerp.bb" "$OUT/gerp_raw.bed" "$OUT/hcne.bed.gz" "$OUT/hcne_d10.bed" \
      "$OUT/d10ToD11.chain.gz" "$OUT/CNEs.raw.bed"

echo "=== Conservation tracks complete ==="
