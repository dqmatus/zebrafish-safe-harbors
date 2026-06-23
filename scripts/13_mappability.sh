#!/bin/bash
# 13_mappability.sh
# Genome K-mer mappability mask (eLife GSH paper criterion + general uniqueness
# backstop). Safe harbors must be UNIQUELY targetable: a region whose K-mers
# recur elsewhere in the genome is a poor landing pad (ambiguous CRISPR cutting,
# unmappable genotyping/sequencing reads, cryptic segmental duplication not
# caught by RepeatMasker).
#
# Uses GenMap (Pockrandt et al. 2020) to score every position by the uniqueness
# of the K-mer starting there: mappability = 1 / (# genome occurrences). Regions
# scoring below MAPPABILITY_THRESHOLD are written to GRCz11_mappability_low.bed,
# which 08b Step 3 then removes (discard-whole).
#
# Parameters (config.sh): GSH_MAPPABILITY_K, GSH_MAPPABILITY_E, MAPPABILITY_THRESHOLD.
# The FM-index is ~4GB and is deleted on completion unless KEEP_GENMAP_INDEX=1.
#
# Output (data BED is gitignored; THIS SCRIPT is the reproducible artifact):
#   annotation/mappability/GRCz11_mappability_low.bed   (chr-prefixed, merged)

set -euo pipefail

FASTA=reference/GRCz11.fa
[ -f "$FASTA" ] || { echo "ERROR: $FASTA not found"; exit 1; }

K="${GSH_MAPPABILITY_K:-50}"
E="${GSH_MAPPABILITY_E:-0}"
THRESH="${MAPPABILITY_THRESHOLD:-0.8}"
IDX="${GENMAP_INDEX_DIR:-/tmp/genmap_GRCz11_idx}"
OUTDIR=annotation/mappability
PREFIX="$OUTDIR/GRCz11_genmap_K${K}_E${E}"
LOW=annotation/mappability/GRCz11_mappability_low.bed

mkdir -p "$OUTDIR"

echo "=== GenMap mappability (K=$K, E=$E, threshold=$THRESH) ==="

# 1. Build the FM-index. GenMap REQUIRES the output dir to NOT already exist.
if [ -f "$IDX/index.txt.size" ] || [ -f "$IDX/index.ids.concat" ]; then
    echo "  Reusing existing index at $IDX"
else
    rm -rf "$IDX"
    echo "  Building FM-index -> $IDX (~7 min, ~4GB)..."
    genmap index -F "$FASTA" -I "$IDX"
fi

# 2. Compute per-base mappability as a bedgraph.
echo "  Computing mappability bedgraph..."
genmap map -K "$K" -E "$E" -I "$IDX" -O "$PREFIX" -bg -T 4

# 3. Threshold: keep runs scoring BELOW THRESH (= low mappability = to exclude),
#    normalize to chr-prefixed, merge.
echo "  Thresholding < $THRESH and merging..."
awk -v t="$THRESH" 'BEGIN{OFS="\t"} $4 < t { if ($1 !~ /^chr/) $1="chr"$1; print $1,$2,$3 }' \
    "${PREFIX}.bedgraph" \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "$LOW"

echo "  Low-mappability regions: $(wc -l < "$LOW") (-> $LOW)"
echo "  Genome bp masked: $(awk '{s+=$3-$2} END{print s+0}' "$LOW")"

# 4. Reclaim disk: drop the index and the large bedgraph unless asked to keep.
if [ "${KEEP_GENMAP_INDEX:-0}" != "1" ]; then
    rm -rf "$IDX"
    echo "  Removed FM-index to reclaim disk."
fi
[ "${KEEP_GENMAP_BEDGRAPH:-0}" = "1" ] || rm -f "${PREFIX}.bedgraph"

echo "=== Mappability mask complete ==="
