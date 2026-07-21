#!/bin/bash
# 26_make_igv_batch.sh
# Generate an IGV batch script that loads the safe-harbor evidence tracks
# (constitutive ATAC + candidate + the ChIP/regulatory/conservation exclusion
# tracks) and snapshots one or more candidate loci — showing an ATAC peak
# present over the candidate with every excluded track blank underneath.
#
# Usage:
#   scripts/26_make_igv_batch.sh                     # top 3 candidates (by composite rank)
#   scripts/26_make_igv_batch.sh chr12:24289888-24291839 chr6:35551173-35552810
#   scripts/26_make_igv_batch.sh cope_chr1           # by locus name in the donor design
#
# Then in IGV:  Tools -> Run Batch Script -> results/igv/safe_harbors.bat
#          or:  igv.sh -b results/igv/safe_harbors.bat   (headless)
# Snapshots land in results/igv/*.png

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$REPO/results/igv"
BAT="$OUTDIR/safe_harbors.bat"
PAD=4000                 # bp of flank added around each candidate (gene-desert context)
mkdir -p "$OUTDIR"

# Excluded-region foils for side-by-side comparison: intergenic ATAC peaks that
# carry active chromatin marks (H3K27ac+H3K4me1+H3K4me3) — accessible like a safe
# harbor but removed by the ChIP screen. "name|chrom:start-end".
FOILS=(
  "EXCLUDED_enhancer_sema4c|chr5:14843192-14844205"    # 195kb from sema4c; H3K27ac+H3K4me1+H3K4me3
  "EXCLUDED_enhancer_SH2D4B|chr13:21319046-21319511"   # 281kb from SH2D4B; all three marks
)

# ---- tracks (order = top->bottom in IGV); only load those that exist ----
TRACKS=(
  "results/consensus_peaks_100pct.bed"                 # ATAC: open in 100% of stages
  "results/safe_harbors_biological_100pct_sized.bed"   # the candidate call
  "annotation/chip/H3K27ac_peaks.bed"
  "annotation/chip/H3K4me1_peaks.bed"
  "annotation/chip/H3K4me3_peaks.bed"
  "annotation/chip/CTCF_peaks.bed"
  "annotation/chip/H3K27me3_peaks.bed"
  "annotation/chip/other_chip_peaks.bed"
  "annotation/regulatory/validated_enhancers.bed"
  "annotation/regulatory/promoters.bed"
  "annotation/conservation/gerp_constrained_elements.bed"
  "annotation/conservation/CNEs.bed"
)

# ---- resolve requested loci -> "name<TAB>chrom:start-end" ----
resolve_loci() {
  if [ "$#" -eq 0 ]; then
    # top 3 candidates on real chromosomes, by composite rank (col 15)
    tail -n +2 "$REPO/results/candidate_feature_table.tsv" \
      | awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-5])$/' \
      | sort -t$'\t' -k15,15 -gr | head -3 \
      | awk -F'\t' '{gsub(/[^A-Za-z0-9]/,"_",$6); print $6"_"$1"\t"$1":"$2"-"$3}'
  else
    for q in "$@"; do
      if [[ "$q" == *:*-* ]]; then
        echo -e "$(echo "$q" | tr ':-' '__')\t$q"
      else
        # locus name from the donor design (e.g. cope_chr1)
        awk -F'\t' -v s="$q" '$1==s{print s"\t"$3":"($7~/:/?"":"")""}' "$REPO/results/knockin_donor_design.tsv" >/dev/null 2>&1
        row=$(awk -F'\t' -v s="$q" '$1==s{print $3"\t"$7}' "$REPO/results/knockin_donor_design.tsv")
        chrom=$(echo "$row" | cut -f1); cut=$(echo "$row" | cut -f2 | cut -d: -f2)
        [ -n "$cut" ] && echo -e "${q}\t${chrom}:$((cut-800))-$((cut+800))" || echo "  (locus $q not found)" >&2
      fi
    done
  fi
}

{
  echo "new"
  echo "genome danRer11"
  echo "snapshotDirectory $OUTDIR"
  echo "maxPanelHeight 2000"
  for t in "${TRACKS[@]}"; do
    [ -f "$REPO/$t" ] && echo "load $REPO/$t"
  done
  # keep the excluded tracks compact so blank rows read cleanly
  for t in H3K27ac_peaks H3K4me1_peaks H3K4me3_peaks CTCF_peaks H3K27me3_peaks other_chip_peaks \
           validated_enhancers promoters gerp_constrained_elements CNEs; do
    echo "squish $t"
  done
  while IFS=$'\t' read -r name loc; do
    [ -z "${loc:-}" ] && continue
    chrom=${loc%%:*}; rest=${loc#*:}; s=${rest%-*}; e=${rest#*-}
    echo "goto ${chrom}:$((s-PAD))-$((e+PAD))"
    echo "snapshot safe_harbor_${name}.png"
  done < <(resolve_loci "$@")
  # excluded-region foils (comparison: same accessible/intergenic context, WITH ChIP marks)
  for fo in "${FOILS[@]}"; do
    fname=${fo%%|*}; floc=${fo#*|}
    fchrom=${floc%%:*}; frest=${floc#*:}; fs=${frest%-*}; fe=${frest#*-}
    echo "goto ${fchrom}:$((fs-PAD))-$((fe+PAD))"
    echo "snapshot ${fname}.png"
  done
  echo "exit"   # close IGV after snapshots (remove this line to keep IGV open for interactive use)
} > "$BAT"

echo "wrote $BAT"
echo "  tracks: $(grep -c '^load ' "$BAT") | snapshots: $(grep -c '^snapshot ' "$BAT")"
echo "  run in IGV: Tools -> Run Batch Script -> $BAT   (or: igv.sh -b $BAT)"
grep -E '^goto|^snapshot' "$BAT" | sed 's/^/    /'
