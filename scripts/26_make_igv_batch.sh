#!/bin/bash
# 26_make_igv_batch.sh
# Generate a STYLED IGV batch that snapshots safe-harbor candidates next to
# "accessible-but-excluded" foils, with clean track names, a color scheme, and a
# uniform window — publication-ready panels showing an ATAC peak present over the
# candidate with every excluded track blank underneath (and filled at the foils).
#
# Usage:
#   scripts/26_make_igv_batch.sh                          # top 3 candidates + foils
#   scripts/26_make_igv_batch.sh chr12:24289888-24291839  # explicit loci
#   scripts/26_make_igv_batch.sh cope_chr1                # donor-design locus name
#
# Run: the command printed at the end (uses your installed IGV.app + bundled JDK).
# Snapshots + styled tracks land in results/igv/.

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$REPO/results/igv"; TDIR="$OUTDIR/tracks"
BAT="$OUTDIR/safe_harbors.bat"
WINDOW=10000             # fixed panel width (bp) for the peak/mark comparison
WIDE=120000              # wide gene-context panel (bp) per safe harbor -> shows the gene desert
mkdir -p "$TDIR"

# Excluded-region foils: intergenic ATAC peaks WITH active marks (accessible like a
# safe harbor, removed by the ChIP screen). "name|chrom:start-end".
FOILS=(
  "EXCLUDED_enhancer_sema4c|chr5:14843192-14844205"
  "EXCLUDED_enhancer_SH2D4B|chr13:21319046-21319511"
)

# track styling: "relpath | display name | R,G,B"  (order = top->bottom)
TRACKSPEC=(
  "results/consensus_peaks_100pct.bed|ATAC 100% consensus|33,102,172"
  "results/safe_harbors_biological_100pct_sized.bed|safe-harbor candidate|35,139,69"
  "annotation/chip/H3K27ac_peaks.bed|H3K27ac|202,0,32"
  "annotation/chip/H3K4me1_peaks.bed|H3K4me1|244,109,67"
  "annotation/chip/H3K4me3_peaks.bed|H3K4me3|178,24,43"
  "annotation/chip/CTCF_peaks.bed|CTCF|118,42,131"
  "annotation/chip/H3K27me3_peaks.bed|H3K27me3 (Polycomb)|128,128,128"
  "annotation/chip/other_chip_peaks.bed|other ChIP|140,81,10"
  "annotation/regulatory/validated_enhancers.bed|validated enhancers|230,97,1"
  "annotation/regulatory/promoters.bed|promoters|166,54,3"
  "annotation/conservation/gerp_constrained_elements.bed|conservation (GERP)|27,158,119"
  "annotation/conservation/CNEs.bed|CNEs|102,194,165"
)

# ---- stage styled copies (UCSC track line: name + color, collapsed) ----
STYLED=()
i=0
for spec in "${TRACKSPEC[@]}"; do
  IFS='|' read -r rel name rgb <<< "$spec"
  [ -f "$REPO/$rel" ] || continue
  out="$TDIR/$(printf '%02d' "$i")_$(basename "$rel")"
  { echo "track name=\"$name\" color=$rgb visibility=1 itemRgb=Off"; cut -f1-4 "$REPO/$rel"; } > "$out"
  STYLED+=("$out"); i=$((i+1))
done

# ---- resolve requested loci -> "name<TAB>chrom:start-end" ----
resolve_loci() {
  if [ "$#" -eq 0 ]; then
    tail -n +2 "$REPO/results/candidate_feature_table.tsv" \
      | awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-5])$/' \
      | sort -t$'\t' -k15,15 -gr | head -3 \
      | awk -F'\t' '{gsub(/[^A-Za-z0-9]/,"_",$6); print $6"_"$1"\t"$1":"$2"-"$3}'
  else
    for q in "$@"; do
      if [[ "$q" == *:*-* ]]; then
        echo -e "$(echo "$q" | tr ':-' '__')\t$q"
      else
        row=$(awk -F'\t' -v s="$q" '$1==s{print $3"\t"$7}' "$REPO/results/knockin_donor_design.tsv")
        chrom=$(echo "$row" | cut -f1); cut=$(echo "$row" | cut -f2 | cut -d: -f2)
        [ -n "$cut" ] && echo -e "${q}\t${chrom}:$((cut-800))-$((cut+800))" || echo "  (locus $q not found)" >&2
      fi
    done
  fi
}

# fixed-width view centered on a locus: win_goto chrom start end [width]
win_goto() {
  local chrom=$1 s=$2 e=$3 w=${4:-$WINDOW} c
  c=$(( (s + e) / 2 ))
  echo "goto ${chrom}:$(( c - w/2 ))-$(( c + w/2 ))"
}

{
  echo "new"
  echo "genome danRer11"
  echo "snapshotDirectory $OUTDIR"
  echo "maxPanelHeight 2000"
  for t in "${STYLED[@]}"; do echo "load $t"; done
  while IFS=$'\t' read -r name loc; do
    [ -z "${loc:-}" ] && continue
    chrom=${loc%%:*}; rest=${loc#*:}; s=${rest%-*}; e=${rest#*-}
    win_goto "$chrom" "$s" "$e"            # 10kb peak/mark comparison
    echo "snapshot safe_harbor_${name}.png"
    win_goto "$chrom" "$s" "$e" "$WIDE"    # wide gene-context (shows the gene desert)
    echo "snapshot safe_harbor_${name}_genecontext.png"
  done < <(resolve_loci "$@")
  for fo in "${FOILS[@]}"; do
    fname=${fo%%|*}; floc=${fo#*|}; fchrom=${floc%%:*}; frest=${floc#*:}
    win_goto "$fchrom" "${frest%-*}" "${frest#*-}"
    echo "snapshot ${fname}.png"
  done
  echo "exit"   # remove to keep IGV open for interactive use
} > "$BAT"

echo "wrote $BAT  (styled tracks in $TDIR)"
echo "  ${#STYLED[@]} tracks | $(grep -c '^snapshot ' "$BAT") snapshots | ${WINDOW}bp uniform window"
grep -E '^snapshot' "$BAT" | sed 's/^/    /'
APP=$(ls -d /Applications/IGV_*.app 2>/dev/null | sort -V | tail -1)
if [ -n "$APP" ]; then
  JV=$(find "$APP/Contents" -maxdepth 3 -path '*/bin/java' -type f 2>/dev/null | head -1)
  echo ""
  echo "  run from terminal:"
  echo "    \"$JV\" @\"$APP/Contents/Java/igv.args\" -Xmx4g --module-path=\"$APP/Contents/Java/lib\" \\"
  echo "      --module=org.igv/org.broad.igv.ui.Main -b \"$BAT\""
fi
