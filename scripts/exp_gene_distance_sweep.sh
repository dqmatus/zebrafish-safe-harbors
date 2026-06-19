#!/bin/bash
# exp_gene_distance_sweep.sh  (branch: relaxed-gene-distance)
#
# Sensitivity analysis: how does the gene-distance buffer affect the inert
# candidate count, and does relaxing it surface validated landing sites
# (tpl102, pIGLET14a/24b)?
#
# Holds the inertness (ChIP/regulatory) and repeat/essential/telomere filters
# fixed; sweeps only the gene flanking buffer. Prints, per buffer:
#   - candidates after gene-distal + size filter (08-equivalent)
#   - candidates after the full regulatory (Phase 3 discard-whole) screen
#   - whether each validated site and chr12/nrxn1a is IN the final list.
#
# Finding (2026-06): relaxing 10kb -> gene-body-only adds ~875 inert candidates
# (+38%), but recovers NONE of the validated sites — they fail on repeat overlap
# (tpl102, 392bp), gene-body insertion (pIGLET24b, in skilb), or active marks.

set -euo pipefail
GS=reference/GRCz11.chrom.sizes
norm(){ awk 'BEGIN{OFS="\t"}{if($1!~/^chr/)$1="chr"$1;print}'; }

# Combined regulatory exclusion: all Phase 3 discard-whole tracks.
cat annotation/chip/H3K27ac_peaks.bed annotation/chip/H3K4me1_peaks.bed \
    annotation/chip/H3K4me3_peaks.bed annotation/chip/H3K27me3_peaks.bed \
    annotation/chip/CTCF_peaks.bed annotation/chip/other_chip_peaks.bed \
    annotation/regulatory/validated_enhancers.bed annotation/regulatory/promoters.bed \
    annotation/regulatory/novel_lncrna.bed \
  | cut -f1-3 | norm | sort -k1,1 -k2,2n | bedtools merge > /tmp/REG.bed

# Fixed (non-gene) exclusions: repeats + essential genes +/-20kb + telomeres.
cat annotation/repeats_transposons.bed annotation/essential_genes_plus_20kb.bed \
    annotation/telomere_regions.bed \
  | cut -f1-3 | norm | sort -k1,1 -k2,2n | bedtools merge > /tmp/FIXED.bed

printf 'chr3\t45364050\t45364508\n'   > /tmp/tpl.bed       # tpl102
printf 'chr14\t32825753\t32826753\n'  > /tmp/p14.bed       # pIGLET14a
printf 'chr24\t26386249\t26387249\n'  > /tmp/p24.bed       # pIGLET24b
printf 'chr12\t24289888\t24291839\n'  > /tmp/c12.bed       # chr12/nrxn1a top hit

printf "%-10s %14s %12s  %-7s %-9s %-9s %-6s\n" \
  "gene_buf" "preP3(gene+sz)" "FINAL_inert" "tpl102" "pIGLET14a" "pIGLET24b" "chr12"
for B in 10000 5000 2000 1000 0; do
  if [ "$B" -eq 0 ]; then gb=annotation/genes.bed
  else bedtools slop -i annotation/genes.bed -g "$GS" -b "$B" > /tmp/gb.bed; gb=/tmp/gb.bed; fi
  cat "$gb" /tmp/FIXED.bed | cut -f1-3 | norm | sort -k1,1 -k2,2n | bedtools merge > /tmp/MASTER.bed
  bedtools subtract -a results/consensus_peaks_100pct.bed -b /tmp/MASTER.bed \
    | awk '($3-$2)>=300' | sort -k1,1 -k2,2n | bedtools merge > /tmp/preP3.bed
  bedtools intersect -a /tmp/preP3.bed -b /tmp/REG.bed -v > /tmp/final.bed
  IN(){ [ "$(bedtools intersect -a "$1" -b /tmp/final.bed -u | wc -l | tr -d ' ')" -gt 0 ] && echo IN || echo out; }
  lbl=$([ "$B" -eq 0 ] && echo "0(body)" || echo "$((B/1000))kb")
  printf "%-10s %14s %12s  %-7s %-9s %-9s %-6s\n" "$lbl" \
    "$(wc -l </tmp/preP3.bed|tr -d ' ')" "$(wc -l </tmp/final.bed|tr -d ' ')" \
    "$(IN /tmp/tpl.bed)" "$(IN /tmp/p14.bed)" "$(IN /tmp/p24.bed)" "$(IN /tmp/c12.bed)"
done
