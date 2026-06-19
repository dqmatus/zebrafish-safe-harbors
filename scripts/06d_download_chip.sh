#!/bin/bash
# 06d_download_chip.sh
# Downloads DANIO-CODE ChIP-seq narrowPeak peak files for the Phase 3 regulatory
# screen: H3K4me1, H3K4me3 (Skarmeta Lab, same series as the existing H3K27ac),
# and CTCF (Horsfield Lab — only "Bud" stage exists in all of DANIO-CODE).
#
# NOTE: the DCC file host migrated. Files are now served from
#   https://danio-code-dcc.genereg.net/files/annotated_files/ChIP-seq/<BS>/<file>
# (the old danio-code.zfin.org/files/... path now 301-redirects to a homepage).
# H3K27ac was downloaded previously and already lives in annotation/chip/.

set -euo pipefail

BASE="https://danio-code-dcc.genereg.net/files/annotated_files/ChIP-seq"

mkdir -p annotation/chip
cd annotation/chip

echo "=== Downloading ChIP-seq peaks (H3K4me1, H3K4me3, CTCF) ==="

# --- H3K4me1 (Skarmeta Lab, USERpanosfirbas, pval1e-07.500K) — 4 stages ---
echo "H3K4me1 (Skarmeta, 4 stages)..."
for u in \
  "DCD006178BS/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000657SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006177BS/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000656SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006175BS/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000655SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006173BS/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000650SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" ; do
    curl -fsSL -O "$BASE/$u"
done

# --- H3K4me3 (Skarmeta Lab, USERpanosfirbas, pval1e-07.500K) — 6 stages, 2 reps ---
echo "H3K4me3 (Skarmeta, 6 stages x reps)..."
for u in \
  "DCD006178BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000647SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006178BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000649SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006177BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000646SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006177BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000643SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006176BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000648SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006176BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000645SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006175BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000642SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006175BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000639SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006174BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000644SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006174BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000641SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006173BS/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000640SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" ; do
    curl -fsSL -O "$BASE/$u"
done

# --- CTCF (Horsfield Lab, USERmickael.dong, control-subtracted) — Bud stage only ---
echo "CTCF (Horsfield, Bud stage only)..."
for u in \
  "DCD007229BS/ChIP-seq_Horsfield_Lab_CTCF_0004AS.DCD002910SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep1.pval0.01.500K.narrowPeak.gz" \
  "DCD007229BS/ChIP-seq_Horsfield_Lab_CTCF_0004AS.DCD002891SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep2.pval0.01.500K.narrowPeak.gz" ; do
    curl -fsSL -O "$BASE/$u"
done

# --- Additional EMBRYONIC-stage ChIP from other DANIO-CODE labs ---------------
# Cairns / Mueller / Schier / Bartel / Zon labs (embryonic stages only; adult
# "90 Days-2 Years" deliberately excluded to match the 0-48 hpf ATAC window).
# Note: Zon's H3K27ac file is named "H3K27Ac" (capital A) — 06c matches case-
# insensitively so it still merges into the H3K27ac union.
echo "Other-lab embryonic H3K4me1/H3K4me3/H3K27ac..."
for u in \
  "DCD007130BS/ChIP-seq_Cairns_Lab_H3K4me1_0002AS.DCD002675SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007130BS/ChIP-seq_Cairns_Lab_H3K4me1_0002AS.DCD002679SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K4me1_0002AS.DCD002714SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K4me1_0002AS.DCD002666SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007025BS/ChIP-seq_Bartel_Lab_H3K4me3_0001AS.DCD002371SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007026BS/ChIP-seq_Bartel_Lab_H3K4me3_0001AS.DCD002380SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD006678BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001524SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006677BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001529SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006675BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001525SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006675BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001540SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006673BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001530SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006673BS/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001517SQ.danRer11.USERdunjanik.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD007150BS/ChIP-seq_Schier_Lab_H3K4me3_0002AS.DCD002755SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006634BS/ChIP-seq_Schier_Lab_H3K4me3_0002AS.DCD002758SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006633BS/ChIP-seq_Schier_Lab_H3K4me3_0002AS.DCD002752SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006562BS/ChIP-seq_Zon_Lab_H3K4me3_0002AS.DCD001005SQ.danRer11.USERsyang.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K27ac_0002AS.DCD002644SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K27ac_0002AS.DCD002715SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD006562BS/ChIP-seq_Zon_Lab_H3K27Ac_0002AS.DCD001011SQ.danRer11.USERsyang.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" ; do
    curl -fsSL -O "$BASE/$u"
done

# --- H3K27me3 (Polycomb-silenced; repressive-mark exclusion) -----------------
# Embryonic only (Skarmeta + Schier; adult excluded). Used to avoid landing in
# silenced chromatin where a transgene could itself be silenced.
echo "H3K27me3 (Skarmeta + Schier, embryonic)..."
for u in \
  "DCD006633BS/ChIP-seq_Schier_Lab_H3K27me3_0002AS.DCD002753SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006178BS/ChIP-seq_Skarmeta_Lab_H3K27me3_0002AS.DCD000631SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006173BS/ChIP-seq_Skarmeta_Lab_H3K27me3_0002AS.DCD000630SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006175BS/ChIP-seq_Skarmeta_Lab_H3K27me3_0002AS.DCD000632SQ.danRer11.USERpanosfirbas.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" ; do
    curl -fsSL -O "$BASE/$u"
done

# --- Other ChIP targets (H3K36me3, H2A.Z/H2AFV, H3K14ac, Pol II, dev TFs) -----
# Embryonic only; 06c merges these into other_chip_peaks.bed (08b Step 1d).
echo "Other active/TF marks (H3K36me3/H2A.Z/H3K14ac/PolII/TFs)..."
for u in \
  "DCD007132BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002652SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep3.pval1e-07.500K.narrowPeak.gz" \
  "DCD007132BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002656SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007132BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002669SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007129BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002641SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007129BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002711SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007128BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002697SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep3.pval1e-07.500K.narrowPeak.gz" \
  "DCD007128BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002704SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007128BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002712SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007127BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002710SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007127BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002699SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep3.pval1e-07.500K.narrowPeak.gz" \
  "DCD007127BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002680SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002676SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002721SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007123BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002684SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007123BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002700SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep3.pval1e-07.500K.narrowPeak.gz" \
  "DCD007123BS/ChIP-seq_Cairns_Lab_H2AFV_0002AS.DCD002657SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K14ac_0002AS.DCD002658SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K14ac_0002AS.DCD002683SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep2.pval1e-07.500K.narrowPeak.gz" \
  "DCD007124BS/ChIP-seq_Cairns_Lab_H3K14ac_0002AS.DCD002702SQ.danRer11.USERdanio-user.R1.merged.nodup_x_ctl_for_rep3.pval1e-07.500K.narrowPeak.gz" \
  "DCD007025BS/ChIP-seq_Bartel_Lab_H3K36me3_0001AS.DCD002372SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD007026BS/ChIP-seq_Bartel_Lab_H3K36me3_0001AS.DCD002381SQ.danRer11.USERmickael.dong.R1.merged.nodup_x_ctl_for_rep1.pval1e-07.500K.narrowPeak.gz" \
  "DCD006633BS/ChIP-seq_Schier_Lab_H3K36me3_0002AS.DCD002754SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006633BS/ChIP-seq_Schier_Lab_Pol_II_4H8_0001AS.DCD002757SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  "DCD006633BS/ChIP-seq_Schier_Lab_Pol_II_8WG16_0001AS.DCD002751SQ.danRer11.USERdanio-user.R1.merged.nodup.pval1e-07.500K.narrowPeak.gz" \
  ; do
    curl -fsSL -O "$BASE/$u"
done

echo ""
echo "=== Download complete ==="
echo "H3K4me1: $(ls *H3K4me1*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "H3K4me3: $(ls *H3K4me3*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "H3K27ac: $(ls *H3K27[Aa]c*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "H3K27me3:$(ls *H3K27me3*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "CTCF:    $(ls *CTCF*.narrowPeak.gz 2>/dev/null | wc -l) files"
