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

echo ""
echo "=== Download complete ==="
echo "H3K4me1: $(ls *H3K4me1*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "H3K4me3: $(ls *H3K4me3*.narrowPeak.gz 2>/dev/null | wc -l) files"
echo "CTCF:    $(ls *CTCF*.narrowPeak.gz 2>/dev/null | wc -l) files"
