#!/bin/bash
# 04e_extract_cope_dope.sh
# Download the DANIO-CODE COPE/DOPE element tracks (Baranasic et al. 2022) for the
# positive "open but inert" cross-reference flagged in the candidate table by 09
# (Step 7b) and 08b Step 6.
#
# COPE = Constitutive Orphan open-chromatin Promoter Elements: regions with an
# ATAC-seq signal but WITHOUT CRE-associated histone marks ("open but inert").
# DOPE = the Developmentally regulated (pooled) counterpart.
#
# Source: DANIO-CODE UCSC track hub (danRer11). NB the zfin.org host serves HTML
# for the hub path; use the genereg.net mirror, which serves the bigBed files.
#
# Requires bigBedToBed on PATH (UCSC; https://hgdownload.soe.ucsc.edu/admin/exe/).
#
# Output (data BEDs gitignored; THIS SCRIPT is the reproducible artifact):
#   annotation/cope_dope/COPE_regions.bed
#   annotation/cope_dope/DOPE_regions.bed

set -euo pipefail
command -v bigBedToBed >/dev/null 2>&1 || { echo "ERROR: bigBedToBed not on PATH"; exit 1; }

OUT=annotation/cope_dope
mkdir -p "$OUT"
HUB="https://danio-code-dcc.genereg.net/trackhub/danRer11"
norm() { awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }'; }

fetch() {  # $1=remote .bb basename  $2=output bed
    curl -sL --fail --max-time 90 "$HUB/$1" -o "$OUT/tmp.bb"
    bigBedToBed "$OUT/tmp.bb" "$OUT/tmp.bed"
    norm < "$OUT/tmp.bed" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$2"
    rm -f "$OUT/tmp.bb" "$OUT/tmp.bed"
    echo "  $2: $(wc -l < "$2") elements, $(awk '{s+=$3-$2}END{print s}' "$2") bp"
}

echo "=== DANIO-CODE COPE/DOPE element tracks (danRer11) ==="
fetch "copes_dr11.bb"     "$OUT/COPE_regions.bed"
fetch "dopes_all_dr11.bb" "$OUT/DOPE_regions.bed"
echo "=== COPE/DOPE tracks complete ==="
