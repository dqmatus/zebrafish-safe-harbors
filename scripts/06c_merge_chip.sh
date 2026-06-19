#!/bin/bash
# 06c_merge_chip.sh
# Build per-mark union ChIP-seq peak BEDs from downloaded DANIO-CODE narrowPeak
# samples. For each histone mark / TF, all replicate/stage samples are merged
# into a single non-overlapping union track used as a regulatory exclusion set
# in 08b_biological_filters.sh.
#
# Input : annotation/chip/ChIP-seq_*<MARK>*.narrowPeak[.gz]  (multiple samples)
# Output: annotation/chip/<MARK>_peaks.bed                   (merged union)
#
# Reproduces the manually-built H3K27ac_peaks.bed and extends the same recipe
# to H3K4me1, H3K4me3, and CTCF once those samples are downloaded.

set -euo pipefail

# Match mark tokens case-insensitively: some DANIO-CODE files use "H3K27Ac"
# (capital A) vs our "H3K27ac" token (e.g. Zon Lab). nullglob so an absent mark
# expands to nothing rather than a literal pattern.
shopt -s nocaseglob nullglob

cd annotation/chip

# ----------------------------------------------------------------------
# Marks to build: "MATCH_TOKEN:OUTPUT_FILE"
# MATCH_TOKEN is matched against the sample filenames (case-sensitive).
# OUTPUT_FILE must match the path expected by config.sh.
# ----------------------------------------------------------------------
MARKS=(
    "H3K27ac:H3K27ac_peaks.bed"
    "H3K4me1:H3K4me1_peaks.bed"
    "H3K4me3:H3K4me3_peaks.bed"
    "H3K27me3:H3K27me3_peaks.bed"
    "CTCF:CTCF_peaks.bed"
)

# ----------------------------------------------------------------------
# Stream a narrowPeak sample (gzipped or plain) to stdout.
# ----------------------------------------------------------------------
read_peaks() {
    local f="$1"
    case "$f" in
        *.gz) gzip -dc "$f" ;;
        *)    cat "$f" ;;
    esac
}

echo "=== Merging ChIP-seq samples into per-mark union tracks ==="
echo ""

for entry in "${MARKS[@]}"; do
    token="${entry%%:*}"
    out="${entry##*:}"

    # Collect samples matching this mark, deduped by DANIO-CODE accession
    # (DCD#######SQ). When both a .gz and a plain copy exist for the same
    # accession, prefer the .gz (this is what the original H3K27ac merge used).
    # Glob .gz first so it wins the dedup. (Plain `seen`/`picked` strings keep
    # this compatible with macOS bash 3.2, which lacks associative arrays.)
    picked=()
    seen=" "
    for f in ChIP-seq_*"${token}"*.narrowPeak.gz ChIP-seq_*"${token}"*.narrowPeak; do
        [ -e "$f" ] || continue
        acc=$(echo "$f" | grep -oE 'DCD[0-9]+SQ' || true)
        [ -n "$acc" ] || acc="$f"
        case "$seen" in
            *" $acc "*)
                echo "  [skip] $f  (accession $acc already taken)"
                continue ;;
        esac
        seen="$seen$acc "
        picked+=("$f")
    done

    n=${#picked[@]}
    if [ "$n" -eq 0 ]; then
        echo "[--] $token: no samples found — skipping (download samples to build $out)"
        echo ""
        continue
    fi

    echo "[OK] $token: merging $n sample(s) -> $out"

    # chrom/start/end only, chr-prefix normalized, genome-sorted, union-merged.
    for f in "${picked[@]}"; do read_peaks "$f"; done \
        | cut -f1-3 \
        | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > "$out"

    peaks=$(wc -l < "$out")
    bp=$(awk '{sum+=$3-$2} END{print sum+0}' "$out")
    echo "     $peaks merged peaks, $bp bp covered"
    echo ""
done

# ----------------------------------------------------------------------
# Catch-all: the OTHER active/transcribed ChIP targets — H3K36me3 (transcribed
# bodies), H2A.Z/H2AFV, H3K14ac, Pol II — merged into one union. They get
# identical discard-whole treatment in 08b, so a single combined track is simplest.
# Excludes (case-insensitively):
#   - the 5 named marks built above (H3K27ac/H3K4me1/H3K4me3/H3K27me3/CTCF)
#   - developmental pioneer TFs (Pou5f3/nanog/Zic3/Sall4/Cdx4/Mxtx2): these bind
#     the early-embryo genome promiscuously (Pou5f3 ~3.8% of the genome) and are
#     NOT discrete regulatory elements, so they are unsuitable as an exclusion mask.
# ----------------------------------------------------------------------
shopt -s nocasematch
other=()
seen=" "
for f in ChIP-seq_*.narrowPeak.gz ChIP-seq_*.narrowPeak; do
    [ -e "$f" ] || continue
    case "$f" in
        *h3k27ac*|*h3k4me1*|*h3k4me3*|*h3k27me3*|*ctcf*) continue ;;
        *pou5f3*|*nanog*|*zic3*|*sall4*|*cdx4*|*mxtx2*) continue ;;
    esac
    acc=$(echo "$f" | grep -oE 'DCD[0-9]+SQ' || true); [ -n "$acc" ] || acc="$f"
    case "$seen" in *" $acc "*) continue ;; esac
    seen="$seen$acc "
    other+=("$f")
done
shopt -u nocasematch
if [ ${#other[@]} -gt 0 ]; then
    echo "[OK] other targets: merging ${#other[@]} sample(s) -> other_chip_peaks.bed"
    for f in "${other[@]}"; do read_peaks "$f"; done \
        | cut -f1-3 \
        | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > other_chip_peaks.bed
    echo "     $(wc -l < other_chip_peaks.bed) merged peaks, $(awk '{s+=$3-$2}END{print s+0}' other_chip_peaks.bed) bp"
    echo ""
fi

echo "=== ChIP-seq merge complete ==="
