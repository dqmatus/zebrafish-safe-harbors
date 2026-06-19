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

echo "=== ChIP-seq merge complete ==="
