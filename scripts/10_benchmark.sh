#!/bin/bash
# 10_benchmark.sh
# Phase 5: Validation harness
#
# Benchmarks the pipeline output against known integration sites:
# 1. pIGLET14a / pIGLET24b — pipeline should rank 24b above 14a
# 2. Sleeping Beauty founder sites — cross-check with sampling-bias caveats

set -euo pipefail

# Load configuration
if [ ! -f "./config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ./config.sh

echo "=== Phase 5: Validation Benchmarks ==="
echo ""

cd results
# Resolve paths before changing directory
PIGLET_14A_BED="../${PIGLET_14A_BED}"
PIGLET_24B_BED="../${PIGLET_24B_BED}"
SB_FOUNDER_SITES_BED="../${SB_FOUNDER_SITES_BED}"

# ======================================================================
# Utility: check if coordinate file exists
# ======================================================================
check_benchmark_data() {
    local file="$1"
    local label="$2"
    if [ -f "$file" ] && [ -s "$file" ]; then
        echo "  [FOUND] $label: $file"
        return 0
    else
        echo "  [SKIP]  $label not found: $file"
        echo "          >> Add to DATA_NEEDED.md or provide coordinates"
        return 1
    fi
}

# ======================================================================
# Benchmark 1: pIGLET Sites
# ======================================================================
echo "--- Benchmark 1: pIGLET Integration Sites ---"
echo ""

PIGLET_AVAILABLE=false

if check_benchmark_data "$PIGLET_14A_BED" "pIGLET14a coordinates" && \
   check_benchmark_data "$PIGLET_24B_BED" "pIGLET24b coordinates"; then
    PIGLET_AVAILABLE=true
fi

if [ "$PIGLET_AVAILABLE" = true ]; then
    echo ""
    echo "  Checking pIGLET sites against consensus tiers..."
    echo ""

    # Normalize chr names
    awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$PIGLET_14A_BED" > piglet14a_chr.bed
    awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$PIGLET_24B_BED" > piglet24b_chr.bed

    # Check each tier
    echo "  pIGLET14a overlaps:"
    for tier in 100 90 80; do
        PEAKS="consensus_peaks_${tier}pct.bed"
        if [ -f "$PEAKS" ]; then
            OVERLAP=$(bedtools intersect -a piglet14a_chr.bed -b "$PEAKS" -u | wc -l)
            echo "    ${tier}% consensus peaks: $OVERLAP hit(s)"
        fi
    done

    echo ""
    echo "  pIGLET24b overlaps:"
    for tier in 100 90 80; do
        PEAKS="consensus_peaks_${tier}pct.bed"
        if [ -f "$PEAKS" ]; then
            OVERLAP=$(bedtools intersect -a piglet24b_chr.bed -b "$PEAKS" -u | wc -l)
            echo "    ${tier}% consensus peaks: $OVERLAP hit(s)"
        fi
    done

    # Check against final candidate list
    echo ""
    echo "  Checking against final candidate feature table..."

    if [ -f "candidate_feature_table.tsv" ]; then
        # Find rank of pIGLET sites in the feature table
        echo ""
        echo "  pIGLET14a rank in final candidates:"
        while IFS=$'\t' read -r chr start end rest; do
            RANK=$(awk -v c="$chr" -v s="$start" -v e="$end" \
                'NR>1 && $1==c && $2<=e && $3>=s {print NR-1; exit}' \
                candidate_feature_table.tsv)
            if [ -n "$RANK" ]; then
                echo "    $chr:$start-$end → Rank #$RANK"
            else
                echo "    $chr:$start-$end → NOT in final candidate list"
            fi
        done < piglet14a_chr.bed

        echo ""
        echo "  pIGLET24b rank in final candidates:"
        while IFS=$'\t' read -r chr start end rest; do
            RANK=$(awk -v c="$chr" -v s="$start" -v e="$end" \
                'NR>1 && $1==c && $2<=e && $3>=s {print NR-1; exit}' \
                candidate_feature_table.tsv)
            if [ -n "$RANK" ]; then
                echo "    $chr:$start-$end → Rank #$RANK"
            else
                echo "    $chr:$start-$end → NOT in final candidate list"
            fi
        done < piglet24b_chr.bed

        # Validation check: 24b should rank above 14a
        echo ""
        RANK_14A=$(while IFS=$'\t' read -r chr start end rest; do
            awk -v c="$chr" -v s="$start" -v e="$end" \
                'NR>1 && $1==c && $2<=e && $3>=s {print NR-1; exit}' \
                candidate_feature_table.tsv
        done < piglet14a_chr.bed)

        RANK_24B=$(while IFS=$'\t' read -r chr start end rest; do
            awk -v c="$chr" -v s="$start" -v e="$end" \
                'NR>1 && $1==c && $2<=e && $3>=s {print NR-1; exit}' \
                candidate_feature_table.tsv
        done < piglet24b_chr.bed)

        echo "  ╔══════════════════════════════════════════════════════╗"
        if [ -n "$RANK_24B" ] && [ -n "$RANK_14A" ]; then
            if [ "$RANK_24B" -lt "$RANK_14A" ]; then
                echo "  ║  ✓ PASS: pIGLET24b (rank $RANK_24B) > pIGLET14a (rank $RANK_14A)  ║"
            else
                echo "  ║  ✗ FAIL: pIGLET24b (rank $RANK_24B) ≤ pIGLET14a (rank $RANK_14A)  ║"
                echo "  ║         Pipeline ranking may need recalibration!           ║"
            fi
        elif [ -n "$RANK_24B" ] && [ -z "$RANK_14A" ]; then
            echo "  ║  ✓ PASS: pIGLET24b ranked (#$RANK_24B), 14a not in list     ║"
        elif [ -z "$RANK_24B" ] && [ -z "$RANK_14A" ]; then
            echo "  ║  ⚠ INCONCLUSIVE: Neither site in final candidate list        ║"
        else
            echo "  ║  ✗ FAIL: pIGLET14a ranked but 24b is not!                    ║"
            echo "  ║         This should not happen — investigate.                 ║"
        fi
        echo "  ╚══════════════════════════════════════════════════════╝"
    else
        echo "  WARNING: candidate_feature_table.tsv not found. Run script 09 first."
    fi

    # Clean up
    rm -f piglet14a_chr.bed piglet24b_chr.bed
else
    echo ""
    echo "  pIGLET benchmark skipped — coordinates not available."
    echo "  >> Add pIGLET14a and pIGLET24b coordinates to DATA_NEEDED.md locations"
fi
echo ""

# ======================================================================
# Benchmark 2: Sleeping Beauty Founder Sites
# ======================================================================
echo "--- Benchmark 2: Sleeping Beauty Founder Cross-check ---"
echo ""

if check_benchmark_data "$SB_FOUNDER_SITES_BED" "SB founder insertion sites"; then
    echo ""

    # Normalize chr names
    awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
        "$SB_FOUNDER_SITES_BED" > sb_sites_chr.bed

    NUM_SB=$(wc -l < sb_sites_chr.bed)
    echo "  Total SB founder sites: $NUM_SB"
    echo ""

    # CAVEAT: SB has TA-dinucleotide + AT-rich/bendability bias
    echo "  ┌─────────────────────────────────────────────────────────────┐"
    echo "  │ SAMPLING BIAS CAVEAT:                                       │"
    echo "  │ SB transposon has an obligate TA dinucleotide target and    │"
    echo "  │ AT-rich/DNA bendability integration bias. Comparisons are   │"
    echo "  │ only valid within SB-reachable sequence space. SB absence   │"
    echo "  │ in a candidate is NON-INFORMATIVE for TA-poor regions.      │"
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""

    # Check overlap with consensus tiers
    echo "  SB site overlaps with consensus peaks:"
    for tier in 100 90 80; do
        PEAKS="consensus_peaks_${tier}pct.bed"
        if [ -f "$PEAKS" ]; then
            OVERLAP=$(bedtools intersect -a sb_sites_chr.bed -b "$PEAKS" -u | wc -l)
            PCT=$(echo "scale=1; $OVERLAP * 100 / $NUM_SB" | bc)
            echo "    ${tier}% consensus: $OVERLAP / $NUM_SB ($PCT%)"
        fi
    done

    # Check overlap with final candidates
    echo ""
    echo "  SB site overlaps with final candidates:"
    if [ -f "safe_harbors_biological_100pct_sized.bed" ]; then
        OVERLAP=$(bedtools intersect -a sb_sites_chr.bed \
            -b safe_harbors_biological_100pct_sized.bed -u | wc -l)
        echo "    Biological filtered (100%): $OVERLAP / $NUM_SB"
    fi

    # TA-content analysis of candidates (for context)
    echo ""
    echo "  TA-dinucleotide context:"
    GENOME_FASTA="../reference/GRCz11.fa"
    if [ -f "$GENOME_FASTA" ] && [ -f "safe_harbors_biological_100pct_sized.bed" ]; then
        # Count candidates containing at least one TA
        bedtools getfasta \
            -fi "$GENOME_FASTA" \
            -bed safe_harbors_biological_100pct_sized.bed \
            -tab \
        | awk '{
            seq = toupper($2)
            if (index(seq, "TA") > 0) ta_count++
            total++
        }
        END {
            printf "    Candidates with ≥1 TA dinucleotide: %d / %d (%.1f%%)\n", ta_count, total, ta_count*100/total
            printf "    (SB can only integrate at TA sites — candidates without TA are SB-unreachable)\n"
        }'
    else
        echo "    [SKIP] Genome FASTA not available for TA analysis"
    fi

    # Clean up
    rm -f sb_sites_chr.bed
else
    echo "  SB benchmark skipped — founder site coordinates not available."
    echo "  >> Add SB founder insertion coordinates to DATA_NEEDED.md location"
fi
echo ""

# ======================================================================
# Summary
# ======================================================================
echo "=== Benchmark Summary ==="
echo ""
echo "  pIGLET benchmark: $([ "$PIGLET_AVAILABLE" = true ] && echo "RAN" || echo "SKIPPED (data needed)")"
echo "  SB founder benchmark: $([ -f "$SB_FOUNDER_SITES_BED" ] && echo "RAN" || echo "SKIPPED (data needed)")"
echo ""
echo "  See DATA_NEEDED.md for instructions on obtaining benchmark coordinates."
echo ""
echo "=== Benchmarking complete ==="