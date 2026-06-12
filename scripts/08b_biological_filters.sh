#!/bin/bash
# 08b_biological_filters.sh
# Phase 3: Apply biological criteria filters to safe harbor candidates
#
# Filters applied (each is a no-op if the required data file is absent):
# 1. Regulatory element exclusion (H3K27ac, H3K4me1, H3K4me3, CTCF)
# 2. Conservation exclusion (phastCons, CNEs)
# 3. Mappability / blacklist / segmental duplication filter
# 4. GSH-style distance criteria (cancer genes, dosage-sensitive, miRNA/ncRNA, UCEs)
# 5. Hi-C TAD boundary flagging (annotation only, not exclusion)
# 6. COPE/DOPE cross-reference (annotation only)

set -euo pipefail

# Load configuration
if [ ! -f "./config.sh" ]; then
    echo "ERROR: config.sh not found in parent directory"
    exit 1
fi
source ./config.sh

echo "=== Phase 3: Biological Criteria Filtering ==="
echo ""

cd results

# ======================================================================
# Utility function: check if a BED file exists and is non-empty
# ======================================================================
check_data() {
    local file="$1"
    local label="$2"
    if [ -f "$file" ] && [ -s "$file" ]; then
        echo "  [FOUND] $label: $file ($(wc -l < "$file") regions)"
        return 0
    else
        echo "  [SKIP]  $label: $file not found or empty"
        return 1
    fi
}

# ======================================================================
# Utility function: subtract BED from candidates with chr-name guard
# ======================================================================
subtract_track() {
    local candidates="$1"
    local track="$2"
    local output="$3"
    local label="$4"

    # Verify shared chromosome names
    local shared
    shared=$(comm -12 \
        <(cut -f1 "$candidates" | sort -u) \
        <(cut -f1 "$track" | sort -u) | wc -l)

    if [ "$shared" -eq 0 ]; then
        echo "    WARNING: Zero shared chr names between candidates and $label"
        echo "    Attempting chr-prefix normalization on track..."
        # Try normalizing the track
        awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$track" \
            > "${track}.chr_normalized.tmp"
        shared=$(comm -12 \
            <(cut -f1 "$candidates" | sort -u) \
            <(cut -f1 "${track}.chr_normalized.tmp" | sort -u) | wc -l)
        if [ "$shared" -eq 0 ]; then
            echo "    ERROR: Still zero shared chr names after normalization. Skipping $label."
            cp "$candidates" "$output"
            rm -f "${track}.chr_normalized.tmp"
            return 1
        fi
        bedtools subtract -a "$candidates" -b "${track}.chr_normalized.tmp" > "$output"
        rm -f "${track}.chr_normalized.tmp"
    else
        bedtools subtract -a "$candidates" -b "$track" > "$output"
    fi

    local before after removed
    before=$(wc -l < "$candidates")
    after=$(wc -l < "$output")
    removed=$((before - after))
    echo "    $label: removed $removed candidates ($before → $after)"
    return 0
}

# ======================================================================
# Input: use the sized candidates from script 08
# ======================================================================
INPUT_100="safe_harbors_100pct_sized.bed"
INPUT_90="safe_harbors_90pct_sized.bed"
INPUT_80="safe_harbors_80pct_sized.bed"

# We'll work primarily with the 100% tier, then apply same filters to others
if [ ! -s "$INPUT_100" ]; then
    echo "ERROR: $INPUT_100 not found or empty. Run 08_find_safe_harbors.sh first."
    exit 1
fi

echo "Starting candidates:"
echo "  100% consensus: $(wc -l < "$INPUT_100")"
echo "  90% consensus: $(wc -l < "$INPUT_90")"
echo "  80% consensus: $(wc -l < "$INPUT_80")"
echo ""

# Track current working file for each tier
CURRENT_100="$INPUT_100"
CURRENT_90="$INPUT_90"
CURRENT_80="$INPUT_80"

# ======================================================================
# Step 1: Regulatory Element Exclusion (ChIP-seq)
# ======================================================================
echo "--- Step 1: Regulatory Element Exclusion ---"
echo ""

CHIP_TRACKS=()
CHIP_LABELS=()

if check_data "$CHIP_H3K27AC_BED" "H3K27ac (active enhancers)"; then
    CHIP_TRACKS+=("$CHIP_H3K27AC_BED")
    CHIP_LABELS+=("H3K27ac")
fi
if check_data "$CHIP_H3K4ME1_BED" "H3K4me1 (poised enhancers)"; then
    CHIP_TRACKS+=("$CHIP_H3K4ME1_BED")
    CHIP_LABELS+=("H3K4me1")
fi
if check_data "$CHIP_H3K4ME3_BED" "H3K4me3 (active promoters)"; then
    CHIP_TRACKS+=("$CHIP_H3K4ME3_BED")
    CHIP_LABELS+=("H3K4me3")
fi
if check_data "$CHIP_CTCF_BED" "CTCF (insulators)"; then
    CHIP_TRACKS+=("$CHIP_CTCF_BED")
    CHIP_LABELS+=("CTCF")
fi

if [ ${#CHIP_TRACKS[@]} -gt 0 ]; then
    echo ""
    echo "  Applying ChIP-seq exclusions..."

    # Combine all available ChIP tracks into one exclusion set
    cat "${CHIP_TRACKS[@]}" \
        | cut -f1-3 \
        | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > regulatory_elements_combined.bed

    echo "  Combined regulatory regions: $(wc -l < regulatory_elements_combined.bed)"

    subtract_track "$CURRENT_100" "regulatory_elements_combined.bed" "filtered_100_chip.bed" "ChIP-seq regulatory"
    subtract_track "$CURRENT_90" "regulatory_elements_combined.bed" "filtered_90_chip.bed" "ChIP-seq regulatory"
    subtract_track "$CURRENT_80" "regulatory_elements_combined.bed" "filtered_80_chip.bed" "ChIP-seq regulatory"

    CURRENT_100="filtered_100_chip.bed"
    CURRENT_90="filtered_90_chip.bed"
    CURRENT_80="filtered_80_chip.bed"
else
    echo "  No ChIP-seq tracks available. Skipping regulatory element exclusion."
    echo "  >> See DATA_NEEDED.md for required files."
fi
echo ""

# ======================================================================
# Step 2: Conservation Exclusion
# ======================================================================
echo "--- Step 2: Conservation Exclusion ---"
echo ""

CONSERVATION_APPLIED=false

if check_data "$CONSERVATION_BED" "phastCons conserved regions"; then
    echo ""
    subtract_track "$CURRENT_100" "$CONSERVATION_BED" "filtered_100_cons.bed" "Conservation (phastCons)"
    subtract_track "$CURRENT_90" "$CONSERVATION_BED" "filtered_90_cons.bed" "Conservation (phastCons)"
    subtract_track "$CURRENT_80" "$CONSERVATION_BED" "filtered_80_cons.bed" "Conservation (phastCons)"
    CURRENT_100="filtered_100_cons.bed"
    CURRENT_90="filtered_90_cons.bed"
    CURRENT_80="filtered_80_cons.bed"
    CONSERVATION_APPLIED=true
fi

if check_data "$CNE_BED" "Conserved non-coding elements (CNEs)"; then
    echo ""
    subtract_track "$CURRENT_100" "$CNE_BED" "filtered_100_cne.bed" "CNEs"
    subtract_track "$CURRENT_90" "$CNE_BED" "filtered_90_cne.bed" "CNEs"
    subtract_track "$CURRENT_80" "$CNE_BED" "filtered_80_cne.bed" "CNEs"
    CURRENT_100="filtered_100_cne.bed"
    CURRENT_90="filtered_90_cne.bed"
    CURRENT_80="filtered_80_cne.bed"
    CONSERVATION_APPLIED=true
fi

if [ "$CONSERVATION_APPLIED" = false ]; then
    echo "  No conservation tracks available. Skipping."
    echo "  >> See DATA_NEEDED.md for required files."
fi
echo ""

# ======================================================================
# Step 3: Mappability / Blacklist / Segmental Duplication
# ======================================================================
echo "--- Step 3: Mappability & Blacklist Filtering ---"
echo ""

MAPPABILITY_APPLIED=false

if check_data "$MAPPABILITY_BED" "Low-mappability regions"; then
    echo ""
    subtract_track "$CURRENT_100" "$MAPPABILITY_BED" "filtered_100_map.bed" "Low mappability"
    subtract_track "$CURRENT_90" "$MAPPABILITY_BED" "filtered_90_map.bed" "Low mappability"
    subtract_track "$CURRENT_80" "$MAPPABILITY_BED" "filtered_80_map.bed" "Low mappability"
    CURRENT_100="filtered_100_map.bed"
    CURRENT_90="filtered_90_map.bed"
    CURRENT_80="filtered_80_map.bed"
    MAPPABILITY_APPLIED=true
fi

if check_data "$SEGDUP_BED" "Segmental duplications"; then
    echo ""
    subtract_track "$CURRENT_100" "$SEGDUP_BED" "filtered_100_segdup.bed" "Segmental duplications"
    subtract_track "$CURRENT_90" "$SEGDUP_BED" "filtered_90_segdup.bed" "Segmental duplications"
    subtract_track "$CURRENT_80" "$SEGDUP_BED" "filtered_80_segdup.bed" "Segmental duplications"
    CURRENT_100="filtered_100_segdup.bed"
    CURRENT_90="filtered_90_segdup.bed"
    CURRENT_80="filtered_80_segdup.bed"
    MAPPABILITY_APPLIED=true
fi

if check_data "$BLACKLIST_BED" "Blacklist regions"; then
    echo ""
    subtract_track "$CURRENT_100" "$BLACKLIST_BED" "filtered_100_bl.bed" "Blacklist"
    subtract_track "$CURRENT_90" "$BLACKLIST_BED" "filtered_90_bl.bed" "Blacklist"
    subtract_track "$CURRENT_80" "$BLACKLIST_BED" "filtered_80_bl.bed" "Blacklist"
    CURRENT_100="filtered_100_bl.bed"
    CURRENT_90="filtered_90_bl.bed"
    CURRENT_80="filtered_80_bl.bed"
    MAPPABILITY_APPLIED=true
fi

if [ "$MAPPABILITY_APPLIED" = false ]; then
    echo "  No mappability/blacklist tracks available. Skipping."
    echo "  >> See DATA_NEEDED.md for required files."
fi
echo ""

# ======================================================================
# Step 4: GSH-style Distance Criteria
# ======================================================================
echo "--- Step 4: GSH Distance Criteria (Papapetrou/Sadelain) ---"
echo ""

GSH_APPLIED=false

# 4a: miRNA/ncRNA buffer
if check_data "$MIRNA_NCRNA_BED" "miRNA/ncRNA coordinates"; then
    echo ""
    echo "  Applying ${GSH_MIRNA_BUFFER_KB}kb buffer around miRNAs/ncRNAs..."
    bedtools slop \
        -i "$MIRNA_NCRNA_BED" \
        -g ../reference/GRCz11.chrom.sizes \
        -b $((GSH_MIRNA_BUFFER_KB * 1000)) \
    | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > mirna_ncrna_buffered.bed

    subtract_track "$CURRENT_100" "mirna_ncrna_buffered.bed" "filtered_100_mirna.bed" "miRNA/ncRNA ±${GSH_MIRNA_BUFFER_KB}kb"
    subtract_track "$CURRENT_90" "mirna_ncrna_buffered.bed" "filtered_90_mirna.bed" "miRNA/ncRNA ±${GSH_MIRNA_BUFFER_KB}kb"
    subtract_track "$CURRENT_80" "mirna_ncrna_buffered.bed" "filtered_80_mirna.bed" "miRNA/ncRNA ±${GSH_MIRNA_BUFFER_KB}kb"
    CURRENT_100="filtered_100_mirna.bed"
    CURRENT_90="filtered_90_mirna.bed"
    CURRENT_80="filtered_80_mirna.bed"
    GSH_APPLIED=true
fi

# 4b: Cancer-associated genes buffer
if check_data "$CANCER_ASSOCIATED_GENES_BED" "Cancer-associated genes"; then
    echo ""
    echo "  Applying ${GSH_CANCER_GENE_BUFFER_KB}kb buffer around cancer genes..."
    bedtools slop \
        -i "$CANCER_ASSOCIATED_GENES_BED" \
        -g ../reference/GRCz11.chrom.sizes \
        -b $((GSH_CANCER_GENE_BUFFER_KB * 1000)) \
    | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > cancer_genes_buffered.bed

    subtract_track "$CURRENT_100" "cancer_genes_buffered.bed" "filtered_100_cancer.bed" "Cancer genes ±${GSH_CANCER_GENE_BUFFER_KB}kb"
    subtract_track "$CURRENT_90" "cancer_genes_buffered.bed" "filtered_90_cancer.bed" "Cancer genes ±${GSH_CANCER_GENE_BUFFER_KB}kb"
    subtract_track "$CURRENT_80" "cancer_genes_buffered.bed" "filtered_80_cancer.bed" "Cancer genes ±${GSH_CANCER_GENE_BUFFER_KB}kb"
    CURRENT_100="filtered_100_cancer.bed"
    CURRENT_90="filtered_90_cancer.bed"
    CURRENT_80="filtered_80_cancer.bed"
    GSH_APPLIED=true
fi

# 4c: Dosage-sensitive genes buffer
if check_data "$DOSAGE_SENSITIVE_GENES_BED" "Dosage-sensitive genes"; then
    echo ""
    echo "  Applying ${GSH_DOSAGE_SENSITIVE_BUFFER_KB}kb buffer around dosage-sensitive genes..."
    bedtools slop \
        -i "$DOSAGE_SENSITIVE_GENES_BED" \
        -g ../reference/GRCz11.chrom.sizes \
        -b $((GSH_DOSAGE_SENSITIVE_BUFFER_KB * 1000)) \
    | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > dosage_sensitive_buffered.bed

    subtract_track "$CURRENT_100" "dosage_sensitive_buffered.bed" "filtered_100_dosage.bed" "Dosage-sensitive ±${GSH_DOSAGE_SENSITIVE_BUFFER_KB}kb"
    subtract_track "$CURRENT_90" "dosage_sensitive_buffered.bed" "filtered_90_dosage.bed" "Dosage-sensitive ±${GSH_DOSAGE_SENSITIVE_BUFFER_KB}kb"
    subtract_track "$CURRENT_80" "dosage_sensitive_buffered.bed" "filtered_80_dosage.bed" "Dosage-sensitive ±${GSH_DOSAGE_SENSITIVE_BUFFER_KB}kb"
    CURRENT_100="filtered_100_dosage.bed"
    CURRENT_90="filtered_90_dosage.bed"
    CURRENT_80="filtered_80_dosage.bed"
    GSH_APPLIED=true
fi

# 4d: Ultraconserved elements
if check_data "$UCE_BED" "Ultraconserved elements (UCEs)"; then
    echo ""
    echo "  Applying ${GSH_UCE_BUFFER_KB}kb buffer around UCEs..."
    bedtools slop \
        -i "$UCE_BED" \
        -g ../reference/GRCz11.chrom.sizes \
        -b $((GSH_UCE_BUFFER_KB * 1000)) \
    | awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > uce_buffered.bed

    subtract_track "$CURRENT_100" "uce_buffered.bed" "filtered_100_uce.bed" "UCEs ±${GSH_UCE_BUFFER_KB}kb"
    subtract_track "$CURRENT_90" "uce_buffered.bed" "filtered_90_uce.bed" "UCEs ±${GSH_UCE_BUFFER_KB}kb"
    subtract_track "$CURRENT_80" "uce_buffered.bed" "filtered_80_uce.bed" "UCEs ±${GSH_UCE_BUFFER_KB}kb"
    CURRENT_100="filtered_100_uce.bed"
    CURRENT_90="filtered_90_uce.bed"
    CURRENT_80="filtered_80_uce.bed"
    GSH_APPLIED=true
fi

if [ "$GSH_APPLIED" = false ]; then
    echo "  No GSH gene lists available. Skipping extended distance criteria."
    echo "  >> See DATA_NEEDED.md for required files."
fi
echo ""

# ======================================================================
# Step 5: Hi-C TAD Boundary Flagging (annotation only — not exclusion)
# ======================================================================
echo "--- Step 5: TAD Boundary Flagging ---"
echo ""

if check_data "$TAD_BOUNDARIES_BED" "Hi-C TAD boundaries"; then
    echo ""
    echo "  Flagging candidates near TAD boundaries (for information only)..."
    bedtools closest \
        -a <(sort -k1,1 -k2,2n "$CURRENT_100") \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$TAD_BOUNDARIES_BED" | sort -k1,1 -k2,2n) \
        -d \
    > tad_boundary_distances_100.txt
    echo "  TAD boundary distances written to tad_boundary_distances_100.txt"
else
    echo "  No TAD boundary data available. Skipping."
    echo "  >> See DATA_NEEDED.md for required files."
fi
echo ""

# ======================================================================
# Step 6: COPE/DOPE Cross-reference (annotation — positive signal)
# ======================================================================
echo "--- Step 6: COPE/DOPE Cross-reference ---"
echo ""

if check_data "$COPE_BED" "COPE regions (open but inert)"; then
    echo ""
    echo "  Flagging candidates overlapping COPE regions (positive: open but lacking active marks)..."
    bedtools intersect \
        -a "$CURRENT_100" \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$COPE_BED") \
        -u \
    > cope_overlapping_100.bed
    echo "  Candidates overlapping COPE: $(wc -l < cope_overlapping_100.bed) / $(wc -l < "$CURRENT_100")"
else
    echo "  No COPE track available. Skipping."
fi

if check_data "$DOPE_BED" "DOPE regions (distal open chromatin)"; then
    echo ""
    echo "  Flagging candidates overlapping DOPE regions..."
    bedtools intersect \
        -a "$CURRENT_100" \
        -b <(awk 'BEGIN{OFS="\t"} { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$DOPE_BED") \
        -u \
    > dope_overlapping_100.bed
    echo "  Candidates overlapping DOPE: $(wc -l < dope_overlapping_100.bed) / $(wc -l < "$CURRENT_100")"
else
    echo "  No DOPE track available. Skipping."
fi
echo ""

# ======================================================================
# Final output: copy filtered candidates to final names
# ======================================================================
echo "--- Final Filtered Output ---"
echo ""

cp "$CURRENT_100" safe_harbors_biological_100pct.bed
cp "$CURRENT_90" safe_harbors_biological_90pct.bed
cp "$CURRENT_80" safe_harbors_biological_80pct.bed

echo "Final candidates after biological filtering:"
echo "  100% consensus: $(wc -l < safe_harbors_biological_100pct.bed)"
echo "  90% consensus: $(wc -l < safe_harbors_biological_90pct.bed)"
echo "  80% consensus: $(wc -l < safe_harbors_biological_80pct.bed)"
echo ""

# Re-apply minimum size filter (subtraction can create fragments)
awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_biological_100pct.bed \
    > safe_harbors_biological_100pct_sized.bed
awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_biological_90pct.bed \
    > safe_harbors_biological_90pct_sized.bed
awk -v min="$MIN_HARBOR_SIZE" '($3 - $2) >= min' safe_harbors_biological_80pct.bed \
    > safe_harbors_biological_80pct_sized.bed

echo "After re-applying minimum size filter (>=${MIN_HARBOR_SIZE}bp):"
echo "  100% consensus: $(wc -l < safe_harbors_biological_100pct_sized.bed)"
echo "  90% consensus: $(wc -l < safe_harbors_biological_90pct_sized.bed)"
echo "  80% consensus: $(wc -l < safe_harbors_biological_80pct_sized.bed)"
echo ""

# Clean up intermediate files
rm -f filtered_*.bed regulatory_elements_combined.bed
rm -f mirna_ncrna_buffered.bed cancer_genes_buffered.bed
rm -f dosage_sensitive_buffered.bed uce_buffered.bed

echo "=== Phase 3: Biological filtering complete ==="