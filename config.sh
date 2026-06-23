#!/bin/bash
# config.sh
# Central configuration for the zebrafish safe-harbors pipeline
# All magic numbers and thresholds are defined here for reproducibility

# ======================================================================
# Consensus Peak Summit Parameters
# ======================================================================
# Summit window size (bp on each side of the narrowPeak col-10 summit)
# Default: 150bp on each side = 300bp total window
SUMMIT_WINDOW_HALF=150

# Consensus tiers (percentage of stages/biosamples required for consistency)
# These now represent: "open in X% of developmental stages or biosamples"
# After collapsing replicates to biosample level, then to stage level
CONSENSUS_MIN_PCT_80=80
CONSENSUS_MIN_PCT_90=90
CONSENSUS_MIN_PCT_100=100

# ======================================================================
# Score Normalization & Aggregation
# ======================================================================
# Method for score aggregation across samples:
# - "quantile": quantile-normalize scores across samples, then median
# - "rank": rank scores, then median
# - "max": use maximum score (conservative, less recommended)
SCORE_NORMALIZATION_METHOD="quantile"

# ======================================================================
# Exclusion & Filtering Parameters
# ======================================================================
# Gene buffer distance (bp on each side of gene body)
GENE_BUFFER_KB=10
# Relaxed gene buffer for the alternate (more-permissive) candidate tier
RELAXED_GENE_BUFFER_KB=2

# Essential gene buffer (bp)
ESSENTIAL_GENE_BUFFER_KB=20

# Telomere distance from chromosome end (bp)
TELOMERE_END_DISTANCE=500000

# Centromere model (currently disabled — zebrafish midpoint is inaccurate)
CENTROMERE_USE=false

# Minimum safe harbor size (bp — must fit landing pad construct)
MIN_HARBOR_SIZE=300

# Restrict candidates to primary chromosomes (chr1-25). Unplaced scaffolds
# (chrUn_*, *_alt, *_random) and chrM make poor landing pads (unreliable
# assembly/coordinates, hard to validate) and — lacking gene annotation — would
# otherwise pass the gene-distal filter and even rank highly, since bedtools
# closest reports gene_distance = -1 there, which 09 treats as maximally distal.
PRIMARY_CHROM_REGEX='^chr([1-9]|1[0-9]|2[0-5])$'

# ======================================================================
# Phase 3: GSH-style Distance Criteria (Papapetrou/Sadelain)
# ======================================================================
# Extended gene buffer for GSH criteria (bp)
GSH_GENE_BUFFER_KB=50

# Distance from cancer-associated genes (bp)
GSH_CANCER_GENE_BUFFER_KB=300

# Distance from dosage-sensitive genes (bp)
GSH_DOSAGE_SENSITIVE_BUFFER_KB=100

# Distance from miRNAs/ncRNAs (bp)
GSH_MIRNA_BUFFER_KB=50

# Distance from lncRNAs (bp). eLife GSH paper uses 100kb (human); scaled down
# for the ~1.4Gb zebrafish genome (100kb would over-mask; 25kb adds protection
# beyond the 10kb gene buffer without decimating the candidate list).
GSH_LNCRNA_BUFFER_KB=25

# Distance from ultraconserved elements (bp)
GSH_UCE_BUFFER_KB=50

# ======================================================================
# Phase 3: Regulatory Element Exclusion Tracks
# ======================================================================
# ChIP-seq peak files (BED format, chr-prefixed)
# Set to empty string or nonexistent path if not available
CHIP_H3K27AC_BED="annotation/chip/H3K27ac_peaks.bed"
CHIP_H3K4ME1_BED="annotation/chip/H3K4me1_peaks.bed"
CHIP_H3K4ME3_BED="annotation/chip/H3K4me3_peaks.bed"
CHIP_CTCF_BED="annotation/chip/CTCF_peaks.bed"
# H3K27me3 = repressive/Polycomb mark; excluded to avoid landing in silenced
# chromatin where a transgene could itself be silenced (distinct from active marks)
CHIP_H3K27ME3_BED="annotation/chip/H3K27me3_peaks.bed"
# Combined union of all OTHER ChIP targets (H3K36me3, H2A.Z/H2AFV, H3K14ac,
# Pol II, and TFs) — active/transcribed/TF-bound chromatin, built by 06c.
CHIP_OTHER_BED="annotation/chip/other_chip_peaks.bed"

# Conservation tracks
CONSERVATION_BED="annotation/conservation/phastCons_conserved.bed"
CNE_BED="annotation/conservation/CNEs.bed"

# Hi-C TAD boundaries (flagging only, not exclusion)
TAD_BOUNDARIES_BED="annotation/hic/TAD_boundaries.bed"

# DANIO-CODE COPE/DOPE tracks
COPE_BED="annotation/cope_dope/COPE_regions.bed"
DOPE_BED="annotation/cope_dope/DOPE_regions.bed"

# Supplement-derived regulatory tracks (Baranasic et al. 2022, danRer11)
# Built by scripts/04b_extract_supplement_annotations.py
VALIDATED_ENHANCERS_BED="annotation/regulatory/validated_enhancers.bed"  # T10+T11, experimentally validated
PROMOTERS_BED="annotation/regulatory/promoters.bed"                      # T4 consensus + T13 core promoters
NOVEL_LNCRNA_BED="annotation/regulatory/novel_lncrna.bed"               # T3 intergenic lncRNA/TUCP bodies

# ======================================================================
# Phase 3: Mappability & Blacklist
# ======================================================================
MAPPABILITY_BED="annotation/mappability/GRCz11_mappability_low.bed"
SEGDUP_BED="annotation/mappability/segmental_duplications.bed"
BLACKLIST_BED="annotation/mappability/GRCz11_blacklist.bed"

# Mappability threshold (regions below this are excluded)
MAPPABILITY_THRESHOLD=0.8
# genmap K-mer length and allowed errors for the uniqueness scan (scripts/13).
# K=50 / E=0 = "is this 50-mer unique in the genome?" (Illumina-read-length,
# perfect-match). With E=0, mappability is 1 for unique k-mers and <0.8 for any
# k-mer occurring >=2x, so MAPPABILITY_THRESHOLD=0.8 excludes all non-unique runs.
GSH_MAPPABILITY_K=50
GSH_MAPPABILITY_E=0

# ======================================================================
# Phase 3: GSH Gene Lists
# ======================================================================
GENE_BED_COMPREHENSIVE="annotation/genes_comprehensive.bed"
MIRNA_NCRNA_BED="annotation/mirna_ncrna_coordinates.bed"
LNCRNA_BED="annotation/lncrna_coordinates.bed"
DOSAGE_SENSITIVE_GENES_BED="annotation/dosage_sensitive_genes.bed"
CANCER_ASSOCIATED_GENES_BED="annotation/cancer_associated_genes.bed"
UCE_BED="annotation/ultraconserved_elements.bed"

# ======================================================================
# Phase 4: sgRNA Targetability
# ======================================================================
SGRNA_BED="annotation/sgrna/candidate_guides.bed"
SGRNA_MIN_SPECIFICITY=0.95
SGRNA_TOOL="CRISPOR"

# ======================================================================
# Phase 5: Validation Benchmarks
# ======================================================================
PIGLET_14A_BED="annotation/benchmarks/pIGLET14a_coords.bed"
PIGLET_24B_BED="annotation/benchmarks/pIGLET24b_coords.bed"
# tpl102: position-effect-resistant attP landing site (Roberts et al. 2014, Development)
TPL102_BED="annotation/benchmarks/tpl102_coords.bed"
SB_FOUNDER_SITES_BED="annotation/benchmarks/SB_founder_sites.bed"

# ======================================================================
# Annotation & Data Provenance
# ======================================================================
# Ensembl release version
ENSEMBL_RELEASE=110

# Genome build
GENOME_BUILD="GRCz11"

# UCSC RepeatMasker release (rmsk.txt source)
RMSK_SOURCE="https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz"

# DANIO-CODE accession metadata (SQ → biosample → stage mapping)
# Format: one line per sample: "SQ_ID\tBiosample_ID\tStage\tLab"
# Example: DCD0031275Q  DCD007335BS 0hpf    Zhang_Lab
# If not present, pipeline assumes each file is a unique biosample/stage
DANIO_CODE_METADATA="annotation/danio_code_accession_map.txt"

# ======================================================================
# Quality Control & Logging
# ======================================================================
# Enable verbose logging
VERBOSE=true

# Log file location (relative to script directory)
LOGFILE="${PWD}/pipeline.log"

# ======================================================================
# Export for use in other scripts
# ======================================================================
export SUMMIT_WINDOW_HALF
export CONSENSUS_MIN_PCT_80 CONSENSUS_MIN_PCT_90 CONSENSUS_MIN_PCT_100
export SCORE_NORMALIZATION_METHOD
export GENE_BUFFER_KB RELAXED_GENE_BUFFER_KB ESSENTIAL_GENE_BUFFER_KB TELOMERE_END_DISTANCE
export CENTROMERE_USE MIN_HARBOR_SIZE PRIMARY_CHROM_REGEX
export GSH_GENE_BUFFER_KB GSH_CANCER_GENE_BUFFER_KB GSH_DOSAGE_SENSITIVE_BUFFER_KB
export GSH_MIRNA_BUFFER_KB GSH_LNCRNA_BUFFER_KB GSH_UCE_BUFFER_KB
export CHIP_H3K27AC_BED CHIP_H3K4ME1_BED CHIP_H3K4ME3_BED CHIP_CTCF_BED CHIP_H3K27ME3_BED CHIP_OTHER_BED
export CONSERVATION_BED CNE_BED TAD_BOUNDARIES_BED
export COPE_BED DOPE_BED
export VALIDATED_ENHANCERS_BED PROMOTERS_BED NOVEL_LNCRNA_BED
export MAPPABILITY_BED SEGDUP_BED BLACKLIST_BED MAPPABILITY_THRESHOLD GSH_MAPPABILITY_K GSH_MAPPABILITY_E
export GENE_BED_COMPREHENSIVE MIRNA_NCRNA_BED LNCRNA_BED
export DOSAGE_SENSITIVE_GENES_BED CANCER_ASSOCIATED_GENES_BED UCE_BED
export SGRNA_BED SGRNA_MIN_SPECIFICITY SGRNA_TOOL
export PIGLET_14A_BED PIGLET_24B_BED TPL102_BED SB_FOUNDER_SITES_BED
export ENSEMBL_RELEASE GENOME_BUILD RMSK_SOURCE DANIO_CODE_METADATA
export VERBOSE LOGFILE