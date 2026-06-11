#!/bin/bash
# config.sh
# Central configuration for the zebrafish safe-harbors pipeline
# All magic numbers and thresholds are defined here for reproducibility

# ============================================================================
# Consensus Peak Summit Parameters
# ============================================================================
# Summit window size (bp on each side of the narrowPeak col-10 summit)
# Default: 150bp on each side = 300bp total window
SUMMIT_WINDOW_HALF=150

# Consensus tiers (percentage of biosamples/stages required)
CONSENSUS_MIN_PCT_80=80
CONSENSUS_MIN_PCT_90=90
CONSENSUS_MIN_PCT_100=100

# ============================================================================
# Exclusion & Filtering Parameters
# ============================================================================
# Gene buffer distance (bp on each side of gene body)
GENE_BUFFER_KB=10

# Essential gene buffer (bp)
ESSENTIAL_GENE_BUFFER_KB=20

# Telomere distance from chromosome end (bp)
TELOMERE_END_DISTANCE=500000

# Centromere model (currently disabled - zebrafish midpoint is inaccurate)
CENTROMERE_USE=false

# Minimum safe harbor size (bp - must fit landing pad construct)
MIN_HARBOR_SIZE=300

# ============================================================================
# Annotation & Data Provenance
# ============================================================================
# Ensembl release version
ENSEMBL_RELEASE=110

# Genome build
GENOME_BUILD="GRCz11"

# UCSC RepeatMasker release (rmsk.txt source)
RMSK_SOURCE="https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz"

# ============================================================================
# Quality Control & Logging
# ============================================================================
# Enable verbose logging
VERBOSE=true

# Log file location (relative to script directory)
LOGFILE="${PWD}/pipeline.log"

# ============================================================================
# Export for use in other scripts
# ============================================================================
export SUMMIT_WINDOW_HALF
export CONSENSUS_MIN_PCT_80 CONSENSUS_MIN_PCT_90 CONSENSUS_MIN_PCT_100
export GENE_BUFFER_KB ESSENTIAL_GENE_BUFFER_KB TELOMERE_END_DISTANCE
export MIN_HARBOR_SIZE ENSEMBL_RELEASE GENOME_BUILD VERBOSE LOGFILE
