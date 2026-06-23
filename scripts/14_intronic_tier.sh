#!/bin/bash
# 14_intronic_tier.sh
# Produce the INTRON-PERMISSIVE candidate list: instead of excluding whole gene
# bodies (±GENE_BUFFER_KB), this tier excludes only EXONS (±INTRONIC_EXON_BUFFER_KB)
# plus essential genes, repeats, and telomeres — so deep-intronic / exon-distal
# sites of non-essential genes become eligible. The Phase 3 inertness screen
# (ChIP/regulatory/mappability/distance, 08b) is IDENTICAL to the strict tier;
# only the gene-distance geometry changes.
#
# Rationale (eLife GSH paper + collaborator review): many validated landing
# sites are intronic, and "≥10kb from ANY gene" discards a large, biologically
# reasonable search space. This tier surfaces it for comparison.
#
# CAVEAT: pIGLET24b is NOT recovered here — its integration point sits ~92bp from
# a skilb exon, so it fails any safe exon buffer. It is an exon-proximal site,
# not the deep-intronic class this tier targets (documented in config.sh).
#
# Reuses the parameterized 08 -> 08b -> 09 chain via env vars, writing
# tier-suffixed outputs; strict/relaxed outputs are left untouched:
#   results/safe_harbors_intronic_*            (08)
#   results/safe_harbors_intronic_biological_* (08b)
#   results/candidate_feature_table_intronic.tsv (09)
#
# Prereq: scripts/05_combine_exclusions.sh has been run (builds
#         annotation/master_exclusions_intronic.bed).

set -euo pipefail

if [ ! -s annotation/master_exclusions_intronic.bed ]; then
    echo "ERROR: annotation/master_exclusions_intronic.bed missing — run scripts/05_combine_exclusions.sh first."
    exit 1
fi

export MASTER_EXCL="../annotation/master_exclusions_intronic.bed"
export HARBOR_PREFIX="safe_harbors_intronic"
export FEATURE_TABLE="candidate_feature_table_intronic.tsv"

echo "=== Intron-permissive candidate list (exclude exons±INTRONIC_EXON_BUFFER_KB, not gene bodies) ==="
bash scripts/08_find_safe_harbors.sh
bash scripts/08b_biological_filters.sh
bash scripts/09_annotate_candidates.sh

echo ""
echo "=== Intron-permissive tier complete ==="
echo "  candidates: results/candidate_feature_table_intronic.tsv"
