#!/bin/bash
# 12_relaxed_tier.sh
# Produce the ALTERNATE, more-permissive candidate list using a relaxed gene
# flanking buffer (config: RELAXED_GENE_BUFFER_KB, default 2kb) instead of the
# strict GENE_BUFFER_KB (10kb). Inertness (Phase 3 ChIP/regulatory screen) is
# identical to the strict tier — only the gene-distance margin is relaxed.
#
# Reuses the parameterized 08 -> 08b -> 09 chain via env vars, writing
# tier-suffixed outputs; the strict outputs are left untouched:
#   results/safe_harbors_relaxed_*            (08)
#   results/safe_harbors_relaxed_biological_* (08b)
#   results/candidate_feature_table_relaxed.tsv (09)
#
# Prereq: scripts/05_combine_exclusions.sh has been run (builds
#         annotation/master_exclusions_relaxed.bed).

set -euo pipefail

if [ ! -s annotation/master_exclusions_relaxed.bed ]; then
    echo "ERROR: annotation/master_exclusions_relaxed.bed missing — run scripts/05_combine_exclusions.sh first."
    exit 1
fi

export MASTER_EXCL="../annotation/master_exclusions_relaxed.bed"
export HARBOR_PREFIX="safe_harbors_relaxed"
export FEATURE_TABLE="candidate_feature_table_relaxed.tsv"

echo "=== Relaxed-tier candidate list (gene buffer = RELAXED_GENE_BUFFER_KB) ==="
bash scripts/08_find_safe_harbors.sh
bash scripts/08b_biological_filters.sh
bash scripts/09_annotate_candidates.sh

echo ""
echo "=== Relaxed tier complete ==="
echo "  candidates: results/candidate_feature_table_relaxed.tsv"
