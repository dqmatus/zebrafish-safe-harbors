# Zebrafish Genomic Safe Harbor Discovery Pipeline

Identifies constitutively accessible, intergenic, gene-distal chromatin regions
in zebrafish (GRCz11) as candidate single-copy integration / landing-pad sites,
using DANIO-CODE ATAC-seq data.

## Pipeline Overview
01_setup_reference.sh       → Download/index GRCz11 genome
02_download_annotations.sh  → Download Ensembl + UCSC annotations
03_process_annotations.sh   → Process into BED format (genes, exons, repeats, TEs)
04_create_exclusions.py     → Generate exclusion region lists
05_combine_exclusions.sh    → Merge all exclusion BEDs (with chr normalization)
06_download_atac.sh         → Download DANIO-CODE ATAC-seq narrowPeak files
06b_normalize_scores.py     → Quantile-normalize scores across samples
07_find_consensus_peaks.sh  → Summit-based consensus with replicate/stage collapsing
08_find_safe_harbors.sh     → Filter consensus peaks against exclusions
08b_biological_filters.sh   → Phase 3: regulatory/conservation/mappability filters
09_annotate_candidates.sh   → Rich feature table with composite ranking
10_benchmark.sh             → Validation against pIGLET/SB known sites

## Methodology

### Consensus Peak Calling (Phase 2 redesign)
1. Extract peak summits (narrowPeak col-10) with quantile-normalized scores
2. Collapse replicates to biosample level (median score)
3. Collapse biosamples to developmental stage level (median score)
4. Create fixed windows (±150bp) around stage-level summits
5. Merge windows and count how many stages have a summit in each region
6. Filter by consensus tiers (80%, 90%, 100% of stages)

### Exclusion Criteria
- Exons (Ensembl 110)
- Essential genes ±20kb
- Transposable elements (DNA, LINE, SINE, LTR, Retroposon, RC classes)
- Gene bodies ±10kb (all Ensembl genes)
- Telomere regions (500kb from chromosome ends)
- Minimum region size: 300bp

### Composite Ranking
Candidates are ranked by a composite score incorporating:
- Stage count (open in X/10 stages) — primary weight
- Region size (larger = better for construct design)
- Distance from nearest gene (penalty for distance=0, bonus for >50kb)
- sgRNA availability (when data provided)
- Conservation overlap (penalty, when data provided)
- Repeat content (penalty)

### Biological Criteria (Phase 3 — requires external data)
When data is provided (see DATA_NEEDED.md), additional filters are applied:
- H3K27ac/H3K4me1/H3K4me3/CTCF exclusion (regulatory elements)
- Conservation/CNE exclusion (functional non-coding)
- Mappability/blacklist/segmental duplication removal
- GSH-style distance criteria (cancer genes, dosage-sensitive, miRNA/ncRNA, UCEs)

## Results

### Current Run Statistics
| Metric | Value |
|--------|-------|
| ATAC-seq samples | 25 (20 in accession map) |
| Biosamples | 11 |
| Developmental stages | 10 (0hpf – 12hpf) |
| Consensus peaks (100%) | 31,617 |
| Consensus peaks (90%) | 60,178 |
| Consensus peaks (80%) | 94,390 |
| Safe harbors after filtering (100%) | 20,075 |
| Safe harbors after filtering (90%) | 39,566 |
| Safe harbors after filtering (80%) | 61,563 |

### Top 5 Candidates (IGV-validated)

| # | Locus | Size | Nearest Gene | Distance | GC% |
|---|-------|------|--------------|----------|-----|
| 1 | chr16:683,375-684,873 | 1,498bp | irx2a | 120kb | 48.3% |
| 2 | chr17:47,977,414-47,979,001 | 1,587bp | slc39a9 | 148kb | 49.7% |
| 3 | chr20:4,566,962-4,568,172 | 1,210bp | ENSDARG00000116556 | 112kb | 50.3% |
| 4 | chr18:23,726,829-23,728,108 | 1,279bp | ENSDARG00000097764 | 102kb | 51.4% |
| 5 | chr13:5,731,152-5,732,028 | 876bp | ENSDARG00000095303 | 137kb | 47.6% |

All candidates are:
- Open in 10/10 developmental stages (constitutive)
- Located in gene deserts >250kb
- >100kb from nearest Ensembl gene
- Zero repeat overlap
- GC content 47-51%

### Validation Benchmark
- **pIGLET24b** (experimentally validated, strong expression): In 100% consensus ✓, ranked in final candidate list ✓
- **pIGLET14a** (weaker expression): Only 80% consensus, NOT in final candidate list ✓
- Pipeline correctly ranks 24b above 14a: **PASS**

## Important Limitations

### "Open ≠ Inert"
Constitutively accessible intergenic chromatin is not guaranteed to be
transcriptionally neutral. These regions may be:
- Distal enhancers (especially between developmental TF loci like irx, meis)
- Insulator/boundary elements
- Non-coding RNA loci not in current annotations

Phase 3 data (ChIP-seq, conservation) is needed to distinguish "open and inert"
from "open and regulatory." See DATA_NEEDED.md.

### Annotation Discrepancies
Ensembl gene models are often more inclusive than RefSeq (longer UTRs, alternative
transcripts). Some candidates may appear "inside" an Ensembl gene while being
clearly intergenic by RefSeq. The pIGLET24b site demonstrates this — it works
experimentally despite being within an Ensembl gene model.

### Pipeline Must Be Re-run
⚠️ **These results reflect the corrected pipeline.** The bug fixes in Phase 1
(chr-naming mismatch, transposon class filter, gene-distal enforcement) and the
Phase 2 methodology changes (summit consensus, replicate collapsing, stage
stratification) fundamentally change what gets filtered. **Do not use candidates
from previous runs.** Design sgRNAs only from the current `candidate_feature_table.tsv`.

## Configuration

All parameters are centralized in `config.sh`:
- Summit window: ±150bp
- Consensus tiers: 80/90/100%
- Gene buffer: 10kb
- Essential gene buffer: 20kb
- Telomere exclusion: 500kb
- Minimum harbor size: 300bp
- Score normalization: quantile

## Data Provenance

| Resource | Version/Date |
|----------|-------------|
| Genome | GRCz11 (danRer11) |
| Gene annotation | Ensembl release 110 |
| Repeat annotation | UCSC RepeatMasker (rmsk.txt) |
| ATAC-seq | DANIO-CODE (Skarmeta Lab + Zhang Lab) |

## Requirements

- bedtools (≥2.26)
- Python 3 with numpy, pandas
- bash with pipefail support
- ~5GB disk space for annotations + peaks

See `DATA_NEEDED.md` for external data required for full Phase 3-5 execution.

## Running the Pipeline

```bash
conda activate zebrafish_atac
cd zebrafish_safe_harbor

# Process annotations (if re-running from scratch)
bash scripts/03_process_annotations.sh
bash scripts/05_combine_exclusions.sh

# Normalize and find consensus
python3 scripts/06b_normalize_scores.py
cd scripts && bash 07_find_consensus_peaks.sh && cd ..

# Find and annotate safe harbors
bash scripts/08_find_safe_harbors.sh
bash scripts/08b_biological_filters.sh
bash scripts/09_annotate_candidates.sh
bash scripts/10_benchmark.sh