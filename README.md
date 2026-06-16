# Zebrafish Genomic Safe Harbor Discovery Pipeline

Identifies constitutively accessible, intergenic, gene-distal chromatin regions
in zebrafish (GRCz11) as candidate single-copy integration / landing-pad sites,
using DANIO-CODE ATAC-seq data.

## Pipeline Overview

```
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
```

## Methodology

### Consensus Peak Calling (Phase 2)
1. Extract peak summits (narrowPeak col-10) with quantile-normalized scores
2. Collapse replicates to biosample level (median score)
3. Collapse biosamples to developmental stage level (median score)
4. Create fixed windows (±150bp) around stage-level summits
5. Merge windows and count how many stages have a summit in each region
6. Filter by consensus tiers (80%, 90%, 100% of stages)

### Exclusion Criteria
- Gene bodies ±10kb (all Ensembl genes) — enforced via master_exclusions.bed
- Exons (Ensembl 110)
- Essential genes ±20kb
- Transposable elements (DNA, LINE, SINE, LTR, Retroposon, RC classes)
- Telomere regions (500kb from chromosome ends)
- Minimum region size: 300bp

### Composite Ranking
Candidates are ranked by a composite score incorporating:
- Stage count (open in X/10 stages) — primary weight
- Region size (larger = better for construct design)
- Distance from nearest gene (+50 for >50kb, +30 for >10kb, -200 for inside gene)
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
| Developmental stages | 10 (0 – 48 hpf) |
| Stage range | Unfertilized egg → Long-pec (2 dpf) |
| Consensus peaks (100%) | ~31,000 |
| Consensus peaks (90%) | ~60,000 |
| Consensus peaks (80%) | ~94,000 |
| Safe harbors after filtering (100%) | 4,503 |

### Developmental Stages
| Stage | Time | Lab | Biosamples |
|-------|------|-----|------------|
| Unfertilized egg | 0 hpf | Zhang | 1 |
| 2-cell | 1 hpf | Zhang | 1 (2 reps) |
| High | 2 hpf | Zhang | 1 (2 reps) |
| Mid blastula | 3 hpf | Zhang | 1 |
| Dome | 4.33 hpf | Skarmeta | 1 (2 reps) |
| Shield | 6 hpf | Skarmeta + Zhang | 2 (4 reps) |
| 75%-epiboly | 8 hpf | Skarmeta | 1 (2 reps) |
| 5-9 somites | 12 hpf | Skarmeta | 1 (2 reps) |
| Prim-5 | 24 hpf | Skarmeta | 1 (2 reps) |
| Long-pec | 48 hpf | Skarmeta | 1 (2 reps) |

Stage assignments verified against Marlétaz et al. 2018 supplemental table
(DOI: 10.1038/s41586-018-0734-6) and DANIO-CODE data portal.

### Top 5 Candidates (IGV-validated, main chromosomes, ≥800bp)

| # | Locus | Size | Nearest Gene | Distance | GC% |
|---|-------|------|--------------|----------|-----|
| 1 | chr16:683,375-684,868 | 1,493bp | irx2a | 120kb | 51.6% |
| 2 | chr7:45,155,247-45,156,838 | 1,591bp | vstm2b | 79kb | 42.2% |
| 3 | chr14:48,494,914-48,496,373 | 1,459bp | 5S_rRNA | 80kb | 46.2% |
| 4 | chr15:47,252,117-47,253,549 | 1,432bp | ENSDARG00000103069 | 59kb | 43.6% |
| 5 | chr10:35,762,801-35,764,213 | 1,412bp | ovca2 | 94kb | 40.2% |

All candidates are:
- Open in 10/10 developmental stages (0 – 48 hpf, constitutive)
- Located in gene deserts >200kb
- >50kb from nearest Ensembl gene
- Zero repeat overlap
- GC content 40-52%

### Validation Benchmark
- **pIGLET24b** (experimentally validated, strong expression): In 100% consensus peaks ✓,
  but NOT in final candidate list (inside Ensembl gene ENSDARG00000079167/skilb)
- **pIGLET14a** (weaker expression): Only 80% consensus, NOT in final candidate list
- **Result: INCONCLUSIVE** — pIGLET24b is excluded by strict gene-distal criteria despite
  being a validated functional site. This highlights that intronic sites can work
  experimentally but don't meet conservative GSH definitions.

## Important Limitations

### "Open ≠ Inert"
Constitutively accessible intergenic chromatin is not guaranteed to be
transcriptionally neutral. These regions may be:
- Distal enhancers (especially between developmental TF loci like irx, meis)
- Insulator/boundary elements
- Non-coding RNA loci not in current annotations

Phase 3 data (ChIP-seq, conservation) is needed to distinguish "open and inert"
from "open and regulatory." The DANIO-CODE COPE/DOPE tracks (open chromatin
lacking active histone marks) would be particularly valuable. See DATA_NEEDED.md.

### Annotation Discrepancies
Ensembl gene models are often more inclusive than RefSeq (longer UTRs, alternative
transcripts). Some candidates may appear "inside" an Ensembl gene while being
clearly intergenic by RefSeq. The pIGLET24b site demonstrates this — it works
experimentally despite being within an Ensembl gene model.

### Phase 3 Not Yet Active
The regulatory element screen (H3K27ac, CTCF, conservation, mappability) is
implemented but has no data loaded. Candidates near developmental TF loci
(e.g., #1 near irx2a) may be distal enhancers. Adding DANIO-CODE ChIP tracks
is the highest-priority next step.

### Pipeline Must Be Re-run After Changes
⚠️ Any changes to exclusion criteria, consensus parameters, or input data
require re-running scripts 07 → 08 → 08b → 09 → 10. Do not design sgRNAs
from stale results.

## Configuration

All parameters are centralized in `config.sh`:
- Summit window: ±150bp
- Consensus tiers: 80/90/100%
- Gene buffer: 10kb (via master_exclusions.bed)
- Essential gene buffer: 20kb
- Telomere exclusion: 500kb
- Minimum harbor size: 300bp
- Score normalization: quantile

## Data Provenance

| Resource | Version/Source |
|----------|---------------|
| Genome | GRCz11 (danRer11) |
| Gene annotation | Ensembl release 110 |
| Repeat annotation | UCSC RepeatMasker (rmsk.txt) |
| ATAC-seq | DANIO-CODE (Skarmeta Lab + Zhang Lab) |
| Stage verification | Marlétaz et al. 2018 (10.1038/s41586-018-0734-6) |

## Requirements

- bedtools (≥2.26)
- Python 3 with numpy, pandas
- bash with pipefail support
- conda environment: `zebrafish_atac`

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

## Repository Structure
zebrafish_safe_harbor/
├── scripts/          # Pipeline scripts (01-10)
├── config.sh         # Central configuration
├── annotation/       # Gene/repeat annotations (not in git - see setup)
│   ├── benchmarks/   # pIGLET coordinates
│   └── danio_code_accession_map.txt
├── peaks/            # ATAC-seq narrowPeak files (not in git)
├── reference/        # GRCz11 genome (not in git)
├── results/          # Pipeline output
│   └── candidate_feature_table.tsv
├── DATA_NEEDED.md    # External data requirements
└── README.md

## Branch
This work is on branch pipeline-review.