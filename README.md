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
04b_extract_supplement_annotations.py → Extract validated enhancers / promoters /
                              novel lncRNA from Baranasic 2022 supplement (danRer11)
05_combine_exclusions.sh    → Merge all exclusion BEDs (with chr normalization)
06_download_atac.sh         → Download DANIO-CODE ATAC-seq narrowPeak files
06b_normalize_scores.py     → Per-sample rank-normalize narrowPeak scores
06c_merge_chip.sh           → Merge per-mark ChIP samples into union peak BEDs
06d_download_chip.sh        → Download DANIO-CODE H3K4me1/H3K4me3/CTCF ChIP-seq
07_find_consensus_peaks.sh  → Summit-based consensus with replicate/stage collapsing
08_find_safe_harbors.sh     → Filter consensus peaks against exclusions
08b_biological_filters.sh   → Phase 3: regulatory/conservation/mappability filters
                              (discard-whole: any candidate touching a track is dropped)
09_annotate_candidates.sh   → Rich feature table with composite ranking
10_benchmark.sh             → Validation against pIGLET/SB known sites
11_design_sgrnas.py         → Enumerate on-target SpCas9 (NGG) guides for top
                              candidates (specificity/efficiency scored in CRISPOR;
                              picks in results/sgrna_shortlist.tsv)
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
| ATAC-seq samples | 31 (25 in accession map) |
| Developmental stages | 11 (0 – 48 hpf) |
| Stage range | Unfertilized egg → Long-pec (2 dpf) |
| Consensus peaks (100% / 11 stages) | 32,601 |
| Consensus peaks (90%) | 99,715 |
| Consensus peaks (80%) | 139,793 |
| Gene-distal + size-filtered (100%) | 4,771 |
| **Final safe harbors after Phase 3 (100%)** | **2,372** |

Phase 3 regulatory screen (discard-whole) removes 4,771 → 2,372 at the 100% tier:
ChIP marks −2,360, validated enhancers −3, promoters −25, novel lncRNA −11.
~49% of gene-distal constitutive peaks carry an active/poised-enhancer or insulator
mark. ChIP = **46 samples across 6 DANIO-CODE labs** (Skarmeta, Cairns, Mueller,
Schier, Bartel, Zon; Horsfield for CTCF), embryonic stages only — H3K27ac (13),
H3K4me1 (8), H3K4me3 (23), CTCF (2).

### Developmental Stages
| Stage | Time | Lab |
|-------|------|-----|
| Unfertilized egg | 0 hpf | Zhang |
| 2-cell | 1 hpf | Zhang |
| High | 2 hpf | Zhang |
| Mid-blastula / 1k-cell | 3 hpf | Zhang |
| Oblong | 3.7 hpf | Zhang |
| Dome | 4.33 – 4.5 hpf | Skarmeta + Zhang |
| Shield | 6 hpf | Skarmeta + Zhang |
| 75%-epiboly | 8 hpf | Skarmeta |
| 5-9 somites | 12 hpf | Skarmeta |
| Prim-5 | 24 hpf | Skarmeta |
| Long-pec | 48 hpf | Skarmeta |

Oblong (3.7 hpf) and a cross-lab Zhang Dome were added 2026-06-16 (peri-ZGA
wild-type ATAC). Pre-ZGA stages (64-/256-cell) were deliberately excluded — a
constitutive landing site need not be accessible before zygotic genome activation.

Stage assignments verified against Marlétaz et al. 2018 supplemental table
(DOI: 10.1038/s41586-018-0734-6) and DANIO-CODE data portal.

### Top 5 Candidates (main chromosomes, post full Phase 3 screen)

| # | Locus | Size | Nearest Gene | Distance | GC% |
|---|-------|------|--------------|----------|-----|
| 1 | chr12:24,289,888-24,291,839 | 1,951bp | nrxn1a | 50kb | 44.2% |
| 2 | chr20:11,445,664-11,447,454 | 1,790bp | ENSDARG00000096520 | 182kb | 45.5% |
| 3 | chr24:3,128,791-3,130,466 | 1,675bp | inhbaa | 66kb | 48.0% |
| 4 | chr6:35,551,173-35,552,810 | 1,637bp | rgs8 | 62kb | 45.5% |
| 5 | chr14:27,882,540-27,884,172 | 1,632bp | ENSDARG00000075666 | 114kb | 46.2% |

All candidates are:
- Open in 11/11 developmental stages (0 – 48 hpf, constitutive)
- >50kb from nearest Ensembl gene; zero repeat overlap
- Free of H3K27ac/H3K4me1/H3K4me3/CTCF ChIP signal, validated enhancers,
  promoters, and novel lncRNA (discard-whole regulatory screen)

**Note — prior top 5 changed after Phase 3:** the earlier top hits were filtered out.
`irx2a` carries H3K4me1 + H3K4me3 (poised enhancer + promoter — the suspected
regulatory role, now confirmed); `ovca2` overlaps a CTCF insulator; `vstm2b` is no
longer open across all 11 stages (lost constitutivity at the added Oblong stage).
`5S_rRNA` and ENSDARG00000103069 survive the full screen.

### Candidate sgRNA Guides (SpCas9, danRer11)

Top candidates were submitted to **CRISPOR** (danRer11, NGG PAM) for guide
specificity (MIT/CFD) + efficiency (Doench '16) scoring; off-targets were
cross-checked for in-gene (exon/intron) hits. Recommended guide per candidate
(full top-3 with all metrics in `results/sgrna_shortlist.tsv`):

| # | Locus | Guide (20nt) | PAM | MIT | Off-tgt | Doench | In-gene off-tgts |
|---|-------|--------------|-----|-----|---------|--------|------------------|
| 1 | chr12 nrxn1a | `TGCGCGGTTAGCGCTACTCG` | TGG | 100 | 2 | 64 | none ✅ |
| 5 | chr14 | `AGTTAACTCGGCGTACTTGG` | GGG | 99 | 3 | 64 | 1 intron |
| 4 | chr6 rgs8 | `GCTTAGCGGCAAGTCGCATT` | CGG | 100 | 4 | 50 | 4 intron, 0 exon |
| 3 | chr24 inhbaa | `AGGTCTGTCGACACGGATTG` | TGG | 99 | 7 | 66 | 2 exon (4mm, low score) |
| 2 | chr20 | `ACTTTCATGGGGACCCCCCT` | CGG | 100 | 4 | 58 | 1 exon; cut near edge |

- **chr12/nrxn1a** is the standout: perfect specificity, zero in-gene off-targets,
  central cut site.
- **chr20 and chr24** are weakest — no top guide fully avoids exonic off-targets
  (all such hits are 4-mismatch with low MIT off-target scores, so low risk).

> **NOTE — off-target tooling:** local `cas-offinder` was unreliable on this
> machine (OpenCL CPU device absent; the GPU backend produced spurious hits that
> failed genome verification). **Use CRISPOR for specificity**, not the local run.
> `scripts/11_design_sgrnas.py` is on-target enumeration only.

### Validation Benchmark
- **pIGLET24b** (experimentally validated, strong expression): In 100% consensus peaks ✓,
  but NOT in final candidate list (inside Ensembl gene ENSDARG00000079167/skilb)
- **pIGLET14a** (weaker expression): In 100% consensus peaks ✓ (was 80% before the
  11-stage expansion), still NOT in final candidate list (intronic)
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

### Phase 3 — Active (regulatory) / Partial (conservation, mappability)
The regulatory element screen is **active**: all four ChIP marks
(H3K27ac, H3K4me1, H3K4me3, CTCF) plus supplement-derived validated enhancers,
promoters, and novel lncRNA are loaded and applied (discard-whole). This confirmed
the suspected enhancer role of the former #1 candidate near irx2a (H3K4me1+H3K4me3).

Still **not loaded** (filters implemented, no-op until data added): phastCons/CNE
conservation, mappability/blacklist/segmental duplications, GSH gene lists
(cancer/dosage-sensitive/miRNA/UCE), Hi-C TAD boundaries, and DANIO-CODE COPE/DOPE.
CTCF is from a single stage (Horsfield "Bud", ~10 hpf) — the only CTCF in DANIO-CODE.

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
python3 scripts/04b_extract_supplement_annotations.py   # supplement-derived tracks
bash scripts/05_combine_exclusions.sh

# Regulatory ChIP tracks (Phase 3)
bash scripts/06d_download_chip.sh    # H3K4me1/H3K4me3/CTCF (H3K27ac pre-existing)
bash scripts/06c_merge_chip.sh       # merge per-mark union BEDs

# Normalize and find consensus
# NOTE: 06b reads UNCOMPRESSED *.narrowPeak in peaks/ — gunzip -k any new .gz first
python3 scripts/06b_normalize_scores.py
cd scripts && bash 07_find_consensus_peaks.sh && cd ..

# Find and annotate safe harbors
bash scripts/08_find_safe_harbors.sh
bash scripts/08b_biological_filters.sh
bash scripts/09_annotate_candidates.sh
bash scripts/10_benchmark.sh

## Repository Structure
zebrafish_safe_harbor/
├── scripts/          # Pipeline scripts (01-11)
├── config.sh         # Central configuration
├── annotation/       # Gene/repeat annotations (not in git - see setup)
│   ├── benchmarks/   # pIGLET coordinates
│   ├── chip/         # ChIP-seq peaks + per-mark union BEDs (not in git)
│   ├── regulatory/   # supplement-derived enhancer/promoter/lncRNA BEDs (not in git)
│   └── danio_code_accession_map.txt
├── peaks/            # ATAC-seq narrowPeak files (not in git)
├── reference/        # GRCz11 genome (not in git)
├── results/          # Pipeline output
│   ├── candidate_feature_table.tsv
│   └── sgrna_shortlist.tsv   # CRISPOR-scored guides for top candidates
├── DATA_NEEDED.md    # External data requirements
└── README.md

## Branch
This work is on branch pipeline-review-v2.