# Zebrafish Genomic Safe Harbor Discovery Pipeline

Identifies constitutively accessible, intergenic, gene-distal chromatin regions
in zebrafish (GRCz11) as candidate single-copy integration / landing-pad sites.
Constitutive accessibility is called from DANIO-CODE ATAC-seq; inertness is
screened with DANIO-CODE ChIP-seq, sequence conservation, and the DANIO-CODE
COPE/DOPE element tracks, and cross-checked for cell-type accessibility against
the Zebrahub single-nucleus atlas.

📊 **[Pipeline schematic & flowchart → docs/PIPELINE.md](docs/PIPELINE.md)**
(figure: [`results/pipeline_schematic.png`](results/pipeline_schematic.png))

## Pipeline Overview

```
01_setup_reference.sh       → Download/index GRCz11 genome
02_download_annotations.sh  → Download Ensembl + UCSC annotations
03_process_annotations.sh   → Process into BED format (genes, exons, repeats, TEs)
04_create_exclusions.py     → Generate exclusion region lists
04b_extract_supplement_annotations.py → Extract validated enhancers / promoters /
                              novel lncRNA from Baranasic 2022 supplement (danRer11)
04c_extract_ncrna_tracks.sh → miRNA / lncRNA tracks (eLife distance criteria)
04d_extract_conservation_tracks.sh → GERP constrained elements + Ancora CNEs
04e_extract_cope_dope.sh    → DANIO-CODE COPE/DOPE "open-but-inert" element tracks
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
11_design_sgrnas.py         → Enumerate on-target SpCas9 (NGG) guides for top candidates
11b_offtarget_bwa.py        → Local, autonomous off-target screen (bwa, exon-aware) +
                              CRISPRscan on-target efficiency (zebrafish-trained) →
                              safety-ranked results/sgrna_shortlist_auto.tsv (exon-clean
                              guides only; replaces manual CRISPOR — 15/16 optional)
11c_consolidate_shortlist.py → Audit designed loci: autonomous top pick vs the
                              injected guide + its rank, per locus (uses 11b's scorers)
12_relaxed_tier.sh          → Alternate candidate list at a relaxed gene-distance
                              buffer (2kb); reuses 08→08b→09 via env vars
13_mappability.sh           → genmap K-mer uniqueness mask (eLife GSH criterion)
14_intronic_tier.sh         → intron-permissive tier (exon-distal, not gene-distal)
15–18 (COPE loci)           → sgRNA / HDR donor / genotyping / GenBank maps for the two COPE targets
19_pipeline_figure.py       → pipeline schematic (docs/PIPELINE.md)
22–24 (oligo master)        → consolidated oligo order + per-locus design summary (all knock-in loci)
25_zebrahub_celltype_accessibility.py → pan-cell-type (snATAC) accessibility cross-check
26_make_igv_batch.sh        → IGV batch: per-candidate evidence snapshots
27–28 (KNIT)                → nickase (KNIT) biotinylated-donor variant (branch: knit-adaptation)
```

Scripts 08/08b/09 are parameterized by env vars (`MASTER_EXCL`, `HARBOR_PREFIX`,
`FEATURE_TABLE`) with defaults that preserve the strict pipeline; `12_relaxed_tier.sh`
sets them to emit the relaxed tier without touching strict outputs.

## Methodology

### Consensus Peak Calling (Phase 2)
1. Extract peak signal with quantile-normalized scores
2. Collapse replicates to biosample level (median score)
3. Collapse biosamples to developmental stage level (median score)
4. Create fixed windows (±150bp) around stage-level summits
5. Merge windows and count how many stages have a summit in each region
6. Filter by consensus tiers (80%, 90%, 100% of stages)

> **Note on scores & summits:** per-sample accessibility scores are quantile-normalized
> and median-collapsed, but the consensus tiers and the composite ranking currently use
> only summit **presence/absence** per stage — constitutive *strength* is not yet used to
> select or rank candidates. The ±150 bp window is centered on the merged-region midpoint
> (the per-sample narrowPeak col-10 summit is not preserved after stage merging).

### Exclusion Criteria
- Gene bodies ±10kb (all Ensembl genes) — enforced via master_exclusions.bed
- Exons (Ensembl 110)
- Essential genes ±20kb
- Transposable elements (DNA, LINE, SINE, LTR, Retroposon, RC classes)
- Telomere regions (500kb from chromosome ends)
- Non-primary chromosomes (primary chr1–25 only; drops unplaced scaffolds `chrUn_*`, alt/random, chrM)
- Minimum region size: 300bp

### Composite Ranking
Candidates are ranked by a composite score incorporating:
- Stage count (open in X/11 stages) — primary weight
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
| **Final safe harbors after Phase 3 + eLife (100%)** | **1,693** (strict) / **2,090** (relaxed) / **4,131** (intron-permissive) |

### Candidate tiers (gene-distance)
Three candidate lists are produced, differing in the gene/exon flanking buffer
(inertness / ChIP / conservation screen is identical):
- **Strict** — `GENE_BUFFER_KB=10` → 1,693 candidates → `candidate_feature_table.tsv`
- **Relaxed** — `RELAXED_GENE_BUFFER_KB=2` → 2,090 → `candidate_feature_table_relaxed.tsv`
  (build with `scripts/12_relaxed_tier.sh`)
- **Intron-permissive** — exon-distal (not whole-gene-distal) → 4,131 →
  `candidate_feature_table_intronic.tsv` (build with `scripts/14_intronic_tier.sh`)

Counts reflect the eLife-criteria regeneration: the ncRNA distance criteria
(miRNA ±50 kb, lncRNA ±25 kb) plus a scaffold fix took the strict list from
2,316 → 1,693.

The relaxed tier surfaces gene-proximal-but-inert sites (in the spirit of the
tpl102 benchmark, which works experimentally despite being ~341bp from a gene).
chr12/nrxn1a is #1 in both tiers. Note: tpl102 itself still doesn't surface even
at 2kb — it also overlaps a repeat (see Validation Benchmark).

The ChIP regulatory screen (discard-whole) takes 4,771 → 2,316 at the 100% tier;
the eLife ncRNA-distance criteria then reduce it to the final **1,693** (strict). ChIP breakdown:
ChIP active marks −2,360, validated enhancers −3, promoters −25, novel lncRNA −11,
H3K27me3 (Polycomb) −2, other active/transcribed ChIP −54. ~49% of gene-distal
constitutive peaks carry an active/poised-enhancer or insulator mark; almost none
are Polycomb-silenced (candidates are open chromatin) — confirming genuine inertness.

ChIP = **74 samples across 7 DANIO-CODE labs** (Skarmeta, Cairns, Mueller, Schier,
Bartel, Zon, Horsfield), embryonic stages only, 9 targets:
H3K27ac (13), H3K4me1 (8), H3K4me3 (23), H3K27me3 (4), CTCF (2), plus a combined
"other active/transcribed" track — H3K36me3 (3), H2A.Z/H2AFV (16), H3K14ac (3),
Pol II (2). **Developmental pioneer TFs (Pou5f3/Nanog/Sall4/Zic3/Cdx4/Mxtx2) were
evaluated and deliberately excluded** — they bind the early genome promiscuously
(Pou5f3 ~3.8% of the genome) and are unsuitable as a regulatory-element mask.

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

### Top 5 Candidates (main chromosomes, post full Phase 3 + eLife screen)

| # | Locus | Size | Nearest Gene | Distance | GC% |
|---|-------|------|--------------|----------|-----|
| 1 | chr12:24,289,888-24,291,839 | 1,951bp | nrxn1a | 50kb | 44.2% |
| 2 | chr24:3,128,791-3,130,466 | 1,675bp | inhbaa | 66kb | 48.0% |
| 3 | chr6:35,551,173-35,552,810 | 1,637bp | rgs8 | 63kb | 45.5% |
| 4 | chr14:27,882,540-27,884,172 | 1,632bp | ENSDARG00000075666 | 115kb | 46.2% |
| 5 | chr16:49,331,441-49,333,192 | 1,751bp | ENSDARG00000062621 | 60kb | 46.8% |

(Ranking refreshed after the eLife regeneration; chr20/ENSDARG00000096520 — previously
#2 — now falls just outside the top 5 but remains among the designed loci, see
`results/locus_design_summary.tsv`.)

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

**Final designed loci** (guide + HDR donor + genotyping/junction primers + WT/INTEGRATED
GenBank maps) are consolidated in `results/locus_design_summary.tsv` and
`results/oligo_master.tsv`. The six pursued sites are: **chr12/nrxn1a**, **chr20**,
**chr24/inhbaa** (top de novo GSH), **cope_chr1** and **cope_chr14** (COPE-anchored
constitutive elements), and **pIGLET24b_v2** (re-guided validated benchmark site).

> **NOTE — off-target tooling:** guide specificity AND on-target efficiency are now
> scored **autonomously and locally** by `scripts/11b_offtarget_bwa.py`: bwa aln → NGG-PAM
> filter → exon-aware off-target tiers, plus CRISPRscan (Moreno-Mateos 2015, zebrafish-
> injection-trained) on-target efficiency. Guides are ranked **safety-first** — an exonic
> off-target is a HARD exclusion, then fewest close off-targets → highest CRISPRscan →
> nearest ATAC-peak center — replacing the manual CRISPOR web round-trip (local
> `cas-offinder` was unreliable here — OpenCL CPU device absent). CRISPOR (MIT/CFD + Doench)
> stays optional via 15/16. `11_design_sgrnas.py` is the on-target enumerator that feeds 11b.

### Validation Benchmark
Three experimentally-validated attP/integration sites are benchmarked. All three are
in the 100% consensus (constitutively open) but excluded from the final candidate
list — and with the full ChIP screen loaded, *why* each is excluded is now explicit:

| Site | Phenotype | Active marks | Gene context | Why excluded |
|------|-----------|--------------|--------------|--------------|
| **pIGLET24b** | strong expression | H3K4me3 1000bp + H3K27ac 683bp + H3K4me1 340bp | gene-associated | **active chromatin** |
| pIGLET14a | weak expression | none | gene-associated | gene-associated |
| **tpl102** (Roberts 2014) | position-effect-*resistant* | minimal (25bp H3K4me3) | **341bp from ube2ia + repeat** | gene-proximal + repeat |

- **pIGLET24b** expresses strongly *because* it sits in active chromatin (full H3K4me3
  promoter + H3K27ac/H3K4me1 enhancer) — correctly rejected as "open but not inert."
- **tpl102** is largely inert (almost no active marks) yet works as a consistent,
  position-effect-resistant landing site **despite being only ~341bp from a gene**
  (`ube2ia`) and overlapping a repeat — so it fails the gene-distal/repeat criteria.

**Interpretation — not misses, by design.** None of the three validated sites lie in
the deep, fully-inert gene deserts this pipeline targets: pIGLET24b is active
chromatin, and pIGLET14a/tpl102 are gene-proximal. The pipeline optimizes *safety via
inertness*; the validated sites optimize *expression / position-effect resistance* by
other means. tpl102 is the most instructive: it shows a gene-proximal site can be a
good landing pad — i.e. our strict >10 kb gene-distance criterion may be conservative.
The trade-off cuts both ways (inert candidates like chr12/nrxn1a may express more
modestly), so the decisive test remains empirical.

## Important Limitations

### "Open ≠ Inert"
Constitutively accessible intergenic chromatin is not guaranteed to be
transcriptionally neutral. These regions may be:
- Distal enhancers (especially between developmental TF loci like irx, meis)
- Insulator/boundary elements
- Non-coding RNA loci not in current annotations

The Phase 3 ChIP screen now **actively distinguishes** "open and inert" from "open
and regulatory" — 74 DANIO-CODE ChIP samples (9 targets, active/insulator/Polycomb)
plus validated enhancers/promoters. pIGLET24b is the proof case: constitutively
open yet covered by active H3K4me3/H3K27ac, so it is correctly removed. The DANIO-CODE
COPE/DOPE "open-but-inert" tracks are now loaded and cross-referenced (19 candidates
overlap a COPE, 52 a DOPE); sequence conservation (GERP constrained elements + Ancora
CNEs) is implemented but data-gated. See DATA_NEEDED.md.

### Annotation Discrepancies
Ensembl gene models are often more inclusive than RefSeq (longer UTRs, alternative
transcripts). Some candidates may appear "inside" an Ensembl gene while being
clearly intergenic by RefSeq. The pIGLET24b site demonstrates this — it works
experimentally despite being within an Ensembl gene model.

### Phase 3 — Active (regulatory + COPE/DOPE) / Partial (conservation, mappability, GSH lists, TAD)
The regulatory element screen is **active**: all four ChIP marks
(H3K27ac, H3K4me1, H3K4me3, CTCF) plus supplement-derived validated enhancers,
promoters, and novel lncRNA are loaded and applied (discard-whole). This confirmed
the suspected enhancer role of the former #1 candidate near irx2a (H3K4me1+H3K4me3).

**DANIO-CODE COPE/DOPE are now loaded** and cross-referenced as a positive
"open-but-inert" corroboration (annotation, not exclusion): of the 1,693 strict
candidates, **19 overlap a COPE and 52 a DOPE** (`in_cope`/`in_dope`). The eLife
ncRNA distance criteria (miRNA ±50 kb, lncRNA ±25 kb) are also applied.

Filter code is implemented (GERP constrained elements + Ancora CNEs for conservation;
genmap K-mer mappability / blacklist / segmental duplications; cancer/dosage-sensitive/UCE
GSH gene lists; Hi-C TAD boundaries, annotation-only) but **remains a no-op until the
respective data tracks are downloaded** (gitignored; see DATA_NEEDED.md and scripts
04c–04e, 13). `conservation_frac` and `tad_distance` columns are present but currently
empty pending their data. CTCF is from a single stage (Horsfield "Bud", ~10 hpf) — the
only CTCF in DANIO-CODE.

### Pipeline Must Be Re-run After Changes
⚠️ Any changes to exclusion criteria, consensus parameters, or input data
require re-running scripts 07 → 08 → 08b → 09 → 10 → 11 → 11b. Do not design
sgRNAs from stale results.

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
| ATAC-seq | DANIO-CODE (Skarmeta + Zhang labs; 31 samples, 11 stages) |
| ChIP-seq | DANIO-CODE (74 samples, 7 labs; H3K27ac/K4me1/K4me3/K27me3/CTCF + others) |
| Conservation | GERP constrained elements + Ancora CNEs (danRer11) |
| Open-but-inert elements | DANIO-CODE COPE/DOPE track hub (Baranasic et al. 2022) |
| Mappability | genmap K-mer uniqueness (GRCz11) |
| Cell-type accessibility | Zebrahub snATAC (Bhattacharya et al. 2024) |
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

# sgRNA design — autonomous (local off-target, no manual CRISPOR)
python3 scripts/11_design_sgrnas.py 5      # enumerate on-target guides for the top 5
python3 scripts/11b_offtarget_bwa.py       # local off-target screen -> sgrna_shortlist_auto.tsv
```

## Repository Structure
```
zebrafish_safe_harbor/
├── scripts/          # Pipeline scripts (01-28)
├── config.sh         # Central configuration
├── annotation/       # Gene/repeat annotations (not in git - see setup)
│   ├── benchmarks/   # pIGLET coordinates
│   ├── chip/         # ChIP-seq peaks + per-mark union BEDs (not in git)
│   ├── regulatory/   # supplement-derived enhancer/promoter/lncRNA BEDs (not in git)
│   ├── conservation/ # GERP constrained elements + Ancora CNEs (not in git)
│   ├── cope_dope/    # DANIO-CODE COPE/DOPE element BEDs (not in git)
│   └── danio_code_accession_map.txt
├── peaks/            # ATAC-seq narrowPeak files (not in git)
├── reference/        # GRCz11 genome (not in git)
├── results/          # Pipeline output
│   ├── candidate_feature_table.tsv
│   └── sgrna_shortlist.tsv   # CRISPOR-scored guides for top candidates
├── DATA_NEEDED.md    # External data requirements
└── README.md
```

## Branch
The pipeline lives on `main`. The KNIT nickase-donor adaptation is staged on the
`knit-adaptation` branch (see `KNIT_ADAPTATION.md`), to be merged once implemented.