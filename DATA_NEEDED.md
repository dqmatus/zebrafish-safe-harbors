# DATA_NEEDED.md

## External Data Required for Full Pipeline Execution

This document tracks data files and resources that are not included in this repository but are required for full pipeline execution or specific phases.

### PHASE 3: Biological Criteria Strengthening

#### ChIP-seq Tracks (DANIO-CODE) — ✅ ACQUIRED
All four marks are downloaded and merged into per-mark union BEDs in
`annotation/chip/` by `scripts/06d_download_chip.sh` + `scripts/06c_merge_chip.sh`,
and applied (discard-whole) in 08b. Host: `danio-code-dcc.genereg.net` (the old
`danio-code.zfin.org/files/...` path now redirects to a homepage).

Embryonic stages from all DANIO-CODE labs (adult "90 Days-2 Years" excluded to
match the 0-48 hpf ATAC window). Downloaded by `06d_download_chip.sh`.
- **H3K27ac** (active enhancer) — ✅ 13 samples (Skarmeta/Cairns/Zon), 34,297 peaks. `CHIP_H3K27AC_BED`
- **H3K4me1** (poised enhancer) — ✅ 8 samples (Skarmeta/Cairns), 36,162 peaks. `CHIP_H3K4ME1_BED`
- **H3K4me3** (promoter) — ✅ 23 samples (Skarmeta/Mueller/Schier/Bartel/Zon/Cairns), 41,028 peaks. `CHIP_H3K4ME3_BED`
- **H3K27me3** (Polycomb-silenced; repressive) — ✅ 4 samples (Skarmeta/Schier), `CHIP_H3K27ME3_BED`.
  Distinct rationale: avoid landing in silenced chromatin (transgene could be silenced).
  Applied as 08b Step 1c. Low impact (−2 at 100% tier) since candidates are open chromatin.
- **CTCF** (insulator) — ✅ Horsfield Lab, 43,679 peaks. `CHIP_CTCF_BED`
- **Other active/transcribed** (combined, `CHIP_OTHER_BED`, 08b Step 1d) — ✅ 24 samples,
  ~22 Mb: H3K36me3 (3, transcribed bodies), H2A.Z/H2AFV (16), H3K14ac (3), Pol II (2).
  - Developmental pioneer TFs (Pou5f3/Nanog/Sall4/Zic3/Cdx4/Mxtx2) are available on
    DANIO-CODE but **deliberately NOT used** — they bind promiscuously (Pou5f3 ~3.8%
    of the genome), so they over-mask and are unsuitable as a regulatory-element filter.
  - ⚠️ CAVEAT: only "Bud" stage (~10 hpf) exists in all of DANIO-CODE; applied across all stages.
  - Note: Zon's H3K27ac file is named `H3K27Ac` (capital A); `06c` matches case-insensitively.

#### Supplement-derived tracks (Baranasic et al. 2022) — ✅ ACQUIRED
Extracted from `annotation/baranasic_2022_MOESM3.xlsx` by
`scripts/04b_extract_supplement_annotations.py` into `annotation/regulatory/`:
- **Validated enhancers** (T10/T11, 350 regions) — `VALIDATED_ENHANCERS_BED`
- **Promoters** (consensus T4 + core T13, 20,990 regions) — `PROMOTERS_BED`
- **Novel intergenic lncRNA/TUCP** (T3, 381 regions) — `NOVEL_LNCRNA_BED`

#### Conservation Tracks — ⚙️ IMPLEMENTED (data-gated)

Filter + extraction implemented (`scripts/04d_extract_conservation_tracks.sh`),
applied as exclusion in 08b Step 2; no-op until the data BEDs are present
(`conservation_frac` is currently empty in the candidate table).

- **GERP constrained elements** — `CONSERVATION_BED`
  - Source: Ensembl/UCSC GERP++ constrained elements (danRer11)
  - Usage: exclude highly constrained non-coding sequence (potential regulatory)
- **Ancora CNEs** — `CNE_BED`
  - Source: Ancora conserved non-coding elements (danRer11)
  - Usage: explicit CNE exclusion list

#### Hi-C / TAD Data

- **TAD Boundaries**
  - Status: NOT INCLUDED
  - Source: Published Hi-C for zebrafish (e.g., from Peng et al. or equivalent)
  - Usage: Flag candidates near TAD boundaries (potential position-effect variegation issues)
  - Config parameter: `TAD_BOUNDARIES_BED`
  - Note: Currently only for flagging, not exclusion

#### DANIO-CODE COPE/DOPE Tracks — ✅ ACQUIRED

Downloaded by `scripts/04e_extract_cope_dope.sh` from the DANIO-CODE track hub
(`danio-code-dcc.genereg.net`; `copes_dr11.bb` / `dopes_all_dr11.bb`) and
cross-referenced in 08b Step 6 / 09 as a positive "open-but-inert" annotation
(not an exclusion). Of the 1,693 strict candidates, 19 overlap a COPE and 52 a DOPE
(`in_cope` / `in_dope`).

- **COPE** — constitutively open chromatin lacking active CRE histone marks
  ("open but inert"). `COPE_BED`
- **DOPE** — the developmentally-regulated counterpart. `DOPE_BED`

### PHASE 3: Mappability & Blacklist

- **Mappability Map (genmap K-mer)** — ⚙️ IMPLEMENTED (data-gated)
  - Status: filter implemented (`scripts/13_mappability.sh`, genmap K-mer uniqueness);
    no-op until the mappability BED is generated/present
  - Source: genmap (danRer11) — the eLife GSH uniqueness criterion
  - Usage: Exclude low-mappability / multi-mapping regions
  - Config parameter: `MAPPABILITY_BED`

- **Segmental Duplication Regions**
  - Status: NOT INCLUDED
  - Source: UCSC Genome Browser or dbSNP
  - Usage: Exclude segmental duplications (potential off-target integration)
  - Config parameter: `SEGDUP_BED`

- **ENCODE Blacklist (adapted for danRer11)**
  - Status: NOT INCLUDED
  - Source: ENCODE blacklist or custom artifact regions
  - Usage: Exclude known problematic regions
  - Config parameter: `BLACKLIST_BED`

### PHASE 3: GSH Criteria Enhancement

- **miRNA/ncRNA Coordinates** — ✅ ACQUIRED
  - Status: extracted by `scripts/04c_extract_ncrna_tracks.sh`; the eLife ncRNA
    distance criteria (miRNA ±50 kb, lncRNA ±25 kb) are applied
  - Source: Ensembl BioMart (gene biotype = miRNA / lncRNA)
  - Config parameter: `MIRNA_NCRNA_BED`, `LNCRNA_BED`

- **Dosage-Sensitive Genes**
  - Status: NOT INCLUDED
  - Source: ClinVar, OMIM, or curated list
  - Usage: Extra buffer around dosage-sensitive genes to avoid haploinsufficiency
  - Config parameter: `DOSAGE_SENSITIVE_GENES_BED`

- **Cancer-Associated Genes**
  - Status: NOT INCLUDED
  - Source: COSMIC database
  - Usage: Extra distance from cancer-associated loci
  - Config parameter: `CANCER_ASSOCIATED_GENES_BED`

- **Ultraconserved Elements (UCEs)**
  - Status: NOT INCLUDED
  - Source: UCSC Genome Browser
  - Usage: Exclude ultraconserved elements (potential regulatory regions)
  - Config parameter: `UCE_BED`

### PHASE 4: sgRNA Targetability

- **sgRNA Pre-computation (CRISPOR)** — ✅ DONE
  - Status: on-target guides enumerated (`scripts/11_design_sgrnas.py`,
    `15_cope_target_guides.py`), scored in CRISPOR (MIT/CFD + Doench '16), curated in
    `results/sgrna_shortlist.tsv` and `results/locus_design_summary.tsv`
  - Config parameter: `SGRNA_BED`
  - Note: local cas-offinder was unreliable — CRISPOR is authoritative for specificity

### PHASE 5: Validation Benchmarks

#### pIGLET Coordinates — ✅ ACQUIRED

Coordinates in `annotation/benchmarks/` (`pIGLET14a_coords.bed`, `pIGLET24b_coords.bed`),
benchmarked by `scripts/10_benchmark.sh`. Both sit in the 100% consensus but are excluded
from the final list (pIGLET24b = active chromatin; pIGLET14a = gene-associated) — a
validation-by-design outcome, not a miss. A re-guided `pIGLET24b_v2` is among the designed
loci (`results/locus_design_summary.tsv`). Config: `PIGLET_14A_BED`, `PIGLET_24B_BED`.

#### Sleeping Beauty Founder Sites

- **SB Founder Insertion Coordinates**
  - Status: NOT INCLUDED
  - Source: Lab SB integration mapping / sequencing data
  - Usage: Cross-check recovered founders with pipeline predictions
  - Format: BED file with insertion sites
  - Config parameter: `SB_FOUNDER_SITES_BED`
  - Caveat: SB has TA-dinucleotide bias; only compare within TA-reachable regions

### PHASE 6: Documentation & Reproducibility

- **bedtools version**
  - Current assumption: bedtools >= 2.26
  - Recommended: Document exact version used (e.g., bedtools v2.31)
  - Impact: Minor differences in edge-case handling

- **DANIO-CODE Manifest with Checksums**
  - Status: NOT INCLUDED
  - Source: DANIO-CODE project metadata
  - Usage: Verify file integrity of downloaded ATAC-seq samples
  - Recommended: Compute and store SHA256 checksums for all ATAC files

### Summary Table

| Data Type | Phase | Status | Priority | Notes |
|-----------|-------|--------|----------|-------|
| H3K27ac, H3K4me1, H3K4me3, CTCF | 3 | ✅ Acquired | High | Regulatory exclusion (applied) |
| Conservation (GERP + Ancora CNE) | 3 | ⚙️ Implemented, data-gated | High | Exclusion in 04d / 08b |
| TAD boundaries | 3 | ⚙️ Implemented, data-gated | Medium | Flagging only (`tad_distance` empty) |
| COPE/DOPE | 3 | ✅ Acquired | Low | Cross-validation (19 COPE / 52 DOPE) |
| Mappability (genmap) | 3 | ⚙️ Implemented, data-gated | Medium | Reduce artifacts (script 13) |
| miRNA/ncRNA annotations | 3 | ✅ Acquired | High | eLife distance criteria applied (04c) |
| Dosage-sensitive / cancer / UCE | 3 | ⚙️ Buffers configured, data-gated | Low | Gene-list data still needed |
| sgRNA (CRISPOR) | 4 | ✅ Done | Medium | sgrna_shortlist / locus_design_summary |
| pIGLET14a/24b coords | 5 | ✅ Acquired | High | annotation/benchmarks/ |
| SB founder sites | 5 | Missing | Medium | Validation benchmark |
| DANIO-CODE checksums | 6 | Missing | Low | Reproducibility |

### How to Add Missing Data

1. **Download**: Obtain the data file (external source or lab data).
2. **Format**: Convert to standard BED format (chrom, start, end, [optional fields]).
3. **Normalize chr prefix**: Ensure consistent `chr`-naming (chr1 vs 1) before use.
4. **Update config.sh**: Add new parameter and file path.
5. **Update scripts**: Reference new parameter in Phase 3+ scripts.
6. **Update this file**: Mark as INCLUDED, note source and version.

### Questions or Data Availability?

For pIGLET, SB founder, or DANIO-CODE data inquiries, contact the lab or check:
- DANIO-CODE: https://www.danio-code.org
- pIGLET publications: Original manuscripts (doi pending)
- SB integration data: Lab records or supplementary materials

---

## OPTIONAL ENHANCEMENT: Zebrahub multiome snATAC (pan-cell-type accessibility)

Adds an orthogonal "open across *cell types*" axis to our bulk "open across
*stages*" criterion (Bhattacharya et al. 2024, bioRxiv 2024.10.18.618987; GRCz11).
10x multiome, timepoints **10/12/14/16/19/24 hpf** — includes the 24 hpf scoring
window with real cell-type diversity.

- Raw: SRA BioProject **PRJNA1164307**
- Processed (download one): snATAC **`.h5ad`**, or a pseudobulk **peak matrix
  (cell-type x peak) + peak BED** — Google Drive folder `1lk3VdHQuEiXTHZufMXp0C84T3pncqdCF`
  (see repo github.com/czbiohub-sf/zebrahub-multiome-analysis, notebooks
  `Fig_peak_umap/` and `EDA_peak_parts_list/`)
- Portal: https://zebrahub.ds.czbiohub.org/

Run (annotates candidate_feature_table with pan_celltype_frac):
  # .h5ad needs anndata -> use the zmarker env
  /opt/homebrew/Caskroom/miniconda/base/envs/zmarker/bin/python \
    scripts/25_zebrahub_celltype_accessibility.py --h5ad <snATAC.h5ad> --timepoint 24hpf
  # or from a matrix:
  python3 scripts/25_zebrahub_celltype_accessibility.py --peaks <peaks.bed> --matrix <celltype_x_peak.csv>

May need a one-line tweak for the object's cell-type obs key (--celltype-key) /
peak var format once the file is inspected.
