# DATA_NEEDED.md

## External Data Required for Full Pipeline Execution

This document tracks data files and resources that are not included in this repository but are required for full pipeline execution or specific phases.

### PHASE 3: Biological Criteria Strengthening

#### ChIP-seq Tracks (DANIO-CODE) — ✅ ACQUIRED
All four marks are downloaded and merged into per-mark union BEDs in
`annotation/chip/` by `scripts/06d_download_chip.sh` + `scripts/06c_merge_chip.sh`,
and applied (discard-whole) in 08b. Host: `danio-code-dcc.genereg.net` (the old
`danio-code.zfin.org/files/...` path now redirects to a homepage).

- **H3K27ac** (active enhancer) — ✅ Skarmeta Lab, 10 samples, 31,051 peaks. `CHIP_H3K27AC_BED`
- **H3K4me1** (poised enhancer) — ✅ Skarmeta Lab, 4 stages, 32,511 peaks. `CHIP_H3K4ME1_BED`
- **H3K4me3** (promoter) — ✅ Skarmeta Lab, 6 stages × reps, 28,108 peaks. `CHIP_H3K4ME3_BED`
- **CTCF** (insulator) — ✅ Horsfield Lab, 43,679 peaks. `CHIP_CTCF_BED`
  - ⚠️ CAVEAT: only "Bud" stage (~10 hpf) exists in all of DANIO-CODE; applied across all stages.

#### Supplement-derived tracks (Baranasic et al. 2022) — ✅ ACQUIRED
Extracted from `annotation/baranasic_2022_MOESM3.xlsx` by
`scripts/04b_extract_supplement_annotations.py` into `annotation/regulatory/`:
- **Validated enhancers** (T10/T11, 350 regions) — `VALIDATED_ENHANCERS_BED`
- **Promoters** (consensus T4 + core T13, 20,990 regions) — `PROMOTERS_BED`
- **Novel intergenic lncRNA/TUCP** (T3, 381 regions) — `NOVEL_LNCRNA_BED`

#### Conservation Tracks

- **Conservation (PhyloCSF / phastCons)**
  - Status: NOT INCLUDED
  - Source: UCSC Genome Browser or PhyloCSF track for danRer11
  - Usage: Exclude highly conserved non-coding elements (potential regulatory regions)
  - Config parameter: `CONSERVATION_BED`
  - Recommendation: Use phastCons score >= 0.8 (highly conserved)

- **Conserved Non-Coding Elements (CNEs)**
  - Status: NOT INCLUDED
  - Source: Genomicus or UCSC Genome Browser
  - Usage: Explicit CNE exclusion list
  - Config parameter: `CNE_BED`

#### Hi-C / TAD Data

- **TAD Boundaries**
  - Status: NOT INCLUDED
  - Source: Published Hi-C for zebrafish (e.g., from Peng et al. or equivalent)
  - Usage: Flag candidates near TAD boundaries (potential position-effect variegation issues)
  - Config parameter: `TAD_BOUNDARIES_BED`
  - Note: Currently only for flagging, not exclusion

#### DANIO-CODE COPE/DOPE Tracks

- **COPE (Chromatin Open but lacking Promoter marks)**
  - Status: NOT INCLUDED
  - Source: DANIO-CODE derived or custom analysis
  - Usage: Identify truly "open but inert" regions (open chromatin without active marks)
  - Config parameter: `COPE_BED`

- **DOPE (Distal Open Chromatin)**
  - Status: NOT INCLUDED
  - Source: DANIO-CODE derived
  - Usage: Cross-check consensus peaks for distal regulatory signals
  - Config parameter: `DOPE_BED`

### PHASE 3: Mappability & Blacklist

- **Mappability Map (50-mer or 36-mer)**
  - Status: NOT INCLUDED
  - Source: GEM mappability (danRer11) or ENCODE blacklist derivative
  - Usage: Exclude low-mappability / multi-mapping regions
  - Config parameter: `MAPPABILITY_BED`
  - Threshold: Recommended min mappability >= 0.8

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

- **miRNA/ncRNA Coordinates**
  - Status: NOT INCLUDED (currently only protein-coding genes in filters)
  - Source: Ensembl BioMart (gene biotype = miRNA or ncRNA)
  - Usage: Gene ± buffer criteria should include small RNAs for complete "gene-distal" definition
  - Required file: `mirna_ncrna_coordinates.bed`
  - Config parameter: `GENE_BED_COMPREHENSIVE` (replace `genes.bed` with version including ncRNAs)

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

- **sgRNA Pre-computation (CRISPOR or FlashFry)**
  - Status: NOT INCLUDED
  - Source: CRISPOR web server or local FlashFry run
  - Usage: Filter candidates by sgRNA availability and off-target scores
  - Output format: BED with sgRNA score in column 5+
  - Config parameter: `SGRNA_BED`
  - Note: Must be pre-computed per candidate; recommend min specificity >= 0.95

### PHASE 5: Validation Benchmarks

#### pIGLET Coordinates

- **pIGLET14a insertion site**
  - Status: NOT INCLUDED
  - Source: Original pIGLET manuscript or lab provenance
  - Usage: Validate that pIGLET14a ranks lower than pIGLET24b
  - Expected: chr?, ?-? bp
  - Config parameter: `PIGLET_14A_BED`

- **pIGLET24b insertion site**
  - Status: NOT INCLUDED
  - Source: pIGLET manuscript / lab data
  - Usage: Validate that pIGLET24b ranks highly (should be in 100% consensus)
  - Expected: chr?, ?-? bp
  - Config parameter: `PIGLET_24B_BED`

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
| H3K27ac, H3K4me1, H3K4me3, CTCF | 3 | Missing | High | Regulatory element exclusion |
| Conservation / CNE | 3 | Missing | High | Distal enhancer exclusion |
| TAD boundaries | 3 | Missing | Medium | Flagging only |
| COPE/DOPE | 3 | Missing | Low | Cross-validation |
| Mappability / Blacklist | 3 | Missing | Medium | Reduce artifacts |
| miRNA/ncRNA annotations | 3 | Missing | High | Complete gene-distal def. |
| Dosage-sensitive / cancer genes | 3 | Missing | Low | Extra safety margin |
| UCEs | 3 | Missing | Low | Conservative filtering |
| sgRNA pre-computation | 4 | Missing | Medium | Testability |
| pIGLET14a/24b coords | 5 | Missing | High | Validation benchmark |
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
