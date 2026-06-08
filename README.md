# Zebrafish Safe Harbor Sites for CRISPR Landing Pads

## Project Overview
Identification of constitutively open chromatin regions in zebrafish (0-6 dpf) for CRISPR landing pad insertion.

## Final Safe Harbor Candidates

| Site | Coordinates (GRCz11) | Distance to Gene | Accessibility Score |
|------|---------------------|------------------|---------------------|
| 1 | chr15:47,250,490-47,253,818 | 56.9 kb | 85 |
| 2 | chr20:54,645,976-54,647,801 | 18.4 kb | 1,672 |
| 3 | chr11:15,324,566-15,328,259 | 17.4 kb | ~100 |

## Methods
- ATAC-seq data from DANIO-CODE (25 wild-type samples)
- Consensus peaks identified (open in 100% of samples)
- Filtered against genes, transposons, essential genes
- gRNAs designed using CRISPOR

## Repository Contents
- `scripts/` - Analysis pipeline scripts
- `results/` - Output files and candidate regions
- `annotation/` - Processed annotation files
- `docs/` - Summary documents

## Landing Pad Constructs
Three recombinase systems designed for RMCE:
1. Cre: loxP / lox2272
2. Flp: FRT / F3
3. Bxb1: attP-GT / attP-GA
