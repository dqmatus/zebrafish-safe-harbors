#!/usr/bin/env python3
"""
04_create_exclusions.py
Creates additional exclusion regions (telomeres, centromeres, essential genes)
"""

import pandas as pd
import os

# Change to annotation directory
os.chdir('annotation')

# -----------------------------------------
# Create chromosome end regions to avoid
# -----------------------------------------
print("Creating chromosome end regions...")

chrom_sizes = pd.read_csv(
    '../reference/GRCz11.chrom.sizes',
    sep='\t',
    header=None,
    names=['chrom', 'size']
)

# Keep only standard chromosomes (1-25)
standard_chroms = [str(i) for i in range(1, 26)]
chrom_sizes = chrom_sizes[chrom_sizes['chrom'].isin(standard_chroms)]

telomere_regions = []
centromere_regions = []

for _, row in chrom_sizes.iterrows():
    chrom = row['chrom']
    size = row['size']
    
    # Telomeres: 500kb from each end
    telomere_regions.append({'chrom': chrom, 'start': 0, 'end': 500000})
    telomere_regions.append({'chrom': chrom, 'start': max(0, size - 500000), 'end': size})
    
    # Centromere: 2Mb around chromosome center (approximate)
    center = size // 2
    centromere_regions.append({
        'chrom': chrom,
        'start': max(0, center - 2000000),
        'end': min(size, center + 2000000)
    })

pd.DataFrame(telomere_regions).to_csv(
    'telomere_regions.bed', sep='\t', header=False, index=False
)

pd.DataFrame(centromere_regions).to_csv(
    'centromere_regions.bed', sep='\t', header=False, index=False
)

print("Created telomere and centromere exclusion regions")

# -----------------------------------------
# Create essential genes list
# -----------------------------------------
print("Creating essential genes list...")

essential_prefixes = [
    'rpl', 'rps',      # Ribosomal proteins
    'eif', 'eef',      # Translation factors
    'atp',             # ATP synthase
    'nduf',            # Mitochondrial complex I
    'cox',             # Cytochrome oxidase
    'polr',            # RNA polymerase
    'mcm',             # DNA replication
    'pcna',            # DNA replication
    'hist',            # Histones
]

# Read genes
genes = pd.read_csv(
    'genes.bed',
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand']
)

def is_essential(gene_id):
    gene_lower = str(gene_id).lower()
    return any(gene_lower.startswith(prefix) for prefix in essential_prefixes)

essential = genes[genes['gene_id'].apply(is_essential)]
essential[['chrom', 'start', 'end', 'gene_id']].to_csv(
    'essential_genes.bed', sep='\t', header=False, index=False
)

print(f"Identified {len(essential)} essential genes")
print("=== Exclusion file creation complete ===")