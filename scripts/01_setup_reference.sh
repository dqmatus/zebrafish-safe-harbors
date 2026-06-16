#!/bin/bash
# 01_setup_reference.sh
# Downloads and indexes the zebrafish reference genome (GRCz11)

set -e

echo "=== Setting up reference genome ==="

# Create directory if needed
mkdir -p reference
cd reference

# Download GRCz11 from Ensembl
echo "Downloading genome..."
curl -O https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz

echo "Decompressing..."
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
mv Danio_rerio.GRCz11.dna.primary_assembly.fa GRCz11.fa

echo "Indexing for samtools..."
samtools faidx GRCz11.fa

echo "Creating chromosome sizes file..."
cut -f1,2 GRCz11.fa.fai > GRCz11.chrom.sizes

echo "Building bowtie2 index (this takes 30-60 minutes)..."
bowtie2-build --threads 8 GRCz11.fa GRCz11

echo "=== Reference setup complete ==="