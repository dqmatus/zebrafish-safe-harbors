#!/bin/bash
# 02_download_annotations.sh
# Downloads gene annotations and repeat annotations

set -e

echo "=== Downloading annotations ==="

# Create directory if needed
mkdir -p annotation
cd annotation

# Gene annotations from Ensembl
echo "Downloading gene annotations..."
curl -O https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz
gunzip Danio_rerio.GRCz11.110.gtf.gz

# Repeat annotations from UCSC
echo "Downloading repeat annotations..."
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz
gunzip rmsk.txt.gz

echo "=== Annotation downloads complete ==="