#!/usr/bin/env python3
"""
04b_extract_supplement_annotations.py

Extract danRer11/GRCz11 coordinate tracks from the Baranasic et al. 2022
(Nat Genet, DANIO-CODE multiomic atlas) supplementary tables into BED files
for the Phase 3 regulatory screen (08b).

Source: annotation/baranasic_2022_MOESM3.xlsx (Supplementary Tables 1-15).
All coordinate tables in that file are already in danRer11.

Tracks produced (in annotation/regulatory/):
  validated_enhancers.bed  <- T10 (transgenic reporter-validated enhancers)
                              + T11 (validated enhancers active at Prim-5)
  promoters.bed            <- T4 (consensus promoters, full regions)
                              + T13 (core-promoter dominant TSS, buffered +/-TSS_PAD)
  novel_lncrna.bed         <- T3 (novel TUCP/lncRNA gene bodies, intergenic)

Each output is chr-normalized, sorted, and merged.
"""
import os
import subprocess
import sys
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUPP = os.path.join(REPO, "annotation", "baranasic_2022_MOESM3.xlsx")
OUTDIR = os.path.join(REPO, "annotation", "regulatory")
TSS_PAD = 250  # bp padding applied to 1-bp core-promoter TSS points (T13)

os.makedirs(OUTDIR, exist_ok=True)
if not os.path.exists(SUPP):
    sys.exit(f"ERROR: supplement not found: {SUPP}")

xl = pd.ExcelFile(SUPP, engine="openpyxl")


def rows_from(sheet, chrom_col, start_col, end_col, header, pad=0):
    """Yield (chrom,start,end) from a sheet, chr-normalized, padded, validated."""
    df = xl.parse(sheet, header=header)
    out = []
    for _, r in df.iterrows():
        try:
            c = str(r.iloc[chrom_col]).strip()
            s = int(float(r.iloc[start_col])) - pad
            e = int(float(r.iloc[end_col])) + pad
        except (ValueError, TypeError, IndexError):
            continue
        if c in ("nan", ""):
            continue
        if not c.startswith("chr"):
            c = "chr" + c
        s = max(0, s)
        if e > s:
            out.append((c, s, e))
    return out


def write_merged(intervals, out_path, label):
    """Write intervals as sorted+merged BED via bedtools."""
    tmp = out_path + ".tmp"
    with open(tmp, "w") as fh:
        for c, s, e in intervals:
            fh.write(f"{c}\t{s}\t{e}\n")
    with open(out_path, "w") as fh:
        p1 = subprocess.Popen(["sort", "-k1,1", "-k2,2n", tmp], stdout=subprocess.PIPE)
        subprocess.run(["bedtools", "merge"], stdin=p1.stdout, stdout=fh, check=True)
        p1.stdout.close()
        p1.wait()
    os.remove(tmp)
    n = sum(1 for _ in open(out_path))
    bp = sum((int(l.split('\t')[2]) - int(l.split('\t')[1])) for l in open(out_path))
    print(f"  {label}: {n} regions, {bp} bp -> {os.path.relpath(out_path, REPO)}")


print("=== Extracting Baranasic 2022 supplement tracks (danRer11) ===")

# Validated enhancers: T10 (chrom,start,end @ cols 0,1,2) + T11 (seqnames,start,end @ 0,1,2)
enh = rows_from("Supplementary_table_10", 0, 1, 2, header=1)
enh += rows_from("Supplementary_table_11", 0, 1, 2, header=1)
write_merged(enh, os.path.join(OUTDIR, "validated_enhancers.bed"), "validated enhancers (T10+T11)")

# Promoters: T4 consensus promoters (Chromosome,Start,End @ 0,1,2 full regions)
#          + T13 core-promoter dominant TSS (chromosome,domTssStart,domTssEnd @ 1,2,3), padded
prom = rows_from("Supplementary_table_4", 0, 1, 2, header=1)
prom += rows_from("Supplementary_table_13", 1, 2, 3, header=1, pad=TSS_PAD)
write_merged(prom, os.path.join(OUTDIR, "promoters.bed"), f"promoters (T4 + T13 TSS+/-{TSS_PAD}bp)")

# Novel lncRNA/TUCP gene bodies: T3 (seqnames @ col2, start @ col6, end @ col7)
lnc = rows_from("Supplementary_table_3", 2, 6, 7, header=1)
write_merged(lnc, os.path.join(OUTDIR, "novel_lncrna.bed"), "novel lncRNA/TUCP (T3)")

print("=== Supplement track extraction complete ===")
