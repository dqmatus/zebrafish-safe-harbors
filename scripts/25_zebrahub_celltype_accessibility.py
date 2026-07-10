#!/usr/bin/env python3
"""
25_zebrahub_celltype_accessibility.py

Annotate safe-harbor candidates with PAN-CELL-TYPE accessibility from the
Zebrahub single-nucleus multiome snATAC (Bhattacharya et al. 2024; GRCz11;
timepoints 10/12/14/16/19/24 hpf).  Our bulk-ATAC criterion means "open across
developmental *stages*"; this adds the orthogonal "open across *cell types*"
axis, so a landing pad can be shown constitutive in both time and cell type
(especially at the 24 hpf scoring window).

INPUT — one of:
  A) --h5ad FILE     snATAC AnnData (cells x peaks). Auto-detects the peak
                     coordinates (var) and a cell-type obs column; optionally
                     subset to one timepoint. A peak is "accessible" in a cell
                     type if detected (count>0) in >= --access-threshold of that
                     type's cells.
  B) --peaks BED --matrix CSV   a peak BED (chrom,start,end) + a cell-type x peak
                     (or peak x cell-type) accessibility matrix; binarized at
                     --access-threshold.

OUTPUT: candidate table + columns  n_celltypes_accessible, n_celltypes_total,
        pan_celltype_frac  (fraction of cell types whose accessible peaks overlap
        the candidate).  Also writes a peaks->celltypes BED for reuse.

NOTE: the exact obs/var field names in the Zebrahub object may need a one-line
tweak once the file is in hand — see AUTODETECT blocks and --celltype-key /
--timepoint-key overrides.
"""
import argparse, subprocess, sys, os, re

def sh(cmd, **kw): return subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)

def norm_chrom(c):
    c = str(c)
    return c if c.startswith("chr") else "chr" + c

def parse_peak_id(s):
    """'chr1:100-200' or 'chr1-100-200' or 'chr1_100_200' -> (chr1,100,200) or None."""
    m = re.match(r"^(chr[^\W_]+|[0-9]+|[XYM]+|MT|KN\S+)[:\-_](\d+)[\-_](\d+)$", str(s))
    if not m: return None
    return norm_chrom(m.group(1)), int(m.group(2)), int(m.group(3))

# ---------------------------------------------------------------- accessibility
def peaks_from_h5ad(path, celltype_key, timepoint_key, timepoint, thresh):
    try:
        import anndata as ad
        import numpy as np
    except ImportError:
        sys.exit("ERROR: anndata/numpy needed for --h5ad. Install (e.g. in the zmarker/"
                 "cell-type env): pip install anndata scipy numpy  — or use --peaks/--matrix.")
    A = ad.read_h5ad(path)
    print(f"  loaded {A.shape[0]} cells x {A.shape[1]} peaks", file=sys.stderr)
    # AUTODETECT cell-type column
    if celltype_key is None:
        for k in ("cell_type", "celltype", "cell_type_broad", "annotation", "cell_states",
                  "cell_state", "leiden", "clusters", "seurat_clusters"):
            if k in A.obs: celltype_key = k; break
    if celltype_key not in A.obs:
        sys.exit(f"ERROR: cell-type column not found; obs has {list(A.obs.columns)}. Use --celltype-key.")
    # optional timepoint subset
    if timepoint:
        if timepoint_key is None:
            for k in ("timepoint", "stage", "age", "hpf", "time", "dev_stage"):
                if k in A.obs: timepoint_key = k; break
        if timepoint_key and timepoint_key in A.obs:
            A = A[A.obs[timepoint_key].astype(str).str.contains(str(timepoint))].copy()
            print(f"  subset to timepoint '{timepoint}': {A.shape[0]} cells", file=sys.stderr)
        else:
            print(f"  WARN: timepoint key not found; using all cells", file=sys.stderr)
    # peak coordinates from var
    coords = []
    for v in A.var_names:
        p = parse_peak_id(v)
        coords.append(p)
    if all(c is None for c in coords):  # maybe in var columns
        cols = {c.lower(): c for c in A.var.columns}
        cc = cols.get("chrom") or cols.get("chr") or cols.get("seqnames")
        sc = cols.get("start"); ec = cols.get("end")
        if cc and sc and ec:
            coords = [(norm_chrom(A.var[cc][i]), int(A.var[sc][i]), int(A.var[ec][i])) for i in range(A.n_vars)]
        else:
            sys.exit("ERROR: could not parse peak coordinates from var_names or var columns.")
    X = A.X
    import scipy.sparse as sp
    ct = A.obs[celltype_key].astype(str).values
    types = sorted(set(ct))
    # per cell type: detection fraction per peak (count>0), threshold -> accessible
    peak_ct = {}   # peak_index -> set(celltypes)
    for t in types:
        mask = ct == t
        n = int(mask.sum())
        if n == 0: continue
        sub = X[mask]
        det = (sub > 0)
        frac = np.asarray(det.sum(axis=0)).ravel() / n
        for j in np.where(frac >= thresh)[0]:
            peak_ct.setdefault(int(j), set()).add(t)
    rows = []
    for j, ts in peak_ct.items():
        c = coords[j]
        if c: rows.append((c[0], c[1], c[2], ",".join(sorted(ts))))
    return rows, types

def peaks_from_matrix(peaks_bed, matrix_csv, thresh):
    import csv
    peaks = [ln.split("\t")[:3] for ln in open(peaks_bed) if ln.strip() and not ln.startswith(("#", "track"))]
    import pandas as pd
    M = pd.read_csv(matrix_csv, index_col=0)
    # orient so rows=peaks, cols=celltypes (heuristic: match n_peaks)
    if M.shape[0] != len(peaks) and M.shape[1] == len(peaks):
        M = M.T
    types = list(M.columns)
    rows = []
    for i, (c, s, e) in enumerate(peaks):
        ts = [t for t in types if float(M.iloc[i][t]) >= thresh]
        if ts: rows.append((norm_chrom(c), int(s), int(e), ",".join(ts)))
    return rows, types

# ------------------------------------------------------------------ annotate
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad"); ap.add_argument("--peaks"); ap.add_argument("--matrix")
    ap.add_argument("--celltype-key"); ap.add_argument("--timepoint-key")
    ap.add_argument("--timepoint", default="24hpf", help="subset to this timepoint substring (h5ad only); '' = all")
    ap.add_argument("--access-threshold", type=float, default=0.10)
    ap.add_argument("--candidates", default="results/candidate_feature_table.tsv")
    ap.add_argument("--peaks-out", default="annotation/zebrahub/snatac_peaks_by_celltype.bed")
    ap.add_argument("--out", default="results/candidate_feature_table_zebrahub.tsv")
    a = ap.parse_args()

    tp = a.timepoint or None
    if a.h5ad:
        rows, types = peaks_from_h5ad(a.h5ad, a.celltype_key, a.timepoint_key, tp, a.access_threshold)
    elif a.peaks and a.matrix:
        rows, types = peaks_from_matrix(a.peaks, a.matrix, a.access_threshold)
    else:
        sys.exit("Provide --h5ad OR (--peaks and --matrix).")
    os.makedirs(os.path.dirname(a.peaks_out), exist_ok=True)
    with open(a.peaks_out, "w") as f:
        for c, s, e, ts in sorted(rows, key=lambda r: (r[0], r[1])):
            f.write(f"{c}\t{s}\t{e}\t{ts}\n")
    ntypes = len(types)
    print(f"  {len(rows)} accessible peaks across {ntypes} cell types -> {a.peaks_out}", file=sys.stderr)

    # candidate BED (cols 1-3), keep header + rows
    lines = open(a.candidates).read().splitlines()
    hdr = lines[0]
    with open("/tmp/_cand.bed", "w") as f:
        for ln in lines[1:]:
            c = ln.split("\t"); f.write(f"{c[0]}\t{c[1]}\t{c[2]}\n")
    # intersect candidate x peaks, collect union of cell types per candidate
    r = sh(f"sort -k1,1 -k2,2n /tmp/_cand.bed > /tmp/_c.bed; "
           f"sort -k1,1 -k2,2n '{a.peaks_out}' > /tmp/_p.bed; "
           f"bedtools intersect -a /tmp/_c.bed -b /tmp/_p.bed -wa -wb 2>/dev/null")
    from collections import defaultdict
    union = defaultdict(set)
    for ln in r.stdout.splitlines():
        p = ln.split("\t")
        key = (p[0], p[1], p[2])
        union[key].update(p[6].split(","))
    with open(a.out, "w") as f:
        f.write(hdr + "\tn_celltypes_accessible\tn_celltypes_total\tpan_celltype_frac\n")
        for ln in lines[1:]:
            c = ln.split("\t"); key = (c[0], c[1], c[2])
            n = len(union.get(key, set()))
            f.write(f"{ln}\t{n}\t{ntypes}\t{n/ntypes:.3f}\n")
    print(f"  wrote {a.out} (+ n_celltypes_accessible, pan_celltype_frac)", file=sys.stderr)

if __name__ == "__main__":
    main()
