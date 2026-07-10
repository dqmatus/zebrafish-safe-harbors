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
  A) --h5ad FILE     Zebrahub master pseudobulk object
                     `peaks_by_ct_tp_master_anno.h5ad` (peaks x conditions, where
                     each condition = cell-type x timepoint; ~640830 x 190).
                     Auto-detects orientation, parses peak coords (obs) and the
                     cell-type/timepoint of each condition (var), drops low-cell
                     conditions (var['n_cells'] < --min-cells), optionally keeps
                     one timepoint, and calls a peak "accessible" in a cell type
                     if its (log-norm) pseudobulk value > --access-threshold.
                     Also handles a plain cells x peaks AnnData if given.
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
def _coords_from_index(names, frame):
    """peak coords from an index of 'chr:start-end' ids, else chrom/start/end columns."""
    coords = [parse_peak_id(v) for v in names]
    if all(c is None for c in coords):
        cols = {c.lower(): c for c in frame.columns}
        cc = cols.get("chrom") or cols.get("chr") or cols.get("seqnames")
        sc = cols.get("start"); ec = cols.get("end")
        if cc and sc and ec:
            coords = [(norm_chrom(frame[cc][i]), int(frame[sc][i]), int(frame[ec][i])) for i in range(len(names))]
    return coords

def parse_condition(name):
    """condition 'celltype x timepoint' -> (celltype, hpf|None). Handles ...24hpf..., ..._24_..."""
    m = re.search(r"(\d+)\s*hpf", str(name), re.I) or re.search(r"[_|\-](10|12|14|16|19|24)([_|\-]|$)", str(name))
    tp = m.group(1) + "hpf" if m else None
    ct = re.sub(r"[_|\-]?\d+\s*hpf[_|\-]?", "_", str(name), flags=re.I).strip("_|- ") if tp else str(name)
    return (ct or str(name)), tp

def peaks_from_h5ad(path, celltype_key, timepoint_key, timepoint, thresh, min_cells):
    try:
        import anndata as ad
        import numpy as np
    except ImportError:
        sys.exit("ERROR: anndata/numpy needed for --h5ad. Use the zmarker env:\n"
                 "  /opt/homebrew/Caskroom/miniconda/base/envs/zmarker/bin/python scripts/25_...py --h5ad ...\n"
                 "or provide --peaks/--matrix instead.")
    A = ad.read_h5ad(path)
    print(f"  loaded AnnData {A.shape}", file=sys.stderr)

    # ---- detect orientation: pseudobulk master (peaks x conditions) vs cells x peaks ----
    obs_is_peaks = sum(parse_peak_id(v) is not None for v in A.obs_names[:20]) >= 10
    var_is_peaks = sum(parse_peak_id(v) is not None for v in A.var_names[:20]) >= 10

    if obs_is_peaks and not var_is_peaks:
        # PSEUDOBULK MASTER: obs=peaks, var=conditions (cell-type x timepoint)
        coords = _coords_from_index(A.obs_names, A.obs)
        M = A.layers["log_norm"] if "log_norm" in A.layers else A.X
        import numpy as np, scipy.sparse as sp
        M = M.toarray() if sp.issparse(M) else np.asarray(M)      # peaks x conditions
        ncells = A.var["n_cells"].values if "n_cells" in A.var else np.full(A.n_vars, min_cells)
        peak_ct = {}
        types = set()
        for j, cond in enumerate(A.var_names):
            if ncells[j] < min_cells: continue
            ct, tp = parse_condition(cond)
            if timepoint and tp and str(timepoint).rstrip("hpf") not in str(tp): continue
            types.add(ct)
            col = M[:, j]
            for i in np.where(col > thresh)[0]:
                peak_ct.setdefault(int(i), set()).add(ct)
        rows = [(coords[i][0], coords[i][1], coords[i][2], ",".join(sorted(cts)))
                for i, cts in peak_ct.items() if coords[i]]
        print(f"  pseudobulk master: {A.n_obs} peaks x {A.n_vars} conditions; "
              f"kept {len(types)} cell types (timepoint={timepoint or 'all'}, n_cells>={min_cells})", file=sys.stderr)
        return rows, sorted(types)

    # ---- else: cells x peaks AnnData (fallback) ----
    coords = _coords_from_index(A.var_names, A.var)
    if celltype_key is None:
        for k in ("cell_type", "celltype", "cell_type_broad", "annotation", "cell_states",
                  "cell_state", "leiden", "clusters"):
            if k in A.obs: celltype_key = k; break
    if celltype_key not in A.obs:
        sys.exit(f"ERROR: cell-type column not found; obs has {list(A.obs.columns)}. Use --celltype-key.")
    if timepoint:
        if timepoint_key is None:
            for k in ("timepoint", "stage", "age", "hpf", "time", "dev_stage"):
                if k in A.obs: timepoint_key = k; break
        if timepoint_key and timepoint_key in A.obs:
            A = A[A.obs[timepoint_key].astype(str).str.contains(str(timepoint))].copy()
    import numpy as np
    ct = A.obs[celltype_key].astype(str).values
    types = sorted(set(ct)); peak_ct = {}
    for t in types:
        mask = ct == t; n = int(mask.sum())
        if n == 0: continue
        frac = np.asarray((A.X[mask] > 0).sum(axis=0)).ravel() / n
        for j in np.where(frac >= thresh)[0]:
            peak_ct.setdefault(int(j), set()).add(t)
    rows = [(coords[j][0], coords[j][1], coords[j][2], ",".join(sorted(cts)))
            for j, cts in peak_ct.items() if coords[j]]
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
    ap.add_argument("--timepoint", default="24hpf", help="keep only this timepoint (e.g. 24hpf); '' = all")
    ap.add_argument("--min-cells", type=int, default=50, help="drop conditions with fewer nuclei (pseudobulk reliability)")
    ap.add_argument("--access-threshold", type=float, default=0.0,
                    help="pseudobulk: log-norm value above which a peak is 'accessible' (>0=detected). cells x peaks: detection fraction.")
    ap.add_argument("--candidates", default="results/candidate_feature_table.tsv")
    ap.add_argument("--peaks-out", default="annotation/zebrahub/snatac_peaks_by_celltype.bed")
    ap.add_argument("--out", default="results/candidate_feature_table_zebrahub.tsv")
    a = ap.parse_args()

    tp = a.timepoint or None
    if a.h5ad:
        rows, types = peaks_from_h5ad(a.h5ad, a.celltype_key, a.timepoint_key, tp, a.access_threshold, a.min_cells)
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
