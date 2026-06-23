#!/usr/bin/env python3
"""
16_curate_crispor_shortlist.py

Curate a specificity-ranked sgRNA shortlist from CRISPOR output for one COPE
target, cross-referencing the on-target guides enumerated by 15_cope_target_guides.py.

CRISPOR's guides/offtargets downloads are binary .xls with a metadata preamble;
the real table header is on row 9 (pandas header=8).

Ranking = SAFETY first (fewest close [<=3 mismatch] off-targets in exons, then
fewest close off-targets overall), then MIT specificity, then Doench '16
efficiency. Restricted to guides whose blunt cut falls inside the COPE.

Usage:
  python3 scripts/16_curate_crispor_shortlist.py LABEL GUIDES.xls OFFTARGETS.xls OUT.tsv
    LABEL = e.g. cope_chr1 (must match a label in cope_guides_chr1_chr14.tsv)
"""
import sys, csv
import pandas as pd

label, guides_xls, offt_xls, out_tsv = sys.argv[1:5]
COPE_GUIDES = "results/cope_guides_chr1_chr14.tsv"

g = pd.read_excel(guides_xls, header=8).rename(
    columns={"#guideId": "guideId", "Doench '16-Score": "doench"})
g["proto20"] = g["targetSeq"].astype(str).str[:20]
for c in ["mitSpecScore", "cfdSpecScore", "offtargetCount", "doench"]:
    g[c] = pd.to_numeric(g[c], errors="coerce")

o = pd.read_excel(offt_xls, header=8)
o["mismatchCount"] = pd.to_numeric(o["mismatchCount"], errors="coerce")
loc = [c for c in o.columns
       if o[c].astype(str).str.contains("intron|exon|intergenic", na=False).any()][0]

ours = {r["proto"]: r for r in csv.DictReader(open(COPE_GUIDES), delimiter="\t")
        if r["label"] == label}
if not ours:
    sys.exit(f"No guides with label {label} in {COPE_GUIDES}")

m = g[g["proto20"].isin(ours)].copy()
m = m[m["proto20"].map(lambda p: ours[p]["in_cope"] == "1")].copy()

def safety(gid):
    s = o[o["guideId"].astype(str) == str(gid)]
    close = s[s["mismatchCount"] <= 3]
    return len(close), int(close[loc].astype(str).str.startswith("exon").sum())

rows = []
for _, r in m.iterrows():
    le3, exle3 = safety(r["guideId"]); u = ours[r["proto20"]]
    rows.append(dict(chrom=u["chrom"], protospacer=r["proto20"], pam=u["pam"], strand=u["strand"],
                     cut=int(u["cut"]), dist_cope_center=int(u["dist_cope_center"]),
                     mit=int(r["mitSpecScore"]), cfd=int(r["cfdSpecScore"]), doench=int(r["doench"]),
                     ot_total=int(r["offtargetCount"]), ot_le3mm=le3, ot_exon_le3mm=exle3,
                     guideId=str(r["guideId"])))
sh = (pd.DataFrame(rows)
        .sort_values(["ot_exon_le3mm", "ot_le3mm", "mit", "doench"],
                     ascending=[True, True, False, False]))
sh.insert(0, "rank", range(1, len(sh) + 1))
sh.to_csv(out_tsv, sep="\t", index=False)
print(f"{label}: {len(sh)} in-COPE guides -> {out_tsv}")
print(sh.head(5).to_string(index=False))
