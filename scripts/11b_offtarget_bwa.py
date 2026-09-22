#!/usr/bin/env python3
"""
11b_offtarget_bwa.py

Local, autonomous genome-wide off-target screen for the sgRNA candidates enumerated by
11_design_sgrnas.py -- replaces the manual CRISPOR web round-trip (cas-offinder's OpenCL
backend is unreliable on this machine). Ported from the zmarker knock-in pipeline
(scripts/31_offtarget_bwa.py); the GRCz11 bwa index already lives in reference/.

For each candidate guide it finds every genomic site matching the 20 nt protospacer with
<=4 mismatches (bwa aln, CRISPOR-style params), keeps those with a valid NGG PAM, and
classifies them by mismatch tier and whether the off-target cut lands in an exon (the
dangerous class -- an off-target that disrupts a gene). Then it ranks the guides per
candidate SAFETY-FIRST, matching the criteria 16_curate_crispor_shortlist.py used on
CRISPOR output:

    sort by  ot_exon_le3mm ASC, ot_le3mm ASC, dist_center ASC

and writes the best guide per candidate (the autonomous shortlist) plus a full per-guide
off-target table. CRISPOR remains available (scripts 15/16) as an OPTIONAL final MIT/CFD/
Doench check on the chosen guides, but is no longer required to produce a shortlist.

Inputs : results/sgrna_candidates.tsv  (from 11_design_sgrnas.py)
Outputs: results/sgrna_offtarget.tsv       (every guide, scored)
         results/sgrna_shortlist_auto.tsv  (best guide per candidate)

Usage: python3 scripts/11b_offtarget_bwa.py [results/sgrna_candidates.tsv]
       (run in the `zebrafish_atac` conda env: needs bwa + samtools on PATH)
"""
import bisect
import csv
import os
import subprocess
import sys
import tempfile

REF = "reference/GRCz11.fa"
GTF = "annotation/Danio_rerio.GRCz11.110.gtf"
_args = [a for a in sys.argv[1:] if not a.startswith("--")]
IN = _args[0] if _args else "results/sgrna_candidates.tsv"
PER_CAND = 20                          # score only the K nearest-center guides per candidate (speed)
if "--per-candidate" in sys.argv:
    PER_CAND = int(sys.argv[sys.argv.index("--per-candidate") + 1])
OUT_ALL = "results/sgrna_offtarget.tsv"
OUT_SHORT = "results/sgrna_shortlist_auto.tsv"

_C = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(_C)[::-1]


def load_exons():
    """chrom -> sorted (start,end) exon intervals (Ensembl-named, matches reference/index)."""
    ex = {}
    for line in open(GTF):
        if "\texon\t" not in line:
            continue
        c = line.split("\t")
        ex.setdefault(c[0], []).append((int(c[3]), int(c[4])))
    for k in ex:
        ex[k].sort()
    return ex


def in_exon(ex, chrom, pos):
    ivs = ex.get(chrom)
    if not ivs:
        return False
    i = bisect.bisect_right(ivs, (pos, 1 << 62)) - 1
    for j in range(max(0, i - 3), min(len(ivs), i + 2)):
        if ivs[j][0] <= pos <= ivs[j][1]:
            return True
    return False


def faidx(chrom, a, b):
    r = subprocess.run(["samtools", "faidx", REF, f"{chrom}:{a}-{b}"],
                       capture_output=True, text=True)
    return "".join(r.stdout.split("\n")[1:]).upper()


def bwa_hits(proto):
    """All genomic hits of the 20 nt protospacer (<=4 mm) as (chrom, pos1, strand, nm)."""
    with tempfile.TemporaryDirectory() as d:
        fq, sai, sam = f"{d}/g.fq", f"{d}/g.sai", f"{d}/g.sam"
        open(fq, "w").write(f"@g\n{proto}\n+\n{'I' * len(proto)}\n")
        with open(sai, "wb") as fh:
            subprocess.run(["bwa", "aln", "-n", "4", "-o", "0", "-l", "20", "-k", "4",
                            "-N", REF, fq], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        with open(sam, "w") as fh:
            subprocess.run(["bwa", "samse", "-n", "100000", REF, sai, fq],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        hits = []
        for line in open(sam):
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "*":
                continue
            flag = int(f[1])
            hits.append((f[2], int(f[3]), "-" if flag & 16 else "+",
                         next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")), 0)))
            for tag in f[11:]:
                if tag.startswith("XA:Z:"):
                    for alt in tag[5:].split(";"):
                        if not alt:
                            continue
                        ch, pos, cig, nm = alt.split(",")
                        st = "-" if pos.startswith("-") else "+"
                        hits.append((ch, abs(int(pos)), st, int(nm)))
        return hits


def score_guide(proto, ex):
    """Off-target counts over valid-NGG-PAM sites, by mismatch tier + exon overlap."""
    ot_le3, ot_exon, tot = 0, 0, 0
    for chrom, pos, strand, nm in bwa_hits(proto):
        if nm == 0:
            continue                                    # on-target
        if strand == "+":
            pam = faidx(chrom, pos + 20, pos + 22)
            cut = pos + 17
        else:
            pam = rc(faidx(chrom, pos - 3, pos - 1))
            cut = pos + 3
        if len(pam) < 3 or pam[1:3] != "GG":
            continue                                    # no NGG PAM -> not a Cas9 off-target
        tot += 1
        if nm <= 3:
            ot_le3 += 1
            if in_exon(ex, chrom, cut):
                ot_exon += 1
    return dict(ot_total=tot, ot_le3mm=ot_le3, ot_exon_le3mm=ot_exon)


def flag_of(s):
    if s["ot_exon_le3mm"] == 0 and s["ot_le3mm"] <= 2:
        return "OK"
    return "CHECK" if s["ot_exon_le3mm"] == 0 else "RISK: exonic OT"


def main():
    if not os.path.exists(IN):
        sys.exit(f"{IN} not found -- run 11_design_sgrnas.py first")
    guides = list(csv.DictReader(open(IN), delimiter="\t"))
    # keep only the K nearest-center guides per candidate (script 11 already ranks by center)
    bycand = {}
    for g in guides:
        bycand.setdefault(g["cand"], []).append(g)
    guides = []
    for cand, gs in bycand.items():
        gs.sort(key=lambda g: int(g["dist_center"]))
        guides.extend(gs[:PER_CAND])
    ex = load_exons()
    print(f"scoring {len(guides)} guides ({len(bycand)} candidates x <= {PER_CAND})\n")
    print(f"{'cand':22} {'strand':>6} {'protospacer':22} {'dist':>5} "
          f"{'OT<=3':>6} {'exonOT':>7} {'OTtot':>6}  flag")
    scored = []
    for g in guides:
        s = score_guide(g["proto"], ex)
        row = dict(g, **s, flag=flag_of(s))
        scored.append(row)
        print(f"{g['cand']:22} {g['strand']:>6} {g['proto']:22} {g['dist_center']:>5} "
              f"{s['ot_le3mm']:>6} {s['ot_exon_le3mm']:>7} {s['ot_total']:>6}  {row['flag']}")

    # rank per candidate: safety first (exonic OT, then close OT), then proximity to center
    for r in scored:
        r["_key"] = (r["ot_exon_le3mm"], r["ot_le3mm"], int(r["dist_center"]))
    scored.sort(key=lambda r: (r["cand"], r["_key"]))

    cols = ["rank", "cand", "gene", "strand", "proto", "pam", "gc", "cut", "dist_center",
            "ot_total", "ot_le3mm", "ot_exon_le3mm", "flag"]
    with open(OUT_ALL, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in scored:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # best guide per candidate = the autonomous shortlist
    best, seen = [], set()
    for r in scored:
        if r["cand"] not in seen:
            best.append(r); seen.add(r["cand"])
    with open(OUT_SHORT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(best, key=lambda r: int(r["rank"])):
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\n[out] {OUT_ALL} ({len(scored)} guides)")
    print(f"[out] {OUT_SHORT} ({len(best)} candidates -> best guide each)")


if __name__ == "__main__":
    main()
