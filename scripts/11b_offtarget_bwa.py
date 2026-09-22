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
dangerous class -- an off-target that disrupts a gene). It also computes on-target
efficiency locally with CRISPRscan (Moreno-Mateos 2015), the zebrafish-injection-trained
model (more appropriate here than Doench Rule Set 2, which is mammalian-cell-trained).
Then it ranks the guides per candidate SAFETY-FIRST, then by efficiency:

    sort by  ot_exon_le3mm ASC, ot_le3mm ASC, CRISPRscan DESC, dist_center ASC

and writes the best guide per candidate (the autonomous shortlist) plus a full per-guide
table. CRISPOR (MIT/CFD + Doench) remains available (scripts 15/16) as an OPTIONAL final
check, but is no longer required to produce a shortlist.

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


def ucsc_to_ensembl(chrom):
    """chr12 -> 12 ; chrM -> MT ; chrUn_KN150617v1 -> KN150617.1 (match reference/GRCz11.fa)."""
    if chrom.startswith("chrUn_"):
        tail = chrom[len("chrUn_"):]
        if "v" in tail:
            acc, ver = tail.rsplit("v", 1)
            return f"{acc}.{ver}"
        return tail
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if c in ("M", "MT") else c


# ---- CRISPRscan on-target efficiency (Moreno-Mateos 2015; zebrafish-trained) ----
# Coefficients verbatim from CRISPOR (maximilianh/crisporWebsite: crisporEffScores.py).
# Input = 35 nt context: 6 nt 5' + 20 nt protospacer + 3 nt PAM + 6 nt 3' (1-based positions).
# Score = int(100 * (intercept + sum of matched position-specific mono/di-nucleotide weights)).
CRISPRSCAN_INTERCEPT = 0.183930943629
CRISPRSCAN_PARAMS = [
    ('AA',18,-0.097377097),('TT',18,-0.094424075),('TT',13,-0.08618771),('CT',26,-0.084264893),
    ('GC',25,-0.073453609),('T',21,-0.068730497),('TG',23,-0.066388075),('AG',23,-0.054338456),
    ('G',30,-0.046315914),('A',4,-0.042153521),('AG',34,-0.041935908),('GA',34,-0.037797707),
    ('A',18,-0.033820432),('C',25,-0.031648353),('C',31,-0.030715556),('G',1,-0.029693709),
    ('C',16,-0.021638609),('A',14,-0.018487229),('A',11,-0.018287292),('T',34,-0.017647692),
    ('AA',10,-0.016905415),('A',19,-0.015576499),('G',34,-0.014167123),('C',30,-0.013182733),
    ('GA',31,-0.01227989),('T',24,-0.011996172),('A',15,-0.010595296),('G',4,-0.005448869),
    ('GG',9,-0.00157799),('T',23,-0.001422243),('C',15,-0.000477727),('C',26,-0.000368973),
    ('T',27,-0.000280845),('A',31,0.00158975),('GT',18,0.002391744),('C',9,0.002449224),
    ('GA',20,0.009740799),('A',25,0.010506405),('A',12,0.011633235),('A',32,0.012435231),
    ('T',22,0.013224035),('C',20,0.015089514),('G',17,0.01549378),('G',18,0.016457816),
    ('T',30,0.017263162),('A',13,0.017628924),('G',19,0.017916844),('A',27,0.019126815),
    ('G',11,0.020929039),('TG',3,0.022949996),('GC',3,0.024681785),('G',14,0.025116714),
    ('GG',10,0.026802158),('G',12,0.027591138),('G',32,0.03071249),('A',22,0.031930909),
    ('G',20,0.033957008),('C',21,0.034262921),('TT',17,0.03492881),('T',13,0.035445171),
    ('G',26,0.036146649),('A',24,0.037466478),('C',22,0.03763162),('G',16,0.037970942),
    ('GG',12,0.041883009),('TG',18,0.045908991),('TG',31,0.048136812),('A',35,0.048596259),
    ('G',15,0.051129717),('C',24,0.052972314),('TG',15,0.053372822),('GT',11,0.053678436),
    ('GC',9,0.054171402),('CA',30,0.057759851),('GT',24,0.060952114),('G',13,0.061360905),
    ('CA',24,0.06221937),('AG',10,0.063717093),('G',10,0.067739182),('C',13,0.069495944),
    ('GT',31,0.07342535),('GG',13,0.074355848),('C',27,0.079933922),('G',27,0.085151052),
    ('CC',21,0.088919601),('CC',23,0.095072286),('G',22,0.10114438),('G',24,0.105488325),
    ('GT',23,0.106718563),('GG',25,0.111559441),('G',9,0.114600681),
]


def calc_crisprscan(seq35):
    """35 nt guide-strand context -> integer efficiency; None if the sequence is invalid."""
    if not seq35 or len(seq35) != 35:      # N bases simply match no motif (as in CRISPOR)
        return None
    score = CRISPRSCAN_INTERCEPT
    for motif, pos, w in CRISPRSCAN_PARAMS:
        if seq35[pos - 1:pos - 1 + len(motif)] == motif:
            score += w
    return int(100 * score)


def guide_context35(cand, gstart, gend, strand):
    """Fetch the 35 nt guide-strand context (6 + 20 + 3 + 6) for a guide, or None at a chrom edge."""
    chrom = ucsc_to_ensembl(cand.split(":")[0])
    a = int(gstart) - 6                       # 0-based inclusive
    b = int(gend) + 6                         # 0-based exclusive
    if a < 0:
        return None
    w = faidx(chrom, a + 1, b)                # samtools faidx is 1-based inclusive
    if len(w) != 35:
        return None
    return w if strand == "+" else rc(w)


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
    assert calc_crisprscan("TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA") == 77, "CRISPRscan self-test failed"
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
          f"{'OT<=3':>6} {'exonOT':>7} {'OTtot':>6} {'eff':>4}  flag")
    scored = []
    for g in guides:
        s = score_guide(g["proto"], ex)
        eff = calc_crisprscan(guide_context35(g["cand"], g["gstart"], g["gend"], g["strand"]))
        row = dict(g, **s, crisprscan=("" if eff is None else eff), flag=flag_of(s))
        scored.append(row)
        print(f"{g['cand']:22} {g['strand']:>6} {g['proto']:22} {g['dist_center']:>5} "
              f"{s['ot_le3mm']:>6} {s['ot_exon_le3mm']:>7} {s['ot_total']:>6} "
              f"{str(row['crisprscan']):>4}  {row['flag']}")

    # rank per candidate: safety first (exonic OT, then close OT), then CRISPRscan efficiency, then center
    for r in scored:
        eff = int(r["crisprscan"]) if r["crisprscan"] != "" else -1
        r["_key"] = (r["ot_exon_le3mm"], r["ot_le3mm"], -eff, int(r["dist_center"]))
    scored.sort(key=lambda r: (r["cand"], r["_key"]))

    cols = ["rank", "cand", "gene", "strand", "proto", "pam", "gc", "cut", "dist_center",
            "ot_total", "ot_le3mm", "ot_exon_le3mm", "crisprscan", "flag"]
    with open(OUT_ALL, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in scored:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # best guide per candidate = the autonomous shortlist.
    # HARD GATE for safe harbors: never ship a guide with an exonic off-target (it could
    # disrupt a gene and confound the line). scored is sorted safety-first, so the first
    # row per candidate is exon-clean if any exists; if not, flag the site for review.
    best, seen = [], set()
    for r in scored:
        if r["cand"] in seen:
            continue
        seen.add(r["cand"])
        if r["ot_exon_le3mm"] > 0:
            r = dict(r, flag="REVIEW: no exon-clean guide")
        best.append(r)
    with open(OUT_SHORT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(best, key=lambda r: int(r["rank"])):
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\n[out] {OUT_ALL} ({len(scored)} guides)")
    print(f"[out] {OUT_SHORT} ({len(best)} candidates -> best guide each)")


if __name__ == "__main__":
    main()
