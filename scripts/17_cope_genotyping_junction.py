#!/usr/bin/env python3
"""
17_cope_genotyping_junction.py

For the two COPE lab targets (chr1, chr14):
  1) update the donor full-primers in knockin_donor_design.tsv to the longer
     plasmid-anneal tails actually ordered (DM416 fwd / DM417 rev);
  2) design Tm-matched genomic genotyping primers flanking the cut;
  3) build junction primers (genomic F/R + universal StayGold primers), using
     the constant cassette offsets measured from the existing loci
     (5' = 653 bp, 3' = 358 bp into the cassette).

Tm is reported by Biopython Tm_NN for INTERNAL matching only; the user finalises
annealing temp in NEB's Tm Calculator (Phusion 2X), which reads ~8-10C higher.

Appends rows to results/{genotyping_primers,junction_primers}.tsv and rewrites
the COPE donor full-primers in results/knockin_donor_design.tsv.
"""
import csv, subprocess
from Bio.SeqUtils import MeltingTemp as mt

REF = "reference/GRCz11.fa"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
DM416 = "CTGCAACTTGTTTGAGGAATTTGAGG"          # landing-pad F anneal (ordered)
DM417 = "AAACTTAGAGGCATCATAGTGGAGATCA"        # landing-pad R anneal (ordered)
SG_R_UNIV = "CGGTGAAGGACTTGCCATTGAT"          # StayGold R (5' junction)
SG_F_UNIV = "CGTGATCACATCATCCAGTCGG"          # StayGold F (3' junction)
OFF5, OFF3 = 653, 358                          # cassette offsets (constant)

def rc(s): return s.translate(COMP)[::-1]
def faidx(ens, a, b):  # 0-based half-open [a,b)
    r = subprocess.run(["samtools", "faidx", REF, f"{ens}:{a+1}-{b}"], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()
def gc(s): return (s.count("G") + s.count("C")) / len(s)
def tm(s): return mt.Tm_NN(s)

LENGTHS = (20, 21, 22, 23, 24, 25, 26)   # allow extra length for AT-rich flanks
                                         # (window stays narrow, so amplicon stays modest)
def good(p):
    return (20 <= len(p) <= 26 and 0.40 <= gc(p) <= 0.60
            and p[-1] in "GC" and not any(b*5 in p for b in "ACGT"))

def pick(cands, target):
    """cands: list of (primer, anchor_pos); valid one with Tm_NN closest to target.
    Returns (tm, primer, pos) or None."""
    v = [(abs(tm(p)-target), tm(p), p, pos) for p, pos in cands if good(p)]
    if not v:
        return None
    _, t, p, pos = min(v)
    return (t, p, pos)

def design_pair(ens, cut0, target=57.0):
    lo, hi = cut0-300, cut0+680
    win = faidx(ens, lo, hi); base = lo
    # forward: F start ~140-280bp upstream; 3' end upstream of cut
    fwd = [(win[s-base:s-base+L], s)
           for s in range(cut0-280, cut0-130)
           for L in LENGTHS if s+L <= cut0-5]
    # reverse: R 5' end ~380-620bp downstream (keeps WT amplicon ~550-850bp)
    rev = [(rc(win[e-L-base:e-base]), e)
           for e in range(cut0+380, cut0+620)
           for L in LENGTHS if e-L >= cut0+5 and e <= hi]
    F = pick(fwd, target)
    R = pick(rev, F[0]) if F else None   # match R to F's actual Tm
    return F, R

ENS = {"chr1": "1", "chr14": "14"}
# (locus, chrom, cut0[0-based break]) ; recorded 1-based insertion = cut0+1
TARGETS = [("cope_chr1", "chr1", 15550248), ("cope_chr14", "chr14", 43942400)]

geno_rows, junc_rows, donor_primers = [], [], {}
for site, chrom, cut0 in TARGETS:
    ens = ENS[chrom]; ins = cut0 + 1
    (_, Fp, Fs), (_, Rp, Re) = design_pair(ens, cut0)
    Fpos1 = Fs + 1                     # 1-based F start
    Rpos1 = Re                         # 1-based R 5' start (rightmost)
    wt = Rpos1 - Fpos1 + 1
    j5 = (ins - Fpos1) + OFF5
    j3 = (Rpos1 - ins) + OFF3
    geno_rows.append([site, chrom, str(ins), Fp, f"{tm(Fp):.1f}", f"{chrom}:{Fpos1}",
                      Rp, f"{tm(Rp):.1f}", f"{chrom}:{Rpos1}", str(wt)])
    junc_rows.append([site, chrom, Fp, SG_R_UNIV, str(j5), SG_F_UNIV, Rp, str(j3)])
    # donor full primers: HA5 + DM416 ; revcomp(HA3) + DM417
    ha5 = faidx(ens, cut0-45, cut0); ha3 = faidx(ens, cut0, cut0+45)
    donor_primers[site] = (ha5 + DM416, rc(ha3) + DM417)
    print(f"{site}: F {Fp} (Tm_NN {tm(Fp):.1f}) / R {Rp} (Tm_NN {tm(Rp):.1f}) | "
          f"WT {wt}bp | 5'junc {j5} | 3'junc {j3}")

def upsert(path, rows):
    """Drop any existing cope_* rows, then append fresh ones (idempotent)."""
    labels = {r[0] for r in rows}
    keep = [ln for ln in open(path).read().splitlines()
            if ln.split("\t")[0] not in labels]
    with open(path, "w") as f:
        f.write("\n".join(keep) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
upsert("results/genotyping_primers.tsv", geno_rows)
upsert("results/junction_primers.tsv", junc_rows)

# rewrite COPE donor full-primers in place
lines = open("results/knockin_donor_design.tsv").read().splitlines()
hdr = lines[0].split("\t"); iF, iR = hdr.index("fwd_primer_full"), hdr.index("rev_primer_full")
out = [lines[0]]
for ln in lines[1:]:
    c = ln.split("\t")
    if c[0] in donor_primers:
        c[iF], c[iR] = donor_primers[c[0]]
    out.append("\t".join(c))
open("results/knockin_donor_design.tsv", "w").write("\n".join(out) + "\n")
print("\nUpdated donor full-primers (DM416/DM417); appended genotyping + junction rows.")
