#!/usr/bin/env python3
"""
23_oligo_master.py
Assemble the single source-of-truth oligo sheet (DM416-DM447) for the knock-in
project: universal anneal + donor ultramers + genotyping + StayGold junction
primers, with ordered/to-order status and an annotation column in the lab's
style (F primer carries WT + Tm + 5' junction size; R primer carries 3' junction).

  DM416/417  universal cassette-anneal (Phusion, plasmid)
  DM418-425  donor ultramers  chr12 / chr20 / chr24 / pIGLET24b(original)   [ORDERED]
  DM426-433  genotyping F/R   chr12 / chr20 / chr24 / pIGLET24b(original)   [ORDERED]
  DM434/435  StayGold junction F / R (universal)                           [ORDERED*]
  DM436-441  donor ultramers  cope_chr1 / cope_chr14 / pIGLET24b_v2         [TO ORDER]
  DM442-447  genotyping F/R   cope_chr1 / cope_chr14 / pIGLET24b_v2         [TO ORDER]

* DM435's ordered sequence is a paste error; correct StayGold-R is set here.
Output: results/oligo_master.tsv
"""
import csv
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
DM416 = "CTGCAACTTGTTTGAGGAATTTGAGG"; DM417 = "AAACTTAGAGGCATCATAGTGGAGATCA"
SG_F = "CGTGATCACATCATCCAGTCGG"; SG_R = "CGGTGAAGGACTTGCCATTGAT"
def mod(seq): return "/5AmMC6/" + f"{seq[0]}*{seq[1]}*{seq[2]}*{seq[3]}" + seq[4:]

d = {r["site"]: r for r in csv.DictReader(open("results/knockin_donor_design.tsv"), delimiter="\t")}
g = {r["locus"]: r for r in csv.DictReader(open("results/genotyping_primers.tsv"), delimiter="\t")}
j = {r["locus"]: r for r in csv.DictReader(open("results/junction_primers.tsv"), delimiter="\t")}

rows = [
    ("DM416", "universal", "cassette-anneal F", "ORDERED", DM416, "plasmid-priming; 3' of every donor_F"),
    ("DM417", "universal", "cassette-anneal R", "ORDERED", DM417, "plasmid-priming; 3' of every donor_R"),
]
DONORS = [("DM418","DM419","chr12_nrxn1a","ORDERED"),("DM420","DM421","chr20_ENSDARG96520","ORDERED"),
          ("DM422","DM423","chr24_inhbaa","ORDERED"),("DM424","DM425","pIGLET24b_locus","ORDERED"),
          ("DM436","DM437","cope_chr1","TO ORDER"),("DM438","DM439","cope_chr14","TO ORDER"),
          ("DM440","DM441","pIGLET24b_v2","TO ORDER")]
GENO = [("DM426","DM427","chr12_nrxn1a","ORDERED"),("DM428","DM429","chr20_ENSDARG96520","ORDERED"),
        ("DM430","DM431","chr24_inhbaa","ORDERED"),("DM432","DM433","pIGLET24b_locus","ORDERED"),
        ("DM442","DM443","cope_chr1","TO ORDER"),("DM444","DM445","cope_chr14","TO ORDER"),
        ("DM446","DM447","pIGLET24b_v2","TO ORDER")]

def sup(s): return " — SUPERSEDED (71% GC), use DM440/441" if s == "pIGLET24b_locus" else (" — replaces DM424/425" if s == "pIGLET24b_v2" else "")
for nf, nr, s, st in DONORS:
    F = mod(d[s]["HA5_45nt"] + DM416); R = mod(rc(d[s]["HA3_45nt"]) + DM417)
    rows.append((nf, s, "donor_F", st, F, "HDR donor ultramer HA5+DM416; /5AmMC6/+3'PS" + sup(s)))
    rows.append((nr, s, "donor_R", st, R, "HDR donor ultramer revcomp(HA3)+DM417; /5AmMC6/+3'PS" + sup(s)))
# StayGold junction primers sit between the donor blocks in the DM series
rows.append(("DM434", "universal", "StayGold junction F", "ORDERED", SG_F, "3' junction (with any geno_R)"))
rows.append(("DM435", "universal", "StayGold junction R", "ORDERED*", SG_R, "5' junction (with any geno_F); *ordered seq was a paste error"))
for nf, nr, s, st in GENO:
    gg, jj = g[s], j[s]
    annF = f"no insert {gg['WT_amplicon_bp']}, Tm(NN) {gg['F_Tm']}, +DM435 (5' junc) = {jj['amplicon_5pjunction_bp']}bp"
    rows.append((nf, s, "genotyping F", st, gg["genomic_F"], annF))
    rows.append((nr, s, "genotyping R", st, gg["genomic_R"], f"+DM434 (3' junc) = {jj['amplicon_3pjunction_bp']}bp"))

# emit in DM-number order
rows.sort(key=lambda r: int(r[0].replace("DM", "").rstrip("*")))
with open("results/oligo_master.tsv", "w") as f:
    f.write("name\tlocus\trole\tstatus\tsequence\tannotation\n")
    for r in rows: f.write("\t".join(r) + "\n")
print(f"wrote results/oligo_master.tsv ({len(rows)} oligos, DM416-DM447)")
print(f"  ORDERED {sum(1 for r in rows if r[3].startswith('ORDERED'))} | TO ORDER {sum(1 for r in rows if r[3]=='TO ORDER')}")
