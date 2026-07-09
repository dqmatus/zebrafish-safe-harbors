#!/usr/bin/env python3
"""
22_oligo_order_summary.py
Combine donor + genotyping + junction designs into one ordering-ready summary
across all active knock-in loci. Junction PCR reuses the genotyping genomic F/R
plus the two universal StayGold primers, so the only per-locus oligos to order
are: donor_F/R ultramers + genotyping_F/R.

Outputs:
  results/oligo_order_summary.tsv   (one row per oligo — paste into a bulk order)
  results/locus_design_summary.tsv  (one row per locus — guide/cut/arms/amplicons)
"""
import csv
gc = lambda s: 100*(s.count("G")+s.count("C"))/len(s) if s else 0

donor = {r["site"]: r for r in csv.DictReader(open("results/knockin_donor_design.tsv"), delimiter="\t")}
geno  = {r["locus"]: r for r in csv.DictReader(open("results/genotyping_primers.tsv"), delimiter="\t")}
junc  = {r["locus"]: r for r in csv.DictReader(open("results/junction_primers.tsv"), delimiter="\t")}

# active loci (original pIGLET24b_locus superseded by pIGLET24b_v2)
ORDER = ["chr12_nrxn1a", "chr20_ENSDARG96520", "chr24_inhbaa", "cope_chr1", "cope_chr14", "pIGLET24b_v2"]
SG_R = next(iter(junc.values()))["StayGold_R_universal"]
SG_F = next(iter(junc.values()))["StayGold_F_universal"]

# ---- per-oligo order sheet ----
oligos = []
for L in ORDER:
    d = donor.get(L); g = geno.get(L)
    if d:
        oligos.append((f"{L}_donor_F", d["fwd_primer_full"], "HDR donor ultramer (fwd): HA5 + DM416"))
        oligos.append((f"{L}_donor_R", d["rev_primer_full"], "HDR donor ultramer (rev): revcomp(HA3) + DM417"))
    if g:
        oligos.append((f"{L}_geno_F", g["genomic_F"], "genotyping / junction genomic F"))
        oligos.append((f"{L}_geno_R", g["genomic_R"], "genotyping / junction genomic R"))
# universal oligos (order once)
oligos += [
    ("DM416_landingpad_F", "CTGCAACTTGTTTGAGGAATTTGAGG", "universal cassette-anneal (fwd) — 3' of every donor_F"),
    ("DM417_landingpad_R", "AAACTTAGAGGCATCATAGTGGAGATCA", "universal cassette-anneal (rev) — 3' of every donor_R"),
    ("StayGold_junction_F", SG_F, "universal 3' junction primer (with geno_R)"),
    ("StayGold_junction_R", SG_R, "universal 5' junction primer (with geno_F)"),
]
with open("results/oligo_order_summary.tsv", "w") as f:
    f.write("oligo_name\tsequence\tlength\tnotes\n")
    for n, s, note in oligos:
        f.write(f"{n}\t{s}\t{len(s)}\t{note}\n")

# ---- per-locus design summary ----
with open("results/locus_design_summary.tsv", "w") as f:
    cols = ["locus", "chrom", "cut_GRCz11", "guide_20nt", "pam", "strand", "mit_spec",
            "HA5_GC%", "HA3_GC%", "donor_len_F", "donor_len_R",
            "geno_F_Tm", "geno_R_Tm", "WT_amplicon_bp", "5pjunc_bp", "3pjunc_bp"]
    f.write("\t".join(cols) + "\n")
    for L in ORDER:
        d = donor.get(L, {}); g = geno.get(L, {}); j = junc.get(L, {})
        row = [L, d.get("chrom", ""), d.get("cut_genomic_GRCz11", ""), d.get("guide_protospacer_20nt", ""),
               d.get("pam", ""), d.get("strand", ""), d.get("mit_spec", ""),
               f"{gc(d.get('HA5_45nt','')):.0f}", f"{gc(d.get('HA3_45nt','')):.0f}",
               str(len(d.get("fwd_primer_full", ""))), str(len(d.get("rev_primer_full", ""))),
               g.get("F_Tm", ""), g.get("R_Tm", ""), g.get("WT_amplicon_bp", ""),
               j.get("amplicon_5pjunction_bp", ""), j.get("amplicon_3pjunction_bp", "")]
        f.write("\t".join(row) + "\n")

print(f"wrote results/oligo_order_summary.tsv  ({len(oligos)} oligos: "
      f"{len(ORDER)*4} locus-specific + 4 universal)")
print("wrote results/locus_design_summary.tsv")
