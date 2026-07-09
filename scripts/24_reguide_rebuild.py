#!/usr/bin/env python3
"""
24_reguide_rebuild.py

Re-guide cope_chr1, cope_chr14, pIGLET24b_v2 to higher on-target (IDT/RuleSet3)
guides and rebuild the full reagent set for each: HA arms, donor ultramers
(HA + DM416/DM417), Tm-matched genotyping primers, junction primers, and
WT + INTEGRATED GenBank maps. MIT/offtarget pulled from the CRISPOR guide tables.

New guides (chosen for on-target while keeping specificity + arm GC + cut in the
accessible/inert region):
  cope_chr1     ATGAGAGGTGACAAGGATGG  (RS3 48, MIT 91, 3 close OT/0 exonic)
  cope_chr14    GTCATAGACTGAGCGACTGG  (RS3 96, MIT 97, 1 close OT/1 exonic)
  pIGLET24b_v2  TTGTCTGTCCAATAGCGCGT  (RS3 83, MIT 98, 0 close OT)
"""
import csv, re, subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

REF = "reference/GRCz11.fa"
CDIR = "/Users/dsm16/Documents/Betzig/Swinburne-Gokul-Betzig constructs/fish CRISPR/landing_pads"
CASSETTE = f"{CDIR}/spacer-attP-GT-unc45b-lck-mSG-SV40-attP-GA-spacer.gb"
OUTDIR = f"{CDIR}/integrated_alleles"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
DM416 = "CTGCAACTTGTTTGAGGAATTTGAGG"; DM417 = "AAACTTAGAGGCATCATAGTGGAGATCA"
SG_R = "CGGTGAAGGACTTGCCATTGAT"; SG_F = "CGTGATCACATCATCCAGTCGG"
OFF5, OFF3 = 653, 358
def rc(s): return s.translate(COMP)[::-1]
def gc(s): return (s.count("G")+s.count("C"))/len(s)
def tm(s): return mt.Tm_NN(s)
def faidx(ens, a, b):
    r = subprocess.run(["samtools", "faidx", REF, f"{ens}:{a+1}-{b}"], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()

# locus: (ens, chrom, window0based(a,b), guide, crispor_guides_xls, crispor_offtargets_xls)
DL="/Users/dsm16/Downloads"
LOCI = [
 ("cope_chr1","1","chr1",(15549794,15551236),"GTGCTGGTCACAACCATGAG",
  f"{DL}/guides_danRer11-chr1-15549795-15551236.xls",f"{DL}/offtargets_danRer11-chr1-15549795-15551236.xls"),
 ("cope_chr14","14","chr14",(43942076,43942572),"GTCATAGACTGAGCGACTGG",
  f"{DL}/guides_danRer11-chr14-43942077-43942572.xls",f"{DL}/offtargets_danRer11-chr14-43942077-43942572.xls"),
 ("pIGLET24b_v2","24","chr24",(26386480,26387050),"TTGTCTGTCCAATAGCGCGT",
  f"{DL}/guides_pIGLET24b_danRer11-chr24-26386249-26387249.xls",f"{DL}/offtargets_pIGLET24b_danRer11-chr24-26386249-26387249.xls"),
]

def find_cut(ens, a, b, guide):
    seq = faidx(ens, a, b)
    i = seq.find(guide)
    if i >= 0: return a+i+17, "+", seq[i+20:i+23]      # 0-based break, strand, PAM
    j = seq.find(rc(guide))
    if j >= 0: return a+j+3, "-", rc(seq[j-3:j])
    raise SystemExit(f"guide {guide} not found in {ens}:{a}-{b}")

def crispor_metrics(xls, off_xls, guide):
    g = pd.read_excel(xls, header=8).rename(columns={"#guideId":"guideId"})
    g["p20"] = g["targetSeq"].astype(str).str[:20]
    row = g[g["p20"] == guide].iloc[0]
    return int(pd.to_numeric(row["mitSpecScore"])), int(pd.to_numeric(row["offtargetCount"]))

# genotyping designer (Tm-matched, longer OK for AT-rich flanks)
LENGTHS = tuple(range(20, 33))
def good(p): return 0.38 <= gc(p) <= 0.62 and p[-1] in "GC" and not any(c*6 in p for c in "ACGT")
def pick(cands, target):
    v = [(abs(tm(p)-target), tm(p), p, pos) for p, pos in cands if good(p)]
    return None if not v else min(v)[1:]
def design_pair(ens, cut0, target=62.0):
    a, b = cut0-300, cut0+980; win = faidx(ens, a, b)
    fwd = [(win[s-a:s-a+L], s) for s in range(cut0-280, cut0-130) for L in LENGTHS if s+L <= cut0-5]
    rev = [(rc(win[e-L-a:e-a]), e) for e in range(cut0+380, cut0+900) for L in LENGTHS if e-L >= cut0+5 and e <= b]
    F = pick(fwd, target); R = pick(rev, F[0]) if F else None
    return F, R

def parse_cassette(path):
    text = open(path).read()
    seq = re.sub(r"[^acgtACGT]", "", text.split("ORIGIN")[1].split("//")[0]).upper()
    feats, cur = [], None
    for ln in text.split("FEATURES", 1)[1].split("ORIGIN")[0].splitlines()[1:]:
        m = re.match(r"^ {5}(\S+)\s+(complement\()?(\d+)\.\.(\d+)\)?\s*$", ln)
        if m:
            if cur: feats.append(cur)
            cur = dict(type=m.group(1), strand=-1 if m.group(2) else 1, start=int(m.group(3))-1, end=int(m.group(4)), label=None)
        elif cur is not None:
            q = re.search(r'/(?:label|note)="?([^"\n]+)', ln)
            if q and cur["label"] is None: cur["label"] = q.group(1).strip().strip('"')
    if cur: feats.append(cur)
    return seq, feats
def feat(s, e, label, typ="misc_feature", strand=1):
    f = SeqFeature(FeatureLocation(int(s), int(e), strand=strand), type=typ); f.qualifiers["label"] = [label]; return f
def find(win, sub):
    i = win.find(sub)
    if i >= 0: return (i, i+len(sub), 1)
    i = win.find(rc(sub))
    return (i, i+len(sub), -1) if i >= 0 else None

cas_seq, cas_feats = parse_cassette(CASSETTE)
donor_rows, geno_rows, junc_rows = {}, {}, {}
DH = "site role chrom guide_protospacer_20nt pam strand cut_genomic_GRCz11 mit_spec offtarget_count HA5_45nt HA3_45nt fwd_primer_tail rev_primer_tail notes fwd_primer_full rev_primer_full".split()
GH = "locus chrom insertion_site genomic_F F_Tm F_pos genomic_R R_Tm R_pos WT_amplicon_bp".split()
JH = "locus chrom genomic_F StayGold_R_universal amplicon_5pjunction_bp StayGold_F_universal genomic_R amplicon_3pjunction_bp".split()

for site, ens, chrom, (a, b), guide, gx, ox in LOCI:
    cut0, strand, pam = find_cut(ens, a, b, guide)
    mit, ot = crispor_metrics(gx, ox, guide)
    ha5 = faidx(ens, cut0-45, cut0); ha3 = faidx(ens, cut0, cut0+45)
    fwd = ha5 + DM416; rev = rc(ha3) + DM417
    (tF, Fseq, Fs), (tR, Rseq, Re) = design_pair(ens, cut0)
    ins = cut0+1; Fp = Fs+1; Rp = Re
    wt = Rp-Fp+1; j5 = (ins-Fp)+OFF5; j3 = (Rp-ins)+OFF3
    donor_rows[site] = dict(site=site, role="re-guide (high on-target)", chrom=chrom, guide_protospacer_20nt=guide,
        pam=pam, strand=strand, cut_genomic_GRCz11=f"{chrom}:{ins}", mit_spec=str(mit), offtarget_count=str(ot),
        HA5_45nt=ha5, HA3_45nt=ha3, fwd_primer_tail=ha5, rev_primer_tail=rc(ha3),
        notes=f"re-guided for on-target (RS3); arms {gc(ha5)*100:.0f}/{gc(ha3)*100:.0f}%",
        fwd_primer_full=fwd, rev_primer_full=rev)
    geno_rows[site] = dict(locus=site, chrom=chrom, insertion_site=str(ins), genomic_F=Fseq, F_Tm=f"{tF:.1f}",
        F_pos=f"{chrom}:{Fp}", genomic_R=Rseq, R_Tm=f"{tR:.1f}", R_pos=f"{chrom}:{Rp}", WT_amplicon_bp=str(wt))
    junc_rows[site] = dict(locus=site, chrom=chrom, genomic_F=Fseq, StayGold_R_universal=SG_R,
        amplicon_5pjunction_bp=str(j5), StayGold_F_universal=SG_F, genomic_R=Rseq, amplicon_3pjunction_bp=str(j3))
    print(f"{site}: {guide} ({strand}) cut {chrom}:{ins} | arms {gc(ha5)*100:.0f}/{gc(ha3)*100:.0f}% | MIT {mit} | "
          f"geno Tm {tF:.0f}/{tR:.0f} WT {wt} | junc {j5}/{j3}")

    # GenBank maps
    w0, w1 = min(Fp-1, cut0)-20, max(Rp, cut0)+20
    win = faidx(ens, w0, w1); cutw = cut0-w0
    rec = SeqRecord(Seq(win), id=site, name=f"WT_{site}", description=f"WT {chrom} {site} re-guide (GRCz11)")
    rec.annotations["molecule_type"] = "DNA"
    fs = [feat(0, len(win), f"genomic {chrom} {site}", "source")]
    for sq, lab in [(Fseq, f"genotyping_F (Tm {tF:.1f})"), (Rseq, f"genotyping_R (Tm {tR:.1f})")]:
        h = find(win, sq); fs.append(feat(h[0], h[1], lab, "primer_bind", h[2])) if h else None
    fs += [feat(cutw-45, cutw, "5' homology arm (45 bp)"), feat(cutw, cutw+45, "3' homology arm (45 bp)")]
    h = find(win, guide)
    if h:
        fs.append(feat(h[0], h[1], f"sgRNA spacer: {guide}", "misc_feature", h[2]))
        fs.append(feat(h[1], h[1]+3, f"PAM: {pam}", "misc_feature", 1) if h[2] == 1 else feat(h[0]-3, h[0], f"PAM: {pam}", "misc_feature", -1))
    fs.append(feat(cutw-1, cutw+1, "Cas9 blunt cut = insertion site"))
    rec.features = sorted(fs, key=lambda f: int(f.location.start))
    SeqIO.write(rec, f"{OUTDIR}/WT_{site}.gb", "genbank")
    left, right = win[:cutw], win[cutw:]; integ = left+cas_seq+right; off = len(left); roff = off+len(cas_seq)
    ir = SeqRecord(Seq(integ), id=site, name=f"INT_{site}", description=f"INTEGRATED attP cassette at {chrom} {site}")
    ir.annotations["molecule_type"] = "DNA"
    ifs = [feat(0, len(integ), f"INTEGRATED allele {chrom} {site} (attP cassette)", "source"),
           feat(off-45, off, "5' homology arm (45 bp)"), feat(roff, roff+45, "3' homology arm (45 bp)"),
           feat(off-1, off+1, "5' junction (genomic|cassette)"), feat(roff-1, roff+1, "3' junction (cassette|genomic)")]
    for cf in cas_feats: ifs.append(feat(cf["start"]+off, cf["end"]+off, cf["label"] or cf["type"], cf["type"], cf["strand"]))
    hF = find(left, Fseq); ifs.append(feat(hF[0], hF[1], f"genotyping_F (Tm {tF:.1f})", "primer_bind", hF[2])) if hF else None
    hR = find(right, Rseq); ifs.append(feat(hR[0]+roff, hR[1]+roff, f"genotyping_R (Tm {tR:.1f})", "primer_bind", hR[2])) if hR else None
    for sq, lab in [(SG_R, "StayGold_R (5' junction)"), (SG_F, "StayGold_F (3' junction)")]:
        h = find(cas_seq, sq); ifs.append(feat(h[0]+off, h[1]+off, lab, "primer_bind", h[2])) if h else None
    ir.features = sorted(ifs, key=lambda f: int(f.location.start))
    SeqIO.write(ir, f"{OUTDIR}/INTEGRATED_{site}_attP.gb", "genbank")

# upsert into the deliverable TSVs
def upsert(path, hdr, rows):
    keys = set(rows)
    lines = open(path).read().splitlines()
    keep = [ln for ln in lines if ln.split("\t")[0] not in keys]
    with open(path, "w") as f:
        f.write("\n".join(keep) + "\n" + "\n".join("\t".join(rows[k][h] for h in hdr) for k in rows) + "\n")
upsert("results/knockin_donor_design.tsv", DH, donor_rows)
upsert("results/genotyping_primers.tsv", GH, geno_rows)
upsert("results/junction_primers.tsv", JH, junc_rows)
print("\nupdated donor/genotyping/junction TSVs + GenBank maps for all 3 re-guided loci")
