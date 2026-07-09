#!/usr/bin/env python3
"""
21_piglet24b_v2_build.py

Build the full knock-in reagent set for the pIGLET24b RE-GUIDE (pIGLET24b_v2),
chosen to move the cut out of the 71%-GC core into GC-balanced arms so the donor
ultramers actually amplify (the original pIGLET24b guide's HA5 was 71% GC).

Guide: GTTGTTGAGAAGCGCGATAA (+, PAM AGG), cut chr24:26386978
CRISPOR (danRer11): MIT 97 / CFD 98, offtargets 23 (0 exonic <=3mm), Doench 47
HA arms ~51% / 38% GC; still inside the accessible/inert peak (26386502-26387026).

Produces, matching the other loci:
  - donor row in results/knockin_donor_design.tsv (HA + DM416/DM417 ultramers)
  - genotyping + junction primers (results/{genotyping,junction}_primers.tsv)
  - WT + INTEGRATED GenBank maps in the constructs folder
"""
import csv, re, subprocess
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
DM416 = "CTGCAACTTGTTTGAGGAATTTGAGG"
DM417 = "AAACTTAGAGGCATCATAGTGGAGATCA"
SG_R = "CGGTGAAGGACTTGCCATTGAT"
SG_F = "CGTGATCACATCATCCAGTCGG"
OFF5, OFF3 = 653, 358

# --- target ---
SITE = "pIGLET24b_v2"
ENS, CHROM = "24", "chr24"
GUIDE, PAM, STRAND = "GTTGTTGAGAAGCGCGATAA", "AGG", "+"
CUT0 = 26386978 - 1          # 0-based break (HA3 start); recorded 1-based = CUT0+1
MIT, OT = 97, 23

def rc(s): return s.translate(COMP)[::-1]
def faidx(a, b):
    r = subprocess.run(["samtools", "faidx", REF, f"{ENS}:{a+1}-{b}"], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()
def gc(s): return (s.count("G")+s.count("C"))/len(s)
def tm(s): return mt.Tm_NN(s)

ha5 = faidx(CUT0-45, CUT0); ha3 = faidx(CUT0, CUT0+45)
fwd_full = ha5 + DM416
rev_full = rc(ha3) + DM417
print(f"{SITE}: cut {CHROM}:{CUT0+1}  HA5({gc(ha5)*100:.0f}%) HA3({gc(ha3)*100:.0f}%)")

# ---------- 1. donor row ----------
DHDR = "site role chrom guide_protospacer_20nt pam strand cut_genomic_GRCz11 mit_spec offtarget_count HA5_45nt HA3_45nt fwd_primer_tail rev_primer_tail notes fwd_primer_full rev_primer_full".split()
drow = {k:"" for k in DHDR}
drow.update(site=SITE, role="pIGLET24b_reguide", chrom=CHROM, guide_protospacer_20nt=GUIDE, pam=PAM,
            strand=STRAND, cut_genomic_GRCz11=f"{CHROM}:{CUT0+1}", mit_spec=str(MIT), offtarget_count=str(OT),
            HA5_45nt=ha5, HA3_45nt=ha3, fwd_primer_tail=ha5, rev_primer_tail=rc(ha3),
            notes="re-guide of pIGLET24b: cut moved +260bp into GC-balanced arms (51/38%) so donor amplifies; still in inert peak",
            fwd_primer_full=fwd_full, rev_primer_full=rev_full)
def upsert(path, hdr, rows, keyidx=0):
    lines = open(path).read().splitlines()
    keys = {r[keyidx] if isinstance(r, list) else r[hdr[keyidx]] for r in rows}
    keep = [ln for ln in lines if ln.split("\t")[keyidx] not in keys]
    with open(path, "w") as f:
        f.write("\n".join(keep) + "\n" + "\n".join("\t".join(r[h] for h in hdr) for r in rows) + "\n")
upsert("results/knockin_donor_design.tsv", DHDR, [drow])

# ---------- 2. genotyping primers (Tm-matched, longer OK for high Tm) ----------
LENGTHS = tuple(range(20, 33))
def good(p): return 0.38<=gc(p)<=0.62 and p[-1] in "GC" and not any(b*6 in p for b in "ACGT")
def pick(cands, target):
    v = [(abs(tm(p)-target), tm(p), p, pos) for p, pos in cands if good(p)]
    return None if not v else min(v)[1:]
def design_pair(cut0, target=62.0):
    lo, hi = cut0-300, cut0+980
    win = faidx(lo, hi); base = lo
    fwd = [(win[s-base:s-base+L], s) for s in range(cut0-280, cut0-130) for L in LENGTHS if s+L<=cut0-5]
    rev = [(rc(win[e-L-base:e-base]), e) for e in range(cut0+380, cut0+900) for L in LENGTHS if e-L>=cut0+5 and e<=hi]
    F = pick(fwd, target); R = pick(rev, F[0]) if F else None
    return F, R
(tF, Fseq, Fs), (tR, Rseq, Re) = design_pair(CUT0)
ins = CUT0+1; Fpos1 = Fs+1; Rpos1 = Re
wt = Rpos1-Fpos1+1; j5 = (ins-Fpos1)+OFF5; j3 = (Rpos1-ins)+OFF3
print(f"  genotyping F {Fseq} (Tm {tF:.1f}) / R {Rseq} (Tm {tR:.1f}) | WT {wt}bp | 5'junc {j5} | 3'junc {j3}")

GHDR = "locus chrom insertion_site genomic_F F_Tm F_pos genomic_R R_Tm R_pos WT_amplicon_bp".split()
grow = dict(locus=SITE, chrom=CHROM, insertion_site=str(ins), genomic_F=Fseq, F_Tm=f"{tF:.1f}",
            F_pos=f"{CHROM}:{Fpos1}", genomic_R=Rseq, R_Tm=f"{tR:.1f}", R_pos=f"{CHROM}:{Rpos1}", WT_amplicon_bp=str(wt))
upsert("results/genotyping_primers.tsv", GHDR, [grow])

JHDR = "locus chrom genomic_F StayGold_R_universal amplicon_5pjunction_bp StayGold_F_universal genomic_R amplicon_3pjunction_bp".split()
jrow = dict(locus=SITE, chrom=CHROM, genomic_F=Fseq, StayGold_R_universal=SG_R, amplicon_5pjunction_bp=str(j5),
            StayGold_F_universal=SG_F, genomic_R=Rseq, amplicon_3pjunction_bp=str(j3))
upsert("results/junction_primers.tsv", JHDR, [jrow])

# ---------- 3. GenBank maps ----------
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
w0, w1 = min(Fpos1-1, CUT0)-20, max(Rpos1, CUT0)+20
win = faidx(w0, w1); cutw = CUT0-w0
# WT
rec = SeqRecord(Seq(win), id=SITE, name=f"WT_{SITE}", description=f"WT {CHROM} pIGLET24b re-guide (GRCz11)")
rec.annotations["molecule_type"] = "DNA"
fs = [feat(0, len(win), f"genomic {CHROM} pIGLET24b re-guide", "source")]
for seq_, lab in [(Fseq, f"genotyping_F (Tm {tF:.1f})"), (Rseq, f"genotyping_R (Tm {tR:.1f})")]:
    h = find(win, seq_)
    if h: fs.append(feat(h[0], h[1], lab, "primer_bind", h[2]))
fs += [feat(cutw-45, cutw, "5' homology arm (45 bp)"), feat(cutw, cutw+45, "3' homology arm (45 bp)")]
h = find(win, GUIDE)
if h:
    fs.append(feat(h[0], h[1], f"sgRNA spacer: {GUIDE}", "misc_feature", h[2]))
    fs.append(feat(h[1], h[1]+3, f"PAM: {PAM}", "misc_feature", 1) if h[2]==1 else feat(h[0]-3, h[0], f"PAM: {PAM}", "misc_feature", -1))
fs.append(feat(cutw-1, cutw+1, "Cas9 blunt cut = insertion site"))
rec.features = sorted(fs, key=lambda f: int(f.location.start))
SeqIO.write(rec, f"{OUTDIR}/WT_{SITE}.gb", "genbank")
# INTEGRATED
left, right = win[:cutw], win[cutw:]
integ = left + cas_seq + right; off = len(left); roff = off+len(cas_seq)
irec = SeqRecord(Seq(integ), id=SITE, name=f"INT_{SITE}", description=f"INTEGRATED attP cassette at {CHROM} pIGLET24b re-guide")
irec.annotations["molecule_type"] = "DNA"
ifs = [feat(0, len(integ), f"INTEGRATED allele {CHROM} pIGLET24b_v2 (attP cassette)", "source"),
       feat(off-45, off, "5' homology arm (45 bp)"), feat(roff, roff+45, "3' homology arm (45 bp)"),
       feat(off-1, off+1, "5' junction (genomic|cassette)"), feat(roff-1, roff+1, "3' junction (cassette|genomic)")]
for cf in cas_feats: ifs.append(feat(cf["start"]+off, cf["end"]+off, cf["label"] or cf["type"], cf["type"], cf["strand"]))
hF = find(left, Fseq)
if hF: ifs.append(feat(hF[0], hF[1], f"genotyping_F (Tm {tF:.1f})", "primer_bind", hF[2]))
hR = find(right, Rseq)
if hR: ifs.append(feat(hR[0]+roff, hR[1]+roff, f"genotyping_R (Tm {tR:.1f})", "primer_bind", hR[2]))
for seq_, lab in [(SG_R, "StayGold_R (5' junction)"), (SG_F, "StayGold_F (3' junction)")]:
    h = find(cas_seq, seq_)
    if h: ifs.append(feat(h[0]+off, h[1]+off, lab, "primer_bind", h[2]))
irec.features = sorted(ifs, key=lambda f: int(f.location.start))
SeqIO.write(irec, f"{OUTDIR}/INTEGRATED_{SITE}_attP.gb", "genbank")
print(f"  wrote WT_{SITE}.gb ({len(win)}bp) + INTEGRATED_{SITE}_attP.gb ({len(integ)}bp)")
print("done.")
