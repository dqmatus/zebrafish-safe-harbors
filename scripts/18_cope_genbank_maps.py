#!/usr/bin/env python3
"""
18_cope_genbank_maps.py

Generate annotated WT + INTEGRATED(attP) GenBank maps for the two COPE lab
targets (chr1, chr14), matching the convention of the original 4 loci.

Reads the committed deliverables so annotations stay consistent:
  results/knockin_donor_design.tsv  (guide, pam, strand, cut, HA5, HA3)
  results/genotyping_primers.tsv    (genomic F/R)
  results/junction_primers.tsv      (universal StayGold primers)
Cassette (Benchling .gb whose LOCUS line breaks Biopython's reader) is parsed
manually for sequence + features; maps are WRITTEN with Biopython.

Outputs to the constructs folder:
  .../landing_pads/integrated_alleles/{WT_cope_chrN.gb, INTEGRATED_cope_chrN_attP.gb}
"""
import csv, re, subprocess, sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

REF = "reference/GRCz11.fa"
CDIR = "/Users/dsm16/Documents/Betzig/Swinburne-Gokul-Betzig constructs/fish CRISPR/landing_pads"
CASSETTE = f"{CDIR}/spacer-attP-GT-unc45b-lck-mSG-SV40-attP-GA-spacer.gb"
OUTDIR = f"{CDIR}/integrated_alleles"
ENS = {"chr1": "1", "chr14": "14"}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(s): return s.translate(COMP)[::-1]

def faidx(ens, a, b):  # 0-based half-open
    r = subprocess.run(["samtools", "faidx", REF, f"{ens}:{a+1}-{b}"], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()

def parse_cassette(path):
    text = open(path).read()
    seq = re.sub(r"[^acgtACGT]", "", text.split("ORIGIN")[1].split("//")[0]).upper()
    feats = []
    cur = None
    for ln in text.split("FEATURES", 1)[1].split("ORIGIN")[0].splitlines()[1:]:
        m = re.match(r"^ {5}(\S+)\s+(complement\()?(\d+)\.\.(\d+)\)?\s*$", ln)
        if m:
            if cur:
                feats.append(cur)
            cur = dict(type=m.group(1), strand=-1 if m.group(2) else 1,
                       start=int(m.group(3)) - 1, end=int(m.group(4)), label=None)
        elif cur is not None:
            q = re.search(r'/(?:label|note)="?([^"\n]+)"?', ln)
            if q and cur["label"] is None:
                cur["label"] = q.group(1).strip().strip('"')
    if cur:
        feats.append(cur)
    return seq, feats

def feat(start, end, label, typ="misc_feature", strand=1):
    f = SeqFeature(FeatureLocation(int(start), int(end), strand=strand), type=typ)
    f.qualifiers["label"] = [label]
    return f

def find(win, sub):
    """Return (start,end,strand) of sub in win (search both strands), or None."""
    i = win.find(sub)
    if i >= 0:
        return (i, i + len(sub), 1)
    i = win.find(rc(sub))
    if i >= 0:
        return (i, i + len(sub), -1)
    return None

# ---- load committed values ----
donor = {r["site"]: r for r in csv.DictReader(open("results/knockin_donor_design.tsv"), delimiter="\t")}
geno = {r["locus"]: r for r in csv.DictReader(open("results/genotyping_primers.tsv"), delimiter="\t")}
junc = {r["locus"]: r for r in csv.DictReader(open("results/junction_primers.tsv"), delimiter="\t")}
SG_R = next(iter(junc.values()))["StayGold_R_universal"]
SG_F = next(iter(junc.values()))["StayGold_F_universal"]
cas_seq, cas_feats = parse_cassette(CASSETTE)
print(f"cassette: {len(cas_seq)}bp, {len(cas_feats)} features parsed")

for site in ("cope_chr1", "cope_chr14"):
    d, g, j = donor[site], geno[site], junc[site]
    chrom = d["chrom"]; ens = ENS[chrom]
    cut0 = int(d["cut_genomic_GRCz11"].split(":")[1]) - 1     # 0-based break
    guide, pam, gstrand = d["guide_protospacer_20nt"], d["pam"], d["strand"]
    ha5, ha3 = d["HA5_45nt"], d["HA3_45nt"]
    Fseq, Rseq = g["genomic_F"], g["genomic_R"]
    Fpos0 = int(g["F_pos"].split(":")[1]) - 1
    Rpos = int(g["R_pos"].split(":")[1])                      # rightmost genomic of R
    w0, w1 = min(Fpos0, cut0) - 20, max(Rpos, cut0) + 20
    win = faidx(ens, w0, w1)
    cutw = cut0 - w0
    gene = "COPE-validated safe harbor"

    # ---------- WT ----------
    rec = SeqRecord(Seq(win), id=site, name=f"WT_{site}",
                    description=f"WT {chrom} {gene} (GRCz11)")
    rec.annotations["molecule_type"] = "DNA"
    fs = [feat(0, len(win), f"genomic {chrom} near {gene}", "source")]
    for seq_, lab, typ in [(Fseq, f"genotyping_F (Tm {g['F_Tm']})", "primer_bind"),
                           (Rseq, f"genotyping_R (Tm {g['R_Tm']})", "primer_bind")]:
        h = find(win, seq_)
        if h: fs.append(feat(h[0], h[1], lab, typ, h[2]))
    fs.append(feat(cutw - 45, cutw, "5' homology arm (45 bp)"))
    fs.append(feat(cutw, cutw + 45, "3' homology arm (45 bp)"))
    h = find(win, guide)
    if h:
        fs.append(feat(h[0], h[1], f"sgRNA spacer: {guide}", "misc_feature", h[2]))
        # PAM sits immediately 3' of the protospacer on the guide strand:
        #  + strand -> just after the match; - strand -> just before it (top-strand CCN)
        if h[2] == 1:
            fs.append(feat(h[1], h[1] + 3, f"PAM: {pam}", "misc_feature", 1))
        else:
            fs.append(feat(h[0] - 3, h[0], f"PAM: {pam}", "misc_feature", -1))
    fs.append(feat(cutw - 1, cutw + 1, "Cas9 blunt cut = insertion site"))
    rec.features = sorted(fs, key=lambda f: int(f.location.start))
    from Bio import SeqIO
    SeqIO.write(rec, f"{OUTDIR}/WT_{site}.gb", "genbank")

    # ---------- INTEGRATED ----------
    left = win[:cutw]; right = win[cutw:]
    integ = left + cas_seq + right
    off = len(left)                       # cassette starts here
    roff = len(left) + len(cas_seq)       # right genomic starts here
    irec = SeqRecord(Seq(integ), id=site, name=f"INT_{site}",
                     description=f"INTEGRATED attP cassette at {chrom} {gene}")
    irec.annotations["molecule_type"] = "DNA"
    ifs = [feat(0, len(integ), f"INTEGRATED allele {chrom} (attP cassette)", "source"),
           feat(off - 45, off, "5' homology arm (45 bp)"),
           feat(roff, roff + 45, "3' homology arm (45 bp)"),
           feat(off - 1, off + 1, "5' junction (genomic|cassette)"),
           feat(roff - 1, roff + 1, "3' junction (cassette|genomic)")]
    for cf in cas_feats:
        ifs.append(feat(cf["start"] + off, cf["end"] + off,
                        cf["label"] or cf["type"], cf["type"], cf["strand"]))
    # primers on the integrated allele
    hF = find(left, Fseq)
    if hF: ifs.append(feat(hF[0], hF[1], f"genotyping_F (Tm {g['F_Tm']})", "primer_bind", hF[2]))
    hR = find(right, Rseq)
    if hR: ifs.append(feat(hR[0] + roff, hR[1] + roff, f"genotyping_R (Tm {g['R_Tm']})", "primer_bind", hR[2]))
    for seq_, lab in [(SG_R, "StayGold_R (5' junction)"), (SG_F, "StayGold_F (3' junction)")]:
        h = find(cas_seq, seq_)
        if h: ifs.append(feat(h[0] + off, h[1] + off, lab, "primer_bind", h[2]))
    irec.features = sorted(ifs, key=lambda f: int(f.location.start))
    SeqIO.write(irec, f"{OUTDIR}/INTEGRATED_{site}_attP.gb", "genbank")

    print(f"{site}: WT {len(win)}bp ({len(rec.features)} feats) + INTEGRATED {len(integ)}bp ({len(irec.features)} feats)")

print("\nWrote 4 maps to", OUTDIR)
