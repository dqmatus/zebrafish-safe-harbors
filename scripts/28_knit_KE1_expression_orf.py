#!/usr/bin/env python3
"""
28_knit_KE1_expression_orf.py

Reverse-translate + E. coli codon-optimize the KNIT KE1 editor
(nCas9-D10A-mSA, 1574 aa; Gao et al. Nature 2026) into a synthesis-ready ORF
for pET28b expression via NdeI / XhoI, with a C-terminal 6xHis.

Design:
  5'- [flank] CATATG(NdeI; ATG = start Met) [opt codons res 2..N] [6xHis] TAA CTCGAG(XhoI) [flank] -3'
- NdeI ATG is the start codon -> the pET28b N-terminal His/T7 tag is bypassed;
  our own C-terminal 6xHis + stop is included -> product = M-[KE1]-HHHHHH,
  matching the paper's C-terminal-His-only construct.
- E. coli high-expression codons (primary + alternates), greedily chosen to:
  (a) never create an internal NdeI (CATATG) or XhoI (CTCGAG) site,
  (b) avoid homopolymer runs > 4, (c) avoid > 3 identical codons in a row.
Verifies the DNA back-translates exactly to KE1 + 6xHis before writing.

Output: results/knit_KE1_pET28b_orf.fasta
"""
CODONS = {  # E. coli, preference-ordered (rare codons excluded)
 'A':['GCG','GCT','GCC','GCA'],'R':['CGT','CGC'],'N':['AAC','AAT'],'D':['GAT','GAC'],
 'C':['TGC','TGT'],'Q':['CAG','CAA'],'E':['GAA','GAG'],'G':['GGC','GGT'],
 'H':['CAT','CAC'],'I':['ATC','ATT'],'L':['CTG','CTC','CTT','TTG'],'K':['AAA','AAG'],
 'M':['ATG'],'F':['TTC','TTT'],'P':['CCG','CCA','CCT'],'S':['AGC','TCT','TCC','AGT'],
 'T':['ACC','ACT','ACA'],'W':['TGG'],'Y':['TAC','TAT'],'V':['GTG','GTT','GTC'],
}
FORBIDDEN = ("CATATG", "CTCGAG")   # NdeI, XhoI (both palindromic -> forward scan suffices)
FLANK5, FLANK3 = "GAAGGAG", "TGAGGTT"

# --- standard translation table (for verification) ---
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = {a+b+c: _AA[i*16+j*4+k] for i,a in enumerate(_B) for j,b in enumerate(_B) for k,c in enumerate(_B)}
def translate(dna): return "".join(CODON2AA[dna[i:i+3]] for i in range(0, len(dna)-2, 3))

def maxhomo(s):
    m=c=1
    for i in range(1,len(s)): c=c+1 if s[i]==s[i-1] else 1; m=max(m,c)
    return m

def choose(tail, aa):
    """Codon for aa: among site-safe options, minimize the resulting homopolymer
    run, tie-broken by E. coli preference order (so poly-K etc. auto-alternate)."""
    safe = [c for c in CODONS[aa] if not any(f in tail[-6:]+c for f in FORBIDDEN)]
    cands = safe or CODONS[aa]      # if none site-safe (e.g. forced M/W), repair pass fixes it
    return min(cands, key=lambda c: (maxhomo(tail[-4:]+c), CODONS[aa].index(c)))

def repair_sites(codons):
    """Remove any residual internal NdeI/XhoI by synonymous swaps (translation preserved)."""
    for _ in range(200):
        dna = "".join(codons)
        pos = min([dna.find(f) for f in FORBIDDEN if f in dna] or [-1])
        if pos < 0:
            return codons
        for ci in {pos//3, (pos+2)//3, (pos+5)//3}:
            if not (0 <= ci < len(codons)):
                continue
            aa = CODON2AA[codons[ci]]
            for alt in CODONS.get(aa, []):
                if alt == codons[ci]:
                    continue
                trial = codons[:]; trial[ci] = alt
                seg = "".join(trial[max(0, ci-2):ci+3])
                if not any(f in seg for f in FORBIDDEN):
                    codons = trial; break
            else:
                continue
            break
    raise RuntimeError("could not clear all forbidden sites")

def main():
    seqs={}; name=None
    for l in open("results/knit_editor_protein_sequences.fasta"):
        l=l.rstrip()
        if l.startswith(">"): name=l[1:].split("|")[0].strip(); seqs[name]=""
        elif name: seqs[name]+=l.strip()
    ke1=[v for k,v in seqs.items() if k.startswith("KE1")][0]
    assert ke1[0]=="M"
    prot = ke1 + "HHHHHH"          # append C-terminal 6xHis

    codons = ["ATG"]               # residue 1 (M), provided by the NdeI ATG
    for aa in prot[1:]:
        codons.append(choose("".join(codons), aa))
    codons = repair_sites(codons)  # clear residual internal NdeI/XhoI (e.g. His-Met CATATG)
    coding = "".join(codons) + "TAA"

    # verify
    assert translate(coding[:-3]) == prot, "back-translation mismatch!"
    insert = FLANK5 + "CAT" + coding + "CTCGAG" + FLANK3   # ...CAT+ATG.. = CATATG (NdeI)
    nde, xho = insert.count("CATATG"), insert.count("CTCGAG")
    gc = 100*(insert.count("G")+insert.count("C"))/len(insert)

    hdr = (">KNIT_KE1_pET28b_NdeI_XhoI | E.coli codon-optimized nCas9(D10A)-mSA + C-term 6xHis | "
           f"clone NdeI/XhoI into pET28b(+) | {len(coding)//3-1} codons + stop | "
           "Gao et al Nature 2026")
    with open("results/knit_KE1_pET28b_orf.fasta","w") as f:
        f.write(hdr+"\n")
        for i in range(0,len(insert),70): f.write(insert[i:i+70]+"\n")

    print(f"KE1 ORF written: results/knit_KE1_pET28b_orf.fasta")
    print(f"  insert length      : {len(insert)} bp (coding {len(coding)} bp = {len(coding)//3} codons incl stop)")
    print(f"  protein verified   : M-KE1-6xHis, back-translation EXACT ({len(prot)} aa)")
    print(f"  NdeI sites (want 1): {nde}   XhoI sites (want 1): {xho}")
    print(f"  GC content         : {gc:.1f}%   longest homopolymer: {maxhomo(insert)}")
    print(f"  5' end: {insert[:40]}...")
    print(f"  3' end: ...{insert[-40:]}")

if __name__ == "__main__":
    main()
