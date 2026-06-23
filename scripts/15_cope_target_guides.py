#!/usr/bin/env python3
"""
15_cope_target_guides.py

On-target SpCas9 sgRNA enumeration for the two COPE-validated lab targets
(chr1, chr14) surfaced by the in_cope flag in candidate_feature_table.tsv.

Same logic as 11_design_sgrnas.py (NGG both strands; GC 0.40-0.70; reject polyT
terminator + >=5 homopolymer; blunt cut 3bp 5' of PAM), but restricted to the
candidate ATAC window (cut stays in open chromatin) and ranked by distance of the
cut from the COPE centre (the DANIO-CODE-validated open-but-inert core).

ON-TARGET ONLY. Submit results/cope_candidates_for_crispor.fa to CRISPOR
(danRer11) for MIT/CFD specificity + Doench efficiency + off-target counts, then
curate as we did for the original loci (sgrna_shortlist.tsv).

Outputs:
  results/cope_guides_chr1_chr14.tsv
  results/cope_candidates_for_crispor.fa
"""
import csv, subprocess, sys

REF = "reference/GRCz11.fa"
OUT = "results/cope_guides_chr1_chr14.tsv"
FA = "results/cope_candidates_for_crispor.fa"
GC_MIN, GC_MAX = 0.40, 0.70
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")

# (label, ucsc_chrom, ensembl_chrom, cand_start, cand_end, cope_start, cope_end)
TARGETS = [
    ("cope_chr1",  "chr1",  "1",  15549795, 15551236, 15550222, 15550860),
    ("cope_chr14", "chr14", "14", 43942077, 43942572, 43941618, 43942707),
]


def revcomp(s): return s.translate(COMP)[::-1]
def gc(s): return (s.count("G") + s.count("C")) / len(s) if s else 0


def fetch(chrom_ens, start, end):
    region = f"{chrom_ens}:{start + 1}-{end}"
    out = subprocess.run(["samtools", "faidx", REF, region], capture_output=True, text=True)
    return None if out.returncode != 0 else "".join(out.stdout.splitlines()[1:]).upper()


def quality_ok(p):
    if "TTTT" in p: return False
    if not (GC_MIN <= gc(p) <= GC_MAX): return False
    return not any(b * 5 in p for b in "ACGT")


def find_guides(seq, gstart):
    guides, n = [], len(seq)
    for i in range(0, n - 23 + 1):
        if seq[i + 21:i + 23] == "GG":
            proto = seq[i:i + 20]
            if "N" not in proto:
                guides.append(dict(proto=proto, pam=seq[i + 20:i + 23], strand="+", cut=gstart + i + 17))
        if seq[i:i + 2] == "CC":
            proto = revcomp(seq[i + 3:i + 23])
            if "N" not in proto:
                guides.append(dict(proto=proto, pam=revcomp(seq[i:i + 3]), strand="-", cut=gstart + i + 6))
    return guides


def main():
    rows, fasta = [], []
    for label, ucsc, ens, cs, ce, cope_s, cope_e in TARGETS:
        cope_center = (cope_s + cope_e) // 2
        seq = fetch(ens, cs, ce)
        if not seq:
            print(f"  WARN: no sequence for {label}", file=sys.stderr); continue
        fasta.append(f">{label}_{ucsc}:{cs}-{ce}\n{seq}")
        kept = [g for g in find_guides(seq, cs) if quality_ok(g["proto"])]
        for g in kept:
            g.update(label=label, chrom=ucsc, gc=round(gc(g["proto"]), 2),
                     dist_cope_center=abs(g["cut"] - cope_center),
                     in_cope=int(cope_s <= g["cut"] < cope_e))
        kept.sort(key=lambda g: (not g["in_cope"], g["dist_cope_center"]))
        rows.extend(kept)
        print(f"{label} {ucsc}:{cs}-{ce}: {len(kept)} guides pass on-target filters "
              f"({sum(g['in_cope'] for g in kept)} cut inside the COPE)")

    cols = ["label", "chrom", "strand", "proto", "pam", "gc", "cut", "in_cope", "dist_cope_center"]
    with open(OUT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for g in rows:
            f.write("\t".join(str(g[c]) for c in cols) + "\n")
    with open(FA, "w") as f:
        f.write("\n".join(fasta) + "\n")
    print(f"\nWrote {len(rows)} guides -> {OUT}")
    print(f"Wrote CRISPOR FASTA -> {FA}")


if __name__ == "__main__":
    main()
