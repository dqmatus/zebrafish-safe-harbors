#!/usr/bin/env python3
"""
11_design_sgrnas.py

Find candidate SpCas9 sgRNA sites within the top safe-harbor candidates.

- Reads top-N candidates from results/candidate_feature_table.tsv (by composite_rank)
- Extracts genomic sequence from reference/GRCz11.fa (Ensembl-named) via samtools,
  translating UCSC-style candidate chrom names (chr12, chrUn_KN..v1) -> Ensembl (12, KN...x)
- Enumerates 20nt protospacers with an NGG PAM on both strands
- Applies on-target quality filters (GC, homopolymer / poly-T terminator)
- Scores preference for cut sites near the region center (good for central knock-in)
- Writes results/sgrna_candidates.tsv

This is an ON-TARGET ENUMERATOR ONLY. It does NOT score genome-wide off-target
specificity (a local cas-offinder run was unreliable on this machine's OpenCL
backend). For specificity (MIT/CFD) + efficiency (Doench) scoring, submit the
candidate regions (results/top5_candidates_for_crispor.fa) to CRISPOR (danRer11);
the curated, specificity-ranked picks live in results/sgrna_shortlist.tsv.

Usage: python3 scripts/11_design_sgrnas.py [TOP_N] [--include-scaffolds]
"""
import csv
import subprocess
import sys

REF = "reference/GRCz11.fa"
TABLE = "results/candidate_feature_table.tsv"
OUT = "results/sgrna_candidates.tsv"

TOP_N = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 5
INCLUDE_SCAFFOLDS = "--include-scaffolds" in sys.argv

GC_MIN, GC_MAX = 0.40, 0.70
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def ucsc_to_ensembl(chrom):
    """chr12 -> 12 ; chrM -> MT ; chrUn_KN150617v1 -> KN150617.1"""
    if chrom.startswith("chrUn_"):
        tail = chrom[len("chrUn_"):]              # KN150617v1
        if "v" in tail:
            acc, ver = tail.rsplit("v", 1)
            return f"{acc}.{ver}"
        return tail
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if c in ("M", "MT") else c


def revcomp(s):
    return s.translate(COMP)[::-1]


def gc(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0


def fetch(chrom_ens, start, end):
    """1-based inclusive samtools faidx region from 0-based half-open BED coords."""
    region = f"{chrom_ens}:{start + 1}-{end}"
    out = subprocess.run(["samtools", "faidx", REF, region],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return None
    return "".join(out.stdout.splitlines()[1:]).upper()


def quality_ok(proto):
    if "TTTT" in proto:                # Pol III terminator
        return False, "polyT"
    if not (GC_MIN <= gc(proto) <= GC_MAX):
        return False, f"GC={gc(proto):.2f}"
    for b in "ACGT":                   # long homopolymer
        if b * 5 in proto:
            return False, f"homopolymer{b}"
    return True, "ok"


def find_guides(seq, gstart):
    """Yield dicts for every NGG-PAM 20mer on both strands of seq.
    gstart = genomic start (0-based) of seq[0]."""
    guides = []
    n = len(seq)
    # + strand: protospacer[0:20] then PAM = seq[i+20:i+23] == .GG ; cut between 17/18 of proto
    for i in range(0, n - 23 + 1):
        pam = seq[i + 20:i + 23]
        if pam[1:3] == "GG":
            proto = seq[i:i + 20]
            if "N" in proto:
                continue
            cut = gstart + i + 17                 # blunt cut 3bp 5' of PAM
            guides.append(dict(proto=proto, pam=pam, strand="+",
                               cut=cut, gstart=gstart + i, gend=gstart + i + 23))
    # - strand: PAM = CCN on the top strand at i..i+3, protospacer = revcomp(seq[i+3:i+23])
    for i in range(0, n - 23 + 1):
        pam_top = seq[i:i + 3]
        if pam_top[0:2] == "CC":
            proto = revcomp(seq[i + 3:i + 23])
            if "N" in proto:
                continue
            cut = gstart + i + 6                  # 3bp into protospacer from PAM, on - strand
            guides.append(dict(proto=proto, pam=revcomp(pam_top), strand="-",
                               cut=cut, gstart=gstart + i, gend=gstart + i + 23))
    return guides


def main():
    rows = list(csv.DictReader(open(TABLE), delimiter="\t"))
    for r in rows:
        r["_cr"] = float(r["composite_rank"]) if r["composite_rank"] not in ("", "NA") else -1
    rows.sort(key=lambda r: -r["_cr"])
    if not INCLUDE_SCAFFOLDS:
        rows = [r for r in rows if "_" not in r["chrom"] and r["chrom"] not in ("chrM",)]
    top = rows[:TOP_N]

    out_rows = []
    for rank, r in enumerate(top, 1):
        chrom, s, e = r["chrom"], int(r["start"]), int(r["end"])
        center = (s + e) // 2
        ens = ucsc_to_ensembl(chrom)
        seq = fetch(ens, s, e)
        if not seq:
            print(f"  WARN: no sequence for {chrom} ({ens}) — skipping", file=sys.stderr)
            continue
        guides = find_guides(seq, s)
        kept = []
        for g in guides:
            ok, reason = quality_ok(g["proto"])
            if not ok:
                continue
            g["rank"] = rank
            g["cand"] = f"{chrom}:{s}-{e}"
            g["gene"] = r["nearest_gene"]
            g["gc"] = round(gc(g["proto"]), 2)
            g["dist_center"] = abs(g["cut"] - center)
            kept.append(g)
        kept.sort(key=lambda g: g["dist_center"])
        out_rows.extend(kept)
        print(f"#{rank} {chrom}:{s}-{e} ({r['nearest_gene']}, {e-s}bp): "
              f"{len(guides)} NGG sites -> {len(kept)} pass on-target filters")

    cols = ["rank", "cand", "gene", "strand", "proto", "pam", "gc",
            "cut", "dist_center", "gstart", "gend"]
    with open(OUT, "w") as f:
        f.write("\t".join(cols) + "\n")
        for g in out_rows:
            f.write("\t".join(str(g[c]) for c in cols) + "\n")
    print(f"\nWrote {len(out_rows)} candidate guides -> {OUT}")


if __name__ == "__main__":
    main()
