#!/usr/bin/env python3
"""
27_knit_biotin_donors.py

Derive KNIT (nickase) biotinylated donor oligos from the existing HDR donor
ultramers in results/oligo_master.tsv.

KNIT editing (Gao et al., Nature 2026; doi:10.1038/s41586-026-10819-7) inserts
DNA at an nCas9(D10A) single-strand nick by tethering the donor to a Cas9-fused
monomeric streptavidin (mSA) via a 5' biotin. Our current donors carry /5AmMC6/
(5' amine) for end protection. The ONLY donor change needed for KNIT is to swap
the 5' modifier /5AmMC6/ -> /5Biosg/ (5' biotin), keeping the 5' phosphorothioate
linkages and the 3'-P. Same homology arm, same cassette-priming tail, same guide,
same locus -> the donor still PCR-amplifies off the StayGold cassette identically.

Superseded donors (annotation contains 'SUPERSEDED') are skipped. New oligos are
numbered continuously from START_DM. Both strands (donor_F and donor_R) are
biotinylated so the donor is tethered regardless of orientation; order single-end
(F-only or R-only) variants later to optimize nick/donor geometry if desired.

Input : results/oligo_master.tsv
Output: results/knit_biotin_donor_oligos.tsv
"""
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "results", "oligo_master.tsv")
OUT = os.path.join(REPO, "results", "knit_biotin_donor_oligos.tsv")
START_DM = 448
OLD_MOD, NEW_MOD = "/5AmMC6/", "/5Biosg/"
FIELDS = ["name", "locus", "role", "status", "sequence", "annotation"]


def main():
    with open(SRC, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    donors = [
        r for r in rows
        if r["role"].startswith("donor_")
        and OLD_MOD in r["sequence"]
        and "SUPERSEDED" not in (r.get("annotation") or "")
    ]

    out, n = [], START_DM
    for r in donors:
        out.append({
            "name": f"DM{n}",
            "locus": r["locus"],
            "role": f"{r['role']} (KNIT biotin)",
            "status": "TO ORDER (KNIT)",
            "sequence": r["sequence"].replace(OLD_MOD, NEW_MOD),
            "annotation": (
                f"KNIT donor = {r['name']} with {NEW_MOD} replacing {OLD_MOD}; "
                f"keep 5' phosphorothioate + 3'-P. Biotin tethers the donor to the "
                f"nCas9(D10A)-mSA fusion at the nick."
            ),
        })
        n += 1

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"Wrote {len(out)} KNIT biotin donor oligos "
          f"({out[0]['name']}-{out[-1]['name']}) -> results/{os.path.basename(OUT)}")
    for r in out:
        print(f"  {r['name']}  {r['locus']:<22}  {r['role']}")


if __name__ == "__main__":
    main()
