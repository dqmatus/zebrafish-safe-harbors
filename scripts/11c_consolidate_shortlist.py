#!/usr/bin/env python3
"""
11c_consolidate_shortlist.py

Audit / consolidate the autonomous sgRNA shortlist across the designed loci. For each locus
in results/locus_design_summary.tsv it enumerates local NGG guides in a window around the
design cut, scores them with the autonomous method (11b: bwa off-target + CRISPRscan
on-target efficiency), ranks them SAFETY-FIRST with the exonic-off-target HARD GATE, and
writes the autonomous top pick alongside the originally-designed/injected guide and its rank
-- so you can see at a glance which injected guides the local method endorses, which it would
swap, and why.

Ranking per locus:  ot_exon_le3mm ASC (hard gate) -> ot_le3mm ASC -> CRISPRscan DESC -> dist_center ASC

Note: guides are pooled from a +/-window around the design cut (uniform across the de novo,
COPE and pIGLET loci) and capped to the nearest-center KCAP for speed, so the top pick is the
best within that local pool. Use --window (larger) and --kcap 0 (no cap) for a full-peak scan
(slower). The injected-guide safety assessment is parameter-stable; the exact top pick is not.

Inputs : results/locus_design_summary.tsv  (locus, cut_GRCz11, guide_20nt)
         reference/GRCz11.fa (+ bwa index), annotation/Danio_rerio.GRCz11.110.gtf
Output : results/sgrna_shortlist_auto.tsv

Usage: python3 scripts/11c_consolidate_shortlist.py [--window 350] [--kcap 12]
       (run in the `zebrafish_atac` conda env: needs bwa + samtools on PATH)
"""
import csv
import importlib.util
import os
import sys

SUMMARY = "results/locus_design_summary.tsv"
OUT = "results/sgrna_shortlist_auto.tsv"
WINDOW = int(sys.argv[sys.argv.index("--window") + 1]) if "--window" in sys.argv else 350
KCAP = int(sys.argv[sys.argv.index("--kcap") + 1]) if "--kcap" in sys.argv else 12   # 0 = no cap


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    s11 = _load("s11", os.path.join(here, "11_design_sgrnas.py"))
    s11b = _load("s11b", os.path.join(here, "11b_offtarget_bwa.py"))
    if not os.path.exists(SUMMARY):
        sys.exit(f"{SUMMARY} not found")
    ex = s11b.load_exons()

    cols = ["locus", "chrom", "top_guide", "pam", "strand", "top_cut", "ot_le3mm",
            "ot_exon_le3mm", "crisprscan", "dist_center", "flag",
            "injected_guide", "injected_rank", "injected_flag"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in csv.DictReader(open(SUMMARY), delimiter="\t"):
            locus = r["locus"]
            chrom_ucsc, cut = r["cut_GRCz11"].split(":")
            cut = int(cut)
            injected = r["guide_20nt"]
            seq = s11.fetch(s11.ucsc_to_ensembl(chrom_ucsc), cut - WINDOW, cut + WINDOW)
            if not seq:
                print(f"{locus}: no sequence", flush=True)
                continue
            guides = [g for g in s11.find_guides(seq, cut - WINDOW) if s11.quality_ok(g["proto"])[0]]
            for g in guides:
                g["dist_center"] = abs(g["cut"] - cut)
                g["cand"] = chrom_ucsc
            guides.sort(key=lambda g: g["dist_center"])
            pool = guides if KCAP <= 0 else guides[:KCAP]
            pool += [g for g in guides if g["proto"] == injected and g not in pool]
            for g in pool:
                s = s11b.score_guide(g["proto"], ex)
                eff = s11b.calc_crisprscan(
                    s11b.guide_context35(g["cand"], g["gstart"], g["gend"], g["strand"]))
                g.update(s, crisprscan=(-1 if eff is None else eff), flag=s11b.flag_of(s))
            pool.sort(key=lambda g: (g["ot_exon_le3mm"], g["ot_le3mm"],
                                     -g["crisprscan"], g["dist_center"]))
            top = pool[0]
            if top["ot_exon_le3mm"] > 0:
                top["flag"] = "REVIEW: no exon-clean guide"
            inj = next((g for g in pool if g["proto"] == injected), None)
            inj_rank = str(pool.index(inj) + 1) if inj else "filtered(on-target)"
            inj_flag = inj["flag"] if inj else "n/a"
            fh.write("\t".join(str(x) for x in [
                locus, chrom_ucsc, top["proto"], top["pam"], top["strand"], top["cut"],
                top["ot_le3mm"], top["ot_exon_le3mm"],
                top["crisprscan"] if top["crisprscan"] >= 0 else "",
                top["dist_center"], top["flag"], injected, inj_rank, inj_flag]) + "\n")
            fh.flush()
            print(f"{locus}: top={top['proto']} eff={top['crisprscan']} flag={top['flag']} "
                  f"| injected rank={inj_rank}", flush=True)
    print(f"\n[out] {OUT}")


if __name__ == "__main__":
    main()
