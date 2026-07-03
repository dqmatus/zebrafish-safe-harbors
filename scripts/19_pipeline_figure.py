#!/usr/bin/env python3
"""
19_pipeline_figure.py
Render schematics of the zebrafish safe-harbor discovery pipeline.
Counts are read LIVE from the current result tables.

Outputs:
  results/pipeline_schematic.png / .svg        (detailed, vertical)
  results/pipeline_schematic_talk.png / .svg   (simplified, landscape — for talks)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

C = dict(input="#E8EAED", access="#CFE3F7", excl="#FBE0D0", rank="#D6ECD3",
         design="#E3D7F0", tier="#FFF4CC", out="#2E5D34", edge="#4A4A4A")

def count(path, header=False):
    try:
        with open(path) as f:
            n = sum(1 for _ in f)
        return f"{(n - 1 if header else n):,}"
    except OSError:
        return "—"

N = dict(
    consensus=count("results/consensus_peaks_100pct.bed"),
    strict=count("results/candidate_feature_table.tsv", header=True),
    relaxed=count("results/candidate_feature_table_relaxed.tsv", header=True),
    intronic=count("results/candidate_feature_table_intronic.tsv", header=True),
)

def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(f"results/{stem}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


# ============================ detailed (vertical) ============================
def render_detailed():
    fig, ax = plt.subplots(figsize=(11, 14))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, title, detail, color, tsize=11, dsize=8.5):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.6,rounding_size=1.5",
                                    fc=color, ec=C["edge"], lw=1.3, zorder=2))
        ax.text(x, y+h/2-2.0, title, ha="center", va="top", fontsize=tsize, fontweight="bold", zorder=3)
        if detail:
            ax.text(x, y+h/2-5.2, detail, ha="center", va="top", fontsize=dsize, zorder=3)

    def arrow(y1, y2):
        ax.add_patch(FancyArrowPatch((50, y1), (50, y2), arrowstyle="-|>", mutation_scale=18,
                                     lw=1.6, color=C["edge"], zorder=1))

    ax.text(50, 98.5, "Zebrafish genomic safe-harbor discovery pipeline",
            ha="center", va="top", fontsize=15, fontweight="bold")
    ax.text(50, 95.5, "GRCz11 / danRer11  ·  DANIO-CODE developmental atlas",
            ha="center", va="top", fontsize=9.5, color="#555")

    inputs = ["Genome\nGRCz11", "DANIO-CODE\nATAC-seq\n(11 stages)", "ChIP-seq\nH3K27ac/K4me1/\nK4me3/CTCF/…",
              "Ensembl\ngenes · exons\nmiRNA · lncRNA", "Conservation\nGERP + CNEs", "COPE/DOPE\n+ validated CREs"]
    for x, t in zip([10, 26, 42, 58, 74, 90], inputs):
        ax.add_patch(FancyBboxPatch((x-7.5, 86.5), 15, 6.5, boxstyle="round,pad=0.4,rounding_size=1",
                                    fc=C["input"], ec=C["edge"], lw=1.0, zorder=2))
        ax.text(x, 89.7, t, ha="center", va="center", fontsize=7.3, zorder=3)
    ax.text(50, 84.5, "inputs", ha="center", fontsize=8, style="italic", color="#555")

    box(50, 77, 62, 8, "① Constitutive accessibility  (06–07)",
        "per-stage ATAC summits → collapse replicate→biosample→stage →\n"
        f"peaks open in 100% of 11 stages   →   {N['consensus']} candidate windows", C["access"])
    box(50, 65, 62, 8.5, "② Genomic-context exclusion  (05, 08)",
        "remove gene bodies ±10 kb · exons · essential genes ±20 kb ·\n"
        "transposons · telomeres · non-primary scaffolds · <300 bp", C["excl"])
    box(50, 50.5, 68, 11.5, "③ Regulatory / inertness screen  (08b — discard-whole)",
        "drop any peak overlapping:\n"
        "• active chromatin  H3K27ac · H3K4me1 · H3K4me3 · CTCF\n"
        "• Polycomb H3K27me3 · other ChIP (H2A.Z, PolII…) · validated CREs\n"
        "• conservation (GERP + CNEs) · low mappability (genmap K50)\n"
        "• buffered miRNA ±50 kb · lncRNA ±25 kb", C["excl"], dsize=8.2)
    box(50, 35, 62, 8.5, "④ Annotate & rank  (09)",
        "nearest gene · GC · conservation · sgRNA · TAD · stage count ·\n"
        "COPE/DOPE flag   →   composite rank", C["rank"])
    box(50, 22.5, 62, 6.5, "⑤ Validation benchmarks  (10)",
        "pIGLET14a/24b · tpl102 (exon-proximal → INCONCLUSIVE by design)", C["rank"], dsize=8.2)
    box(50, 10, 68, 9.5, "⑥ Knock-in reagent design  (11–18)",
        "sgRNA enumeration → CRISPOR specificity → PCR donor (45 bp HA + cassette)\n"
        "→ genotyping + junction primers → WT & INTEGRATED GenBank maps", C["design"], dsize=8.2)

    ax.add_patch(FancyArrowPatch((50, 84.4), (50, 81.1), arrowstyle="-|>", mutation_scale=18,
                                 lw=1.6, color=C["edge"], zorder=1))
    for y1, y2 in [(72.8, 69.4), (60.6, 56.4), (44.6, 39.4), (30.7, 25.9), (19.1, 14.9)]:
        arrow(y1, y2)

    ax.add_patch(FancyBboxPatch((84, 31), 30, 12, boxstyle="round,pad=0.4,rounding_size=1",
                                fc=C["tier"], ec=C["edge"], lw=1.1, ls="--", zorder=2))
    ax.text(99, 41.4, "three candidate tiers", ha="center", fontsize=8.6, fontweight="bold")
    ax.text(99, 38.6, f"strict (gene-distal)  {N['strict']:>5}\n"
                      f"relaxed (±2 kb)       {N['relaxed']:>5}\n"
                      f"intron-permissive    {N['intronic']:>5}",
            ha="center", va="center", fontsize=8, family="monospace")
    ax.add_patch(FancyArrowPatch((81, 35), (84, 36), arrowstyle="-|>", mutation_scale=13,
                                 lw=1.3, color=C["edge"], ls="--", zorder=1))

    ax.add_patch(FancyBboxPatch((18, 1.2), 64, 3.4, boxstyle="round,pad=0.3,rounding_size=1",
                                fc=C["out"], ec="none", zorder=2))
    ax.text(50, 2.9, "OUTPUT:  ranked, inert, constitutively-open safe-harbor loci  +  ready-to-order knock-in reagents",
            ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.add_patch(FancyArrowPatch((50, 5.2), (50, 4.7), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.4, color=C["edge"], zorder=1))
    save(fig, "pipeline_schematic")


# ============================ talk (landscape) ============================
def render_talk():
    fig, ax = plt.subplots(figsize=(16, 6.5))
    ax.set_xlim(0, 160); ax.set_ylim(0, 65); ax.axis("off")

    ax.text(80, 62, "Finding safe-harbor loci in the zebrafish genome",
            ha="center", va="top", fontsize=19, fontweight="bold")
    # four selection properties as a banner
    ax.text(80, 55.5, "constitutively OPEN   →   INERT   →   gene-distal & SAFE   →   uniquely TARGETABLE",
            ha="center", va="top", fontsize=12, color="#2E5D34", fontweight="bold")

    steps = [
        ("DANIO-CODE\natlas", "ATAC-seq (11 stages)\n+ ChIP-seq marks", C["input"]),
        ("Open in\nevery stage", f"{N['consensus']} constitutively\naccessible windows", C["access"]),
        ("Strip out the\nnon-inert & unsafe", "genes · enhancers/promoters\nPolycomb · conserved (CNE)\nnon-unique · near ncRNA", C["excl"]),
        ("Rank &\nvalidate", f"{N['strict']} safe-harbor loci\n(gene-distal tier)", C["rank"]),
        ("Knock-in\nreagents", "guides · CRISPOR · donors\ngenotyping · GenBank maps", C["design"]),
    ]
    n = len(steps)
    bw, bh, gap = 26, 24, 5.0
    total = n*bw + (n-1)*gap
    x0 = (160-total)/2 + bw/2
    centers = [x0 + i*(bw+gap) for i in range(n)]
    yc = 30
    for (title, detail, color), xc in zip(steps, centers):
        ax.add_patch(FancyBboxPatch((xc-bw/2, yc-bh/2), bw, bh, boxstyle="round,pad=0.8,rounding_size=2.5",
                                    fc=color, ec=C["edge"], lw=1.6, zorder=2))
        ax.text(xc, yc+bh/2-3.2, title, ha="center", va="top", fontsize=12.5, fontweight="bold", zorder=3)
        ax.text(xc, yc-1.5, detail, ha="center", va="center", fontsize=9.2, zorder=3)
    for i in range(n-1):
        ax.add_patch(FancyArrowPatch((centers[i]+bw/2+0.5, yc), (centers[i+1]-bw/2-0.5, yc),
                                     arrowstyle="-|>", mutation_scale=26, lw=2.4, color=C["edge"], zorder=1))

    ax.text(80, 9, f"alternate tiers:  strict {N['strict']}  ·  relaxed {N['relaxed']}  ·  intron-permissive {N['intronic']}",
            ha="center", fontsize=10, style="italic", color="#555")
    ax.text(80, 5, "GRCz11 / danRer11", ha="center", fontsize=9, color="#888")
    save(fig, "pipeline_schematic_talk")


render_detailed()
render_talk()
print(f"counts — consensus {N['consensus']}, strict {N['strict']}, "
      f"relaxed {N['relaxed']}, intron-permissive {N['intronic']}")
print("wrote results/pipeline_schematic{,_talk}.{png,svg}")
