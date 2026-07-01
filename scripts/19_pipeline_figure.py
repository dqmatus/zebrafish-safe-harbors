#!/usr/bin/env python3
"""
19_pipeline_figure.py
Render a schematic of the zebrafish safe-harbor discovery pipeline.
Outputs results/pipeline_schematic.png (+ .svg).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# palette
C = dict(input="#E8EAED", access="#CFE3F7", excl="#FBE0D0", rank="#D6ECD3",
         design="#E3D7F0", tier="#FFF4CC", edge="#4A4A4A")

fig, ax = plt.subplots(figsize=(11, 14))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, title, detail, color, tsize=11, dsize=8.5):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.6,rounding_size=1.5",
                                fc=color, ec=C["edge"], lw=1.3, zorder=2))
    ax.text(x, y+h/2-2.0, title, ha="center", va="top", fontsize=tsize, fontweight="bold", zorder=3)
    if detail:
        ax.text(x, y+h/2-5.2, detail, ha="center", va="top", fontsize=dsize, zorder=3, wrap=True)
    return (x, y)

def arrow(p1, p2, label=None):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=18,
                                 lw=1.6, color=C["edge"], shrinkA=2, shrinkB=2, zorder=1))
    if label:
        ax.text((p1[0]+p2[0])/2+13, (p1[1]+p2[1])/2, label, ha="left", va="center",
                fontsize=8.5, style="italic", color="#8B0000")

ax.text(50, 98.5, "Zebrafish genomic safe-harbor discovery pipeline",
        ha="center", va="top", fontsize=15, fontweight="bold")
ax.text(50, 95.5, "GRCz11 / danRer11  ·  DANIO-CODE developmental atlas", ha="center", va="top",
        fontsize=9.5, color="#555")

# ---- inputs row ----
inputs = ["Genome\nGRCz11", "DANIO-CODE\nATAC-seq\n(11 stages)", "ChIP-seq\nH3K27ac/K4me1/\nK4me3/CTCF/…",
          "Ensembl\ngenes · exons\nmiRNA · lncRNA", "Conservation\nGERP + CNEs", "COPE/DOPE\n+ validated CREs"]
xs = [10, 26, 42, 58, 74, 90]
for x, t in zip(xs, inputs):
    ax.add_patch(FancyBboxPatch((x-7.5, 86.5), 15, 6.5, boxstyle="round,pad=0.4,rounding_size=1",
                                fc=C["input"], ec=C["edge"], lw=1.0, zorder=2))
    ax.text(x, 89.7, t, ha="center", va="center", fontsize=7.3, zorder=3)
ax.text(50, 84.5, "inputs", ha="center", fontsize=8, style="italic", color="#555")

# ---- main vertical flow ----
p2 = box(50, 77, 62, 8, "① Constitutive accessibility  (06–07)",
         "Per-stage ATAC summits → collapse replicate→biosample→stage →\n"
         "peaks open in 100% of 11 stages   →   32,601 candidate windows",
         C["access"])
p3a = box(50, 65, 62, 8.5, "② Genomic-context exclusion  (05, 08)",
          "Remove gene bodies ±10 kb · exons · essential genes ±20 kb ·\n"
          "transposons · telomeres · non-primary scaffolds · <300 bp",
          C["excl"])
p3b = box(50, 50.5, 68, 11.5, "③ Regulatory / inertness screen  (08b — discard-whole)",
          "Drop any peak overlapping:\n"
          "• active chromatin  H3K27ac · H3K4me1 · H3K4me3 · CTCF\n"
          "• Polycomb H3K27me3 · other ChIP (H2A.Z, PolII…) · validated CREs\n"
          "• conservation (GERP + CNEs) · low mappability (genmap K50)\n"
          "• buffered miRNA ±50 kb · lncRNA ±25 kb",
          C["excl"], dsize=8.2)
p4 = box(50, 35, 62, 8.5, "④ Annotate & rank  (09)",
         "Nearest gene · GC · conservation · sgRNA · TAD · stage count ·\n"
         "COPE/DOPE flag   →   composite rank",
         C["rank"])
p5 = box(50, 22.5, 62, 6.5, "⑤ Validation benchmarks  (10)",
         "pIGLET14a/24b · tpl102 (exon-proximal → INCONCLUSIVE by design)",
         C["rank"], dsize=8.2)
p6 = box(50, 10, 68, 9.5, "⑥ Knock-in reagent design  (11–18)",
         "sgRNA enumeration → CRISPOR specificity → PCR donor (45 bp HA + cassette)\n"
         "→ genotyping + junction primers → WT & INTEGRATED GenBank maps",
         C["design"], dsize=8.2)

# input arrow into phase 1
ax.add_patch(FancyArrowPatch((50, 84.4), (50, 81.1), arrowstyle="-|>", mutation_scale=18,
                             lw=1.6, color=C["edge"], zorder=1))
arrow((50, 72.8), (50, 69.4))
arrow((50, 60.6), (50, 56.4))
arrow((50, 44.6), (50, 39.4))
arrow((50, 30.7), (50, 25.9))
arrow((50, 19.1), (50, 14.9))

# ---- tiers callout (right of phase 4) ----
ax.add_patch(FancyBboxPatch((84, 31), 30, 12, boxstyle="round,pad=0.4,rounding_size=1",
                            fc=C["tier"], ec=C["edge"], lw=1.1, ls="--", zorder=2))
ax.text(99, 41.4, "three candidate tiers", ha="center", fontsize=8.6, fontweight="bold")
ax.text(99, 38.6, "strict (gene-distal)   1,693\n"
                  "relaxed (±2 kb)        2,090\n"
                  "intron-permissive     4,131", ha="center", va="center", fontsize=8, family="monospace")
ax.add_patch(FancyArrowPatch((81, 35), (84, 36), arrowstyle="-|>", mutation_scale=13,
                             lw=1.3, color=C["edge"], ls="--", zorder=1))

# output banner
ax.add_patch(FancyBboxPatch((18, 1.2), 64, 3.4, boxstyle="round,pad=0.3,rounding_size=1",
                            fc="#2E5D34", ec="none", zorder=2))
ax.text(50, 2.9, "OUTPUT:  ranked, inert, constitutively-open safe-harbor loci  +  ready-to-order knock-in reagents",
        ha="center", va="center", fontsize=9, color="white", fontweight="bold")
ax.add_patch(FancyArrowPatch((50, 5.2), (50, 4.7), arrowstyle="-|>", mutation_scale=14,
                             lw=1.4, color=C["edge"], zorder=1))

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
for ext in ("png", "svg"):
    fig.savefig(f"results/pipeline_schematic.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
print("wrote results/pipeline_schematic.png + .svg")
