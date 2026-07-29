# Adapting the safe-harbor knock-in pipeline to KNIT (nickase) editing

Drop-in guide for running **KNIT editing** (Gao et al., *Nature* 2026; nCas9-D10A–mSA + biotinylated donor; doi:10.1038/s41586-026-10819-7) on the reagents this repo already designs. Scope: **large-cassette knock-in at inert safe harbors, delivered as RNP** (this lab's current format), to enable **direct F0 integration of actionable cargo — skipping the landing-pad → F1 → RMCE loop.**

## TL;DR — this pipeline is already ~90% KNIT-ready

KNIT changes **two things** and leaves the rest of the pipeline untouched:

1. **Editor:** swap SpCas9 protein → **nCas9(D10A)–mSA fusion protein** in the RNP. Same delivery format you already use.
2. **Donor 5′ end:** re-order each donor ultramer with **`/5Biosg/` (5′ biotin) in place of `/5AmMC6/`**. Biotin is what the fused mSA grabs to tether the donor at the nick. Everything else about the donor is identical.

Your **already-ORDERED `/5AmMC6/` donors (DM418–423) still run your current Cas9-RNP protocol unchanged** — so the biotin variants are an incremental add for a KNIT arm run *in parallel*, not a re-spend. In fact that gives you a clean, matched head-to-head (same guide, same cassette, same locus) between your DSB protocol and KNIT.

## What stays identical (no redesign)

| Component | Repo artifact | Why it carries over |
|---|---|---|
| Safe-harbor site selection | scripts 01–14, `candidate_feature_table*.tsv` | KNIT doesn't change locus choice. |
| Single SpCas9 NGG guides | `locus_design_summary.tsv`, `sgrna_shortlist.tsv`, CRISPOR-scored | **A nickase uses the identical sgRNA** — nCas9(D10A) nicks instead of double-cuts at the same protospacer. No guide redesign. On-target (RuleSet3/Doench) scores remain a fine proxy for nickase engagement. |
| 45 bp homology arms | `knockin_donor_design.tsv` (`HA5_45nt`/`HA3_45nt`) | KNIT validated arms down to **25 bp**; your 45 bp arms are comfortably in range (same short-homology regime as GeneWeld/pIGLET). |
| StayGold cassette + universal priming | DM416/DM417, donor ultramers | Cassette and PCR-amplification scheme unchanged. |
| Genotyping / junction primers | DM426–435, DM442–447, `genotyping_primers.tsv`, `junction_primers.tsv` | Integration junctions are the same, so 5′/3′ junction PCR + StayGold junction primers read out KNIT integrants unchanged. |

## Change 1 — the editor (RNP route)

You need **purified nCas9(D10A)–mSA–His protein** instead of your current SpCas9. The paper cloned **SpCas9D10A-NLS-mSA-NLS-His in pET28b** and purified it by Ni-NTA — a standard prep. Acquisition options, easiest first:

- **Request the protein or the pET28b expression plasmid from the Wang lab** (KE plasmids are not on Addgene; data page: https://www.wanglab-thu.org.cn/knit).
- **Rebuild:** install D10A in an SpCas9-NLS-His backbone and fuse mSA (the paper's fusion is `SpCas9(D10A)-NLS-mSA-NLS-His`); express/purify in-house or via a contract service.
- Note: commercial "Cas9 D10A nickase" protein exists, but **without the mSA fusion it will not tether the donor** — the fusion is essential. (KE2/SunTag is a later efficiency upgrade; start with KE1 = the single mSA fusion.)

Assemble exactly as your current RNP, then **pre-mix the biotinylated donor with the assembled RNP before injection** so the mSA–biotin tether pre-forms — the configuration KNIT relies on and the one most likely to matter in the rapid 1-cell embryo.

## Change 2 — the donor (biotin instead of amine)

Re-order each donor ultramer identically except the 5′ modifier: **`/5Biosg/` replaces `/5AmMC6/`**; keep the 3× 5′ phosphorothioate and the 3′-P. Rationale: mSA binds biotin (it will not bind your current amine), and biotin is itself a bulky 5′ block — so you keep the end-protection/anti-concatemer benefit your `/5AmMC6/` was providing, plus bound mSA further caps the end.

Ready-to-order table: **`knit_biotin_donor_oligos.tsv`** (DM448–DM459 — biotin variants of DM418/419, 420/421, 422/423, 436/437, 438/439, 440/441; the superseded DM424/425 pIGLET24b_locus set is omitted).

**Biotin strand/geometry.** nCas9(**D10A**) leaves HNH active and nicks the **target (sgRNA-complementary) strand**; the productive tethered-end geometry for single-nick HDR is best settled empirically and it's cheap to do so. Two sensible starting choices:
- **Both ends biotinylated** (order both donor_F and donor_R biotin, per locus — as in the table): most robust, donor is tethered regardless of orientation. Recommended default for the first pilot.
- **Single-end variants** (biotin on F only, or R only): use these to *optimize* geometry once KNIT is working, and to probe which end drives productive repair. (mSA is monomeric = one biotin engaged per editor; KE2/SunTag can engage multiple.)

## Injection, controls, readout (fits your current workflow)

- **KNIT arm:** nCas9(D10A)–mSA RNP (your sgRNA) **pre-mixed with the biotinylated donor** (DM448+), injected 1-cell as usual. Include a **±AZD-7648 sub-arm** — but see the caution below; don't assume your DSB-protocol boost transfers.
- **Benchmark arm (the honest bar):** your **current SpCas9-RNP + `/5AmMC6/` donor (DM418–423) *with* AZD-7648** — your real, optimized protocol — at the *same* loci. This tells you whether KNIT's DSB-free accuracy/large-cargo advantage actually beats what you already do.
- **Readout:** StayGold expression (live imaging) + your existing 5′/3′ junction PCR (DM426–435 / DM442–447). Unchanged.

### Caution on AZD-7648
Your AZD boost comes from suppressing NHEJ at a **double-strand break**. KNIT makes a **single-strand nick**, and the Wang lab tested AZD-7648 (1 µM) specifically with KNIT and saw **no consistent improvement** at single nicks (Extended Data Figs. 1b, 5a,b). Keep a ±AZD arm because the rapidly-replicating embryo converts nicks to one-ended DSBs at forks (so it *might* help more in fish than in HEK293T), and if you want a pharmacological lever, the paper points to **AZD + a Polθ inhibitor (novobiocin/ART558)** rather than AZD alone. But the native KNIT efficiency lever is **KE2 (SunTag)**, not inhibitors.

## First experiment (minimal, decisive)

**Locus:** `chr12_nrxn1a` — your standout guide (MIT 100, 0 in-gene off-targets, central cut). **Editor:** nCas9(D10A)–mSA RNP + your guide. **Donor:** DM448/DM449 (StayGold, biotin, both ends), pre-mixed with the RNP. **Arms:** KNIT ±AZD vs your current SpCas9+AmC6 donor +AZD. **Score:** F0 StayGold + junction PCR; then raise founders.

De-risk option: if you want to test the tethered single-nick mechanism *before* embryos, reconstitute the same RNP + biotin donor in a zebrafish cell line (PAC2/ZF4) at 28 °C first — a one-week flow readout that de-risks the biggest unknown (does single-nick mSA-tethered HDR work at fish temperature at all).

## Strategic fit and one caveat from your own pipeline

KNIT's payoff for you is **direct F0 knock-in of large, actionable cassettes at inert harbors**, skipping landing-pad → RMCE. If you **keep the lox sites** on the cassette, each KNIT integrant also stays RMCE-compatible as a docking-site fallback. Caveat, straight from your README's own trade-off note: deeply inert sites (e.g. chr12/nrxn1a) may express *modestly*, while pIGLET24b expresses strongly *because* it sits in active chromatin — and **KNIT does not change that; expression level is locus-driven, not method-driven.** So pair site choice with your expression expectation, and keep **pIGLET24b_v2 (DM458/459) as the strong-expression positive control** anchoring each KNIT run.
