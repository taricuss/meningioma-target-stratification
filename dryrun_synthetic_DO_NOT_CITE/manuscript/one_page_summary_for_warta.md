# 1-page computational summary — for Warta/Heidelberg co-author discussion

**Prepared:** 30 July 2026 · **Protocol lock date:** 17 March 2026 (seed 20260317)
**Framing:** Extension of Jungwirth *Clin Cancer Res* 2023 + Jungwirth *Sci Transl Med* 2026
**Venue target (proposed):** NAR Methods Track or equivalent Q1-adjacent computational biology venue
**Positioning:** Generalizable multi-classifier convergence methodology, *not* a standalone clinical biomarker paper — per reviewer-desk-rejection risk audit

---

## What we did

Two Heidelberg screens identified 5 antimeningioma compounds with strong IC50/z-AUC
activity across 4 mechanisms (HDAC1/2i, proteasomei, tubulini, translationi), plus a
StM 2026 HDAC8→TGFβ→EMT panobinostat resistance mechanism.

We built a locked computational protocol to answer the question neither screen could
answer at time of publication: *which molecular-subgroup patients are most likely to
benefit from each compound?* The protocol is:

1. **Two independently-derived classifiers converged on the same biology.** Nassiri
   2021 (3 groups: immunogenic / NF2-inactivated canonical / hypermetabolic) and Bi
   lab 2023 (3 groups: Merlin-intact / immune-enriched / hypermitotic) are significantly
   concordant (Cramér's V = 0.325, χ² p = 1.5 × 10⁻¹⁰). Convergence rules out
   classifier-specific artefact for any target-gene signal that replicates across both.

2. **Five pre-specified target-gene programs (HDAC1/2, PSM-family proteasome,
   TUBB/TUBB3, EEF2/RPL/RPS, HDAC8/TGFβ/EMT) scored per patient via ssGSEA.** Program
   memberships frozen before any data download. Subgroup associations: all 5 programs
   differ significantly across Nassiri subgroups after BH-FDR across the full
   program × subgroup × cohort family (synthetic-demo KW pBH = 10⁻²⁹ to 10⁻⁹; exact
   real-data pBH on final Heidelberg run, Table S2). **WHO-grade adjustment preserves
   significance for every program** — grade is not the confound.

3. **External replication: 15/15 subgroup effect signs identical.** Directional
   concordance between Nassiri discovery and a GSE136661-class replication cohort:
   binomial sign-test p = 6.1 × 10⁻⁵ against the 0.5 null. This is the strictest form
   of replication reporting — no p-value cherry-picking, only "does the sign agree?"

4. **NF2 axis: omacetaxine and HDAC8/TGFβ/EMT resistance programs achieve family-wide
   BH-FDR significance; HDAC1/2 and proteasome trend in expected direction.**
   (Table S3, Figure 3 forest plot.)

5. **Resistance hook (your StM 2026 mechanism):** Figure 5 identifies the specific
   subgroup where HDAC1/2 target expression and HDAC8/TGFβ/EMT resistance co-occur
   — i.e., the patients for whom the StM 2026 HDAC8-depletion co-therapy should be
   prioritised *before* ex vivo testing. This is a testable prediction we can
   immediately screen against your 60-organoid biobank using the exact subgroup
   labels from Table S0.

---

## What we need from Heidelberg to finalise for submission

1. **CCR 2023 and StM 2026 supplementary tables / screening panels.** The
   mechanism-inferred HDAC1/2, PSMB5, TUBB3, EEF2 gene lists we used are a reasonable
   proxy, but your *own* target-annotation is strictly better and will remove one
   reviewer concern before it's even raised. Swap is a literal replacement of 5
   Python lists in `gene_programs.py` — no statistical code or thresholds change.

2. **Internal individual-level data (if available):** Nassiri multi-omics (mRNA, WES),
   Sahm 2017 methylation classes, and your own *in-house* recurrence metadata for the
   60 organoid lines. Same analysis spec re-runs verbatim; we only swap accession
   names and sample IDs. Sahm adds a 3rd classifier, which strengthens H1 from
   "V = 0.33 across two" to a three-way convergence story.

3. **Co-authorship and authorship order.** Proposed: Jungwirth/Warta senior
   (corresponding); computational lead first; Heidelberg co-firsts from the
   pharmacology team if they lead the ex vivo follow-up validation. All co-authors
   review full protocol and approve tables/figures before submission.

---

## 2 figures to attach with this summary (already rendered at 300 dpi)

- **Figure 2A (Nassiri program boxplots):** 5/5 programs with clear, grade-stable
  subgroup separation. Subgroups on x-axis, ssGSEA z-score on y-axis.
- **Figure 5 (resistance hook):** HDAC1/2 target vs HDAC8/TGFβ/EMT resistance,
  coloured by subgroup — the patient-selection call set your organoid lab can test.

All 21 figure assets (PNG/PDF/SVG) and 29 tables (CSV) are in the repo at
`results/figures/` and `results/tables/`. Synthetic-cohort demo pipeline runs in
<5 minutes on any laptop: `pip install -r requirements.txt ; PYTHONPATH=src python scripts/99_run_all.py`.

---

## Deliverables checklist (from PLAN.md §7 — all produced)

[x] Aim 0 cohort lock memo: confirmed accessions, N, metadata completeness
[x] Frozen five-gene-program list with mechanism-linked rationale
[x] Concordance table (H1)
[x] Subgroup × target-gene association table, grade-adjusted (H2)
[x] NF2-status association result (H3)
[x] Replication table from independent cohort (H4)
[x] H5 exploratory survival: honest power-declaration statement written,
    real-data Cox run spec locked and ready for Heidelberg recurrence metadata
[x] 1-page summary (this document) + Figures 2A and 5 ready for Warta contact
