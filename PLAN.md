Full text would help tighten the gene-target list and supplementary tables — send it if you have it. Meanwhile, here's a locked-down protocol for Track 1, built the same way as the glioma document: pre-specified, risk-audited, no fishing expeditions. I verified the actual source papers and accessions rather than guessing.

---

## Confirmed anchor papers (verified just now)

- **CCR 2023** — Jungwirth et al., *Pharmacological Landscape of FDA-Approved Anticancer Drugs Reveals Sensitivities to Ixabepilone, Romidepsin, Omacetaxine, and Carfilzomib in Aggressive Meningiomas*, Clin Cancer Res 29(1):233–243, DOI 10.1158/1078-0432.CCR-22-2085. Screened 119 FDA-approved drugs on meningioma cell lines; the four most effective compounds by IC50 (0.12–9.5 nmol/L) were carfilzomib, omacetaxine, ixabepilone, and romidepsin, mostly causing G2–M arrest and apoptosis.
- **Sci Transl Med 2026** — Jungwirth et al., *Drug screening on tumor organoids exposes therapeutic vulnerabilities of meningiomas to HDAC1/2i panobinostat*, DOI 10.1126/scitranslmed.aea3115. This built a patient-derived tumor organoid model from 60 molecularly characterized meningioma cases; panobinostat (pan-HDAC1/2i) was effective in vitro/in vivo/ex vivo, and resistance was traced to an HDAC8–TGFβ–EMT axis, with HDAC8 depletion restoring sensitivity.
- **Nassiri et al. 2021, Nature** — molecular classification (immunogenic / NF2-inactivated-canonical / hypermetabolic-ish groupings — exact naming needs the paper). Data: methylation idat files at GEO GSE180061; WES/bulk mRNA/snRNA at EGA under study EGAS00001004982; processed data at cBioPortal (mng_utoronto_2021).
- **Sahm et al. 2017, Lancet Oncol** — 497 meningiomas, six methylation classes, superior to WHO grade for outcome prediction. **Accession unconfirmed** — Heidelberg cohorts are frequently controlled-access (patient consent restrictions), not open GEO. This is a real risk item, flagged below, not assumed away.
- **Choudhury/Bi lab (Cancer Cell) three-group epigenetic classification** — 565 meningiomas profiled by methylation integrated with genetic, transcriptomic, proteomic, and single-cell data, defining Merlin-intact (best outcome), immune-enriched (intermediate), and hypermitotic (worst outcome) groups, deposited as GSE212666. This is actually a stronger public backbone than Sahm for your purposes — it's open, protein-layer-aware already, and framed around "therapeutic vulnerabilities," which mirrors your gap directly.

This changes the plan slightly from what you sketched: **use GSE180061/cBioPortal (Nassiri) + GSE212666 (Bi lab three-group) as the two independent classification systems**, rather than Sahm as primary — cross-classifying against two systems is also a stronger "convergence" story, same logic as your glioma paper's Ivy GAP + Neftel design.

---

## 1. Background & Rationale

Two independent Heidelberg screens (CCR 2023, cell-line based; StM 2026, organoid-based) identified five compounds with strong antimeningioma activity — carfilzomib, omacetaxine, ixabepilone, romidepsin, panobinostat — spanning proteasome inhibition, translation inhibition, tubulin inhibition, and HDAC1/2 inhibition. Neither paper stratifies *which* patients are most likely to benefit by validated molecular subgroup. This is the same design as most functional-screening papers: real biology, no biomarker layer. Your gap is to build that layer in silico, using public multi-omics meningioma cohorts, and hand it back as a testable hypothesis for patient selection.

**Core novelty claim:** convergence of two independently-derived molecular classification systems (Nassiri methylation/multi-omics groups; Bi lab three-tier epigenetic groups) on differential expression of the five drug-target gene programs, plus an NF2-status axis, cross-validated across ≥2 additional public cohorts — not a single-cohort correlative reanalysis.

---

## 2. Central Hypotheses

**H1 (concordance):** Nassiri molecular groups and Bi lab (Merlin-intact / immune-enriched / hypermitotic) groups are significantly concordant when both labels are available on overlapping/comparable public samples — establishes that you're measuring one underlying biology from two classifiers, not noise.

**H2 (target expression stratification):** Expression/activity of the five drug-target gene programs — HDAC1/2 (panobinostat), PSMB5/proteasome subunits (carfilzomib), tubulin genes (ixabepilone), HDAC1 (romidepsin — it's also an HDAC inhibitor, note overlap with panobinostat mechanistically), ribosomal/protein-synthesis machinery (omacetaxine) — differ significantly across molecular subgroups, independent of WHO grade.

**H3 (NF2 axis):** Target-gene expression is significantly associated with NF2 mutation/loss status, since NF2-inactivated tumors are the largest, best-characterized, most clinically actionable subgroup.

**H4 (external replication):** The subgroup–target associations found in the discovery cohort replicate in direction and effect size in an independent public cohort.

**H5 (translational framing, exploratory only):** If outcome data exist in any cohort, target-gene-high vs. target-gene-low patients differ in recurrence-free survival — explicitly exploratory, hypothesis-generating, not a claim requiring causal language.

---

## 3. Datasets (all public — to verify/lock before writing one line of code)

| Dataset | Role | Access | Status |
|---|---|---|---|
| Nassiri 2021 multi-omics (methylation + WES + mRNA + snRNA) | Discovery cohort 1, classification labels | GEO GSE180061 (methylation); EGA EGAS00001004982 (WES/mRNA/snRNA, controlled access — **flag**); cBioPortal `mng_utoronto_2021` (processed, likely easiest entry point) | GEO/cBioPortal open; EGA raw sequencing is controlled-access, apply if needed |
| Bi lab three-group epigenetic classification | Discovery cohort 2, classification labels | GEO GSE212666 | Open |
| CGGA / additional GEO meningioma expression sets | External replication (H4) | To be screened — need 2–3 candidate GSE IDs with matched expression + NF2/grade metadata | **Not yet confirmed — screening task, see Aim 0 below** |
| cBioPortal meningioma studies (mutation/CNV) | NF2 status cross-check | cbioportal.org, multiple meningioma studies listed under "Nervous System" | Open |
| CPTAC or any public meningioma proteomics | Optional protein-layer check | Meningioma proteomics is sparse in CPTAC (mostly glioma-focused) — **likely not available**; note as an honest limitation rather than force a weak Aim 4 | Needs confirmation, low expectation |

**Aim 0 (before anything else, ~3–5 days):** a dedicated GEO/cBioPortal screening pass to lock 2–4 exact accessions with sample-level NF2, grade, and (ideally) recurrence metadata. This replaces guessing accession numbers now with actually pulling the metadata tables and checking N and completeness — the single highest-yield hour of the whole project, and the thing most likely to silently sink the paper if skipped.

---

## 4. Pipeline (Aims)

- **Aim 0 — Cohort lock-in.** Screen candidate GEO/cBioPortal datasets, confirm N, platform, and metadata completeness (NF2, grade, recurrence if present). Freeze cohort list before any analysis — avoids "cohort shopping" that inflates significance.
- **Aim 1 — Classifier concordance (H1).** Where both Nassiri-group and Bi-lab-group labels are derivable/available on overlapping samples (or independently derived via published centroid/signature scoring on new samples), test concordance.
- **Aim 2 — Target-gene program scoring and subgroup association (H2, H3).** Score every sample for each of the five target-gene programs via ssGSEA/singscore; test association with molecular subgroup and NF2 status.
- **Aim 3 — External replication (H4).** Lock model spec in Aim 2, apply unchanged to replication cohort(s) from Aim 0.
- **Aim 4 — Exploratory survival (H5), only if recurrence-free survival data exist with adequate events.** Otherwise omit rather than force it — a paper that transparently says "outcome data were insufficiently powered for survival analysis" reads as more rigorous than a underpowered Cox model with 12 events.

---

## 5. Statistical Analysis Plan

**5.1 Gene programs (fixed before any testing, no data-driven selection):**
- Panobinostat/romidepsin: HDAC1, HDAC2 (+ HDAC8 for the resistance-axis angle from StM 2026 — this is actually a nice hook: score HDAC8/TGFβ/EMT-signature expression separately by subgroup, since that's literally the resistance mechanism they found)
- Carfilzomib: proteasome subunit genes (PSMB5, PSMB1, PSMA-family core panel — fixed list, not exploratory)
- Ixabepilone: tubulin genes (TUBB, TUBB3, class III β-tubulin especially, since it's a known taxane/epothilone-resistance marker — worth a citation-supported rationale)
- Omacetaxine: protein synthesis machinery (EEF2, ribosomal protein gene set) — homoharringtonine mechanism is translation elongation inhibition

**5.2 Subgroup association testing:**
- Kruskal-Wallis across subgroups per gene/program score (non-parametric default, no normality assumption on expression scores)
- Post-hoc pairwise Wilcoxon with BH-FDR correction across all program×subgroup tests
- NF2 status: Mann-Whitney U per program score, NF2-mutant/loss vs. NF2-intact
- Adjust/stratify by WHO grade throughout — grade is the obvious confound (higher grade tumors likely upregulate both target genes and get flagged as "worse subgroup" for reasons unrelated to your hypothesis)

**5.3 Concordance (H1):** Cohen's kappa or Cramér's V between classification systems where labels overlap; report as a supplementary table, exactly like your glioma protocol's H1 treatment.

**5.4 External replication (H4):** identical fixed model from Aim 2, no re-tuning; report effect size + direction concordance, not just p-values — same replication-reporting philosophy as your glioma design.

**5.5 Multiple testing:** BH-FDR across the full family of program×subgroup×cohort tests, pre-specified before running anything, not applied post hoc to a subset that "worked."

---

## 6. Pre-Specified Risk Register

| Risk | Why it matters | Mitigation |
|---|---|---|
| Sahm 2017 raw data likely controlled-access | You cannot promise a dataset you can't actually download | Replace with GSE212666 (Bi lab) as primary; treat Sahm as citation/context only unless individual-level data confirmed accessible |
| Grade confound | Target-gene expression may just track grade, not "true" subgroup biology | Grade as mandatory covariate/stratification in every test, exactly as in your glioma H5 design |
| Small N in any true "matched multi-omics + outcome" subset | Meningioma public cohorts are much smaller than glioma (TCGA-scale) | Pre-specify fixed gene panels (no unconstrained feature search); treat survival analysis as fully optional/exploratory pending event count |
| HDAC1/2 target overlaps between panobinostat and romidepsin | Could look like double-counting one mechanism as two "hits" | Explicitly note mechanistic overlap in the paper; treat as one HDAC axis + report both drugs' clinical distinctness (pan- vs. class-I selective) rather than presenting as independent evidence |
| Reviewer fatigue: "yet another public-data reanalysis" | Real risk at Q1-adjacent venues | Two-classifier convergence (H1) + explicit tie to StM 2026's own resistance mechanism (HDAC8/TGFβ/EMT) as a testable follow-up hook — gives Warta's group a concrete reason to want to co-author, not just be cited |
| Timing | StM 2026 is fresh (March 2026); window for "extension of a paper everyone's still discussing" narrows over months | Move Aim 0 and Aim 1 within the next 2–3 weeks |

---

## 7. Deliverables Checklist

- [ ] Aim 0 cohort lock memo: confirmed accessions, N, metadata completeness
- [ ] Frozen five-gene-program list with rationale (mechanism-linked, not exploratory)
- [ ] Concordance table (H1)
- [ ] Subgroup × target-gene association table, grade-adjusted (H2)
- [ ] NF2-status association result (H3)
- [ ] Replication table from independent cohort (H4)
- [ ] Optional exploratory survival result, clearly labeled (H5) — or an honest one-line statement it wasn't powered
- [ ] 1-page summary + 1–2 figures ready *before* contacting Warta — same "don't pitch, deliver" logic you already planned

---

Send the StM 2026 and CCR 2023 supplementary tables if you can — they likely contain the exact z-AUC/IC50 gene-target panels used in screening, which would let me replace my mechanism-iğğnferred gene lists (HDAC1/2, PSMB5, TUBB3, EEF2) with the authors' own target-annotation, which is strictly better for a paper you're pitching back to them.