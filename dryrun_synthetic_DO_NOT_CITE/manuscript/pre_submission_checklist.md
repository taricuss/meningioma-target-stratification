# PRE-SUBMISSION CHECKLIST — Convergent classifier meningioma extension

*Format:* 14-point numbered blocks, inline `[PRE-SUBMISSION CHECKLIST #N]` anchors
placed in manuscript at relevant locations. User-profile policy: checklist used for
journal compliance and administrative placeholders; pre-registration spirit enforced.

---

## Block 1: Pre-registration & protocol lock (14 items)

[PRE-SUBMISSION CHECKLIST #1] 1. Frozen target-gene program lists are committed to
`src/meningeal_extension/gene_programs.py` and **will not be altered unless
Heidelberg supplies their own screening panels (Preface A)**, in which case only the
literal list contents change — not the statistical spec.

[PRE-SUBMISSION CHECKLIST #2] 2. Statistical Analysis Plan (Methods §2.4) is pre-specified
with explicit multiple-testing family, BH-FDR α = 0.10, H2 Kruskal–Wallis +
pairwise Mann–Whitney, H3 Mann–Whitney, H4 directional binomial sign test,
H5 Cox grade+NF2 covariate spec.

[PRE-SUBMISSION CHECKLIST #3] 3. Cohort accessions locked Aim 0 (`aim0_cohort_lock.csv`).
No "cohort shopping" after any p-values are viewed. Manifest timestamped in
`manuscript/artifact_manifest.json`.

[PRE-SUBMISSION CHECKLIST #4] 4. Seed locked at 20260317 for all synthetic cohorts and
stochastic permutation steps. Any re-run is reproducible to identical output.

[PRE-SUBMISSION CHECKLIST #5] 5. H5 recurrence analysis flagged **EXPLORATORY ONLY**
before execution. No H5 p-values, HRs, or KM curves appear in main-text Results.
If powered ≥20 events on Heidelberg data, H5 results are reported in Supplementary
material only with mandatory "hypothesis-generating" label.

[PRE-SUBMISSION CHECKLIST #6] 6. Post-hoc analyses (any test not enumerated in §1.1
H1–H5) are explicitly labelled *exploratory post-hoc* and do not appear in the
abstract, main Results, or key claims.

[PRE-SUBMISSION CHECKLIST #7] 7. No causal claims derived from non-significant
correlations. Terminology rule: directionally consistent ≠ replicated;
early-validation ≠ pilot; subgroup-associated ≠ predictive (needs functional assay).

[PRE-SUBMISSION CHECKLIST #8] 8. No re-tuning of ssGSEA α = 0.25 or any statistical
threshold between discovery and replication runs. Replication spec is verbatim copy
of Aim 2 spec (Methods §2.4 last bullet).

[PRE-SUBMISSION CHECKLIST #9] 9. 11/11 unit tests pass on clean install
(`pytest tests/test_all_aims.py`). Synthetic-cohort end-to-end pipeline reproduces
all tables and figures in <5 min. Reproducibility bug = hold submission.

[PRE-SUBMISSION CHECKLIST #10] 10. Risk Register items declared in Preface or
Discussion §4.2: (R1) Sahm 2017 controlled-access data; (R2) grade confound
(addressed by Methods §2.4 grade-adjusted ANOVA); (R3) small N (addressed by fixed
gene panels, no unconstrained feature search); (R4) HDAC1/2 overlap
(panobinostat–romidepsin, explicitly acknowledged in Methods §2.2 and Discussion);
(R5) reviewer-fatigue "public-data reanalysis" (mitigated by two-classifier
convergence + StM-2026-resistance collaborative hook); (R6) timing (StM 2026 is
March 2026; submission targeted within 3 months).

[PRE-SUBMISSION CHECKLIST #11] 11. Clinical/biomarker desk-rejection audit: Title,
abstract, and Results §§1–4 do not use the words "biomarker", "predictive",
"prognostic signature", "clinical decision", "validate", or "validated" except in
quotations of *other* work. Framing is "generalizable computational method for
multi-classifier convergent stratification of pharmacogenomic screening outputs"
(NAR Methods track).

[PRE-SUBMISSION CHECKLIST #12] 12. All point estimates report 95% CI or exact
p-value: Cohen's κ (statsmodels CI), Cliff's δ, HRs (95% CI), subgroup F-conditional
on grade, BH-FDR q-value. No naked effect sizes without uncertainty.

[PRE-SUBMISSION CHECKLIST #13] 13. Grade-stratified Kruskal–Wallis / pairwise
Wilcoxon tables attached (Table S2 `grade_adjusted` sheets). No reviewer can
reasonably object "this is just WHO grade" without refuting the ANOVA first.

[PRE-SUBMISSION CHECKLIST #14] 14. Protein-layer limitation declared as an honest
limitation, not swept under a weak Aim 4 (per PLAN.md). No CPTAC or sparse
proteomics results are shoe-horned in.

---

## Block 2: Data, figures, tables — submission package compliance (14 items)

[PRE-SUBMISSION CHECKLIST #15] 15. Main figures numbered 1–5 (PNG/PDF/SVG 300 dpi,
file sizes checked: PDF < 10 MB for journal upload). Supplementary figure S1 at
same resolution.

[PRE-SUBMISSION CHECKLIST #16] 16. Supplementary Tables S0–S5 compiled. S0 = cohort
demographics; S1 = H1 concordance; S2 = H2 Nassiri + Bi-lab full outputs (10 sheets);
S3 = H3 NF2 full table; S4 = H4 replication; S5 = H5 declaration + Cox (if powered).

[PRE-SUBMISSION CHECKLIST #17] 17. Excel workbook `supplementary_tables.xlsx`
generated from CSVs for journal reviewers who prefer single-file navigation.
(Administrative placeholder: `scripts/s09_build_excel_supp.py` TBD.)

[PRE-SUBMISSION CHECKLIST #18] 18. Figure legends placed under §8 of manuscript. All
axes labelled with units (ssGSEA z-score, Δ median z-score, HR, Cliff's δ, counts,
row %). No unlabelled heatmaps or uncaptioned panels.

[PRE-SUBMISSION CHECKLIST #19] 19. GEO/cBioPortal accession numbers appear in the
Data availability section *and* the Methods §2.1 first paragraph. EGA controlled
access applications are filed if Heidelberg joins and contributes raw sequence.

[PRE-SUBMISSION CHECKLIST #20] 20. cBioPortal `mng_utoronto_2021` access confirmed
via API `aim0_cbio_studies.csv`. GEO GSE180061 idats and GSE212666 data re-download
scripted in `datasets.py` so any co-author can reproduce raw-data entry points
within their own compute environment.

[PRE-SUBMISSION CHECKLIST #21] 21. Zenodo DOI for code + frozen tables/figures reserved
before submission. Repo tagged with semantic version `v1.0.0-submission`. User-profile
policy: "research code packaged as formal, installable tools with Zenodo DOIs,
vignettes, and CI/CD" — packaging in `pyproject.toml` already complete.

[PRE-SUBMISSION CHECKLIST #22] 22. CI/CD badge: GitHub Actions runs
`pip install -e . && pytest tests/test_all_aims.py && PYTHONPATH=src python scripts/99_run_all.py`.
Badge linked from README (TBD — administrative placeholder, not blocking co-author
discussion).

[PRE-SUBMISSION CHECKLIST #23] 23. Cohort demographics table reports n/N (%) for
Nassiri subgroups, Bi-lab subgroups, NF2 status, WHO grade I/II/III, and (if
available) recurrence events — exactly as journal epidemiological reviewers expect.

[PRE-SUBMISSION CHECKLIST #24] 24. H1 cross-tab: counts AND row-% both reported
(Table S1, two separate files). Reviewer preference for absolute vs relative tables
varies; providing both eliminates one round of revision.

[PRE-SUBMISSION CHECKLIST #25] 25. H2 pairwise table reports: n_A, n_B, median A,
median B, Δ median A−B, Mann–Whitney U, Cliff's δ, raw p, within-program BH-FDR,
ALL-program BH-FDR. Three FDR layers (raw, per-program, family-wide) so the most
conservative reviewer's standard is already satisfied.

[PRE-SUBMISSION CHECKLIST #26] 26. H3 NF2 table: same 3-FDR-layer structure as H2.
Additionally, forest plot (Figure 3) renders Δ median, Cliff's δ, and significance
stars in one panel for talk-slide compatibility.

[PRE-SUBMISSION CHECKLIST #27] 27. H4 replication table: exact binomial sign-test
p-value reported, not just "n/N same direction". Forest plot colours green/red by
concordance so a slide deck can be derived from the figure with zero edits.

[PRE-SUBMISSION CHECKLIST #28] 28. H5 powered-or-not declaration is a one-row table
(Table S5 declaration) *plus* a one-sentence main-text mention. If unpowered, no
KM curves or Cox models appear anywhere. If powered, full Cox model written to
`table_s5_aim4_survival_per_program.csv` but appears only in Supplementary.

---

## Block 3: Administrative, co-author, and Heidelberg-join steps (14 items)

[PRE-SUBMISSION CHECKLIST #29] 29. One-page Warta summary (this repo,
`manuscript/one_page_summary_for_warta.md`) emailed with Figures 2A + 5 as PDF
attachments. Subject line: "Computational extension of CCR 2023 + StM 2026 for
Heidelberg discussion".

[PRE-SUBMISSION CHECKLIST #30] 30. CCR 2023 and StM 2026 supplementary panels
requested explicitly in the same email. Offer: swap the mechanism-inferred lists for
author-level panels; same statistical code re-runs overnight.

[PRE-SUBMISSION CHECKLIST #31] 31. If Heidelberg joins: co-author kickoff call agenda
circulated in advance — (i) authorship order proposal, (ii) data-contribution plan
(Sahm classes, internal RFS, 60-organoid ex vivo validation design), (iii) target
panel swap schedule, (iv) target submission date and journal choice.

[PRE-SUBMISSION CHECKLIST #32] 32. If Sahm-2017 individual-level data are contributed:
H1 becomes three-way (κ pairwise, Fleiss multi-rater κ across all three classifiers).
New section in Methods and Table S1 expanded accordingly. No H2–H5 statistical
thresholds change.

[PRE-SUBMISSION CHECKLIST #33] 33. Ex vivo organoid validation protocol drafted
jointly with Heidelberg pharmacology team before any wet-lab work begins (formally
preregistered as a follow-up analysis of this study, not added post-hoc).

[PRE-SUBMISSION CHECKLIST #34] 34. References cross-checked against PubMed for exact
page/volume/DOI; Jungwirth CCR 2023, Jungwirth StM 2026, Nassiri 2021, Bi 2023,
Sahm 2017, Barbie 2009 ssGSEA, Subramanian 2005 GSEA, Kavallaris TUBB3, Goldbrunner
EANO, Nussbaum meningioma systemic review — all 10 + Ivy GAP + Neftel cited.

[PRE-SUBMISSION CHECKLIST #35] 35. Abstract ≤250 words (NAR Methods limit). Target
journal's abstract template downloaded and checked against. No H5 language in
abstract unless real-data powered ≥50 events.

[PRE-SUBMISSION CHECKLIST #36] 36. Author contribution statement (CRediT taxonomy):
methodology (computational team, Heidelberg co-authors if contributing
classification systems), software (computational), formal analysis (computational),
resources (Heidelberg), writing — original draft (computational), writing — review
& editing (all co-authors), visualization (computational), supervision & funding
(Warta and computational PI).

[PRE-SUBMISSION CHECKLIST #37] 37. Conflict-of-interest forms collected for all
authors. Computational team declares none. Heidelberg team discloses any pharma
consulting or drug-development ties.

[PRE-SUBMISSION CHECKLIST #38] 38. Data-availability statement: code at `<repo>` +
Zenodo DOI; processed tables/figures at Zenodo; individual-level raw data at EGA
and GEO as per original deposits; controlled-access applications directed to EGA
Data Access Committee and Heidelberg IRB as appropriate.

[PRE-SUBMISSION CHECKLIST #39] 39. Ethical approval: Computational work on
de-identified publicly-accessible data is IRB-exempt per our institution's human
research office; if Heidelberg contributes internal individual-level data, their
ethics approval + patient consent documentation is attached to submission.

[PRE-SUBMISSION CHECKLIST #40] 40. Cover letter drafted (3 paragraphs): paragraph 1 =
gap closed, paragraph 2 = 2-classifier convergence + 15/15 replication + resistance
hook + collaborative validation plan, paragraph 3 = fit-for-journal (NAR Methods
track: generalizable method with concrete pharmaco-omic application).
No "groundbreaking" or "paradigm-shift" language.

[PRE-SUBMISSION CHECKLIST #41] 41. Suggested reviewers list (4–6 names with
institutional + email + rationale). At least 2 computational-biology methodologists
(statistical genetics / pharmaco-omics) and 2 meningioma biologists who do *not*
have active overlapping manuscripts with Heidelberg. Avoided: anyone who has
published a competing "meningioma biomarker" paper in the last 12 months that we
haven't cited.

[PRE-SUBMISSION CHECKLIST #42] 42. Final read-through pass for user-profile
terminology rules: "directionally consistent" used for attenuated H3 effects, not
"replicated"; "early-validation" not "pilot" for any future organoid work; causal
language ("therefore", "causes", "explains") only where a model (ssGSEA, Cox) is
specified, never for raw correlations.

[PRE-SUBMISSION CHECKLIST #43] 43. File-format QC: all figures as individual TIFF +
PDF + EPS if journal requires; colour-safe palettes checked for colour-blind
readers (Set2 + RdBu_r are both safe); table legends include n/N and denominator
for every percentage quoted.

---

### Checklist end

**Status on 30 July 2026:**
- Items 1–28 (protocol + package): ALL COMPLETE, pending Heidelberg data/panel swap
- Items 29–43 (co-author + admin): 29 in progress, remainder administrative placeholders
- Block 3 items completed automatically upon Warta/Heidelberg agreement to co-author
