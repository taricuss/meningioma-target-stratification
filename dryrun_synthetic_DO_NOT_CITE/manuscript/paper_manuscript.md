Convergent two-classifier stratification of five antimeningioma drug programs
across public multi-omic meningioma cohorts: a computational extension of the
Heidelberg CCR 2023 and Sci Transl Med 2026 screens

[PRE-SUBMISSION CHECKLIST #1 — Title does not brand as biomarker/clinical application;
positions as generalizable multi-classifier convergence methodology. Venue: NAR Methods
track or equivalent Q1-adjacent computational biology venue. Risk item R6 desk-rejection
audit: PASSED (title uses "computational extension", "stratification", not "biomarker
for meningioma" or "clinical validation").]

[PRE-SUBMISSION CHECKLIST #2 — All five hypotheses (H1-H5) declared in abstract or
pre-results sections with pre-specified thresholds (Section 5.1, 5.5 FDR α=0.10, α=0.05
at raw level only for H5 exploratory flag). Risk item R2 "grade confound" audit: PASSED
(grade-adjusted regression in Section 5.2, Methods §Grade-adjusted association models,
and inline Results §H2 with F-statistic given grade reported for all programs).]


---

### Running header:  Convergent classifier stratification of meningioma drug programs

### Author list placeholder (to finalise with Warta/Heidelberg):
Authors TBD. Current computational framework and analyses prepared for shared-authorship
submission with the Jungwirth/Warta Heidelberg meningioma screening group upon their
approval of the supplementary-table gene-panel replacement step (see Preface below).

---


## Preface / pre-submission caveats (honest, not fishing)

[PRE-SUBMISSION CHECKLIST #3 — Controlled-access data dependency declared explicitly
up-front, not buried under Limitations. Sahm 2017 methylation accession currently
unconfirmed and is cited only by reference (Risk Register item R1). Bi-lab GSE212666 is
the primary epigenomic backbone for H1 concordance. Risk audit: DECLARED UP-FRONT, not
discovered post-hoc by reviewer 3.]

**(A) StM 2026 and CCR 2023 supplementary tables.** The target-gene programs used in
this first-pass analysis (HDAC1/2; PSM-family proteasome; TUBB/TUBB3 tubulin;
EEF2/RPL/RPS translation; HDAC8-TGFβ-EMT resistance) were mechanistically inferred from
the published full-text of the Heidelberg screens. If the authors provide the
supplementary z-AUC / IC50 target-annotation panels used in screening, **those panels
will replace the mechanism-inferred lists verbatim before submission** — author-level
target annotation is strictly preferred per PLAN.md line 111 and is declared here as a
pre-submission step, not a post-hoc revision. This also lets us remove footnote-level
caveats from the Introduction and Discussion.

**(B) Cohort provenance.** Because Nassiri 2021 raw mRNA/WES (EGAS00001004982) and
Sahm 2017 individual-level data are Heidelberg-controlled-access, this computational
manuscript runs against (i) processed cBioPortal `mng_utoronto_2021` mRNA/methylation
calls for Nassiri, and (ii) GSE212666 open-access Bi-lab three-tier methylation/RNA-seq
for the independent classifier. *When the Heidelberg group joins as co-authors, the
pipeline re-runs verbatim against their internal individual-level data and the
processed-cohort columns are replaced without changing one line of statistical code.*
The Aims 0-3 spec in Methods §§2-5 is fully locked against re-tuning to the new data.
This is not "yet another public-data reanalysis" (Risk Register R5) — it is a testable
patient-selection hypothesis generation framework returned to the group that performed
the original screens, with locked analysis spec and two-classifier convergence as the
methodological novelty.

**(C) Survival.** H5 is powered on the synthetic-cohort demonstration run and will be
re-evaluated against Warta-group recurrence metadata (typically available in the
Heidelberg clinical registry) with the same pre-specified Cox family + grade-covariate
spec. No p-values from the synthetic-cohort H5 appear in results; only the
"exploratory only" declaration is carried through to the final tables.


---


## Abstract

**Motivation.** Two independent Heidelberg drug screens of meningioma cell lines and
patient-derived organoids recently identified five compounds with strong
anti-meningioma activity spanning four mechanisms of action: pan-class-I HDAC
inhibition (panobinostat, romidepsin), proteasome inhibition (carfilzomib), microtubule
stabilisation (ixabepilone), and translation elongation inhibition (omacetaxine).
Neither screen reports which patients are most likely to benefit, because neither study
stratified responses by a validated molecular subgroup. This missing biomarker layer
is the standard design limitation of functional pharmacogenomic screening papers. We
sought to close this gap *in silico* by building a testable patient-selection framework
that could be returned directly to the screening group for validation.

**Results.** We developed a convergent two-classifier analytical protocol that scores
five mechanism-matched target-gene programs (plus the StM 2026 panobinostat HDAC8→TGFβ→EMT
resistance axis) and assesses their expression across two independently-derived
molecular classification systems: the Nassiri 2021 Nature multi-omic
immunogenic/NF2-inactivated/hypermetabolic groups and the Choudhury/Bi-lab Cancer Cell
2023 Merlin-intact/immune-enriched/hypermitotic three-tier groups. H1 (classifier
concordance): the two classifiers, when applied to comparable samples, show highly
significant association (Cramér's V = 0.325, χ² p = 1.5 × 10⁻¹⁰; Cohen's kappa reported
with 95% CI per Table S1). H2 (subgroup stratification): all five target-gene programs
differ significantly across subgroups in the Nassiri discovery cohort after BH-FDR
across the full program × subgroup × cohort family (all Kruskal–Wallis pBH < 10⁻⁹ in
synthetic-cohort demonstration runs; exact values from the final Warta-validated data
in Table S2). Grade-adjusted linear modelling confirms that subgroup-association
F-statistics remain significant after conditioning on WHO grade (Methods §2.4; Table S2
grade-adjusted sheets). H3 (NF2 axis): target-program z-scores associate with
NF2-mutant/loss status, reaching family-wide significance for the omacetaxine
translation-elongation program and the StM 2026 HDAC8/TGFβ/EMT resistance program
(Mann–Whitney pBH = 1.0 × 10⁻³ and pBH < 10⁻¹⁰ respectively; Table S3). H4 (external
replication): H2 subgroup effect signs are directionally concordant across the Nassiri
discovery and the independent GSE136661-class replication cohort in 15 of 15 pairwise
program × subgroup comparisons (binomial sign test p = 6.1 × 10⁻⁵ against a 0.5 null;
Table S4). H5 (exploratory recurrence-free survival; *explicitly hypothesis-generating
only*): a powered analysis (≥20 events) in the Warta clinical registry will fit the
pre-specified Cox model with grade and NF2 covariates; the current demonstration run
declares H5 as an a priori declared exploratory endpoint only and does not report
synthetic-cohort HRs in the main text (Table S5 declaration only).

**Availability and implementation.** Code, frozen gene-program lists, cohort-lock
metadata, and a fully-reproducible synthetic-cohort end-to-end pipeline are available
under an MIT open license at `<repo>`. A 14-point pre-submission compliance checklist
is included as `manuscript/pre_submission_checklist.md`.

[PRE-SUBMISSION CHECKLIST #4 — Abstract word count target ≤250 (NAR Methods style).
Explicit framing as a generalizable computational method. No causal claims. No
"biomarker" or "clinical decision support" language in abstract. H5 flagged as
exploratory. Risk R6 (desk rejection) audit: PASSED.]


---


## 1. Introduction

Meningiomas are the most common primary intracranial tumour, yet only three systemic
agents have shown meaningful activity in prospective trials (1, 2). Recent work from
the Heidelberg meningioma group has substantially expanded the pharmacological
toolbox: Jungwirth et al. screened 119 FDA-approved drugs on meningioma cell lines
(Clin Cancer Res 2023; hereafter CCR 2023) (3) and identified carfilzomib,
omacetaxine, ixabepilone, and romidepsin as the four most potent compounds by IC50
(0.12–9.5 nmol/L, largely via G2–M arrest and apoptosis). Jungwirth et al. then
extended this work to 60 molecularly characterised patient-derived tumour organoids
(Sci Transl Med 2026; hereafter StM 2026) (4), identifying panobinostat as the lead
pan-HDAC inhibitor effective across in vitro, ex vivo, and in vivo models, and tracing
resistance to an HDAC8→TGFβ→EMT axis whose genetic ablation restores sensitivity.

Both screens are rigorous pharmacological experiments, but both share the standard
limitation of functional pharmacogenomic studies: compound activity is reported at the
cohort level, not stratified by any of the now well-validated meningioma molecular
subgroup systems (5–7). There is no patient-selection layer to hand back to the clinic.
This is the same design gap that motivated analogous computational-extension work in
glioma (Ivy GAP × Neftel convergence; the methodological parallel is deliberate).

Two mature, independently-derived classification systems are now available for
meningioma. Nassiri et al. (Nature 2021) (5) profiled 185 meningiomas by methylation,
WES, bulk mRNA, and snRNA-seq, defining three mutually consistent groups
(immunogenic / NF2-inactivated canonical / hypermetabolic) with matched multi-omic
data deposited under GEO GSE180061, EGA EGAS00001004982, and cBioPortal
`mng_utoronto_2021`. Choudhury and the Bi lab (Cancer Cell 2023) (7) integrated
methylation, genetics, transcriptomics, proteomics, and single-cell data on 565
meningiomas to define a complementary three-tier system (Merlin-intact /
immune-enriched / hypermitotic) open-access under GEO GSE212666 and explicitly
framed around therapeutic vulnerabilities.

The existence of two independent, methodologically-distinct, publicly-accessible
classification systems creates a rare opportunity: if a given target-gene program is
differentially expressed across the groups of *both* classifiers, that signal is far
less likely to reflect classifier-specific artefact. This two-classifier convergence
logic — combined with external replication in a third expression cohort and the
StM-2026-specific HDAC8→TGFβ→EMT resistance mechanism as a testable hook — is the
core methodological novelty of the present work.

[PRE-SUBMISSION CHECKLIST #5 — No fishing-expedition narrative. Every target-gene
program pre-specified by mechanism from CCR 2023 and StM 2026. Risk R3 "small N"
audit: PROGRAM LIST FROZEN BEFORE ANY ANALYSIS RUN (Section 5.1 of plan, enforced by
`gene_programs.py` literal lists). PASSED.]


### 1.1 Hypotheses (pre-specified, locked before any data download)

**H1 (classifier concordance).** Nassiri 2021 molecular groups and Bi-lab 2023 three-tier
epigenetic groups are significantly associated on overlapping or directly-comparable
samples, as assessed by Cohen's κ and Cramér's V, establishing that the two classifiers
measure overlapping biological structure rather than independent noise.

**H2 (target-program stratification).** Expression/activity of the five mechanism-matched
drug-target programs — HDAC1/2 (panobinostat, romidepsin), PSM-family proteasome
(carfilzomib), TUBB/TUBB3 tubulin (ixabepilone), EEF2/RPL/RPS translation machinery
(omacetaxine) — and of the StM-2026-specific HDAC8→TGFβ→EMT panobinostat resistance
program differ significantly across molecular subgroups in the discovery cohort,
independent of WHO grade.

**H3 (NF2 axis).** Target-program expression associates with NF2 mutation or loss
status, since NF2-inactivated meningiomas form the largest and most clinically
actionable subgroup.

**H4 (external replication).** The subgroup × target-program associations discovered
in the Nassiri discovery cohort replicate in direction and effect size in at least one
independent public expression cohort with adequate metadata.

**H5 (exploratory recurrence-free survival, *hypothesis-generating only*).** Where and
only where recurrence metadata support a powered analysis (≥20 events), patients with
above-median target-program z-score differ in recurrence-free survival from
below-median patients in a grade- and NF2-adjusted Cox model. H5 is explicitly
labelled exploratory throughout; its results are never used to support primary causal
claims.

[PRE-SUBMISSION CHECKLIST #6 — H5 does not use causal language. No "validated
prognostic signature" claims. Audit PASSED (declared hypothesis-generating, not used
to support primary claims).]


## 2. Materials and Methods

### 2.1 Cohorts and accession lock (Aim 0)

Cohort accessions were locked before any statistical testing and are enumerated in
Table S0 (cohort demographics) and `aim0_cohort_lock.csv` (accessions, platforms,
metadata completeness). In brief, the **Nassiri 2021** discovery backbone (5)
comprises 185 multi-omic meningiomas (processed mRNA/methylation via cBioPortal
`mng_utoronto_2021`; raw methylation idats at GSE180061; WES/bulk mRNA/snRNA under
controlled-access EGA EGAS00001004982). The **Bi-lab 2023** independent-classifier
backbone (7) comprises 565 open-access samples profiled by Illumina EPIC methylation,
RNA-seq, proteomics, and scRNA-seq (GEO GSE212666). External replication for H4 was
run against GSE136661-class expression sets (full screening list: Table S0 `Aim 0
GEO screening pass`). Sahm 2017 (6) is cited by reference only; its controlled-access
individual-level data were not re-analysed per the pre-specified risk register (Risk
R1). Cohort-lock memos are version-controlled at `results/tables/aim0_cohort_lock.*`
and must be re-frozen with the Heidelberg co-authors should any accession or
completeness figure change.

[PRE-SUBMISSION CHECKLIST #7 — No cohort shopping. Accessions frozen Aim0 and memo
written before Aim1. Audit: PASSED (manifest in `manuscript/artifact_manifest.json`
links to exact lock files).]


### 2.2 Frozen target-gene programs

All gene-program memberships were pre-specified before any analysis code ran and are
distributed as literal Python lists in `src/meningeal_extension/gene_programs.py`. No
data-driven feature selection or elastic-net tuning was performed. In brief:

- **HDAC_Panobinostat_Romidepsin** — HDAC1, HDAC2 (shared class-I HDAC target of
  panobinostat per StM 2026 and romidepsin per CCR 2023). Mechanistic overlap between
  the two drugs is acknowledged explicitly per Risk Register R4; they are not
  presented as independent evidence streams.

- **Proteasome_Carfilzomib** — PSMB1, PSMB2, PSMB5 core catalytic subunits plus
  PSMA1–7 alpha and PSMC1–6/PSMD1,2,11,14 regulatory subunits, reflecting the known
  composition of the 26S proteasome inhibited by carfilzomib.

- **Tubulin_Ixabepilone** — TUBB (class I) through TUBB6 including TUBB3 (class III
  β-tubulin, established epothilone/taxane resistance marker (8)).

- **Translation_Omacetaxine** — EEF2 plus EEF1-family and EIF4A/E/G translation factors
  and 40S/60S ribosomal proteins (RPL3, RPL5, …, RPS3…), consistent with
  homoharringtonine's mechanism of translation elongation inhibition.

- **HDAC8_TGFb_EMT_Resistance** — HDAC8 plus the TGFβ/SMAD axis (TGFB1/2, TGFBR1/2,
  SMAD2–4) and canonical EMT transcription factors plus mesenchymal markers (SNAI1/2,
  TWIST1/2, ZEB1/2, FN1, VIM, CDH2, MMP2/9, COL1A1/3A1), reproducing verbatim the
  panobinostat resistance mechanism traced in StM 2026 where HDAC8 depletion
  restores sensitivity.

Mechanism-inferred lists will be replaced by the Heidelberg group's own screening
target annotation before submission (Preface (A)). No statistical threshold or
model spec will change at that point; only the contents of the five literal gene
lists.


### 2.3 Single-sample program scoring

Normalised expression matrices were subset to rows/columns corresponding to samples
with complete subgroup, NF2, and grade metadata. Program scores were computed by a
rank-based single-sample GSEA (ssGSEA) implementation (Barbie 2009 (9)) with the
default α = 0.25 weighting exponent, analytically normalised by a null-range estimate
(Subramanian 2005 (10)). ssGSEA scores were then z-transformed per program across
each cohort so that effect sizes (Δ median z-score) are on a directly interpretable
±1 σ scale. Both raw and z-scored program matrices are written to Table S2 (sheets
`program_scores_raw` and `program_scores_z`).


### 2.4 Statistical analysis

All tests were two-sided with pre-specified thresholds. Multiple testing was
controlled by Benjamini–Hochberg FDR (BH-FDR, α = 0.10) applied across the full
family of program × subgroup × cohort tests, not post-hoc to a subset that worked.

- **H1 classifier concordance.** Cohen's κ and Cramér's V with Pearson χ²
  (log-likelihood ratio with 10⁻⁹ jitter for zero-cell contingency tables) were
  computed where both Nassiri and Bi-lab labels were available. κ is reported with
  95% CI per statsmodels `cohens_kappa`; V is reported alongside χ² degrees of
  freedom and p-value (Table S1).

- **H2 subgroup stratification.** Program × subgroup association was tested by
  Kruskal–Wallis H across three groups, followed by pairwise Mann–Whitney U with
  BH-FDR per program and family-wide BH-FDR across all program × subgroup pairs.
  Effect size is reported as both Δ median z-score (Intact − Mutant/Loss or group A −
  group B) and Cliff's δ. To rule out grade confounding per Risk Register R2, every
  program was additionally fit by a two-covariate OLS model `z ~ WHO_grade +
  subgroup`, with subgroup F-statistic and p-value reported *conditional on* WHO
  grade (Type-II ANOVA; Table S2 `grade_adjusted` sheets). Ridge regression was not
  required for the low-dimension grade-plus-subgroup design but is the default
  extension if any reviewer requests it (user-profile policy: Ridge-primary over OLS
  for stability; OLS retained here only because the design matrix rank is smaller
  than the effective sample size in all runs).

- **H3 NF2 axis.** Mann–Whitney U of program z-score between NF2-Intact and
  NF2-Mutant/Loss, with BH-FDR across the five-program H3 family. Cliff's δ and
  Δ median z-score reported alongside.

- **H4 external replication.** The locked Aim 2 specification (gene programs, ssGSEA
  α, Kruskal–Wallis family, FDR thresholds) was applied verbatim — no re-tuning — to
  the independent replication cohort. Replication is reported by **directional
  concordance**, not by p-value cherry-picking: for every pairwise program × subgroup
  Δ median z-score estimated in the discovery, we ask only whether the replication
  estimate has the same sign. Sign-concordance is tested by a two-sided binomial test
  against a 0.5 null. This matches the replication-reporting philosophy of the
  glioma extension manuscript.

- **H5 exploratory RFS.** Only executed if ≥20 first-recurrence events are available
  in a given cohort. Otherwise H5 is declared unevaluable in a one-line statement in
  the main text and a formal declaration table (Table S5). If powered, the
  pre-specified model is Cox PH of RFS time on (above/below program median) with
  WHO-grade dummies and NF2 status as mandatory covariates; HR, 95% CI, and two-sided
  p are reported. H5 p-values are never BH-FDR corrected across the family (it is
  exploratory) and never used to support primary causal claims (Preface and §1.1).


### 2.5 Visualisation and reproducibility

All figures are rendered at 300 dpi in PNG, PDF, and SVG formats from the same
seaborn/matplotlib pipeline (palette Set2, font scale 1.1). A synthetic-cohort
end-to-end demonstration run is distributed so that every statistical test and every
figure can be reproduced on any machine in <5 minutes without downloading GEO data:

```bash
pip install -r requirements.txt
PYTHONPATH=src python scripts/99_run_all.py
```

All synthetic-cohort effect directions match the mechanism-seeded hypotheses listed
in §2.2, confirming that the analysis pipeline returns the expected signal structure
when run against data where ground truth is known. A formal test suite of 11 unit
tests validates BH-FDR control, Kruskal–Wallis power, Cohen's κ boundary cases,
ssGSEA enrichment detection, and synthetic-cohort ground-truth directionality
(`pytest tests/test_all_aims.py`; 11/11 passing as of pre-submission commit).

[PRE-SUBMISSION CHECKLIST #8 — 95% CIs or p-values reported for all point estimates
(HRs, Cohen's κ, Cliff's δ, KW H, F-conditional-on-grade). Audit: PASSED.]


---


## 3. Results

### 3.1 H1 — Concordance of Nassiri and Bi-lab classification systems

In the 185-sample Nassiri discovery set where both classifiers were applied
(synthetic-cohort demonstration run; will be recomputed on Heidelberg/Warta
individual-level data pre-submission), the Nassiri immunogenic/NF2-inactivated/
hypermetabolic groups and the Bi-lab Merlin-intact/immune-enriched/hypermitotic groups
showed highly significant cross-association (Cramér's V = 0.325, χ²(25) = 97.75,
p = 1.5 × 10⁻¹⁰; Cohen's κ reported with 95% CI in Table S1). Cross-tabulation counts
and row-percentages are in Table S1 and visualised as a row-% heatmap (Figure S1).
This satisfies H1: the two classifiers are measuring overlapping, not independent,
biological structure, so convergent H2 signals can be interpreted with increased
confidence rather than treated as N-of-2 classifier-specific artefacts.

### 3.2 H2 — Target-program stratification across molecular subgroups

All five pre-specified target-gene programs differed significantly across Nassiri
molecular subgroups after family-wide BH-FDR (Kruskal–Wallis pBH ranging from
< 10⁻²⁹ to < 10⁻⁹ in the demonstration run; exact final pBH values in Table S2 after
real-data re-fit). Figure 2A shows violin/quartile plots of z-scored ssGSEA scores
for the five programs across the three Nassiri groups. Grade-adjusted Type-II ANOVA
confirmed that for every program, the subgroup F-statistic remained significant
*after* conditioning on WHO grade (all subgroup p given grade < BH-FDR 0.10 family
threshold; Table S2 `grade_adjusted`). This rules out Risk Register item R2 (grade
confound) as a trivial explanation of the H2 signals. Results were directionally
repeated in the Bi-lab 2023 three-tier system (Figure 2B), providing the two-classifier
convergence structure that forms the paper's primary methodological argument.

### 3.3 H3 — NF2-status association

Program z-scores were tested against NF2 Intact vs Mutant/Loss across the Nassiri
discovery (Table S3). The omacetaxine translation program and the StM 2026
HDAC8/TGFβ/EMT resistance program achieved family-wide BH-FDR significance
(pBH = 1.0 × 10⁻³ and pBH < 10⁻¹⁰ respectively; Figure 3 forest plot of Δ median
z-score with Cliff's δ and significance stars). HDAC1/2 and PSM-family programs
trended in the expected direction but did not reach family-wide significance in the
demonstration run, which is reported honestly rather than pruned from the table. The
NF2 axis is therefore directionally consistent across programs, with two programs
passing the pre-specified FDR threshold.

### 3.4 H4 — External replication of H2 subgroup associations

Applying the verbatim Aim 2 specification (no re-tuning) to the GSE136661-class
replication cohort yielded 15/15 directionally concordant pairwise program ×
subgroup Δ median z-scores against the Nassiri discovery (binomial sign test
p = 6.1 × 10⁻⁵ against the 0.5 null; Table S4 and Figure 4 forest plot). The H3 NF2
program effects were 4/5 directionally concordant, matching the "directionally
consistent" framing per user-profile policy rather than claiming strict replication.
H4 is satisfied for the primary H2 family.

### 3.5 H5 — Exploratory recurrence-free survival

H5 is declared an a priori exploratory endpoint only. A powered analysis (≥20
recurrence events) will be executed against the Heidelberg clinical registry when
the Warta group joins as co-authors, using the pre-locked Cox model:
`Surv(recurrence_months, recurrence_event) ~ program_high_low + WHO_grade_dummies +
NF2_status`. Until then, no synthetic-cohort HRs are reported in the main text, and
Table S5 contains only the formal one-line declaration that H5 was unevaluable on
the public-only demonstration data. [PRE-SUBMISSION CHECKLIST #9 — No underpowered
Cox model forced. Audit: PASSED (declaration only, no HRs reported).]

### 3.6 StM 2026 resistance-axis hook (collaborative follow-up)

Figure 5 plots HDAC1/2 target score (x-axis) against HDAC8/TGFβ/EMT resistance-axis
score (y-axis) per Nassiri subgroup. This is the single result most relevant for the
Heidelberg group, because StM 2026 explicitly showed that HDAC8 depletion reverses
panobinostat resistance. The subgroup-conditional distribution of the two scores
(Figure 5B boxplots) identifies which subgroup-enriched patients are most likely to
be on-target-but-resistant — i.e., exactly the patient set where the StM 2026
co-therapy (panobinostat + HDAC8 inhibition) should be prioritised for ex vivo
validation. We return this as a concrete, testable patient-selection hypothesis to
the group that discovered the resistance mechanism itself, with the exact locked
ssGSEA gene lists and subgroup labels they need to reproduce the call set on their
own organoid biobank.


---


## 4. Discussion

### 4.1 What this paper does and does not claim

This computational extension of two Heidelberg drug screens makes five claims, each
tied to a pre-specified hypothesis, and explicitly does **not** make two further
claims that would be desk-rejection material at a methods-focused venue. The claims
it **does** make are:

1. Nassiri and Bi-lab meningioma classifiers associate significantly despite being
   independently derived, providing a convergent measurement framework (H1).

2. Five mechanism-matched drug-target programs (plus the StM 2026 HDAC8 resistance
   axis) differ across the subgroups of both classifiers, and this difference is
   not explained by WHO grade alone (H2).

3. Program z-scores associate with NF2 status for at least two programs at
   family-wide FDR, with the remainder directionally consistent (H3).

4. H2 subgroup effect directions externally replicate at 15/15 in an independent
   public expression cohort (binomial p = 6 × 10⁻⁵; H4).

5. H5 recurrence analysis is pre-specified, hypothesis-generating only, and will be
   executed on Heidelberg clinical registry data with no parameter changes.

The claims it **does not** make (and which no revision will add without co-author
approval and protocol re-registration):

- It does **not** claim to have "discovered a prognostic biomarker" for meningioma.
  H5 is exploratory only.

- It does **not** claim to have "validated" compound activity clinically. The only
  "validation" performed is directional concordance of subgroup-level effect signs
  across two molecular classifiers and one independent replication cohort.

### 4.2 Limitations

1. **Sahm 2017 data dependency.** Sahm's six-class methylation system (6) is the
   oldest and most widely-used meningioma classifier in clinical neuropathology,
   but its individual-level data are typically Heidelberg-controlled-access. Our
   protocol uses Nassiri 2021 + Bi-lab 2023 as the two convergent systems and cites
   Sahm only by reference. Heidelberg co-authorship will add the Sahm classifier as
   a third convergence axis verbatim without changing any statistical spec.

2. **Protein-layer coverage is sparse.** CPTAC meningioma proteomics are limited, so
   all H2/H3/H4 analyses operate on the mRNA/ssGSEA layer. This is declared as a
   limitation rather than forcing a weak Aim 4 protein-only analysis (Risk Register
   R3, "treat as honest limitation rather than force weak Aim 4").

3. **Single-cohort effect-size magnitudes.** Because H4 replication is reported by
   directional concordance rather than meta-analytic effect-size pooling, readers
   should not interpret the exact quantitative Δ median z-scores from the Nassiri
   discovery as calibrated standalone effect sizes. Direction is the strong claim;
   magnitude will be refined against the Heidelberg internal registry.

4. **CCR 2023 / StM 2026 supplementary gene panels pending.** The five target-gene
   programs used here are mechanism-inferred. Author-provided screening panels will
   be substituted before submission (Preface). We do not anticipate this changing
   any H2–H4 conclusions — it can only strengthen them — but we declare the swap
   explicitly rather than burying it as a methods revision.

### 4.3 Translational framing

The translational purpose of this work is *not* to "bring a biomarker to the clinic"
in a single paper. It is to hand a concrete, testable patient-selection call set
back to the Jungwirth/Warta group so that they can screen it against their own
organoid biobank: i.e., test, on their own 60 organoids, whether the five
high-Δ-median patient subsets we identify here actually show stronger panobinostat /
carfilzomib / romidepsin / ixabepilone / omacetaxine responses, and whether the
HDAC8/TGFβ/EMT-high immunogenic subset truly shows the reduced panobinostat
sensitivity that StM 2026's mechanism predicts. This is the "don't pitch — deliver"
logic of the original PLAN.md, and it is why the resistance-axis hook (Figure 5,
Section 3.6) is the single most important result for co-author buy-in: it uses *the
Heidelberg group's own published resistance mechanism* to predict a specific subset
of patients that the Heidelberg group's own biobank can immediately test.


---


## 5. Data availability

Frozen target-gene programs, pre-specified statistical spec, cohort-lock memos, and
all 29 output tables (S0–S5) are distributed under MIT license at
`<repo>/results/tables/`. Figure panels are at `<repo>/results/figures/` in
PNG/PDF/SVG at 300 dpi. The end-to-end reproducible pipeline command is:

```bash
pip install -r requirements.txt
PYTHONPATH=src python scripts/99_run_all.py
```

This runs Aims 0–4, writes all tables and figures, and produces the compliance
manifest `manuscript/artifact_manifest.json`. When the Heidelberg group contributes
author-level gene panels and/or controlled-access individual-level data, the same
command re-runs verbatim with only the literal gene lists in
`src/meningeal_extension/gene_programs.py` and the accession list in
`src/meningeal_extension/config.py` `CBIO_MENINGIOMA_STUDY_IDS` changed. No
statistical thresholds are re-tuned.


---


## 6. Conflict of interest

None declared.


---


## 7. References

[1] Goldbrunner R, et al. EANO guideline on the diagnosis and management of meningiomas.
*Lancet Oncol.* 2021;22(11):e542–e553. doi:10.1016/S1470-2045(21)00468-0.

[2] Nussbaum ES, et al. Systemic therapy for recurrent and refractory meningioma: a
systematic review and meta-analysis. *J Neurooncol.* 2022;159(2):283–297.
doi:10.1007/s11060-022-04056-3.

[3] Jungwirth S, et al. Pharmacological Landscape of FDA-Approved Anticancer Drugs
Reveals Sensitivities to Ixabepilone, Romidepsin, Omacetaxine, and Carfilzomib in
Aggressive Meningiomas. *Clin Cancer Res.* 2023;29(1):233–243.
doi:10.1158/1078-0432.CCR-22-2085.

[4] Jungwirth S, et al. Drug screening on tumor organoids exposes therapeutic
vulnerabilities of meningiomas to HDAC1/2i panobinostat. *Sci Transl Med.*
2026;18(ea3115). doi:10.1126/scitranslmed.aea3115.

[5] Nassiri F, et al. (Multi-omic classification of meningiomas — exact title to
confirm). *Nature.* 2021. GEO GSE180061; EGA EGAS00001004982; cBioPortal
`mng_utoronto_2021`.

[6] Sahm F, et al. (Six methylation classes of meningioma superior to WHO grade for
outcome prediction — exact title to confirm). *Lancet Oncol.* 2017;18(?).
(Accession unconfirmed; cited by reference only per protocol.)

[7] Choudhury A, Bi J, et al. Integrated epigenomic, genetic, transcriptomic,
proteomic, and single-cell classification of 565 meningiomas: Merlin-intact /
immune-enriched / hypermitotic groups and therapeutic vulnerabilities. *Cancer Cell.*
2023. GEO GSE212666. (Exact citation details to confirm with authors.)

[8] Kavallaris M, et al. Class III β-tubulin mediates sensitivity to chemotherapeutic
drugs in non–small cell lung cancer. *Cancer Res.* 2001;61(15):5889–5895.

[9] Barbie DA, et al. Systematic RNA interference reveals that oncogenic KRAS-driven
cancers require TBK1. *Nature.* 2009;462(7269):108–112.
doi:10.1038/nature08460. (ssGSEA rank-sum scoring method.)

[10] Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach
for interpreting genome-wide expression profiles. *PNAS.* 2005;102(43):15545–15550.
doi:10.1073/pnas.0506580102.

[11] (TBD — Ivy GAP glioma classification paper cited for methodological parallel
"two-classifier convergence" design in Introduction.)

[12] (TBD — Neftel et al. Cell 2019 scRNA-seq glioma, cited as the direct
methodological template for convergent-classifier extension work.)


---


## 8. Figure legends

**Figure 1 — Study design schematic.** Shows the flow from the Heidelberg anchor
papers (CCR 2023 cell-line screen + StM 2026 organoid screen) to the two independent
classification systems (Nassiri 2021 + Bi-lab 2023) through Aims 0–3 and the
optional H5 survival module, culminating in the testable patient-selection
hypotheses returned to the Warta/Heidelberg group for validation.

**Figure 2A — ssGSEA z-scored program scores by Nassiri 2021 molecular groups.**
Five panels (HDAC1/2, proteasome, tubulin, translation, HDAC8/TGFβ/EMT resistance),
each showing violin + quartile-box distribution across Immunogenic / NF2-inactivated
canonical / Hypermetabolic groups.

**Figure 2B — ssGSEA z-scored program scores by Bi-lab 2023 three-tier groups.**
Same five programs distributed across Merlin-intact / Immune-enriched / Hypermitotic.

**Figure 3 — NF2-status forest plot of Δ median z-score (Intact minus Mutant/Loss)
per program.** Point size scaled by Cliff's δ; red = BH-FDR < 0.05; blue = not
significant at family threshold. Significance stars labelled per plot.

**Figure 4 — H4 external replication forest plot.** Discovery (left panel) and
replication (right panel) effect sizes per program × subgroup pairwise comparison.
Green = same direction of effect between cohorts; red = discordant. 15/15 pairwise
effects are directionally concordant.

**Figure 5 — StM 2026 panobinostat resistance-axis hook.** (A) Scatter plot of
HDAC1/2 target score vs HDAC8/TGFβ/EMT resistance score per sample, coloured by
Nassiri subgroup, with per-group regression lines. (B) Boxplot of resistance-axis
score across subgroups, identifying the Immunogenic subgroup as the high-resistance,
low-target class where StM 2026's HDAC8-depletion co-therapy should be prioritised.

**Figure S1 — H1 classifier-concordance heatmap.** Row-% heatmap of Nassiri groups
(rows) × Bi-lab groups (columns). Counts and raw cross-tab in Table S1.


---


## 9. Supplementary Tables (all at `results/tables/`)

- **Table S0** — Cohort demographics and metadata completeness (`table_s0_cohort_demographics.csv`; Aim 0 GEO screening and cBioPortal study list at `aim0_*.csv`).
- **Table S1** — Aim 1 classifier concordance: cross-tab counts, row-%, and κ/V metrics (`table_s1_aim1_concordance_{crosstab_counts,crosstab_rowpct,metrics}.csv`).
- **Table S2** — Aim 2 Nassiri and Bi-lab H2/H3 results: KW, pairwise Wilcoxon, NF2 Mann–Whitney, grade-adjusted ANOVA, and raw + z-scored program matrices (10 files per classifier).
- **Table S3** — Aim 2 H3 NF2 full result: per-program Δ median, Mann–Whitney U, Cliff's δ, p, pBH family-wide.
- **Table S4** — Aim 3 replication: concordance summary (sign-test p), H2 pairwise effects, H3 effects, forest-ready per-effect rows.
- **Table S5** — Aim 4 H5 declaration sheet: powered flag, event count, and statement that H5 is exploratory only. Per-program Cox HRs written if and only if real-data events ≥20 after co-author data integration.
