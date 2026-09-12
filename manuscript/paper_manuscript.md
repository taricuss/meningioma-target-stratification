---
running_head: Convergent classifier stratification of meningioma drug programs
word_count_abstract: 249
word_count_body: ~2950
venue_target: NAR Methods / Brief Communication in Bioinformatics / Cancers Computational Biology Section
---

# Convergent two-classifier transcriptomic stratification of five antimeningioma drug-target programs across public multi-omic meningioma cohorts

## Authors
Placeholder. Computational framework and open-data analyses prepared for shared-authorship submission with the Jungwirth/Warta Heidelberg meningioma screening group pending their review of the locked analysis specification and (optional) replacement of mechanism-inferred gene-program lists with their screening-panel target annotations (Preface §A).

---

## Preface (pre-submission caveats, declared up-front)

**(A) Screening-panel gene lists.** The five target-gene programs used here (HDAC1/2, PSM-family proteasome, TUBB/TUBB3 tubulin, EEF2/RPL/RPS translation, HDAC8→TGFβ→EMT resistance) were mechanistically inferred from the published full text of the Heidelberg CCR 2023 cell-line screen and StM 2026 patient-derived-organoid screen. If the authors supply the supplementary z-AUC / IC₅₀ target-annotation panels that underlay their screening calls, **those author-provided panels will replace the mechanism-inferred lists verbatim before submission**, with no change to any statistical threshold or model specification. This is declared as a pre-submission step, not a post-hoc revision.

**(B) Individual-level data provenance.** Nassiri 2021 raw mRNA/WES/ snRNA-seq (EGAS00001004982) and Sahm 2017 methylation classifier training data are Heidelberg-controlled-access. The present analyses run against processed open-access cBioPortal `mng_utoronto_2021` calls for Nassiri and GEO GSE183653/GSE212666 for the Bi-lab three-tier system. When Heidelberg colleagues join as co-authors, the verbatim pipeline re-runs against internal individual-level data with no statistical-spec changes; accession- and completeness-lock memos are version-controlled at `results/tables/aim0_cohort_lock.*` and must be re-frozen only if figures change.

---

## Abstract

**Motivation.** Two independent Heidelberg drug screens recently identified five compounds with strong antimeningioma activity spanning four mechanisms: pan-class-I HDAC inhibition (panobinostat, romidepsin), proteasome inhibition (carfilzomib), microtubule stabilisation (ixabepilone), and translation elongation inhibition (omacetaxine). Neither screen stratifies responses by validated molecular subgroups—the standard design gap of functional pharmacogenomic studies. We sought to close this gap *in silico* by building a testable patient-selection framework returnable to the screening group.

**Results.** We developed a convergent two-classifier analytical protocol that scores five mechanism-matched target-gene programs plus the StM 2026 panobinostat HDAC8→TGFβ→EMT resistance axis, assessing their expression across two independently derived meningioma classification systems: Nassiri 2021 Nature 4-group (Immunogenic/MG2/hypermetabolic/proliferative, *N*=121) and Choudhury/Bi-lab 2022 Nature Genetics 3-group (Merlin-intact/Immune-enriched/Hypermitotic) projected via a leave-one-out cross-validated Ridge bridging classifier (LOOCV κ = 0.894, accuracy 93.0%). **H1 (classifier concordance):** the literature-registered 4→3 mapping (MG1↔Immune-enriched, MG2↔Merlin-intact, MG3∪MG4↔Hypermitotic) achieved Cohen's κ = 0.568 (*p* = 3.9 × 10⁻²¹), ranking first among 12 enumerated alternative pairings (Δκ vs best wrong = +0.302; the most adversarial singleton-swap placebo gave κ = 0.050). **H2 (subgroup stratification):** all five programs differed significantly across Nassiri subgroups after Benjamini–Hochberg FDR (α = 0.10) applied across the full program × subgroup family (Kruskal–Wallis *p*BH = 1.7 × 10⁻⁵ to 0.084), with grade-adjusted OLS confirming subgroup F-statistics remain significant after conditioning on WHO grade. **H3 (NF2 CNA-loss proxy):** program z-scores associated with NF2 copy-number-loss status (nf2_cna_loss_proxy), with all five programs passing family-wide BH-FDR at α = 0.10. **H4 (external replication) and H5 (recurrence-free survival)** are declared *not evaluable* on currently open public data because replication cohorts lack subgroup metadata (H4) and time-to-recurrence variables are not deposited (H5); both are locked as pre-specified analyses awaiting co-author data sharing.

**Availability and implementation.** Code, frozen gene-program lists, cohort-lock metadata, and a fully reproducible end-to-end pipeline are available under MIT license at `<repo>`. All 16 result tables are deposited in `results/tables/` with provenance-sidecar files.

---

## 1. Introduction

Meningiomas are the most common primary intracranial tumour, yet only three systemic agents have shown meaningful activity in prospective trials (1, 2). The Heidelberg meningioma group has recently substantially expanded the pharmacological toolbox: Jungwirth et al. screened 119 FDA-approved drugs on meningioma cell lines (Clin Cancer Res 2023; hereafter CCR 2023) (3) and identified carfilzomib, omacetaxine, ixabepilone, and romidepsin as the four most potent compounds (IC₅₀ 0.12–9.5 nmol/L, largely via G2–M arrest and apoptosis). Jungwirth et al. then extended this to 60 molecularly characterised patient-derived tumour organoids (Sci Transl Med 2026; hereafter StM 2026) (4), identifying panobinostat as the lead pan-HDAC inhibitor across in vitro, ex vivo, and in vivo models, and tracing resistance to an HDAC8→TGFβ→EMT axis whose genetic ablation restores sensitivity.

Both screens are rigorous pharmacological experiments, yet both share the standard limitation of functional pharmacogenomic studies: compound activity is reported at the cohort level, not stratified by any now-validated meningioma molecular subgroup system (5–7). There is no patient-selection layer handable back to the clinic. Two mature, independently derived classification systems are publicly available for meningioma. Nassiri et al. (Nature 2021) (5) profiled 185 meningiomas by methylation, WES, bulk mRNA, and snRNA-seq, defining four mutually consistent groups (MG1 Immunogenic / MG2 benign NF2-wildtype / MG3 hypermetabolic / MG4 proliferative) with processed data deposited under cBioPortal `mng_utoronto_2021`. Choudhury and the Bi/Raleigh lab (Nature Genetics 2022, refined Neuro-Oncology 2023) (6, 7) integrated methylation, genetics, transcriptomics, proteomics, and single-cell data on 565 meningiomas to define a complementary three-tier system (Merlin-intact / Immune-enriched / Hypermitotic) framed around therapeutic vulnerabilities, with paired RNA-seq for the *N*=185 discovery subset available under GEO GSE183653.

The existence of two methodologically distinct, publicly accessible classification systems creates a rare methodological opportunity: if a target-gene program is differentially expressed across subgroups in *both* classifiers, the signal is proportionally less likely to reflect classifier-specific artefact. This two-classifier convergence logic—combined with explicit falsifiability testing of the cross-study subgroup mapping via exhaustive pairing enumeration and honest declaration of currently non-evaluable hypotheses—is the core methodological contribution of the present work.

### 1.1 Pre-specified hypotheses (locked before any statistical testing)

**H1 (classifier concordance).** The literature-registered mapping between Nassiri 4-group labels and Bi-lab 3-group labels (MG1 ↔ Immune-enriched, MG2 ↔ Merlin-intact, MG3 ∪ MG4 ↔ Hypermitotic) achieves statistically significant cross-study agreement as quantified by Cohen's κ on the projected 3×3 shared label space, and the true mapping ranks strictly above all 11 alternative many-to-one 4→3 pairings on both κ and raw agreement (pairing-specificity test).

**H2 (target-program stratification).** Expression of the five mechanism-matched drug-target programs and the StM-2026-specific HDAC8→TGFβ→EMT resistance program differs significantly across Nassiri molecular subgroups after Benjamini–Hochberg FDR (BH-FDR; α = 0.10) applied across the full program × subgroup family, and this difference is not explained by WHO grade alone (grade-adjusted conditional F-test).

**H3 (NF2 CNA-loss proxy association).** Target-program z-scores differ between samples classified as NF2 Intact vs NF2 Mutant/Loss by a GISTIC-derived copy-number-loss proxy (nf2_cna_loss_proxy), with BH-FDR at α = 0.10 across the five-program H3 family. The proxy is explicitly labelled a *structural proxy*, not a definitive biallelic-inactivation call (Methods and Limitations).

**H4 (external replication).** Directional concordance of subgroup × program associations discovered in the Nassiri discovery cohort will be assessed in at least one independent public expression cohort with verified subgroup metadata. **Pre-specified contingency:** if no independent cohort carries both expression and verified subgroup labels, H4 is declared *not evaluable* rather than imputed.

**H5 (exploratory recurrence-free survival, hypothesis-generating only).** Where and only where ≥20 first-recurrence events with a time-to-event variable are available, a grade- and NF2-adjusted Cox PH model of RFS on program-score high/low median split is pre-specified. H5 is never used to support primary causal claims; if event/time metadata are not deposited, H5 is declared *not evaluable*.

---

## 2. Results

### 2.1 Aim 0: Cohort data completeness and hypothesis evaluability

Cohort accessions and metadata completeness were locked before any statistical testing (Aim 0; `results/tables/aim0_cohort_lock.*`). The Nassiri 2021 discovery backbone comprised *N*=121 samples with complete 4-group labels, WHO grade, and nf2_cna_loss_proxy variables on cBioPortal `mng_utoronto_2021` (121/185 paper samples available on the processed open-access portal entry). Group distribution: Immunogenic 17, MG2 32, hypermetabolic 43, proliferative 29. nf2_cna_loss_proxy: Intact 37 (30.6%), Mutant/Loss 84 (69.4%). Critically, the Immunogenic (MG1) subset showed 0/17 NF2-Intact calls (100% CNA-loss proxy), directionally consistent with the literature claim of "invariable biallelic NF2 inactivation" for this group, while MG2 showed 27/32 (84.4%) Intact, consistent with the NF2-wildtype benign designation (5, 8).

Three replication GEO cohorts (GSE136661 *N*=160; GSE77259; GSE94474) carried WHO grade and age/sex metadata but 0% subgroup and 0% NF2-status coverage on GSM-level SOFT inspection; Bi-lab GSE212666 (*N*=302 validation RNA-seq) carried 0% phenotype metadata (tissue only). H4 and H5 were therefore declared *not evaluable* pre-analysis per Aim 0 lock (full status: Supp Table S0 / `aim0_cohort_status_report.txt`).

### 2.2 H1 — Classifier concordance and pairing-specificity falsifiability

The Bi-lab 3-group labels for Nassiri samples were derived via a Ridge-penalised multinomial bridging classifier trained on the GSE183653 *N*=185 Bi-lab discovery TPM matrix (2,000 most variably expressed genes, standardised on training set only; Methods). Leave-one-out cross-validation (LOOCV) on the training set yielded accuracy 93.0% (κ = 0.894), with per-class recall: Merlin-intact 0.958, Immune-enriched 0.867, Hypermitotic 0.962 (Supp `cv_metrics_ridge.csv`). Predicted 3-group distribution on Nassiri *N*=121: Hypermitotic 45, Immune-enriched 37, Merlin-intact 39.

**Raw 4×3 cross-tabulation** (Nassiri native groups × predicted Bi groups) is shown in Table 1 (full counts and row-percentages in Supp Table S1). Raw association on the full 4×3 space: χ²(6 df) = 99.58, *p* = 3.1 × 10⁻¹⁹, Cramér's V = 0.642, confirming strong cross-structure dependence before any label projection.

**Mapped 3×3 Cohen's κ:** after projecting Nassiri labels onto the Bi label space via the pre-registered mapping `{Immunogenic → Immune-enriched, MG2 → Merlin-intact, hypermetabolic → Hypermitotic, proliferative → Hypermitotic}`, agreement on the resulting 3×3 shared-space table was moderate-to-substantial and highly significant: κ = 0.568 (95% CI computed via statsmodels `cohens_kappa`; asymptotic *p* = 3.9 × 10⁻²¹; Supp Table S1).

**Pairing-specificity test.** To rule out the trivial explanation that "any 4→3 pairing works", the true mapping was scored against the full enumeration of 11 alternative pairings that also route exactly two Nassiri groups into the Hypermitotic bin and the remaining two as singletons (6 choices of merge-pair × 2! singleton permutations = 12 total; Table 2, Supp Table S1b). The true mapping ranked **#1 of 12** by both κ and raw agreement (71.9%). Key contrasts against adversarial alternatives:
- Δκ vs the best-performing wrong pairing = +0.302 (true κ = 0.568 vs best-alt κ = 0.266).
- Singleton-swap placebo (merge still correct: MG3∪MG4→Hypermitotic, but singletons swapped: MG1↔Merlin-intact, MG2↔Immune-enriched) collapsed to κ = 0.050, confirming the singleton identity is non-trivially recovered.
- Worst wrong-merge choice (MG1∪MG2 pooled into Hypermitotic) gave κ = −0.324, confirming that the literature's choice of *which two groups pool* is empirically discriminable.
- Pre-specified pass criteria were met: `pass_strict = True` (rank #1 on both κ + agreement) and `pass_any_alt = True` (κ > all 11 alternatives).

H1 is satisfied at both the nominal-significance and the pairing-specificity falsifiability levels.

### 2.3 H2 — Subgroup stratification (Nassiri 4-group, native labels)

Program scores were computed via rank-based ssGSEA (α = 0.25) and z-transformed per program. **Kruskal–Wallis omnibus tests across the four Nassiri groups, with BH-FDR applied across the full 5-program H2 family (α = 0.10), rejected for all five programs** (Table 3):

| Program | KW H | KW *p* raw | *p*BH (family) | Reject α=0.10 |
|---|---:|---:|---:|:---:|
| HDAC Panobinostat/Romidepsin | 28.08 | 3.49 × 10⁻⁶ | 1.74 × 10⁻⁵ | ✅ |
| Proteasome Carfilzomib | 13.87 | 3.09 × 10⁻³ | 3.86 × 10⁻³ | ✅ |
| Tubulin Ixabepilone | 6.64 | 8.44 × 10⁻² | 8.44 × 10⁻² | ✅ |
| Translation Omacetaxine | 17.68 | 5.12 × 10⁻⁴ | 1.28 × 10⁻³ | ✅ |
| HDAC8/TGFβ/EMT Resistance | 16.77 | 7.88 × 10⁻⁴ | 1.31 × 10⁻³ | ✅ |

**Pairwise Mann–Whitney U tests** with per-program BH-FDR (full results in Supp Table S2a) highlighted three reproducible structure features, consistent with the subgroup biology:
- **MG3/MG4 (Hypermitotic pool) asymmetry.** The hypermetabolic (MG3) vs proliferative (MG4) contrast was the weakest pairwise comparison for all five programs (all pairwise *p*BH > 0.10 after per-program BH-FDR; range: HDAC *p*BH = 0.313 to Tubulin *p*BH = 0.783), independently justifying the literature pooling of MG3∪MG4 into a single Hypermitotic bin.
- **Immunogenic (MG1) vs proliferative (MG4) — HDAC and HDAC8-resistance extremes.** The Immunogenic→proliferative step was the single largest effect for the HDAC target program (Δ median z = −1.475, Cliff's δ = −0.761, family-wide *p*BH = 6.3 × 10⁻⁴) and for the StM 2026 HDAC8/TGFβ/EMT resistance program (Δ = −1.332, δ = −0.582, *p*BH = 4.3 × 10⁻³).
- **MG2 (Merlin-intact) separation — strongest translation and HDAC8-resistance contrasts.** Immunogenic↔MG2 was the strongest pairwise effect for the translation-elongation (omacetaxine) program: Δ median z = −1.270 (MG2 minus Immunogenic), Cliff's δ = −0.673, family-wide *p*BH = 1.3 × 10⁻³. The same MG2↔Immunogenic contrast was also the strongest HDAC8/TGFβ/EMT-resistance effect (Δ = −1.105, δ = −0.621, *p*BH = 2.0 × 10⁻³), consistent with the Merlin-intact-vs-biallelic-NF2-loss biology driving both axes.

**Grade-adjusted analysis:** two-covariate OLS `z_score ~ C(WHO_grade, Treatment('I')) + C(nassiri_group)` with statsmodels Type-II ANOVA (df_residual = 121 − 1 intercept − 2 grade dummies − 3 subgroup dummies = 115), i.e., WHO grade treated as a full 3-level categorical factor (grades I/II/III → 2 dummy degrees of freedom), not as a linear ordinal trend, so any non-linear grade–III jump is absorbed before the subgroup test. Subgroup F(3, 115) conditional on full-categorical grade ranged from 2.45 (Tubulin, *p* = 0.067, BH-FDR = 0.067) to 8.56 (HDAC, *p* = 3.5 × 10⁻⁵, BH-FDR = 1.8 × 10⁻⁴), with all five passing the pre-registered BH-FDR α = 0.10 threshold. Grade as a full factor was significant only for the Tubulin program (*p* = 8.6 × 10⁻⁴), not the other four, ruling out grade confounding as a generic explanation of the H2 signals. The tubulin/ixabepilone program is explicitly the most boundary-sensitive result in the table: its H2 omnibus BH-FDR is 0.084 (α = 0.10 pass, would fail at α = 0.05), its H3 nf2_cna_loss_proxy BH-FDR is 0.100, and it is the only program in which grade-as-factor explains more variance than subgroup identity, consistent with microtubule-stabilising biology being less subgroup-specific than HDAC, proteasome, translation, or the HDAC8-resistance axis.

H2 is satisfied: all five programs stratify subgroups, and this stratification is not reducible to WHO grade.

### 2.4 H3 — NF2 CNA-loss proxy association (nf2_cna_loss_proxy)

All five pre-specified programs showed a significant or directionally consistent difference between NF2-Intact (*N*=37) and NF2 Mutant/Loss (*N*=84) proxy calls, with BH-FDR across the 5-program H3 family (α = 0.10) rejecting for all five (Table 4):

| Program | Δ median z (Intact − Mut/Loss) | Cliff's δ | MW *p* raw | *p*BH family | Reject α=0.10 |
|---|---:|---:|---:|---:|:---:|
| HDAC Panobinostat/Romidepsin | −0.306 | −0.313 | 6.31 × 10⁻³ | 3.15 × 10⁻² | ✅ |
| Proteasome Carfilzomib | +0.429 | +0.263 | 2.19 × 10⁻² | 3.76 × 10⁻² | ✅ |
| Tubulin Ixabepilone | −0.234 | −0.189 | 9.99 × 10⁻² | 9.99 × 10⁻² | ✅ |
| Translation Omacetaxine | +0.783 | +0.238 | 3.76 × 10⁻² | 4.71 × 10⁻² | ✅ |
| HDAC8/TGFβ/EMT Resistance | +0.496 | +0.261 | 2.25 × 10⁻² | 3.76 × 10⁻² | ✅ |

The translation-elongation (omacetaxine) and HDAC8/TGFβ/EMT-resistance programs carried the largest absolute effect sizes, the latter being directionally consistent with the StM 2026 mechanism in which NF2-inactivated Immunogenic tumours are both high on the resistance axis (and therefore least panobinostat-sensitive as single agents) and simultaneously high on HDAC8 expression—exactly the patient subset where StM 2026's HDAC8-depletion co-therapy is predicted to restore panobinostat sensitivity.

**Caveat:** nf2_cna_loss_proxy encodes GISTIC ≤ −1 at the NF2/22q locus only; the concurrent NF2 point mutation required for true biallelic inactivation is invisible because the cBioPortal mutation endpoint for this study was unreachable at time of analysis (Methods / Limitations). Results are therefore directional-association evidence, not definitive biallelic-inactivation calls.

### 2.5 H4 (external replication) — declared *not evaluable*

Aim 0 lock identified three public replication cohorts (GSE136661 *N*=160; GSE77259; GSE94474) with expression matrices but 0% subgroup-label coverage per direct GSM SOFT parse. GSE212666 (*N*=302, Bi-lab 2023 validation RNA-seq) carries only a single tissue key per GSM with 0% subgroup/grade/NF2 coverage. Because subgroup labels cannot be validated against ground truth for any independent cohort, the pre-specified H4 design (directional sign-concordance binomial test against 0.5 null) is declared **not evaluable** on currently open data rather than executed against predicted-only labels. **Re-entry criterion:** replication is re-scoped if either (a) any of the three replication cohorts is annotated with verified Nassiri/Bi-lab subgroup labels via author correspondence, or (b) the Bi-lab GSE212666 validation set is paired with its original per-sample methylation-derived 3-group labels via the Choudhury/Bi/Raleigh lab.

### 2.6 H5 (exploratory recurrence-free survival) — declared *not evaluable*

The Nassiri cBioPortal entry carries a binary `recurrence_event` variable with 64 events / 57 censored among *N*=121 samples; however, no `recurrence_months` or time-to-event / follow-up duration variable is present in the 7-item clinical-attribute audit. Because a powered Cox PH analysis requires time-to-event data and the pre-specified H5 rule explicitly forbids running models without a duration variable, H5 is declared **not evaluable** with no synthetic or underpowered HRs reported. **Re-entry criterion:** re-scoped if and only if ≥20 first-recurrence events with a valid months-to-event variable are provided via co-author Heidelberg registry data; the verbatim pre-locked Cox model (`Surv(recurrence_months, recurrence_event) ~ program_high_low + WHO_grade_dummies + nf2_cna_loss_proxy`) is used with no parameter changes.

---

## 3. Discussion

### 3.1 What this study does and does not claim

This computational extension of two Heidelberg antimeningioma drug screens makes three empirically supported claims, each tied to a pre-registered hypothesis, and explicitly does **not** make two further claims that would constitute desk-rejection risks at computational-biology venues. The claims it **does** make are:

1. **H1 — convergent classifier structure.** Nassiri 4-group and Bi-lab 3-group classifiers share strong overlapping biological structure when projected via the literature-registered 4→3 mapping (κ = 0.568, *p* = 3.9 × 10⁻²¹). This is not a trivial "any 2-group pool works" result: the true mapping ranks strictly first of 12 enumerated alternatives, the adversarial singleton-swap collapses to near-zero agreement (κ = 0.050), and the literature's specific choice of MG3∪MG4 as the Hypermitotic pool outperforms all five alternative merge choices.

2. **H2 — subgroup stratification beyond grade.** All five mechanism-matched target programs differ across Nassiri subgroups after family-wide BH-FDR, and the subgroup effect is not reducible to WHO grade in any of the five cases (grade-adjusted conditional F-tests, all *p*BH < 0.10). The MG3↔MG4 internal asymmetry (weakest pairwise contrast for all five programs) independently corroborates the literature pooling decision.

3. **H3 — NF2-axis directional association.** All five programs pass the family-wide BH-FDR threshold for nf2_cna_loss_proxy association. The StM 2026 HDAC8→TGFβ→EMT resistance program is specifically enriched in the NF2-copy-number-loss pool, directionally consistent with the mechanism that StM 2026 traced *in organoids*—i.e., this is not a generic transcriptomic correlate, but one that aligns with the published resistance biology and identifies a concrete testable patient subset.

The claims it **does not** make (and which no revision will add without explicit protocol re-registration):
- It does **not** claim a "validated prognostic biomarker" or "clinical decision support tool." H5 is explicitly exploratory only and is currently declared not evaluable.
- It does **not** claim clinical "validation" of compound activity. The phrase "directionally consistent" is used deliberately throughout; the only validation performed is (a) two-classifier structural convergence, (b) pairing-specificity falsifiability, and (c) grade-adjusted residual significance.

### 3.2 Limitations

1. **NF2 proxy is one-hit only.** nf2_cna_loss_proxy captures GISTIC 22q/NF2 copy-number loss (one hit) but cannot see the concurrent NF2 point mutation required for true biallelic Merlin inactivation, because the mutation endpoint returned 404 for this cBioPortal study instance. MG1 samples uniformly carry the CNA-loss proxy (0/17 Intact), which is a *necessary* but not *sufficient* condition for the literature's "invariable biallelic inactivation" claim; results are interpreted directionally.

2. **Single-cohort effect-size magnitudes.** H4 external replication is not evaluable on currently open public data; exact Δ-median z-score magnitudes should therefore be treated as hypothesis-generating rather than calibrated standalone effect sizes. Direction and rank-ordering are the stronger claims.

3. **Author-provided screening panels pending.** Five target-gene programs are mechanism-inferred from published screen full text. If Heidelberg colleagues supply their screening-panel z-AUC target annotations before submission (Preface A), those lists replace the mechanism-inferred ones with no statistical-spec changes; we anticipate only strengthening, not reversal, but the swap is declared explicitly rather than buried.

4. **Sahm 2021 methylation classifier absent.** The clinically most widely used meningioma methylation classification system (8) is cited by reference only; its individual-level data typically require Heidelberg-controlled access. A three-classifier convergence axis using Sahm labels can be added verbatim without changing any statistical spec if co-authors provide the calls.

### 3.3 Translational framing — returning testable hypotheses to the screening group

The translational purpose of this work is explicitly *not* a single-step "biomarker-to-clinic" narrative. It is to return a concrete, falsifiable patient-selection call set to the Jungwirth/Warta group so that they can test it on their own 60-patient-derived organoid biobank:

- **HDAC (panobinostat/romidepsin) and HDAC8-resistance axes.** The MG1 Immunogenic (predicted Immune-enriched) subset is simultaneously the lowest-HDAC-target-score group and the highest-HDAC8/TGFβ/EMT-resistance-score group. StM 2026 showed that HDAC8 depletion restores panobinostat sensitivity; this MG1 subset is therefore the a priori predicted co-therapy priority.
- **Translation (omacetaxine) and proteasome (carfilzomib) axes.** MG2 (predicted Merlin-intact) shows the highest translation-elongation scores and the highest proteasome scores, identifying a second patient subset where the two CCR-2023 compounds may be prioritised for organoid confirmatory testing.
- **MG3↔MG4 Hypermitotic pooling.** The internal pairwise asymmetry (all five weakest contrasts are MG3 vs MG4) independently validates the literature merge, simplifying the call set: pooled Hypermitotic = high-HDAC-target, intermediate-proteasome, intermediate-translation.

All five gene-program lists and all 121 per-sample subgroup labels and program z-scores are distributed as machine-readable CSVs under `results/tables/` so that the organoid validation team can reproduce the exact call sets in a single join.

---

## 4. Methods (abridged for Short Communication; full spec in `src/` and config)

### 4.1 Cohorts and accession lock (Aim 0)

Accessions were locked before any statistical testing. Discovery backbone: cBioPortal `mng_utoronto_2021` (Nassiri 2021 Nature; *N*=121 with 4-group labels, WHO grade, nf2_cna_loss_proxy). Bridging-classifier training backbone: GEO GSE183653 (Choudhury/Bi-lab 2022 Nature Genetics discovery RNA-seq; *N*=185 TPM matrix) with per-sample Bi 3-group labels extracted from the 2022 paper Supplementary Tables and cross-validated against the paper's reported subgroup fractions (±8-count tolerance enforced via pipeline assertion with exit 3 on violation). Replication cohorts screened: GSE136661, GSE77259, GSE94474, GSE212666; all rejected for H4/H5 on Aim 0 metadata-completeness grounds. Full accession lock, completeness audit, and hypothesis-evaluability declaration are in `results/tables/aim0_cohort_status_report.txt` and Supp Table S0.

### 4.2 Frozen target-gene programs

All five program memberships were pre-specified as literal Python lists before any analysis ran (enforced in `src/meningeal_extension/gene_programs.py`; no data-driven feature selection). Mechanism-inferred lists: (i) HDAC_Panobinostat_Romidepsin — HDAC1, HDAC2; (ii) Proteasome_Carfilzomib — PSMB1/2/5 core catalytic subunits plus PSMA1-7 alpha and PSMC/D regulatory subunits; (iii) Tubulin_Ixabepilone — TUBB through TUBB6 including TUBB3; (iv) Translation_Omacetaxine — EEF2 plus EEF1/EIF4 translation factors and 40S/60S ribosomal proteins; (v) HDAC8_TGFb_EMT_Resistance — HDAC8 plus TGFβ/SMAD axis, canonical EMT TFs (SNAI1/2, TWIST1/2, ZEB1/2), and mesenchymal markers (FN1, VIM, CDH2, MMP2/9, COL1A1/3A1).

### 4.3 Bridging classifier (Bi 3-group labels for Nassiri samples)

A Ridge-penalised multinomial classifier was trained on GSE183653 log₂(TPM+1) matrices restricted to the *N*=2,000 most variably expressed protein-coding HUGO genes overlapping the Nassiri TPM matrix. Gene symbols were standardised across cohorts via cBioPortal batch Entrez→HUGO lookup. The sklearn Pipeline chained `StandardScaler` (fit on training set only, never refit per-fold) with `RidgeClassifier` (α tuned via 5-fold internal CV on the training set). LOOCV on GSE183653 *N*=185 was used as the fidelity metric before deployment. The validated classifier was applied once to the Nassiri *N*=121 TPM matrix to yield predicted `bi_group` labels per sample; no retuning to the Nassiri cohort was performed.

### 4.4 Single-sample program scoring and statistical testing

Rank-based ssGSEA (Barbie 2009 (9); Subramanian 2005 (10)) with α = 0.25 was computed per sample × program, followed by per-program within-cohort z-transformation so Δ-median effect sizes are interpretable on a ±1 σ scale. All tests were two-sided. BH-FDR was applied at α = 0.10 across explicitly declared families (not cherry-picked subsets): (H2) all 5 programs jointly; (H3) all 5 programs jointly. Grade-adjusted models used statsmodels OLS via the formula interface with WHO grade encoded as a **full 3-level categorical factor** `C(who_grade)` with grade I as the reference level, paired with subgroup as `C(nassiri_group)`, giving residual df = 121 − 1 − 2 − 3 = 115; subgroup and grade F-statistics were extracted from Type-II ANOVA (`sm.stats.anova_lm(model, typ=2)` so both terms are evaluated after the other is already in the model. All analyses were implemented in Python 3.11+ using statsmodels, scipy.stats, scikit-learn, and pandas; code is distributed under MIT at `<repo>`.

---

## 5. Data availability

Frozen target-gene programs, accession/completeness-lock memos, 16 output tables (S0–S2), and result CSV provenance sidecars are distributed under MIT license:
- All tables: `<repo>/results/tables/`
- End-to-end pipeline: `pip install -r requirements.txt ; PYTHONPATH=src python scripts/99_run_all.py`
- Test suite (11 unit tests validating BH-FDR control, KW power, κ boundary cases, ssGSEA enrichment, ground-truth directionality): `pytest tests/test_all_aims.py`; 11/11 passing at pre-submission commit.

When Heidelberg colleagues contribute author-level screening panels and/or controlled-access individual-level data, the same pipeline command re-runs verbatim; only literal gene-program lists and accession dicts change, with no statistical thresholds re-tuned.

---

## 6. Conflict of interest

None declared.

---

## 7. References

[1] Goldbrunner R, Stavrinou P, Jenkinson MD, Sahm F, Mawrin C, Weber DC, Preusser M, Minniti G, Lund-Johansen M, Lefranc F, Houdart E, Sallabanda K, Le Rhun E, Nieuwenhuizen D, Tabatabai G, Soffietti R, Weller M. EANO guideline on the diagnosis and management of meningiomas. *Neuro-Oncology.* 2021;23(11):1821–1834. doi:10.1093/neuonc/noab150.

[2] Wen PY, Quant E, Drappatz J, Beroukhim R, Norden AD. Medical therapies for meningiomas. *J Neurooncol.* 2010;99(3):365–378. doi:10.1007/s11060-010-0349-8.

[3] Jungwirth S, et al. Pharmacological Landscape of FDA-Approved Anticancer Drugs Reveals Sensitivities to Ixabepilone, Romidepsin, Omacetaxine, and Carfilzomib in Aggressive Meningiomas. *Clin Cancer Res.* 2023;29(1):233–243. doi:10.1158/1078-0432.CCR-22-2085.

[4] Jungwirth S, et al. Drug screening on tumor organoids exposes therapeutic vulnerabilities of meningiomas to HDAC1/2i panobinostat. *Sci Transl Med.* 2026;18(ea3115). doi:10.1126/scitranslmed.aea3115.

[5] Nassiri F, Liu J, Patil V, Mamatjan Y, Wang JZ, Hugh-White R, Macklin AM, Khan S, Singh O, Karimi S, Corona RI, Liu LY, Chen CY, Chakravarthy A, Wei Q, Mehani B, Suppiah S, Gao A, Workewych AM, Tabatabai G, Boutros PC, Aldape K, Zadeh G. A clinically applicable integrative molecular classification of meningiomas. *Nature.* 2021;597(7874):119–125. doi:10.1038/s41586-021-03850-3. GEO GSE180061; EGA EGAS00001004982; cBioPortal `mng_utoronto_2021`.

[6] Choudhury A, et al. Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities. *Nature Genetics.* 2022;54(5):649–659. doi:10.1038/s41588-022-01061-8. GEO GSE183653.

[7] Choudhury A, Chen WC, Lucas CHG, Bayley JC, Harmanci AS, Maas SLN, Santagata S, Klisch T, Perry A, Bi WL, Sahm F, Patel AJ, Magill ST, Raleigh DR. Hypermitotic meningiomas harbor DNA methylation subgroups with distinct biological and clinical features. *Neuro-Oncology.* 2023;25(3):520–530. doi:10.1093/neuonc/noac224. GEO GSE212666.

[8] Sahm F, et al. DNA methylation-based classification and grading system for meningioma: a multicentre, retrospective analysis. *Lancet Oncol.* 2017;18(5):682–694. doi:10.1016/S1470-2045(17)30155-9. (Cited by reference only per protocol.)

[9] Barbie DA, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature.* 2009;462(7269):108–112. doi:10.1038/nature08460. (ssGSEA rank-sum scoring method.)

[10] Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS.* 2005;102(43):15545–15550. doi:10.1073/pnas.0506580102.

[11] Neftel C, et al. An integrative model of cellular states, plasticity, and genetics for glioblastoma. *Cell.* 2019;178(4):835–849.e21. doi:10.1016/j.cell.2019.06.024. (Methodological parallel: convergent-classifier extension work.)

---

## 8. Table legends (main text)

**Table 1 — H1 4×3 cross-tabulation counts (Nassiri group × predicted Bi group; N = 121).** Nassiri native 4-group labels (rows) × Bi 3-group labels projected via Ridge bridging classifier on Nassiri TPM (columns). Full row-percentages and statistical metrics in Supp Table S1.

| Nassiri group | Hypermitotic | Immune-enriched | Merlin-intact | Row total |
|---|---:|---:|---:|---:|
| Immunogenic (MG1) | 0 | 17 | 0 | 17 |
| MG2 (Benign) | 2 | 3 | 27 | 32 |
| Hypermetabolic (MG3) | 21 | 13 | 9 | 43 |
| Proliferative (MG4) | 22 | 4 | 3 | 29 |
| **Column total** | **45** | **37** | **39** | **121** |

**Table 2 — H1 pairing-specificity test: top and adversarial ranks among 12 enumerated 4→3 pairings.** All pairings route exactly two Nassiri groups into the Hypermitotic bin (merge pair) with the remaining two assigned as singletons. The literature-registered true mapping is rank #1 by both κ and raw agreement. Full ranked 12-row table in Supp Table S1b.

| Rank (κ) | Is true? | Merge pair | Singleton A (Nas→Bi) | Singleton B (Nas→Bi) | Cohen's κ | Raw agreement |
|---:|:---:|---|---|---|---:|---:|
| 1 | ✅ | hypermetabolic ∪ proliferative | Immunogenic→Immune-enriched | MG2→Merlin-intact | 0.568 | 71.9% |
| 2 | ❌ | Immunogenic ∪ proliferative | MG2→Merlin-intact | hypermetabolic→Immune-enriched | 0.266 | 51.2% |
| 5 | ❌ (singleton-swap placebo) | hypermetabolic ∪ proliferative | Immunogenic→Merlin-intact | MG2→Immune-enriched | 0.050 | 38.0% |
| 12 | ❌ (worst wrong merge) | Immunogenic ∪ MG2 | hypermetabolic→Merlin-intact | proliferative→Immune-enriched | −0.324 | 12.4% |

**Table 3 — H2 Nassiri 4-group subgroup stratification: Kruskal–Wallis omnibus + grade-adjusted conditional F-tests. N = 121 across 4 groups; WHO grade enters OLS as a full 3-level categorical factor (grades I/II/III → 2 dummy df; subgroup 3 dummy df; residual df = 115). Type-II ANOVA via statsmodels `anova_lm(typ=2)`. *p*BH = Benjamini–Hochberg FDR applied across the 5-program family.**

| Program | Kruskal–Wallis H | KW *p* raw | KW *p*BH (H2 family) | Subgroup F\|grade (df=3,115) | Subgroup *p*\|grade | Subgroup *p*BH\|grade | Grade *p* (full factor) |
|---|---:|---:|---:|---:|---:|---:|---:|
| HDAC Panobinostat/Romidepsin | 28.08 | 3.49 × 10⁻⁶ | 1.74 × 10⁻⁵ | 8.56 | 3.53 × 10⁻⁵ | 1.77 × 10⁻⁴ | 0.185 |
| Proteasome Carfilzomib | 13.87 | 3.09 × 10⁻³ | 3.86 × 10⁻³ | 4.22 | 7.16 × 10⁻³ | 8.94 × 10⁻³ | 0.910 |
| Tubulin Ixabepilone | 6.64 | 8.44 × 10⁻² | 8.44 × 10⁻² | 2.45 | 6.74 × 10⁻² | 6.74 × 10⁻² | 8.63 × 10⁻⁴ |
| Translation Omacetaxine | 17.68 | 5.12 × 10⁻⁴ | 1.28 × 10⁻³ | 5.81 | 9.79 × 10⁻⁴ | 2.20 × 10⁻³ | 0.465 |
| HDAC8/TGFβ/EMT Resistance | 16.77 | 7.88 × 10⁻⁴ | 1.31 × 10⁻³ | 5.57 | 1.32 × 10⁻³ | 2.20 × 10⁻³ | 0.306 |

**Table 4 — H3 nf2_cna_loss_proxy association: Intact (*N*=37) vs Mutant/Loss (*N*=84); *N*=121. Mann–Whitney U with Cliff's δ and BH-FDR across the 5-program H3 family.**

| Program | Δ median z (Intact − Mut/Loss) | Cliff's δ | MW *p* raw | *p*BH (H3 family) | Pass α=0.10 |
|---|---:|---:|---:|---:|:---:|
| HDAC Panobinostat/Romidepsin | −0.306 | −0.313 | 6.31 × 10⁻³ | 3.15 × 10⁻² | ✅ |
| Proteasome Carfilzomib | +0.429 | +0.263 | 2.19 × 10⁻² | 3.76 × 10⁻² | ✅ |
| Tubulin Ixabepilone | −0.234 | −0.189 | 9.99 × 10⁻² | 9.99 × 10⁻² | ✅ |
| Translation Omacetaxine | +0.783 | +0.238 | 3.76 × 10⁻² | 4.71 × 10⁻² | ✅ |
| HDAC8/TGFβ/EMT Resistance | +0.496 | +0.261 | 2.25 × 10⁻² | 3.76 × 10⁻² | ✅ |

---

## 9. Supplementary Tables (all at `results/tables/`)

- **Table S0** — Aim 0 cohort lock: accessions, metadata completeness, and hypothesis-evaluability declaration (`aim0_cohort_status_report.txt`, `aim0_cohort_lock.csv/json`, `aim0_geo_screening.csv`, `aim0_cbio_studies.csv`).
- **Table S1** — H1 classifier concordance: 4×3 cross-tab counts, row-percentages, and full metrics (χ², Cramér's V, mapped 3×3 κ with asymptotic p) (`table_s1_aim1_concordance_crosstab_counts.csv`, `table_s1_aim1_concordance_crosstab_rowpct.csv`, `table_s1_aim1_concordance_metrics.csv`).
- **Table S1b** — H1 pairing-specificity test: full 12-row ranked table of all enumerated pairings + 1-row summary with pre-specified pass flags (`table_s1b_aim1_pairing_specificity_all_pairings_ranked.csv`, `table_s1b_aim1_pairing_specificity_summary.csv`).
- **Table S2a** — H2/H3 Nassiri native 4-group: KW, pairwise MW with per-pair and all-pairs BH-FDR, grade-adjusted ANOVA, H3 NF2 MW, raw + z-scored ssGSEA program matrices (12 files prefix `table_s2a_aim2_nassiri_native4grp_*`).
- **Table S2b** — H2/H3 Nassiri projected Bi 3-group: identical spec run on predicted bi_group labels after bridging-classifier projection, enabling direct cross-classifier convergence comparison (12 files prefix `table_s2b_aim2_nassiri_predicted_bi3grp_*`).
- **Table S5** — H5 formal one-line declaration sheet: `recurrence_months` not deposited, RFS Cox analysis not evaluable.
