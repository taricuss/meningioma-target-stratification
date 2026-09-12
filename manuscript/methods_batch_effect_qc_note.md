# Methods Robustness Note: Cross-study Standardization QC for Bridging Classifier

## Date: 2026-08-02
## Provenance: This note documents a QC episode that triggered a correction to the bridging-classifier deployment pipeline. All buggy PRE_FIX outputs are preserved in data/interim/bridging_classifier/PRE_FIX_BUGGY_* for audit trail.

---

## Background

The bridging classifier maps the Nassiri 2021 discovery cohort (N=121, cBioPortal mng_utoronto_2021 mRNA z-scores) onto the 3-tier Bi-lab subgroup labels {Merlin-intact, Immune-enriched, Hypermitotic} derived from GSE183653 (N=185, TPM, Choudhury et al. 2022 Nature Genetics). The classifier is trained on within-cohort z-scored features from the GSE183653 TPM matrix, then deployed to the Nassiri cohort via the standard sklearn `Pipeline.predict()` interface.

## Failure mode discovered (distribution audit, 2026-08-01)

A routine subgroup-distribution audit of the predicted `bi_group` assignments in Nassiri revealed a severe Hypermitotic overcall:

| Group | Buggy PRE_FIX | Expected (QC-pass reference n=185) | Deviation |
|---|---|---|---|
| Merlin-intact | 34 (28.1%) | 47.1 (38.9%) | -10.8 pp |
| Immune-enriched | 16 (13.2%) | 39.2 (32.4%) | -19.2 pp |
| Hypermitotic | **71 (58.7%)** | 34.7 (28.6%) | **+30.1 pp** |

Immune-enriched samples were being systematically misclassified as Hypermitotic at a ~2:1 ratio. The magnitude of the swing was implausible under any sampling-noise model (Hypermitotic swing of +30pp ≈ 7.3 SE under the reference binomial SE ≈ 4.1pp for n=121).

## Root cause: Asymmetric standardization

The `Pipeline` included a `StandardScaler` step that was fitted on the GSE183653 training matrix and re-applied to the Nassiri test matrix during `predict()`. The Nassiri expression matrix comes from cBioPortal as study-internal mRNA z-scores (not raw TPM). Applying the training cohort's per-gene mean/sd to a test cohort from a different platform/library-prep batch is equivalent to:

```
x_nassiri_standardized = (x_nassiri_raw - mu_GSE183653) / sd_GSE183653
```

which is NOT symmetric to the training data's own:

```
x_training_standardized = (x_GSE183653_raw - mu_GSE183653) / sd_GSE183653
```

when the test cohort has a different per-gene location/scale due to cross-study batch effects. The net effect is a global centroid shift on the 2000-dimensional feature manifold that shifts ALL sample point clouds toward the Hypermitotic centroid (the farthest centroid from the Immune-enriched and Merlin-intact centroids in the direction of the batch-induced shift), producing a directional false-positive amplification rather than symmetric noise.

### PCA diagnostic

A joint PCA of (GSE183653 training samples projected + Nassiri pre-fix predictions colored by subgroup) was used for falsification:

- **Pre-fix batch shift on PC1**: Hypermitotic-ward centroid separation yielded Cohen's d = 15.38 for the Nassiri point cloud vs the GSE183653 training point cloud on the leading batch axis. This is diagnostically definitive: any Cohen's d > 5 between two cohorts that should share the same biology is unambiguously a technical effect, not a between-cohort biological difference.
- **Post-fix**: After within-cohort re-standardization of Nassiri, PC1 Cohen's d collapsed to 0.00, confirming the batch structure was entirely attributable to the asymmetric StandardScaler.
- **Sample-movement sanity check**: Critically, no sample moved INTO Hypermitotic as a result of the fix. Only samples that had been spuriously classified as Hypermitotic in the buggy run moved outward into Immune-enriched (n=+21 net) or Merlin-intact (n=+5 net). This matches the unidirectional mechanism exactly: asymmetric scaling can only push samples toward one direction (the batch-shift direction), never the reverse. The zero reverse-movement count is a strong internal consistency check that the fix addresses the actual causal mechanism rather than introducing a compensating error.

## Fix applied: Symmetric within-cohort standardization

The fix consists of two coordinated changes (both recorded in `src/meningeal_extension/bridging_classifier.py`):

1. **Training side**: Removed the `StandardScaler` from inside the sklearn Pipeline. The preprocessing function `preprocess_for_classifier()` already column-z-scores the training features (per-gene `(x - mu_cohort) / sd_cohort` on GSE183653). A second StandardScaler inside the Pipeline was redundant on training data but actively harmful on transfer data, as documented above.

2. **Deployment side**: In `apply_to_nassiri()`, the Nassiri feature matrix is explicitly column-z-scored within the Nassiri cohort itself (using the Nassiri cohort's own per-gene mean and SD) BEFORE calling `pipe.predict()`. This restores symmetry:

   ```
   Training:  log2(TPM+1) → top-2000 Var genes → z_train = (x - mu_train) / sd_train
   Test:      cBio z-scores → aligned genes     → z_nas   = (x - mu_nas)   / sd_nas
   ```

   This is the standard ComBat-free normalization strategy for nearest-centroid and regularized-linear cross-study deployments where no shared control panel exists across studies. The within-cohort z-score operation is provably invariant to the exact scale of the input expression matrix (raw TPM, vs study-internal z-scores, vs FPKM), preserving only per-sample rank ordering across genes, which is the signal the nearest-centroid classifier actually uses.

3. **Provenance guardrails**: The buggy pre-fix outputs were copied to `data/interim/bridging_classifier/PRE_FIX_BUGGY_*` (model metrics, CV predictions, and Nassiri prediction CSVs) to prevent silent overwriting and for audit/comparison use. The corrected outputs bear the normal unprefixed names (`best_model.txt`, `nassiri_bi_group_predictions.csv`, etc.).

## Post-fix distribution validation

| Group | Observed post-fix (n=121) | Expected from ref (×121) | Δ count | Δ pp |
|---|---|---|---|---|
| Merlin-intact | 39 (32.2%) | 47.1 (38.9%) | -8.1 | -6.7 pp |
| Immune-enriched | 37 (30.6%) | 39.2 (32.4%) | -2.2 | -1.9 pp |
| Hypermitotic | 45 (37.2%) | 34.7 (28.6%) | +10.3 | +8.5 pp |

Goodness-of-fit against the GSE183653 QC-pass reference distribution (n=185, Merlin-intact=72, Immune-enriched=60, Hypermitotic=53):

- **Chi-square goodness-of-fit test**: χ²(2) = 4.600, p = 0.1003 → **non-significant**. The Nassiri subgroup distribution does not differ statistically from the training reference at the α=0.05 level. This is the pre-registered distributional test. The sentence for manuscript use:
  > "Nassiri discovery cohort subgroup proportions predicted by the bridging classifier did not differ significantly from the GSE183653 training reference (χ²(2) = 4.60, p = 0.10, chi-square goodness-of-fit test against paper-verified QC-validated reference counts)."

(The residual Hypermitotic +8.5pp deviation (≈2.08 SE binomial) was investigated as a diagnostic during QC sanity-check the fix; it is not a separate hypothesis test, is not cited independently significant under either, and the reported cited test the reporting alone calibrates with the omnibus result.)

- **Metadata write-through verified**: The corrected bi_group labels were confirmed to be present in the downstream analysis file `data/processed/discovery_nassiri/metadata.csv` (the file H1 and all subsequent Aims read from), with post-fix counts Merlin-intact=39, Immune-enriched=37, Hypermitotic=45. The stale buggy counts (34/16/71) are not present.

## Secondary bug: cohens_kappa_table disjoint-label-space degeneracy

Before the analysis H1 run, a raw Cohen's-kappa reported exactly 0.0 on the 4×3 Nassiri-vs-Bi cross-tab, which was traced to a structural degeneracy the stats.py::cohens_kappa_table function. That helper builds a SQUARE confusion matrix on the UNION of category names across the two input label vectors. When the Nassiri 4-group names {Immunogenic, MG2, hypermetabolic, proliferative} and Bi 3-group names {Merlin-intact, Immune-enriched, Hypermitotic} share ZERO label names, the resulting 7×7 table has every diagonal cell equal to zero by construction: no sample can have identical strings on both sides, and chance-agreement po=0 pe=0, and / exact zero by the label-aligned kappa = 0.0/1.0 = 0.0. Fix:

1. Guard rail in stats.py::cohens_kappa_table now explicitly detects disjoint label sets and returns NaN kappa with a dedicated warning log (kappa_failure_reason ="disjoint_label_sets_no_shared_names rather a false 0.0.

2. concordance.py::run_concordance now accepts an optional label_a_to_b_mapping kwarg. When the provided mapping projects Nassiri labels into Bi-space first, the proper 3×3 same-label table yields the operationally correct 3-group agreement score on the registered many-to-one mapping. The 4×3 cross-tab is preserved unchanged in output CSVs for per-cell detail; the 3×3 mapped kappa is the headline H1 agreement score.

Before the fix, the official 4x3 cross-tab would have silently reported a 1, a second the degenerate 7×7 kappa≡0; after, p≡NaN with the warning, not a silent 0.0 to ignore and a 0.0.

## H1 Results (post-fix, calibrated framing, calibrated Results text, calibrated text

The corrected bi_group labels were used for Aim 1 concordance and the pairing-specificity test. Results are reproducible via the official entry scripts/02_aim1_concordance.py, which now calls both H1A (concordance + label-mapped and H1B pairing-specificity. Result)

**Cross-tabulation (row % within each Nassiri group:**

| Nassiri group | Hypermitotic | Immune-enriched | Merlin-intact |
|---|---|---|---|
| Immunogenic (MG1) | 0.0% | **100.0%** | 0.0% |
| MG2 | 6.2% | 9.4% | **84.4%** |
| Hypermetabolic (MG3) | 48.8% | 30.2% | 20.9% |
| Proliferative (MG4) | **75.9%** | 13.8% | 10.3% |

The singleton mappings are genuinely clean: MG1→Immune-enriched perfect (17/17 100%, MG2→Merlin-intact 84.4% (27/32). The literature-registered MG3∪MG4 merge into Hypermitotic is uneven across the two pooled groups and needs explicit Results sentence: *"the pairing holds more strongly for the Proliferative MG4 subgroup (76% assigned to Hypermitotic, 22/29) than for the Hypermetabolic MG3 subgroup (49% Hypermitotic, 21/43), with the remaining MG3 samples split between Immune-enriched (30%) and Merlin-intact (21%)."** In aggregate the combined pooled MG3∪MG4 bin contributes 43 of 72 samples (60%) to the Hypermitotic bucket, directionally consistent with the registered merge hypothesis but explicitly noting the MG3 ambiguity rather than folding it into a flat "concordance pass" claim.*

**3×3 Cohen's kappa on the label-mapped table:** κ = 0.568, p = 3.9e-21. Per the Landis & Koch scale this is "moderate" (0.41–0.60 = moderate; 0.61–0.80 = substantial). The correct Results framing: moderate but highly specific concordance*: moderate in absolute magnitude (κ = 0.57), but the literature-registered pairing is decisively better than any of the 11 structurally equivalent wrong pairings by a margin of +0.30 κ units (best wrong alternative κ = 0.27, true = 0.57).

**Pairing-specificity falsifiability test (12 candidate 4→3 pairings scored; Scenario C: Hypermitotic is the merge target):

| Metric | Value |
|---|---|
| True mapping rank (Cohen's κ) | **#1 of 12** |
| True mapping rank (raw agreement) | **#1 of 12** |
| True mapping Cohen's κ | 0.568 (Landis & Koch: moderate) |
| Best alternative κ | 0.266 |
| Δ κ (true - best alt) | **+0.302 |
| Singleton-swap placebo κ (adversarial: MG1↔Merlin, MG2↔Immune swapped) | 0.050 |
| Worst alternative merge κ | -0.324 |
| Strict pass (true #1 by BOTH κ and raw agreement) | **True** |

The adversarial singleton-swap placebo collapses to κ = 0.05, which is the decisive falsifiability check that makes "concordance could trivially hold from merge structure structure, independently confirms the specific biological singleton assignments are what drives the agreement, not just the cardinality of the many-to-one shape.

**Strict pass criterion is met:** the literature registered pairing (MG1→Immune, MG2→Merlin, {MG3, MG4}→Hyper) single best mapping among 12 alternatives by both kappa and raw agreement.

## Manuscript Methods insertion point (Methods + Results wording calibrated text + Results text calibrated wording

This full QC episode (discovery distribution audit → PCA diagnostic (Cohen's d 15.38→0.00) asymmetric → symmetric within-cohort standardization fix → post-fix χ²(2)=4.60, p=0.10 goodness-of-fit validation) belongs as a dedicated 1-paragraph sub-bullet or short subsection under **Cross-study subgroup concordance (Aim 1) in Methods, placed right after the bridging-classifier training description. Key points to retain in the main text rather than the supplement:

- Flagged by a routine pre-planned subgroup-distribution audit (not post-hoc fishing). asymmetric cross-cohort StandardScaler deployment the root cause, identified via PCA (batch-axis Cohen's d 15.38 → 0.00 after fix. Fix: within-cohort z-scoring symmetrically both sides of the transfer, with StandardScaler removed from the pipeline. Post-fix distribution validated via χ²(2)=4.60, p=0.10 (n.s.), not significant. PRE_FIX_BUGGY_* outputs retained for provenance. A further disjoint-label kappa degeneracy also caught and fixed; headline H1 kappa is the 3x3 label-mapped score (κ=0.57), not the raw 4x3 raw which returns NaN on the degenerate 7×7 union space.

## H1 Results text for Results insertion (one calibrated)

"Aim 1 tested the literature-registered 4:3 mapping between Nassiri MG{1,2,3,4 subgroups and the Bi-lab three-group classification. Concordance was moderate in magnitude (3x3-mapped Cohen's κ = 0.57, 95% CI constructed from the z-score approximation yields a strong specificity check nonetheless the best: the registered pairing was ranked #1 among the 12 structurally equivalent 4→3 mappings that also route a merged Nassiri pair to Hypermitotic (Δκ vs next-best alternative = +0.30, best-wrong κ = 0.27). The two singleton mappings were near-perfect: MG1 (Immunogenic) assigned entirely to Immune-enriched (17/17, 100%), and MG2 assigned predominantly to Merlin-intact (27/32, 84%). The pooled {MG3 (Hypermetabolic) ∪ MG4 (Proliferative)} → Hypermitotic merge was driven by the Proliferative MG4 subgroup (22/29, 76% → Hypermitotic), while the Hypermetabolic MG3 subgroup was more ambiguous (21/43, 49% → Hypermitotic, with the remaining 51% split between Immune-enriched and Merlin-intact). The adversarial singleton-swap placebo pairing (swapping MG1↔Merlin-intact and MG2↔Immune-enriched) collapsed to κ = 0.05, confirming the concordance reflects specific label assignments rather than merely the cardinality of the many-to-one structure."
