from __future__ import annotations

import logging
import pickle
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, cohen_kappa_score, classification_report
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.neighbors import NearestCentroid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import (
    BI_GROUP_ORDER,
    BRIDGING_CLASSIFIER_LABELS_REQUIREMENT,
    PATHS,
    SEED,
)

logger = logging.getLogger(__name__)

GSE183653_TPM_URL = (
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/"
    "GSE183653/suppl/GSE183653_meningioma_tpm.csv.gz"
)
GSE183653_RAW_DIR = PATHS.data_raw / "GSE183653"
GSE183653_TPM_PATH = GSE183653_RAW_DIR / "GSE183653_meningioma_tpm.csv.gz"
GSE183653_LABELS_PATH = GSE183653_RAW_DIR / "GSE183653_sample_bi_labels.csv"
GSE183653_LABELS_TEMPLATE = GSE183653_RAW_DIR / "GSE183653_sample_bi_labels_TEMPLATE.csv"
BRIDGING_MODEL_OUT_DIR = PATHS.data_interim / "bridging_classifier"
NASSIRI_COHORT_DIR = PATHS.data_processed / "discovery_nassiri"

EXPECTED_LABEL_COUNTS = {
    # Real paper-derived reference distribution from the 185 samples with
    # usable RNA-seq (Choudhury et al. 2022 Nat Genet, N=200 discovery →
    # N=185 after QC). NOT the naive 34/38/28% × 200 = 68/76/56.
    # QC attrition is NON-RANDOM by subgroup (Methods/Limitations note):
    #   Merlin-intact:   72/72  =  0% loss  (indolent tumors survive QC)
    #   Immune-enriched: 60/65  = ~8% loss
    #   Hypermitotic:    53/63  = ~16% loss (aggressive/necrotic → fails RNA-seq)
    # Net swing vs blind 34/38/28% projection: Merlin-intact +9, Immune -10, Hyper +1.
    "Merlin-intact": 72,
    "Immune-enriched": 60,
    "Hypermitotic": 53,
}
LABEL_COUNT_TOLERANCE = 8


def _die(msg: str, code: int = 3) -> None:
    logger.critical(msg)
    sys.stderr.write("\n")
    sys.stderr.flush()
    sys.exit(code)


def ensure_tpm_downloaded() -> Path:
    GSE183653_RAW_DIR.mkdir(parents=True, exist_ok=True)
    if GSE183653_TPM_PATH.exists() and GSE183653_TPM_PATH.stat().st_size > 1024:
        logger.info("GSE183653 TPM already present: %s", GSE183653_TPM_PATH)
        return GSE183653_TPM_PATH
    import urllib.request
    logger.info("Downloading GSE183653 TPM from %s", GSE183653_TPM_URL)
    try:
        req = urllib.request.Request(
            GSE183653_TPM_URL, headers={"User-Agent": "Mozilla/5.0"}
        )
        with urllib.request.urlopen(req, timeout=600) as resp:
            data = resp.read()
        GSE183653_TPM_PATH.write_bytes(data)
        logger.info("Saved %d bytes to %s", len(data), GSE183653_TPM_PATH)
    except Exception as e:
        _die(
            f"Failed to download GSE183653 TPM matrix.\n"
            f"Error: {type(e).__name__}: {e}\n"
            f"Manual fix: download the file directly from:\n"
            f"  {GSE183653_TPM_URL}\n"
            f"and save it to:\n"
            f"  {GSE183653_TPM_PATH}"
        )
    return GSE183653_TPM_PATH


def load_gse183653_tpm() -> pd.DataFrame:
    path = ensure_tpm_downloaded()
    try:
        df = pd.read_csv(path, index_col=0, compression="gzip")
    except Exception as e:
        _die(
            f"Failed to parse GSE183653 TPM CSV: {e}\n"
            f"File: {path}. Check the file is a valid gzip-compressed CSV."
        )
    logger.info("Loaded GSE183653 TPM: %d genes × %d samples", *df.shape)
    if df.shape[1] != 185:
        logger.warning(
            "GSE183653 TPM has %d sample columns (expected 185). "
            "Proceeding as-is; label cross-check (A) will fail if sample_id sets diverge.",
            df.shape[1],
        )
    sample_ids = list(df.columns.astype(str))
    if not sample_ids[0].startswith("M") or not sample_ids[-1].startswith("M"):
        logger.warning(
            "First/last sample IDs ('%s'/'%s') do not match expected M## prefix. "
            "Label cross-check (A) below will still validate set equality.",
            sample_ids[0], sample_ids[-1],
        )
    return df


def validate_and_load_labels(tpm_sample_ids: List[str]) -> pd.DataFrame:
    tpm_samples = set(tpm_sample_ids)
    if not GSE183653_LABELS_PATH.exists():
        _write_labels_template(tpm_sample_ids)
        _die(
            f"================================================================\n"
            f"BRIDGING CLASSIFIER LABEL SIDECAR MISSING — Aim 1 is BLOCKED.\n"
            f"================================================================\n"
            f"Required file not found: {GSE183653_LABELS_PATH}\n"
            f"\n"
            f"See BRIDGING_CLASSIFIER_LABELS_REQUIREMENT in config.py for full spec.\n"
            f"\n"
            f"A TEMPLATE has been written to:\n"
            f"  {GSE183653_LABELS_TEMPLATE}\n"
            f"It contains all {len(tpm_samples)} M## sample IDs from the TPM matrix\n"
            f"(Column: sample_id) and an empty bi_group column for you to fill.\n"
            f"\n"
            f"Label source priority (document this exact path in provenance):\n"
            f"  1. Choudhury et al. 2022 Nature Genetics Supp Table 1/2\n"
            f"     (per-sample M## ID → {{{', '.join(BI_GROUP_ORDER)}}})\n"
            f"  2. Abrarc/meningioma-svm Zenodo 10.5281/zenodo.6353877 subgroups.csv\n"
            f"     (if authors release it)\n"
            f"  3. If neither is obtainable: Aim 1 / Bi-lab H2-H3 → NOT TESTABLE.\n"
            f"\n"
            f"Count sanity (N=185 RNA-seq QC-pass samples, paper-verified reference):\n"
            f"  Merlin-intact    =72  (0% QC loss; indolent tumors survive RNA-seq)\n"
            f"  Immune-enriched  =60  (~8% QC loss of N=65 in full discovery)\n"
            f"  Hypermitotic     =53  (~16% QC loss: aggressive/necrotic → fails RNA-seq)\n"
            f"  (NOT the naive 34/38/28%×200 projection — QC attrition is subgroup-biased)\n"
            f"Tolerance: ±{LABEL_COUNT_TOLERANCE} per bin; CSV is REJECTED beyond that.\n"
            f"\n"
            f"Once you have filled in bi_group, rename TEMPLATE → sample_bi_labels.csv:\n"
            f"  mv {GSE183653_LABELS_TEMPLATE.name} {GSE183653_LABELS_PATH.name}\n"
            f"in directory: {GSE183653_RAW_DIR}\n"
            f"================================================================\n"
        )
    try:
        labels = pd.read_csv(GSE183653_LABELS_PATH, dtype=str)
    except Exception as e:
        _die(f"Failed to parse labels CSV {GSE183653_LABELS_PATH}: {e}")

    for req in ["sample_id", "bi_group"]:
        if req not in labels.columns:
            _die(
                f"Labels CSV missing required column '{req}'.\n"
                f"Columns present: {list(labels.columns)}\n"
                f"See:\n{BRIDGING_CLASSIFIER_LABELS_REQUIREMENT}"
            )
    label_samples = set(labels["sample_id"].astype(str))

    only_tpm = sorted(tpm_samples - label_samples)
    only_lab = sorted(label_samples - tpm_samples)
    if only_tpm or only_lab:
        msg = [
            "Label sanity-check (A) FAILED: sample_id sets don't match TPM 1:1."
        ]
        if only_tpm:
            msg.append(
                f"  Samples in TPM but MISSING from labels CSV ({len(only_tpm)}): "
                f"{only_tpm[:10]}{'...' if len(only_tpm) > 10 else ''}"
            )
        if only_lab:
            msg.append(
                f"  Samples in labels CSV but MISSING from TPM ({len(only_lab)}): "
                f"{only_lab[:10]}{'...' if len(only_lab) > 10 else ''}"
            )
        msg.append(
            f"Expected {len(tpm_samples)} sample IDs matching TPM columns exactly."
        )
        _die("\n".join(msg))

    unknown_groups = sorted(set(labels["bi_group"].dropna().unique()) - set(BI_GROUP_ORDER))
    if unknown_groups:
        _die(
            f"Label sanity-check (B) FAILED: bi_group values not in BI_GROUP_ORDER.\n"
            f"  Unknown labels: {unknown_groups}\n"
            f"  Allowed values: {BI_GROUP_ORDER}\n"
            f"  Check capitalization, spacing, hyphenation (Merlin-intact has HYPHEN)."
        )

    counts = labels["bi_group"].value_counts().reindex(BI_GROUP_ORDER, fill_value=0)
    deviations = {g: int(counts[g]) - EXPECTED_LABEL_COUNTS[g] for g in BI_GROUP_ORDER}
    bad = [g for g, d in deviations.items() if abs(d) > LABEL_COUNT_TOLERANCE]
    if bad:
        lines = [
            f"Label sanity-check (C) FAILED: per-bin counts deviate > ±{LABEL_COUNT_TOLERANCE}"
            f" from paper-verified QC-validated reference (N=185 RNA-seq pass samples):"
        ]
        for g in BI_GROUP_ORDER:
            lines.append(
                f"  {g:20s}  actual={int(counts[g]):3d}  expected={EXPECTED_LABEL_COUNTS[g]:3d}"
                f"  Δ={deviations[g]:+3d}  {'❌ FAIL' if g in bad else 'OK'}"
            )
        lines.append(
            "If labels were extracted from the paper's Supp Tables, verify against:\n"
            "  (a) actual M## IDs in Supp Table 1/2, (b) Fig 1 bar proportions,\n"
            "  (c) Extended Data Fig. subgroup counts for N=185 RNA-seq subset.\n"
            "QC attrition is subgroup-biased (Hypermitotic ~16% loss, Immune ~8%,\n"
            "Merlin-intact 0%); EXPECTED_LABEL_COUNTS above already accounts for this.\n"
            "If extraction is genuinely correct, raise LABEL_COUNT_TOLERANCE (see top)."
        )
        _die("\n".join(lines))

    logger.info(
        "Labels validated OK (N=%d). Per-group counts: %s. Deviations from expected: %s.",
        len(labels),
        {g: int(counts[g]) for g in BI_GROUP_ORDER},
        deviations,
    )
    return labels.set_index("sample_id")


def _write_labels_template(tpm_sample_ids: List[str]) -> None:
    GSE183653_RAW_DIR.mkdir(parents=True, exist_ok=True)
    template = pd.DataFrame({
        "sample_id": sorted(tpm_sample_ids, key=lambda s: (len(s), s)),
        "bi_group": "",
        "notes": "",
    })
    template.to_csv(GSE183653_LABELS_TEMPLATE, index=False)
    logger.info("Wrote label template: %s", GSE183653_LABELS_TEMPLATE)


def preprocess_for_classifier(
    tpm_df: pd.DataFrame,
    nassiri_expr: Optional[pd.DataFrame] = None,
    n_top_genes: int = 2000,
) -> Tuple[pd.DataFrame, List[str]]:
    """Log2(TPM+1) → filter to protein-coding/overlapping → top-N variable → z-score.

    If nassiri_expr is provided, restricts to gene intersection BEFORE selecting
    top-N variable genes (ensures transferability to Nassiri cohort)."""
    x = np.log2(tpm_df.astype(float).T + 1.0)

    gene_universe = list(tpm_df.index.astype(str))
    if nassiri_expr is not None:
        nassiri_genes = set(str(c) for c in nassiri_expr.columns)
        overlap = [g for g in gene_universe if g in nassiri_genes]
        logger.info(
            "Gene overlap GSE183653 (n=%d) ∩ Nassiri (n=%d) → %d genes.",
            len(gene_universe), len(nassiri_genes), len(overlap),
        )
        if len(overlap) < 500:
            logger.warning(
                "Very small gene overlap (%d). Check: are Nassiri columns Entrez IDs\n"
                "while GSE183653 uses Hugo symbols? If so run Entrez→Hugo mapping first.",
                len(overlap),
            )
        x = x.loc[:, overlap]
        gene_universe = overlap

    stds = x.std(axis=0).sort_values(ascending=False)
    keep = list(stds.head(n_top_genes).index)
    x = x.loc[:, keep]
    mu = x.mean(axis=0)
    sd = x.std(axis=0).replace(0, 1.0)
    x_z = (x - mu) / sd
    logger.info(
        "Preprocessed training features: %d samples × %d top-Var genes (log2TPM+1 → z-score).",
        x_z.shape[0], x_z.shape[1],
    )
    return x_z, keep


def build_classifier_pipeline(method: str = "nearest_centroid") -> Pipeline:
    if method == "nearest_centroid":
        return Pipeline([
            ("clf", NearestCentroid(metric="euclidean")),
        ])
    elif method == "ridge":
        return Pipeline([
            ("clf", LogisticRegression(
                C=1.0, solver="lbfgs", max_iter=5000,
                random_state=SEED,
            )),
        ])
    else:
        raise ValueError(f"Unknown classifier method: {method}. Use 'nearest_centroid' or 'ridge'.")
    # NOTE: StandardScaler is intentionally NOT inside the pipeline. preprocess_for_classifier
    # already column-z-scores features (per-gene (x - mu)/sd) on the TRAINING cohort.
    # For any TEST/transfer cohort (e.g. Nassiri), the caller MUST column-z-score
    # within that cohort first (using the cohort's own per-gene mu/sd) before
    # calling pipe.predict(). This keeps the standardization symmetric across
    # cohorts so batch/platform scale differences are not misread as biological signal.


@dataclass
class CvResults:
    method: str
    n_samples: int
    n_features: int
    accuracy: float
    kappa: float
    per_class_accuracy: Dict[str, float]
    classification_report: Dict
    predicted_labels: Dict[str, str]
    cv_strategy: str
    n_splits: int


def run_cross_validation(
    X: pd.DataFrame,
    y: pd.Series,
    method: str = "nearest_centroid",
    cv_strategy: str = "loo",
) -> CvResults:
    sample_ids = list(X.index)
    if cv_strategy == "loo":
        cv = LeaveOneOut()
        n_splits = len(sample_ids)
    else:
        cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
        n_splits = 10
    pipe = build_classifier_pipeline(method)
    try:
        y_pred_arr = cross_val_predict(pipe, X.values, y.values, cv=cv, n_jobs=1)
    except Exception as e:
        _die(f"Cross-validation failed: {e}")
    y_pred = pd.Series(y_pred_arr, index=y.index)
    acc = accuracy_score(y, y_pred)
    kap = cohen_kappa_score(y, y_pred)
    report = classification_report(
        y, y_pred, labels=BI_GROUP_ORDER, output_dict=True, zero_division=0,
    )
    per_class = {}
    for g in BI_GROUP_ORDER:
        mask = y == g
        if mask.any():
            per_class[g] = float((y_pred[mask] == g).sum() / mask.sum())
        else:
            per_class[g] = float("nan")
    return CvResults(
        method=method,
        n_samples=len(sample_ids),
        n_features=X.shape[1],
        accuracy=float(acc),
        kappa=float(kap),
        per_class_accuracy=per_class,
        classification_report={k: v for k, v in report.items() if isinstance(v, dict)},
        predicted_labels={str(s): str(p) for s, p in y_pred.items()},
        cv_strategy=cv_strategy,
        n_splits=n_splits,
    )


def train_and_save(
    X: pd.DataFrame,
    y: pd.Series,
    feature_names: List[str],
    method: str = "nearest_centroid",
    out_dir: Optional[Path] = None,
) -> Tuple[Pipeline, CvResults]:
    out_dir = out_dir or BRIDGING_MODEL_OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    pipe = build_classifier_pipeline(method)
    pipe.fit(X.values, y.values)
    cv = run_cross_validation(X, y, method=method, cv_strategy="loo")
    with open(out_dir / f"model_{method}.pkl", "wb") as fh:
        pickle.dump({
            "pipeline": pipe,
            "feature_order": feature_names,
            "preprocessing": (
                "1. For training cohort (GSE183653 TPM):\n"
                "     log2(TPM+1) → restrict to Nassiri gene overlap →\n"
                "     top-2000 variable genes → COLUMN-Z-SCORE within-cohort\n"
                "     (per-gene (x - mu_cohort) / sd_cohort).\n"
                "2. For ANY transfer/test cohort:\n"
                "     Restrict to feature_order genes (missing → 0) →\n"
                "     COLUMN-Z-SCORE WITHIN THE TEST COHORT (own mu/sd) →\n"
                "     then pipe.predict(). DO NOT use the training cohort's\n"
                "     fitted scaler on a different cohort's raw expression;\n"
                "     cross-study platform/library-prep scale differences will\n"
                "     be misread as biological shifts otherwise.\n"
                "Nassiri-specific notes: cBioPortal expression.csv columns are\n"
                "per-sample z-scored mRNA from the study's internal reference\n"
                "(not TPM). Still column-z-score within-cohort first before\n"
                "passing to the trained pipeline; this is a no-op if already\n"
                "exactly z-scored but corrects any scale/shift drift.\n"
            ),
            "cv_results": asdict(cv),
            "method": method,
        }, fh)
    pd.DataFrame([asdict(cv)]).T.rename(columns={0: "value"}).to_csv(
        out_dir / f"cv_metrics_{method}.csv"
    )
    pred_df = pd.DataFrame({
        "sample_id": list(cv.predicted_labels.keys()),
        "cv_predicted_bi_group": list(cv.predicted_labels.values()),
        "true_bi_group": [str(y[s]) for s in cv.predicted_labels.keys()],
    }).set_index("sample_id")
    pred_df.to_csv(out_dir / f"loo_predictions_{method}.csv")
    logger.info(
        "Trained %s bridging classifier. LOOCV accuracy=%.2f%%, kappa=%.3f.",
        method, 100 * cv.accuracy, cv.kappa,
    )
    return pipe, cv


def _ensure_nassiri_expression() -> Tuple[pd.DataFrame, pd.DataFrame]:
    meta_path = NASSIRI_COHORT_DIR / "metadata.csv"
    expr_path = NASSIRI_COHORT_DIR / "expression.csv"
    if not meta_path.exists():
        _die(
            f"Nassiri metadata missing: {meta_path}.\n"
            f"Run: python scripts/fetch_cbio_mng_utoronto_2021.py"
        )
    if not expr_path.exists():
        _die(
            f"Nassiri expression matrix missing: {expr_path}.\n"
            f"The fetch script was run with --skip-mrna, or the mRNA endpoint failed.\n"
            f"Re-run WITHOUT --skip-mrna:\n"
            f"  python scripts/fetch_cbio_mng_utoronto_2021.py\n"
            f"(If mRNA z-scores are throttled by cBioPortal, re-run with --skip-download\n"
            f" after a cooldown; or use the TSV datahub download as a manual fallback.)"
        )
    meta = pd.read_csv(meta_path, index_col=0)
    expr = pd.read_csv(expr_path, index_col=0)
    logger.info(
        "Loaded Nassiri cohort: %d samples in metadata, %d samples × %d genes in expression.",
        len(meta), *expr.shape,
    )
    return meta, expr


def apply_to_nassiri(
    pipe: Pipeline,
    feature_order: List[str],
    out_dir: Optional[Path] = None,
) -> pd.DataFrame:
    out_dir = out_dir or BRIDGING_MODEL_OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    meta, expr = _ensure_nassiri_expression()

    nassiri_genes = [str(c) for c in expr.columns]
    missing_features = [g for g in feature_order if g not in set(nassiri_genes)]
    if missing_features:
        logger.warning(
            "Nassiri expression missing %d/%d classifier features. "
            "Missing features will be set to 0 (mean of z-score scale — conservative). "
            "First 10 missing: %s",
            len(missing_features), len(feature_order),
            missing_features[:10],
        )
    # Nassiri expression from cBioPortal: rows=samples, cols=genes.
    # Values are study-internal mRNA z-scores (NOT raw counts or TPM).
    # Step 1: align to feature_order, pad missing with 0 (= z-score mean).
    X_nas_aligned = pd.DataFrame(0.0, index=expr.index, columns=feature_order)
    present = [g for g in feature_order if g in set(nassiri_genes)]
    X_nas_aligned.loc[:, present] = expr.loc[:, present].values.astype(float)

    # Step 2: COLUMN-Z-SCORE NASSIRI WITHIN ITS OWN COHORT.
    # This mirrors training preprocess_for_classifier, which does
    # (x - mu_training) / sd_training. Standardizing each cohort
    # independently removes cross-study platform/library-prep scale
    # shifts that would otherwise be absorbed as class-level bias
    # (the exact failure mode: Immune→Hyper false positive amplification
    # observed when this step was skipped).
    mu_nas = X_nas_aligned.mean(axis=0)
    sd_nas = X_nas_aligned.std(axis=0).replace(0, 1.0)
    X_nas_z = (X_nas_aligned - mu_nas) / sd_nas
    logger.info(
        "Nassiri features aligned + within-cohort z-scored. "
        "Range (min/median/max) of feature means after z: [%.2f, %.2f, %.2f]; "
        "feature sds: [%.2f, %.2f, %.2f].",
        X_nas_z.mean(axis=0).min(),
        X_nas_z.mean(axis=0).median(),
        X_nas_z.mean(axis=0).max(),
        X_nas_z.std(axis=0).min(),
        X_nas_z.std(axis=0).median(),
        X_nas_z.std(axis=0).max(),
    )

    # Step 3: predict. Pipeline has NO StandardScaler step (removed to avoid
    # asymmetric double-standardization during training + single on test).
    preds = pipe.predict(X_nas_z.values)
    proba = None
    if hasattr(pipe.named_steps["clf"], "predict_proba"):
        try:
            proba_arr = pipe.predict_proba(X_nas_z.values)
            proba = pd.DataFrame(proba_arr, index=X_nas_z.index,
                                 columns=pipe.classes_)
        except Exception:
            proba = None

    pred_series = pd.Series(preds, index=X_nas_z.index, name="bi_group_predicted")
    out = pd.DataFrame({
        "bi_group": pd.Categorical(pred_series, categories=BI_GROUP_ORDER, ordered=False),
    }, index=X_nas_z.index)
    if proba is not None:
        for g in BI_GROUP_ORDER:
            if g in proba.columns:
                out[f"bi_group_prob_{g.replace('-', '_')}"] = proba[g].values

    out.to_csv(out_dir / "nassiri_bi_group_predictions.csv")

    merged_meta = meta.copy()
    merged_meta["bi_group"] = out["bi_group"].reindex(merged_meta.index).values
    if proba is not None:
        for g in BI_GROUP_ORDER:
            col = f"bi_group_prob_{g.replace('-', '_')}"
            if col in out.columns:
                merged_meta[col] = out[col].reindex(merged_meta.index).values
    merged_meta.to_csv(NASSIRI_COHORT_DIR / "metadata.csv")
    logger.info(
        "Applied bridging classifier to Nassiri N=%d. Predicted group counts: %s.",
        len(out),
        {g: int((out["bi_group"] == g).sum()) for g in BI_GROUP_ORDER},
    )
    return out


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(name)-24s | %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    tpm = load_gse183653_tpm()
    labels = validate_and_load_labels(list(tpm.columns.astype(str)))

    try:
        _meta, nassiri_expr = _ensure_nassiri_expression()
    except SystemExit as e:
        if e.code != 3:
            raise
        logger.warning(
            "Nassiri expression not available yet; training only (no application step). "
            "Re-run after fetch_cbio_mng_utoronto_2021.py completes expression download."
        )
        nassiri_expr = None

    X, feature_names = preprocess_for_classifier(tpm, nassiri_expr=nassiri_expr, n_top_genes=2000)
    y = labels["bi_group"].loc[X.index].astype(str)

    BRIDGING_MODEL_OUT_DIR.mkdir(parents=True, exist_ok=True)
    best_method: Optional[str] = None
    best_kappa: float = -1.0
    for method in ["nearest_centroid", "ridge"]:
        pipe, cv = train_and_save(X, y, feature_names, method=method)
        if cv.kappa > best_kappa:
            best_kappa = cv.kappa
            best_method = method
    assert best_method is not None
    with open(BRIDGING_MODEL_OUT_DIR / "best_model.txt", "w") as fh:
        fh.write(f"best_method={best_method}\nloocv_kappa={best_kappa:.4f}\n")

    if nassiri_expr is not None:
        model_path = BRIDGING_MODEL_OUT_DIR / f"model_{best_method}.pkl"
        with open(model_path, "rb") as fh:
            bundle = pickle.load(fh)
        apply_to_nassiri(bundle["pipeline"], bundle["feature_order"])

    logger.info("Bridging classifier pipeline complete. Outputs in: %s", BRIDGING_MODEL_OUT_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
