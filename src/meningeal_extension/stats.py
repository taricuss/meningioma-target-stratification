from __future__ import annotations

import logging
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import (
    chi2_contingency,
    cramervonmises,
    fisher_exact,
    kruskal,
    mannwhitneyu,
    wilcoxon,
)
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.inter_rater import cohens_kappa

from .config import FDR_ALPHA, SIGNIFICANCE_THRESHOLD

logger = logging.getLogger(__name__)


def bh_fdr(pvalues: Sequence[float], alpha: float = FDR_ALPHA) -> pd.DataFrame:
    arr = np.asarray(pvalues, dtype=float)
    mask_valid = ~np.isnan(arr)
    reject = np.full(arr.shape, False)
    pvals_corrected = np.full(arr.shape, np.nan)
    if mask_valid.any():
        rej, pv_corr, _, _ = multipletests(arr[mask_valid], alpha=alpha, method="fdr_bh")
        reject[mask_valid] = rej
        pvals_corrected[mask_valid] = pv_corr
    return pd.DataFrame(
        {"p_raw": arr, "p_bh_fdr": pvals_corrected, "reject_bh_fdr": reject}
    )


def kruskal_wallis_across_groups(
    values: Sequence[float],
    groups: Sequence,
    nan_policy: str = "omit",
) -> Tuple[float, float, int, int]:
    val_arr = np.asarray(values, dtype=float)
    grp_arr = np.asarray(groups)
    if nan_policy == "omit":
        mask = ~np.isnan(val_arr)
        val_arr = val_arr[mask]
        grp_arr = grp_arr[mask]
    unique_groups = list(np.unique(grp_arr))
    n_groups = len(unique_groups)
    samples = [val_arr[grp_arr == g] for g in unique_groups]
    samples = [s for s in samples if len(s) > 0]
    if len(samples) < 2:
        return np.nan, np.nan, n_groups, 0
    try:
        stat, p = kruskal(*samples, nan_policy=nan_policy)
    except ValueError:
        return np.nan, np.nan, n_groups, 0
    return float(stat), float(p), n_groups, int(sum(len(s) for s in samples))


def pairwise_wilcoxon(
    values: Sequence[float],
    groups: Sequence,
    groups_order: Optional[Iterable] = None,
    alternative: str = "two-sided",
) -> pd.DataFrame:
    val_arr = np.asarray(values, dtype=float)
    grp_arr = np.asarray(groups)
    mask = ~np.isnan(val_arr)
    val_arr = val_arr[mask]
    grp_arr = grp_arr[mask]

    if groups_order is None:
        unique = list(np.unique(grp_arr))
    else:
        unique = [g for g in groups_order if g in np.unique(grp_arr)]
    rows = []
    for i, a in enumerate(unique):
        for b in unique[i + 1 :]:
            sa = val_arr[grp_arr == a]
            sb = val_arr[grp_arr == b]
            if len(sa) == 0 or len(sb) == 0:
                continue
            try:
                stat, p = mannwhitneyu(sa, sb, alternative=alternative)
                med_a = float(np.median(sa))
                med_b = float(np.median(sb))
                delta_median = med_a - med_b
                cd_value = _cliffs_delta(sa, sb)
            except ValueError:
                stat, p, med_a, med_b, delta_median, cd_value = (
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                )
            rows.append(
                {
                    "group_a": a,
                    "group_b": b,
                    "n_a": int(len(sa)),
                    "n_b": int(len(sb)),
                    "median_a": med_a,
                    "median_b": med_b,
                    "delta_median_a_minus_b": delta_median,
                    "mann_whitney_u": float(stat) if not np.isnan(stat) else stat,
                    "p_raw": p,
                    "cliffs_delta": cd_value,
                }
            )
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df
    fdr = bh_fdr(df["p_raw"].values)
    df["p_bh_fdr"] = fdr["p_bh_fdr"].values
    df["reject_bh_fdr"] = fdr["reject_bh_fdr"].values
    return df


def _cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    diff = x[:, None] - y[None, :]
    greater = np.sum(diff > 0)
    less = np.sum(diff < 0)
    return float((greater - less) / (len(x) * len(y)))


def mann_whitney_two_group(
    values: Sequence[float],
    groups: Sequence,
    group_labels: Tuple[str, str] = ("Intact", "Mutant/Loss"),
    alternative: str = "two-sided",
) -> pd.DataFrame:
    val_arr = np.asarray(values, dtype=float)
    grp_arr = np.asarray(groups)
    mask = ~np.isnan(val_arr)
    val_arr = val_arr[mask]
    grp_arr = grp_arr[mask]
    a_label, b_label = group_labels
    sa = val_arr[grp_arr == a_label]
    sb = val_arr[grp_arr == b_label]
    if len(sa) == 0 or len(sb) == 0:
        return pd.DataFrame(
            {
                "group_a": [a_label],
                "group_b": [b_label],
                "n_a": [int(len(sa))],
                "n_b": [int(len(sb))],
                "median_a": [np.nan],
                "median_b": [np.nan],
                "delta_median": [np.nan],
                "mann_whitney_u": [np.nan],
                "cliffs_delta": [np.nan],
                "p_raw": [np.nan],
                "p_bh_fdr": [np.nan],
                "reject_bh_fdr": [False],
            }
        )
    stat, p = mannwhitneyu(sa, sb, alternative=alternative)
    med_a = float(np.median(sa))
    med_b = float(np.median(sb))
    fdr = bh_fdr([p])
    return pd.DataFrame(
        {
            "group_a": [a_label],
            "group_b": [b_label],
            "n_a": [int(len(sa))],
            "n_b": [int(len(sb))],
            "median_a": [med_a],
            "median_b": [med_b],
            "delta_median_a_minus_b": [med_a - med_b],
            "mann_whitney_u": [float(stat)],
            "cliffs_delta": [_cliffs_delta(sa, sb)],
            "p_raw": [p],
            "p_bh_fdr": [fdr["p_bh_fdr"].values[0]],
            "reject_bh_fdr": [bool(fdr["reject_bh_fdr"].values[0])],
        }
    )


def cohens_kappa_table(labels_a: Sequence, labels_b: Sequence) -> pd.DataFrame:
    a = np.asarray(labels_a)
    b = np.asarray(labels_b)
    mask = ~(pd.isna(a) | pd.isna(b))
    a, b = a[mask], b[mask]
    if len(a) == 0:
        return pd.DataFrame({"n_paired": [0], "kappa": [np.nan], "p_value": [np.nan]})
    set_a = set(a.tolist())
    set_b = set(b.tolist())
    shared = set_a & set_b
    if len(shared) == 0:
        logger.warning(
            "cohens_kappa_table: label sets are DISJOINT (a=%s, b=%s). "
            "Building a square confusion matrix on the UNION of names will "
            "produce a degenerate kappa = 0.0 by construction (all diagonal "
            "cells are 0). The caller MUST map labels into a shared label "
            "space BEFORE calling this function. Returning kappa=NaN with "
            "n_categories=0 to signal this failure mode rather than "
            "silently returning a misleading 0.0.",
            sorted(set_a), sorted(set_b),
        )
        return pd.DataFrame(
            {
                "n_paired": [int(len(a))],
                "n_categories": [0],
                "cohens_kappa": [np.nan],
                "p_cohens_kappa": [np.nan],
                "chi2": [np.nan],
                "chi2_df": [np.nan],
                "p_chi2": [np.nan],
                "cramers_v": [np.nan],
                "kappa_failure_reason": ["disjoint_label_sets_no_shared_names"],
            }
        )
    if len(shared) < min(len(set_a), len(set_b)):
        logger.warning(
            "cohens_kappa_table: partial label overlap only (shared=%d of "
            "|a|=%d, |b|=%d). Kappa will only reflect the shared-label "
            "subset. Caller is responsible for label-space alignment if "
            "this is unintended. Shared names: %s. A-only: %s. B-only: %s.",
            len(shared), len(set_a), len(set_b),
            sorted(shared), sorted(set_a - shared), sorted(set_b - shared),
        )
    categories = sorted(set_a | set_b)
    n_cat = len(categories)
    idx = {c: i for i, c in enumerate(categories)}
    table = np.zeros((n_cat, n_cat), dtype=float)
    for ai, bi in zip(a, b):
        table[idx[ai], idx[bi]] += 1.0
    diag_sum = float(np.trace(table))
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = cohens_kappa(table.astype(int))
        kappa = float(res.kappa)
        if hasattr(res, "pvalue_one_sided"):
            p_val = float(res.pvalue_one_sided) if np.isfinite(res.pvalue_one_sided) else np.nan
        else:
            p_val = np.nan
        if not np.isfinite(kappa):
            kappa = np.nan
    except Exception:
        kappa, p_val = np.nan, np.nan
    # HARD SAFETY GUARD against the κ = 0.0 union-table artifact:
    # If the confusion matrix has ZERO diagonal entries (no exact label-name
    # matches at all) AND kappa reports exactly 0.0, this is always the
    # degenerate statsmodels output caused by building a square table on the
    # UNION of two disjoint label name spaces. There is NO legitimate
    # statistical scenario with N >= 2 samples, 2+ categories, and a non-null
    # association that produces exactly zero diagonal agreement — the caller
    # has simply forgotten to project both raters into a shared label space.
    # Fail by returning NaN with a loud failure reason rather than silently
    # leaking the misleading 0.0.
    if diag_sum == 0.0 and kappa == 0.0 and n_cat > 2 and len(a) >= 2:
        logger.error(
            "cohens_kappa_table: KAPPA-ARTIFACT TRAP TRIGGERED — diagonal_sum=%.0f, "
            "n_cat=%d, kappa=%.4f exactly. This is ALWAYS the degenerate "
            "statsmodels output on a union-based square confusion matrix of "
            "two DISJOINT label name spaces. The caller MUST project both "
            "label vectors into a shared name space BEFORE calling this "
            "function (see concordance.py::run_concordance parameter "
            "label_a_to_b_mapping=). Categories=%s. Returning NaN with "
            "explicit failure reason to kill any downstream headline numbers.",
            diag_sum, n_cat, kappa, categories,
        )
        return pd.DataFrame(
            {
                "n_paired": [int(len(a))],
                "n_categories": [0],
                "cohens_kappa": [np.nan],
                "p_cohens_kappa": [np.nan],
                "chi2": [np.nan],
                "chi2_df": [np.nan],
                "p_chi2": [np.nan],
                "cramers_v": [np.nan],
                "kappa_failure_reason": ["artifact_zero_diagonal_union_table_detected"],
            }
        )
    # Cramér's V with robust contingency (add 1e-9 jitter and use log-likelihood ratio)
    table_safe = table + 1e-9
    try:
        chi2, p_chi2, dof, _ = chi2_contingency(table_safe, lambda_="log-likelihood")
    except Exception:
        try:
            chi2, p_chi2, dof, _ = chi2_contingency(table_safe)
        except Exception:
            chi2, p_chi2, dof = np.nan, np.nan, np.nan
    n_total = int(table.sum())
    min_dim = min(table.shape) - 1
    if np.isfinite(chi2) and n_total * min_dim > 0:
        cramer_v = float(np.sqrt(chi2 / (n_total * min_dim)))
    else:
        cramer_v = np.nan
    return pd.DataFrame(
        {
            "n_paired": [n_total],
            "n_categories": [n_cat],
            "cohens_kappa": [kappa],
            "p_cohens_kappa": [p_val],
            "chi2": [chi2],
            "chi2_df": [dof],
            "p_chi2": [p_chi2],
            "cramers_v": [cramer_v],
        }
    )


def direction_concordance(
    discovery_effects: pd.Series,
    replication_effects: pd.Series,
) -> pd.DataFrame:
    joined = pd.concat([discovery_effects.rename("discovery"), replication_effects.rename("replication")], axis=1, join="inner")
    if len(joined) == 0:
        return pd.DataFrame({"n_tested": [0], "n_same_direction": [0], "fraction_same_direction": [np.nan]})
    same_dir = (np.sign(joined["discovery"]) == np.sign(joined["replication"])) | (
        (joined["discovery"] == 0) & (joined["replication"] == 0)
    )
    n_same = int(same_dir.sum())
    n_tested = int(len(joined))
    # binomial sign test p-value against 0.5 null
    p_binom = _binomial_twosided(n_same, n_tested, 0.5)
    return pd.DataFrame(
        {
            "n_tested": [n_tested],
            "n_same_direction": [n_same],
            "fraction_same_direction": [n_same / n_tested],
            "p_binomial_sign_test": [p_binom],
        }
    )


def _binomial_twosided(k: int, n: int, p: float) -> float:
    from math import comb
    if n == 0:
        return np.nan
    exp_lo = min(int(n * p), k)
    exp_hi = max(int(n * p), k)
    prob_k = comb(n, k) * (p ** k) * ((1 - p) ** (n - k))
    total = 0.0
    for i in range(n + 1):
        prob_i = comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
        if prob_i <= prob_k + 1e-15:
            total += prob_i
    return float(min(total, 1.0))


def stars_for_p(p: float) -> str:
    if np.isnan(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < SIGNIFICANCE_THRESHOLD:
        return "*"
    return "ns"


cliffs_delta = _cliffs_delta
