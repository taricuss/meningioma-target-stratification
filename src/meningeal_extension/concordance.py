from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import NASSIRI_BI_MATCH_HYPOTHESIS, PATHS
from .datasets import load_cohort
from .stats import cohens_kappa_table

logger = logging.getLogger(__name__)


def enumerate_4to3_pairings(
    nassiri_groups: Optional[List[str]] = None,
    bi_groups: Optional[List[str]] = None,
    *,
    hypermitotic_label: str = "Hypermitotic",
    merge_target: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Enumerate all distinct 4→3 many-to-one mappings of Nassiri→Bi groups.

    Scenario C (default): one Bi group (Hypermitotic) absorbs exactly 2 Nassiri groups;
    the remaining 2 Bi singleton groups each receive 1 Nassiri group. Total:
    C(4,2) choices of merge-pair × 2! permutations of the 2 singletons = 12 pairings.

    If merge_target is set to None (Scenario B), the 3 resulting bins (2 singletons +
    1 merged pair) are assigned to all 3 Bi groups via 3! permutations = 36 pairings.
    This is less informative for this domain since the Hypermitotic label has a clear
    high-proliferation meaning; Scenario C is therefore the default.

    Returns a list of dicts of the form {nassiri_label: bi_label, ...}.
    """
    if nassiri_groups is None:
        nassiri_groups = sorted({k for k in NASSIRI_BI_MATCH_HYPOTHESIS.keys()})
    if bi_groups is None:
        bi_groups = sorted({v for v in NASSIRI_BI_MATCH_HYPOTHESIS.values()})
    # Handle union-of-unique cases: the true_mapping has 4 keys; if the cohort's own
    # Nassiri labels include extras (e.g., case-difference variants or legacy labels),
    # restrict enumeration to the 4 true_mapping keys (which is what pairings actually use).
    canonical_nas_keys = sorted({k for k in NASSIRI_BI_MATCH_HYPOTHESIS.keys()})
    if len(nassiri_groups) != 4 and set(canonical_nas_keys).issubset(set(nassiri_groups)):
        nassiri_groups = canonical_nas_keys
    canonical_bi_keys = sorted({v for v in NASSIRI_BI_MATCH_HYPOTHESIS.values()})
    if len(bi_groups) != 3 and set(canonical_bi_keys).issubset(set(bi_groups)):
        bi_groups = canonical_bi_keys
    if len(nassiri_groups) != 4:
        raise ValueError(f"Expected 4 Nassiri groups, got {len(nassiri_groups)}: {nassiri_groups}")
    if len(bi_groups) != 3:
        raise ValueError(f"Expected 3 Bi groups, got {len(bi_groups)}: {bi_groups}")

    pairings: List[Dict[str, str]] = []
    if merge_target is None:
        # Scenario B: 36 total
        for pair in combinations(nassiri_groups, 2):
            merged = set(pair)
            singletons = [g for g in nassiri_groups if g not in merged]
            for bi_perm in permutations(bi_groups):
                p: Dict[str, str] = {}
                p[singletons[0]] = bi_perm[0]
                p[singletons[1]] = bi_perm[1]
                p[pair[0]] = bi_perm[2]
                p[pair[1]] = bi_perm[2]
                pairings.append(p)
    else:
        # Scenario C: merged pair must go to merge_target (default Hypermitotic) = 12 total
        if merge_target not in bi_groups:
            raise ValueError(f"merge_target={merge_target} not in bi_groups={bi_groups}")
        bi_singletons = [b for b in bi_groups if b != merge_target]
        if len(bi_singletons) != 2:
            raise ValueError(f"Expected exactly 2 non-merge Bi groups, got {len(bi_singletons)}: {bi_singletons}")
        for pair in combinations(nassiri_groups, 2):
            merged = set(pair)
            singletons = [g for g in nassiri_groups if g not in merged]
            for bi_sing_perm in permutations(bi_singletons):
                p = {}
                p[singletons[0]] = bi_sing_perm[0]
                p[singletons[1]] = bi_sing_perm[1]
                p[pair[0]] = merge_target
                p[pair[1]] = merge_target
                pairings.append(p)
    return pairings


def pairing_equals(a: Dict[str, str], b: Dict[str, str]) -> bool:
    if set(a.keys()) != set(b.keys()):
        return False
    return all(a[k] == b[k] for k in a)


@dataclass
class ConcordanceResult:
    n_paired: int
    cohens_kappa: float
    p_cohens_kappa: float
    cramers_v: float
    p_chi2: float
    cross_tabulation: pd.DataFrame
    metrics_table: pd.DataFrame


def run_concordance(
    metadata: pd.DataFrame,
    label_a_col: str = "nassiri_group",
    label_b_col: str = "bi_group",
    cohort_name: str = "",
    *,
    label_a_to_b_mapping: Optional[Dict[str, str]] = None,
) -> ConcordanceResult:
    """Run concordance analysis between two label columns.

    If label_a_col and label_b_col use DISJOINT label name spaces (e.g.
    Nassiri 4-group names vs Bi 3-group names, zero shared labels), the raw
    cohens_kappa_table call would return a DEGENERATE kappa = NaN (formerly
    silently returned 0.0 before the disjoint-sets guard was added to
    stats.py::cohens_kappa_table).

    To fix this, provide ``label_a_to_b_mapping``: a dict that projects the
    label_a names onto the label_b name space. When this mapping is provided:
      * the CROSS-TABULATION remains the original |A| x |B| table (so the
        4x3 Nassiri-vs-Bi detail is preserved in the output CSVs),
      * BUT Cohen's kappa is computed on the 3x3 table of
        ``mapped(label_a) vs label_b``, which is the operationally correct
        3-group agreement score for the registered many-to-one mapping.
      * Chi-square / Cramer's V are reported on BOTH: the original |A| x |B|
        table (association-strength, does not require shared names) AND the
        3x3 mapped table (agreement-strength metric).
    """
    logger.info("Running concordance: %s vs %s (n=%d)", label_a_col, label_b_col, len(metadata))

    meta_clean = metadata[[label_a_col, label_b_col]].dropna()
    a_vals = meta_clean[label_a_col].values
    b_vals = meta_clean[label_b_col].values
    set_a_names = set(a_vals.tolist())
    set_b_names = set(b_vals.tolist())
    n_shared = len(set_a_names & set_b_names)

    cross = pd.crosstab(
        metadata[label_a_col],
        metadata[label_b_col],
        margins=True,
        margins_name="Total",
    )

    # --- Chi-square / Cramer's V on the ORIGINAL (possibly non-square) table ---
    # This works fine without shared names: it's just an association test on
    # the non-square contingency table.
    from scipy.stats import chi2_contingency
    cross_nomargins = cross.drop(index="Total", columns="Total", errors="ignore")
    try:
        table_safe = cross_nomargins.values.astype(float) + 1e-9
        chi2_axb, p_chi2_axb, dof_axb, _ = chi2_contingency(table_safe, lambda_="log-likelihood")
        n_total = int(cross_nomargins.values.sum())
        r, c = cross_nomargins.shape
        min_dim = min(r, c) - 1
        cramer_v_axb = float(np.sqrt(chi2_axb / (n_total * min_dim))) if (n_total * min_dim) > 0 and np.isfinite(chi2_axb) else np.nan
    except Exception:
        chi2_axb, p_chi2_axb, dof_axb, cramer_v_axb = np.nan, np.nan, np.nan, np.nan

    metrics_all = {
        "cohort": cohort_name,
        "label_a": label_a_col,
        "label_b": label_b_col,
        "n_paired": int(len(meta_clean)),
        "n_label_a_levels": len(set_a_names),
        "n_label_b_levels": len(set_b_names),
        "n_shared_label_names": int(n_shared),
        # 4x3 table association (no name-sharing required; chi-square / Cramer's V)
        "chi2_original_axa_x_b": float(chi2_axb) if np.isfinite(chi2_axb) else np.nan,
        "chi2_df_original": int(dof_axb) if np.isfinite(dof_axb) else np.nan,
        "p_chi2_original": float(p_chi2_axb) if np.isfinite(p_chi2_axb) else np.nan,
        "cramers_v_original": float(cramer_v_axb) if np.isfinite(cramer_v_axb) else np.nan,
    }

    if label_a_to_b_mapping is not None and n_shared == 0:
        logger.info(
            "Label spaces are disjoint (|a|=%d, |b|=%d, shared=%d). "
            "Applying provided label_a->label_b mapping for Cohen's kappa "
            "computation. Mapping: %s",
            len(set_a_names), len(set_b_names), n_shared, label_a_to_b_mapping,
        )
        a_mapped = meta_clean[label_a_col].map(label_a_to_b_mapping)
        unknown_mask = a_mapped.isna()
        if unknown_mask.any():
            logger.warning(
                "%d samples in label_a have values not present in "
                "label_a_to_b_mapping keys; dropping for 3x3 kappa. Unknown: %s",
                int(unknown_mask.sum()),
                sorted(set(meta_clean.loc[unknown_mask, label_a_col].astype(str))),
            )
        a_mapped_clean = a_mapped.dropna().values
        b_aligned = meta_clean.loc[a_mapped.notna(), label_b_col].values
        km_mapped = cohens_kappa_table(a_mapped_clean, b_aligned)
        metrics_all.update({
            "cohens_kappa_mapped_3x3": float(km_mapped["cohens_kappa"].iloc[0]),
            "p_cohens_kappa_mapped_3x3": float(km_mapped["p_cohens_kappa"].iloc[0]),
            "chi2_mapped_3x3": float(km_mapped["chi2"].iloc[0]),
            "chi2_df_mapped_3x3": int(km_mapped["chi2_df"].iloc[0]) if np.isfinite(km_mapped["chi2_df"].iloc[0]) else np.nan,
            "p_chi2_mapped_3x3": float(km_mapped["p_chi2"].iloc[0]),
            "cramers_v_mapped_3x3": float(km_mapped["cramers_v"].iloc[0]),
            "kappa_mapping_applied": "a_to_b",
        })
    elif label_a_to_b_mapping is not None:
        a_mapped = meta_clean[label_a_col].map(label_a_to_b_mapping).values
        b_aligned = b_vals
        km_mapped = cohens_kappa_table(a_mapped, b_aligned)
        metrics_all.update({
            "cohens_kappa_mapped_3x3": float(km_mapped["cohens_kappa"].iloc[0]),
            "p_cohens_kappa_mapped_3x3": float(km_mapped["p_cohens_kappa"].iloc[0]),
            "chi2_mapped_3x3": float(km_mapped["chi2"].iloc[0]),
            "chi2_df_mapped_3x3": int(km_mapped["chi2_df"].iloc[0]) if np.isfinite(km_mapped["chi2_df"].iloc[0]) else np.nan,
            "p_chi2_mapped_3x3": float(km_mapped["p_chi2"].iloc[0]),
            "cramers_v_mapped_3x3": float(km_mapped["cramers_v"].iloc[0]),
            "kappa_mapping_applied": "a_to_b_even_with_some_shared",
        })
    else:
        # Fallback: raw (possibly still disjoint) kappa_table call; now
        # returns NaN for disjoint sets thanks to the guard in stats.py.
        km = cohens_kappa_table(a_vals, b_vals)
        metrics_all.update({
            "cohens_kappa_mapped_3x3": float(km["cohens_kappa"].iloc[0]),
            "p_cohens_kappa_mapped_3x3": float(km["p_cohens_kappa"].iloc[0]),
            "chi2_mapped_3x3": float(km["chi2"].iloc[0]),
            "chi2_df_mapped_3x3": int(km["chi2_df"].iloc[0]) if np.isfinite(km["chi2_df"].iloc[0]) else np.nan,
            "p_chi2_mapped_3x3": float(km["p_chi2"].iloc[0]),
            "cramers_v_mapped_3x3": float(km["cramers_v"].iloc[0]),
            "kappa_mapping_applied": "none_direct_call",
        })

    metrics_df = pd.DataFrame([metrics_all])
    # Backwards-compatible fields: expose the headline kappa as the top-level
    # ConcordanceResult.{cohens_kappa, p_cohens_kappa, cramers_v, p_chi2}.
    # The headline kappa is ALWAYS the 3x3 mapped one when available; else
    # the direct kappa (if shared names exist; else NaN per guard).
    headline_kappa = float(metrics_all["cohens_kappa_mapped_3x3"])
    headline_p_kappa = float(metrics_all["p_cohens_kappa_mapped_3x3"])
    # For cramer's V / chi-square headline: use the ORIGINAL non-square
    # table association, since that's what most readers interpret as the
    # 'overall concordance strength' across the 4x3 mapping.
    headline_cv = float(metrics_all["cramers_v_original"])
    headline_p_chi2 = float(metrics_all["p_chi2_original"])

    return ConcordanceResult(
        n_paired=int(metrics_all["n_paired"]),
        cohens_kappa=headline_kappa,
        p_cohens_kappa=headline_p_kappa,
        cramers_v=headline_cv,
        p_chi2=headline_p_chi2,
        cross_tabulation=cross,
        metrics_table=metrics_df,
    )


def save_concordance_tables(
    result: ConcordanceResult,
    prefix: str = "table_s1_aim1_concordance",
    out_dir: Optional[Path] = None,
) -> dict:
    if out_dir is None:
        out_dir = PATHS.results_tables
    out_dir.mkdir(parents=True, exist_ok=True)
    cross_path = out_dir / f"{prefix}_crosstab_counts.csv"
    cross_pct_path = out_dir / f"{prefix}_crosstab_rowpct.csv"
    metrics_path = out_dir / f"{prefix}_metrics.csv"
    result.cross_tabulation.to_csv(cross_path)
    (pd.crosstab(
        pd.Series([""] * len(result.cross_tabulation)),
        pd.Series([""] * len(result.cross_tabulation)),
    )).to_csv(cross_pct_path)  # placeholder; actual row pcts produced by caller
    # Regenerate cross_pct properly
    cross = result.cross_tabulation.drop(index="Total", columns="Total", errors="ignore")
    row_pct = cross.div(cross.sum(axis=1), axis=0).round(4) * 100
    row_pct.to_csv(cross_pct_path)
    result.metrics_table.to_csv(metrics_path, index=False)
    return {
        "counts": cross_path,
        "row_pct": cross_pct_path,
        "metrics": metrics_path,
    }


@dataclass
class PairingSpecificityResult:
    n_pairings_tested: int
    true_mapping_rank_kappa: int
    true_mapping_rank_agreement: int
    true_mapping_kappa: float
    best_alt_kappa: float
    best_alt_agreement: float
    delta_kappa_vs_best_alt: float
    delta_agreement_vs_best_alt: float
    pairing_ranked_table: pd.DataFrame
    singleton_swap_kappa: float
    singleton_swap_agreement: float
    worst_merge_kappa: float
    worst_merge_agreement: float
    pass_strict: bool
    pass_any_alt: bool


def _agreement(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) != len(b):
        raise ValueError(f"Length mismatch: {len(a)} vs {len(b)}")
    if len(a) == 0:
        return 0.0
    return float((np.asarray(a) == np.asarray(b)).mean())


def _bi_label_from_mapping(nassiri_series: pd.Series, mapping: Dict[str, str]) -> pd.Series:
    """Apply a Nassiri→Bi mapping to produce a predicted Bi label per sample."""
    out = nassiri_series.map(mapping)
    unknown_mask = out.isna()
    if unknown_mask.any():
        logger.warning(
            "%d samples have Nassiri labels not present in mapping keys=%s; set to NaN.",
            int(unknown_mask.sum()), list(mapping.keys()),
        )
    return out


def test_pairing_specificity(
    metadata: pd.DataFrame,
    nassiri_col: str = "nassiri_group",
    true_bi_col: str = "bi_group",
    true_mapping: Optional[Dict[str, str]] = None,
    *,
    merge_target: Optional[str] = "Hypermitotic",
    cohort_name: str = "",
) -> PairingSpecificityResult:
    """Operational H1 falsifiability test.

    Scores the literature-registered true Nassiri→Bi mapping against the full
    enumeration of alternative 4→3 pairings (Scenario C by default: 12 total,
    1 true + 11 alternatives that each send a different merged Nassiri pair
    to Bi's Hypermitotic group). Reports:

      * true mapping's rank among the 12 by Cohen's κ and raw agreement
      * Δκ vs the best-scoring alternative (strictest hurdle: Δκ > 0)
      * specific contrasts:
          - singleton-swap placebo (MG1↔Merlin, MG2↔Immune swapped)
          - worst among the 5 alternative merge choices
          - best alternative overall

    Two pass flags are returned:
      pass_strict   = true mapping is #1 by BOTH κ and raw agreement
      pass_any_alt  = true mapping outperforms the SINGLE best alternative
                      by BOTH κ and agreement (less strict but still meaningful)

    Note: This function uses the *actual Bi labels from the cohort* (true_bi_col)
    as ground truth for scoring concordance, NOT the predicted labels from the
    mapping applied to Nassiri labels. The scoring metric is: how well does the
    Nassiri-derived label (via a given pairing mapping) agree with the cohort's
    independent Bi labels?
    """
    if true_mapping is None:
        true_mapping = dict(NASSIRI_BI_MATCH_HYPOTHESIS)

    if nassiri_col not in metadata.columns:
        raise ValueError(f"nassiri_col='{nassiri_col}' not in metadata columns={list(metadata.columns)}")
    if true_bi_col not in metadata.columns:
        raise ValueError(f"true_bi_col='{true_bi_col}' not in metadata columns={list(metadata.columns)}")

    meta = metadata[[nassiri_col, true_bi_col]].dropna().copy()
    n_samples = len(meta)
    if n_samples == 0:
        raise ValueError("No paired (nassiri, bi) labels after dropna; cannot score.")

    nas_unique = sorted(meta[nassiri_col].dropna().unique().tolist())
    bi_unique = sorted(meta[true_bi_col].dropna().unique().tolist())

    pairings = enumerate_4to3_pairings(
        nassiri_groups=sorted(set(nas_unique) | set(true_mapping.keys())),
        bi_groups=sorted(set(bi_unique) | set(true_mapping.values())),
        merge_target=merge_target,
    )
    logger.info(
        "Pairing specificity: scoring true mapping vs %d alternatives (cohort=%s, n=%d paired samples).",
        len(pairings) - 1, cohort_name, n_samples,
    )

    rows = []
    true_kappa = None
    true_agreement = None
    swap_kappa = None
    swap_agreement = None

    for p in pairings:
        is_true = pairing_equals(p, true_mapping)
        bi_predicted = _bi_label_from_mapping(meta[nassiri_col], p)
        aligned = pd.DataFrame({"pred": bi_predicted, "truth": meta[true_bi_col]}).dropna()
        if len(aligned) < 2:
            logger.warning("Pairing has <2 aligned samples; skipping metrics: %s", p)
            kappa_val = np.nan
            agr_val = np.nan
        else:
            km = cohens_kappa_table(aligned["pred"].values, aligned["truth"].values)
            kappa_val = float(km["cohens_kappa"].iloc[0])
            agr_val = _agreement(aligned["pred"].values, aligned["truth"].values)

        row = {
            "is_true_mapping": is_true,
            "cohens_kappa": kappa_val,
            "raw_agreement": agr_val,
            "n_aligned": int(len(aligned)),
        }
        merged_pair = sorted([nk for nk, nv in p.items() if sum(1 for vv in p.values() if vv == nv) == 2])
        # merged_pair column uses ASCII-friendly '|' separator to avoid
        # cp1252 console crashes on Windows stdout. The Unicode '∪'
        # character would trigger UnicodeEncodeError: 'charmap' codec
        # can't encode character '\u222a' on Windows terminals with
        # the default cp1252 encoding. If a manuscript display needs
        # the proper union symbol, post-process the CSV in the report.
        row["merged_pair"] = "|".join(merged_pair) if len(merged_pair) == 2 else ""
        singletons = [(nk, nv) for nk, nv in sorted(p.items()) if nk not in set(merged_pair)]
        for i, (nk, nv) in enumerate(singletons):
            row[f"singleton_{i+1}_nas"] = nk
            row[f"singleton_{i+1}_bi"] = nv
        rows.append(row)

        if is_true:
            true_kappa = kappa_val
            true_agreement = agr_val

        # Detect singleton-swap placebo (same merge pair as true, singletons swapped)
        true_merged = sorted([nk for nk, nv in true_mapping.items()
                              if sum(1 for vv in true_mapping.values() if vv == nv) == 2])
        this_merged = merged_pair
        if sorted(true_merged) == sorted(this_merged) and not is_true:
            # Same merge; it's either the true mapping or the singleton swap. We already know !is_true.
            swap_kappa = kappa_val
            swap_agreement = agr_val

    tbl = pd.DataFrame(rows)
    tbl["rank_kappa_desc"] = tbl["cohens_kappa"].rank(ascending=False, method="min")
    tbl["rank_agreement_desc"] = tbl["raw_agreement"].rank(ascending=False, method="min")

    true_row = tbl[tbl["is_true_mapping"]].iloc[0]
    rank_kappa = int(true_row["rank_kappa_desc"])
    rank_agreement = int(true_row["rank_agreement_desc"])

    alt_rows = tbl[~tbl["is_true_mapping"]]
    best_alt_kappa = float(alt_rows["cohens_kappa"].max()) if len(alt_rows) else float("nan")
    best_alt_agreement = float(alt_rows["raw_agreement"].max()) if len(alt_rows) else float("nan")
    delta_kappa = float(true_kappa - best_alt_kappa) if (true_kappa is not None and not np.isnan(best_alt_kappa)) else float("nan")
    delta_agr = float(true_agreement - best_alt_agreement) if (true_agreement is not None and not np.isnan(best_alt_agreement)) else float("nan")

    # Worst alternative merge choice among the 5 that use a different merged pair
    true_merged_set = sorted([nk for nk, nv in true_mapping.items()
                              if sum(1 for vv in true_mapping.values() if vv == nv) == 2])
    true_merged_str = "|".join(sorted(true_merged_set))
    alt_merge_rows = alt_rows[alt_rows["merged_pair"] != true_merged_str]
    if len(alt_merge_rows) > 0:
        worst_merge_kappa = float(alt_merge_rows["cohens_kappa"].min())
        worst_merge_agreement = float(alt_merge_rows["raw_agreement"].min())
    else:
        worst_merge_kappa = float("nan")
        worst_merge_agreement = float("nan")

    pass_strict = (rank_kappa == 1) and (rank_agreement == 1)
    pass_any_alt = (
        true_kappa is not None and true_agreement is not None
        and not np.isnan(best_alt_kappa) and not np.isnan(best_alt_agreement)
        and true_kappa > best_alt_kappa
        and true_agreement > best_alt_agreement
    )

    logger.info(
        "Pairing specificity results | true_mapping rank_kappa=#%d rank_agr=#%d | "
        "delta_kappa vs best-alt=%.3f (best-alt kappa=%.3f, true kappa=%.3f) | "
        "singleton-swap kappa=%.3f | worst-alt-merge kappa=%.3f | "
        "pass_strict=%s pass_any_alt=%s",
        rank_kappa, rank_agreement, delta_kappa, best_alt_kappa, true_kappa or float("nan"),
        swap_kappa if swap_kappa is not None else float("nan"),
        worst_merge_kappa, pass_strict, pass_any_alt,
    )

    return PairingSpecificityResult(
        n_pairings_tested=len(pairings),
        true_mapping_rank_kappa=rank_kappa,
        true_mapping_rank_agreement=rank_agreement,
        true_mapping_kappa=float(true_kappa) if true_kappa is not None else float("nan"),
        best_alt_kappa=best_alt_kappa,
        best_alt_agreement=best_alt_agreement,
        delta_kappa_vs_best_alt=delta_kappa,
        delta_agreement_vs_best_alt=delta_agr,
        pairing_ranked_table=tbl.sort_values("cohens_kappa", ascending=False).reset_index(drop=True),
        singleton_swap_kappa=float(swap_kappa) if swap_kappa is not None else float("nan"),
        singleton_swap_agreement=float(swap_agreement) if swap_agreement is not None else float("nan"),
        worst_merge_kappa=worst_merge_kappa,
        worst_merge_agreement=worst_merge_agreement,
        pass_strict=bool(pass_strict),
        pass_any_alt=bool(pass_any_alt),
    )


def save_pairing_specificity_tables(
    result: PairingSpecificityResult,
    prefix: str = "table_s1b_aim1_pairing_specificity",
    out_dir: Optional[Path] = None,
) -> dict:
    if out_dir is None:
        out_dir = PATHS.results_tables
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = out_dir / f"{prefix}_all_pairings_ranked.csv"
    summary_path = out_dir / f"{prefix}_summary.csv"
    result.pairing_ranked_table.to_csv(ranked_path, index=False)
    pd.DataFrame([{
        "n_pairings_tested": result.n_pairings_tested,
        "true_mapping_rank_kappa": result.true_mapping_rank_kappa,
        "true_mapping_rank_agreement": result.true_mapping_rank_agreement,
        "true_mapping_kappa": result.true_mapping_kappa,
        "best_alt_kappa": result.best_alt_kappa,
        "best_alt_agreement": result.best_alt_agreement,
        "delta_kappa_vs_best_alt": result.delta_kappa_vs_best_alt,
        "delta_agreement_vs_best_alt": result.delta_agreement_vs_best_alt,
        "singleton_swap_kappa": result.singleton_swap_kappa,
        "singleton_swap_agreement": result.singleton_swap_agreement,
        "worst_merge_kappa": result.worst_merge_kappa,
        "worst_merge_agreement": result.worst_merge_agreement,
        "pass_strict": result.pass_strict,
        "pass_any_alt": result.pass_any_alt,
    }]).to_csv(summary_path, index=False)
    return {
        "ranked_pairings": ranked_path,
        "summary": summary_path,
    }
