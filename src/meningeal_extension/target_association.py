from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols, wls

from .config import (
    FDR_ALPHA,
    NF2_EXPORT_COL,
    NF2_INTERNAL_COL,
    PATHS,
)


from .gene_programs import ALL_PROGRAMS, PROGRAM_DESCRIPTIONS, PROGRAM_TO_DRUGS
from .ssgsea import score_all_programs, zscore_per_program
from .stats import (
    bh_fdr,
    kruskal_wallis_across_groups,
    mann_whitney_two_group,
    pairwise_wilcoxon,
)

logger = logging.getLogger(__name__)


def _resolve_nf2_col(meta: pd.DataFrame, requested: Optional[str] = None) -> Optional[str]:
    if requested is not None and requested in meta.columns:
        return requested
    if NF2_INTERNAL_COL in meta.columns:
        return NF2_INTERNAL_COL
    if NF2_EXPORT_COL in meta.columns:
        return NF2_EXPORT_COL
    return None


def run_aim2(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    subgroup_col: str = "nassiri_group",
    nf2_col: str = "nf2_status",
    grade_col: str = "who_grade",
    cohort_name: str = "",
) -> Dict[str, pd.DataFrame]:
    logger.info("Aim2: scoring %d programs for cohort %s (n=%d)", len(ALL_PROGRAMS), cohort_name, len(expression))
    program_scores = score_all_programs(expression)
    program_scores_z = zscore_per_program(program_scores)

    # Align indices
    common = program_scores_z.index.intersection(metadata.index)
    if len(common) < len(program_scores_z):
        logger.warning("Dropping %d samples with no metadata", len(program_scores_z) - len(common))
    program_scores_z = program_scores_z.loc[common]
    meta = metadata.loc[common]

    nf2_col_resolved = _resolve_nf2_col(meta, nf2_col)
    if nf2_col_resolved is None:
        logger.warning("No NF2 column found in metadata (tried %s, %s); H3 will be skipped.",
                       NF2_INTERNAL_COL, NF2_EXPORT_COL)
    elif nf2_col_resolved != nf2_col:
        logger.info("Resolved NF2 column: requested %s not present; using %s.",
                    nf2_col, nf2_col_resolved)
        nf2_col = nf2_col_resolved

    # --- H2: subgroup association (Kruskal-Wallis + pairwise Wilcoxon with BH-FDR)
    h2_rows = []
    pairwise_parts = []
    for program_name in program_scores_z.columns:
        stat, p, n_groups, n_samples = kruskal_wallis_across_groups(
            program_scores_z[program_name].values,
            meta[subgroup_col].values,
        )
        h2_rows.append(
            {
                "program": program_name,
                "drugs": "; ".join(PROGRAM_TO_DRUGS.get(program_name, [])),
                "subgroup_col": subgroup_col,
                "n_groups": n_groups,
                "n_samples": n_samples,
                "kruskal_wallis_H": stat,
                "kruskal_wallis_p_raw": p,
            }
        )
        pw = pairwise_wilcoxon(
            program_scores_z[program_name].values,
            meta[subgroup_col].values,
            groups_order=list(meta[subgroup_col].cat.categories) if hasattr(meta[subgroup_col], "cat") else None,
        )
        if not pw.empty:
            pw.insert(0, "program", program_name)
            pw.insert(1, "subgroup_col", subgroup_col)
            pairwise_parts.append(pw)

    h2_df = pd.DataFrame(h2_rows)
    # BH-FDR across program×subgroup (family of all KW tests)
    if len(h2_df) > 0:
        fdr = bh_fdr(h2_df["kruskal_wallis_p_raw"].values)
        h2_df["kruskal_wallis_p_bh_fdr"] = fdr["p_bh_fdr"].values
        h2_df["kruskal_wallis_reject_bh_fdr"] = fdr["reject_bh_fdr"].values
        h2_df.insert(0, "cohort", cohort_name)
        h2_df["program_description"] = h2_df["program"].map(PROGRAM_DESCRIPTIONS)

    pairwise_df = pd.concat(pairwise_parts, ignore_index=True) if pairwise_parts else pd.DataFrame()
    if not pairwise_df.empty:
        pairwise_df.insert(0, "cohort", cohort_name)
        # Family-wide BH-FDR across all pairwise tests
        fdr_all = bh_fdr(pairwise_df["p_raw"].values)
        pairwise_df["p_bh_fdr_all_program_pairs"] = fdr_all["p_bh_fdr"].values
        pairwise_df["reject_bh_fdr_all_program_pairs"] = fdr_all["reject_bh_fdr"].values

    # --- H3: NF2 association (Mann-Whitney U)
    h3_df = pd.DataFrame()
    if nf2_col_resolved is not None:
        h3_rows = []
        for program_name in program_scores_z.columns:
            mw = mann_whitney_two_group(
                program_scores_z[program_name].values,
                meta[nf2_col].values,
                group_labels=("Intact", "Mutant/Loss"),
            )
            mw.insert(0, "program", program_name)
            mw.insert(0, "cohort", cohort_name)
            h3_rows.append(mw)
        h3_df = pd.concat(h3_rows, ignore_index=True) if h3_rows else pd.DataFrame()
        if len(h3_df) > 0:
            fdr_h3 = bh_fdr(h3_df["p_raw"].values)
            h3_df["p_bh_fdr_family"] = fdr_h3["p_bh_fdr"].values
            h3_df["reject_bh_fdr_family"] = fdr_h3["reject_bh_fdr"].values
            h3_df["program_description"] = h3_df["program"].map(PROGRAM_DESCRIPTIONS)
            h3_df["drugs"] = h3_df["program"].map(lambda p: "; ".join(PROGRAM_TO_DRUGS.get(p, [])))

    # --- Grade-adjusted association via Ridge regression or linear model with grade covariate
    grade_adj = grade_adjusted_association(
        program_scores_z=program_scores_z,
        meta=meta,
        subgroup_col=subgroup_col,
        grade_col=grade_col,
        cohort_name=cohort_name,
    )

    return {
        "program_scores_raw": program_scores,
        "program_scores_z": program_scores_z,
        "h2_subgroup_kw": h2_df,
        "h2_subgroup_pairwise": pairwise_df,
        "h3_nf2": h3_df,
        "grade_adjusted": grade_adj,
    }


def grade_adjusted_association(
    program_scores_z: pd.DataFrame,
    meta: pd.DataFrame,
    subgroup_col: str,
    grade_col: str,
    cohort_name: str,
) -> pd.DataFrame:
    rows = []
    for program in program_scores_z.columns:
        df = pd.DataFrame({"y": program_scores_z[program]})
        df[subgroup_col] = meta[subgroup_col].astype(str)
        df[grade_col] = meta[grade_col].astype(str)
        df = df.dropna()
        if len(df) < 10 or subgroup_col not in df.columns or grade_col not in df.columns:
            continue
        try:
            model = ols("y ~ C(" + grade_col + ") + C(" + subgroup_col + ")", data=df).fit()
            anova = sm.stats.anova_lm(model, typ=2)
            subgroup_f = anova.loc[f"C({subgroup_col})", "F"]
            subgroup_p = anova.loc[f"C({subgroup_col})", "PR(>F)"]
            grade_f = anova.loc[f"C({grade_col})", "F"]
            grade_p = anova.loc[f"C({grade_col})", "PR(>F)"]
            rows.append(
                {
                    "cohort": cohort_name,
                    "program": program,
                    "subgroup_col": subgroup_col,
                    "grade_col": grade_col,
                    "n_samples": int(len(df)),
                    "subgroup_F_given_grade": subgroup_f,
                    "subgroup_p_given_grade": subgroup_p,
                    "grade_F": grade_f,
                    "grade_p": grade_p,
                    "model_r2": model.rsquared,
                }
            )
        except Exception as e:
            logger.warning("LM failed for %s: %s", program, e)
    out = pd.DataFrame(rows)
    if len(out) > 0:
        fdr = bh_fdr(out["subgroup_p_given_grade"].values)
        out["subgroup_p_bh_fdr_given_grade"] = fdr["p_bh_fdr"].values
        out["subgroup_reject_bh_fdr_given_grade"] = fdr["reject_bh_fdr"].values
    return out


def save_aim2_tables(
    results: Dict[str, pd.DataFrame],
    prefix: str = "table_s2_aim2",
    out_dir: Optional[Path] = None,
) -> Dict[str, Path]:
    if out_dir is None:
        out_dir = PATHS.results_tables
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, df in results.items():
        if isinstance(df, pd.DataFrame):
            p = out_dir / f"{prefix}_{key}.csv"
            df.to_csv(p, index=(df.index.name is not None))
            paths[key] = p
    return paths
