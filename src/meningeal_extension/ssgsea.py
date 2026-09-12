from __future__ import annotations

import logging
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .gene_programs import ALL_PROGRAMS

logger = logging.getLogger(__name__)


def rank_normalize(expr: pd.DataFrame) -> pd.DataFrame:
    ranked = expr.rank(axis=1, method="average", na_option="keep").astype(float)
    n_cols = expr.shape[1] - expr.isna().sum(axis=1).values[:, None]
    n_cols_safe = np.where(n_cols == 0, 1, n_cols)
    return (ranked - 0.5) / n_cols_safe


def ssgsea_score(
    expr: pd.DataFrame,
    gene_set: Sequence[str],
    alpha: float = 0.25,
) -> pd.Series:
    """Single-sample GSEA using a rank-based, analytic-normalized enrichment score.

    Parameters
    ----------
    expr : pd.DataFrame
        Rows = samples, columns = gene symbols. Values are expression (any scale).
    gene_set : list of str
        Gene symbols in the signature.
    alpha : float
        Exponent for the weighting curve (0 = Kolmogorov-like, 0.25 default per Barbie 2009).

    Returns
    -------
    pd.Series
        ssGSEA score per sample, named by the input index.
    """
    expr_aligned = expr.copy()
    genes_present = [g for g in gene_set if g in expr_aligned.columns]
    missing = [g for g in gene_set if g not in expr_aligned.columns]
    if len(genes_present) == 0:
        logger.warning("No genes from signature present in expression matrix")
        return pd.Series(np.nan, index=expr_aligned.index, dtype=float)
    if len(missing) > 0:
        logger.debug("Missing %d/%d genes: %s", len(missing), len(gene_set), missing[:5])

    # Step 1: ranks per sample
    R = expr_aligned.rank(axis=1, method="average", na_option="keep").astype(float)
    # Step 2: order descending
    sorted_idx = np.argsort(-R.values, axis=1)
    n_samples, n_genes = R.shape
    genes_arr = np.asarray(expr_aligned.columns)
    gene_set_mask = np.isin(genes_arr, genes_present)
    gene_set_ind = np.where(gene_set_mask)[0]

    scores = np.full(n_samples, np.nan, dtype=float)
    for i in range(n_samples):
        order = sorted_idx[i]
        # vectorized step
        N = n_genes
        N_H = len(genes_present)
        if N_H == 0 or N <= N_H:
            continue
        P_GW = np.zeros(N, dtype=float)
        P_NOG = np.zeros(N, dtype=float)
        in_set = np.isin(order, gene_set_ind)
        weights = (np.arange(N, dtype=float) + 1.0) ** alpha
        cum_weight_in = np.cumsum(weights * in_set)
        total_weight_in = (np.arange(1, N_H + 1, dtype=float) ** alpha).sum() if N_H > 0 else 1.0
        P_GW = cum_weight_in / max(total_weight_in, 1e-12)
        P_NOG = np.cumsum(~in_set) / max(N - N_H, 1)
        es = np.max(P_GW - P_NOG) + np.min(P_GW - P_NOG)
        # analytic normalisation by randomisation-based range estimate (Subramanian 2005):
        # approximate null max ES ≈ sum_{i=1..N_H} w_i / sum w - (N_H / N)
        max_null = (weights[:N_H].sum() / max(total_weight_in, 1e-12)) - (N_H / N)
        min_null = -max_null
        denom = (max_null - min_null) if (max_null - min_null) != 0 else 1.0
        nes = es / denom
        scores[i] = float(nes)
    return pd.Series(scores, index=expr_aligned.index, dtype=float)


def score_all_programs(
    expr: pd.DataFrame,
    programs: Optional[Dict[str, List[str]]] = None,
) -> pd.DataFrame:
    if programs is None:
        programs = ALL_PROGRAMS
    out = {}
    for program_name, gene_list in programs.items():
        out[program_name] = ssgsea_score(expr, gene_list)
    return pd.DataFrame(out)


def zscore_per_program(program_scores: pd.DataFrame) -> pd.DataFrame:
    return (program_scores - program_scores.mean()) / program_scores.std(ddof=0)
