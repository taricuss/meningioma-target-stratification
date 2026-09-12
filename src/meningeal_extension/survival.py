from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .config import NF2_EXPORT_COL, NF2_INTERNAL_COL

logger = logging.getLogger(__name__)


def _resolve_nf2_col(meta: pd.DataFrame) -> Optional[str]:
    if NF2_INTERNAL_COL in meta.columns:
        return NF2_INTERNAL_COL
    if NF2_EXPORT_COL in meta.columns:
        return NF2_EXPORT_COL
    return None

try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    LIFELINES_OK = True
except ImportError:
    LIFELINES_OK = False
    logger.warning("lifelines not installed — survival module will emit placeholder only.")


def run_survival_analysis(
    metadata: pd.DataFrame,
    program_scores_z: pd.DataFrame,
    time_col: str = "recurrence_months",
    event_col: str = "recurrence_event",
    min_events: int = 20,
) -> Dict[str, object]:
    """
    H5 exploratory survival. Strictly powered only if event count >= min_events.
    Otherwise returns a 'not powered' declaration which the manuscript will include verbatim.
    """
    if not LIFELINES_OK:
        return _powered_off("lifelines package is not installed; skipping Cox model fitting.")

    common = metadata.index.intersection(program_scores_z.index)
    meta = metadata.loc[common].copy()
    t = meta[time_col].dropna()
    ev = meta[event_col].dropna()
    n_events = int(ev.sum()) if ev.dtype.kind in "biuf" else int(ev.astype(bool).sum())

    if n_events < min_events:
        msg = (
            f"Survival analysis not sufficiently powered: {n_events} events observed "
            f"vs minimum {min_events} pre-specified; treating H5 as unevaluable and "
            "presented only as an a priori declared exploratory hypothesis."
        )
        return _powered_off(msg, n_events=n_events, min_events=min_events)

    nf2_col = _resolve_nf2_col(meta)

    rows = []
    km_by_program = {}
    for program in program_scores_z.columns:
        df_dict = {
            "T": meta[time_col],
            "E": meta[event_col].astype(int),
            "score_high": (program_scores_z[program] > program_scores_z[program].median()).astype(int),
            "grade": meta.get("who_grade").astype(str),
        }
        if nf2_col is not None:
            df_dict["nf2"] = meta[nf2_col].astype(str)
        else:
            df_dict["nf2"] = np.nan
        df = pd.DataFrame(df_dict).dropna()
        if int(df["E"].sum()) < min_events // 2:
            continue
        # Log-rank high vs low split
        kmf = KaplanMeierFitter()
        mask_high = df["score_high"] == 1
        T1, E1 = df.loc[mask_high, "T"], df.loc[mask_high, "E"]
        T2, E2 = df.loc[~mask_high, "T"], df.loc[~mask_high, "E"]
        lr = logrank_test(T1, T2, E1, E2)
        kmf.fit(T1, E1, label="High")
        kmf.fit(T2, E2, label="Low")
        km_by_program[program] = kmf

        try:
            cph = CoxPHFitter()
            cph_df = df[["T", "E", "score_high"]].copy()
            if df["grade"].nunique() > 1:
                for lv in sorted(df["grade"].unique())[1:]:
                    cph_df[f"grade_{lv}"] = (df["grade"] == lv).astype(int)
            if df["nf2"].nunique() > 1:
                cph_df["nf2_loss"] = (df["nf2"] == "Mutant/Loss").astype(int)
            cph.fit(cph_df, duration_col="T", event_col="E")
            # Robust across lifelines versions
            if hasattr(cph, "summary") and "coef" in cph.summary.columns:
                coef_sh = float(cph.summary.loc["score_high", "coef"])
            elif hasattr(cph, "params_"):
                coef_sh = float(cph.params_["score_high"])
            elif hasattr(cph, "params"):
                coef_sh = float(cph.params["score_high"])
            else:
                coef_sh = np.nan
            hr_score = float(np.exp(coef_sh)) if np.isfinite(coef_sh) else np.nan
            ci = cph.confidence_intervals_
            ci_low = float(np.exp(ci.loc["score_high", ci.columns[0]]))
            ci_high = float(np.exp(ci.loc["score_high", ci.columns[-1]]))
            if hasattr(cph.summary, "loc") and "p" in cph.summary.columns:
                p_cox = float(cph.summary.loc["score_high", "p"])
            else:
                p_cox = np.nan
        except Exception as e:
            logger.warning("Cox failed for %s: %s", program, e)
            hr_score, ci_low, ci_high, p_cox = np.nan, np.nan, np.nan, np.nan

        rows.append(
            {
                "program": program,
                "n_total": int(len(df)),
                "n_events": int(df["E"].sum()),
                "logrank_p": float(lr.p_value),
                "cox_HR_score_high_vs_low": hr_score,
                "cox_HR_95CI_low": ci_low,
                "cox_HR_95CI_high": ci_high,
                "cox_p": p_cox,
                "exploratory_only_flag": True,
            }
        )

    res_df = pd.DataFrame(rows)
    powered = True
    message = (
        "H5 is declared EXPLORATORY per pre-specified protocol (Section 2 H5). "
        "Survival associations are hypothesis-generating only and not used to support "
        "primary causal claims."
    )
    return {
        "powered": powered,
        "n_events": n_events,
        "min_events": min_events,
        "message": message,
        "per_program": res_df,
        "km_objects": km_by_program,
    }


def _powered_off(message: str, n_events: int = 0, min_events: int = 20) -> Dict[str, object]:
    return {
        "powered": False,
        "n_events": n_events,
        "min_events": min_events,
        "message": message,
        "per_program": pd.DataFrame(),
        "km_objects": {},
    }


def save_survival_tables(surv, out_dir: Optional[Path] = None) -> Dict[str, Path]:
    from .config import PATHS

    if out_dir is None:
        out_dir = PATHS.results_tables
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    summary = pd.DataFrame(
        [
            {
                "powered": surv["powered"],
                "n_events": surv["n_events"],
                "min_events_threshold": surv["min_events"],
                "declaration": surv["message"],
            }
        ]
    )
    p = out_dir / "table_s5_aim4_survival_declaration.csv"
    summary.to_csv(p, index=False)
    paths["declaration"] = p
    if isinstance(surv.get("per_program"), pd.DataFrame) and not surv["per_program"].empty:
        p2 = out_dir / "table_s5_aim4_survival_per_program.csv"
        surv["per_program"].to_csv(p2, index=False)
        paths["per_program"] = p2
    return paths
