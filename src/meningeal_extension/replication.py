from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .config import PATHS
from .stats import direction_concordance, bh_fdr

logger = logging.getLogger(__name__)


def run_replication(
    discovery_result: Dict[str, pd.DataFrame],
    replication_result: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    H4 external replication. Takes aim2 outputs for discovery and replication,
    reports effect size/direction concordance for:
      - H2 subgroup KW effect sign (medians)
      - H3 NF2 effect sign (delta medians)
    """
    h2_disc = discovery_result.get("h2_subgroup_kw", pd.DataFrame())
    h2_rep = replication_result.get("h2_subgroup_kw", pd.DataFrame())
    h3_disc = discovery_result.get("h3_nf2", pd.DataFrame())
    h3_rep = replication_result.get("h3_nf2", pd.DataFrame())

    # H2 direction: sign of H statistic (we compare signed mean z-score differences
    # vs the global null of zero difference)
    # We'll use the signed delta median from the pairwise table between the two
    # extreme groups as the canonical effect direction.
    pw_disc = discovery_result.get("h2_subgroup_pairwise", pd.DataFrame())
    pw_rep = replication_result.get("h2_subgroup_pairwise", pd.DataFrame())

    h2_effects = _collect_pairwise_effects(pw_disc, pw_rep)
    h2_conc = pd.DataFrame()
    if not h2_effects.empty:
        h2_conc = direction_concordance(
            h2_effects.set_index(["program", "group_a", "group_b"])["delta_discovery"],
            h2_effects.set_index(["program", "group_a", "group_b"])["delta_replication"],
        )
        h2_conc.insert(0, "hypothesis", "H2_subgroup_pairwise_delta_median")

    h3_effects = pd.merge(
        h3_disc[["program", "delta_median_a_minus_b", "p_raw"]].rename(
            columns={
                "delta_median_a_minus_b": "delta_discovery",
                "p_raw": "p_discovery",
            }
        ),
        h3_rep[["program", "delta_median_a_minus_b", "p_raw"]].rename(
            columns={
                "delta_median_a_minus_b": "delta_replication",
                "p_raw": "p_replication",
            }
        ),
        on="program",
        how="inner",
    )
    h3_conc = pd.DataFrame()
    if not h3_effects.empty:
        h3_conc = direction_concordance(
            h3_effects.set_index("program")["delta_discovery"],
            h3_effects.set_index("program")["delta_replication"],
        )
        h3_conc.insert(0, "hypothesis", "H3_NF2_delta_median")

    conc = pd.concat([h2_conc, h3_conc], ignore_index=True) if (len(h2_conc) or len(h3_conc)) else pd.DataFrame()

    # Build combined replication table: per-program forest-plot-ready row
    forest_rows = []
    for _, row in h3_effects.iterrows():
        forest_rows.append(
            {
                "analysis": "H3_NF2",
                "program": row["program"],
                "delta_discovery": row["delta_discovery"],
                "p_discovery": row["p_discovery"],
                "delta_replication": row["delta_replication"],
                "p_replication": row["p_replication"],
                "same_direction": (
                    np.sign(row["delta_discovery"]) == np.sign(row["delta_replication"])
                    or row["delta_discovery"] == row["delta_replication"] == 0
                ),
            }
        )
    for _, row in h2_effects.iterrows():
        forest_rows.append(
            {
                "analysis": f"H2_pairwise::{row['group_a']}-{row['group_b']}",
                "program": row["program"],
                "delta_discovery": row["delta_discovery"],
                "p_discovery": row["p_discovery"],
                "delta_replication": row["delta_replication"],
                "p_replication": row["p_replication"],
                "same_direction": (
                    np.sign(row["delta_discovery"]) == np.sign(row["delta_replication"])
                    or row["delta_discovery"] == row["delta_replication"] == 0
                ),
            }
        )
    forest_df = pd.DataFrame(forest_rows)
    return {
        "concordance_summary": conc,
        "per_effect_forest": forest_df,
        "h2_pairwise_effects": h2_effects,
        "h3_effects": h3_effects,
    }


def _collect_pairwise_effects(pw_disc, pw_rep):
    if pw_disc.empty or pw_rep.empty:
        return pd.DataFrame()
    cols = ["program", "group_a", "group_b"]
    merged = pd.merge(
        pw_disc[cols + ["delta_median_a_minus_b", "p_raw"]].rename(
            columns={
                "delta_median_a_minus_b": "delta_discovery",
                "p_raw": "p_discovery",
            }
        ),
        pw_rep[cols + ["delta_median_a_minus_b", "p_raw"]].rename(
            columns={
                "delta_median_a_minus_b": "delta_replication",
                "p_raw": "p_replication",
            }
        ),
        on=cols,
        how="inner",
    )
    return merged


def save_replication_tables(replication, out_dir: Optional[Path] = None) -> Dict[str, Path]:
    if out_dir is None:
        out_dir = PATHS.results_tables
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, df in replication.items():
        if isinstance(df, pd.DataFrame):
            p = out_dir / f"table_s4_aim3_replication_{key}.csv"
            df.to_csv(p, index=False)
            paths[key] = p
    return paths
