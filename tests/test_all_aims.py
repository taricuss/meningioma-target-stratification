from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from meningeal_extension.gene_programs import (
    ALL_PROGRAMS,
    GENE_PROGRAMS,
    HDAC_RESISTANCE_PROGRAM,
    PROGRAM_DESCRIPTIONS,
)
from meningeal_extension.stats import (
    bh_fdr,
    cliffs_delta as _cliffs_delta,
    cohens_kappa_table,
    direction_concordance,
    kruskal_wallis_across_groups,
    mann_whitney_two_group,
    pairwise_wilcoxon,
)
from meningeal_extension.ssgsea import ssgsea_score, score_all_programs, zscore_per_program
from meningeal_extension.datasets import build_synthetic_test_cohort


def test_gene_programs_frozen():
    assert "HDAC_Panobinostat_Romidepsin" in GENE_PROGRAMS
    assert "HDAC1" in GENE_PROGRAMS["HDAC_Panobinostat_Romidepsin"]
    assert "HDAC2" in GENE_PROGRAMS["HDAC_Panobinostat_Romidepsin"]
    assert "Proteasome_Carfilzomib" in GENE_PROGRAMS
    assert "PSMB5" in GENE_PROGRAMS["Proteasome_Carfilzomib"]
    assert "Tubulin_Ixabepilone" in GENE_PROGRAMS
    assert "TUBB3" in GENE_PROGRAMS["Tubulin_Ixabepilone"]
    assert "Translation_Omacetaxine" in GENE_PROGRAMS
    assert "EEF2" in GENE_PROGRAMS["Translation_Omacetaxine"]
    assert "HDAC8_TGFb_EMT_Resistance" in HDAC_RESISTANCE_PROGRAM
    assert "HDAC8" in HDAC_RESISTANCE_PROGRAM["HDAC8_TGFb_EMT_Resistance"]
    for k in ALL_PROGRAMS:
        assert k in PROGRAM_DESCRIPTIONS


def test_bh_fdr_uniform_and_all_null():
    pvals = np.random.default_rng(0).uniform(0, 1, size=200)
    out = bh_fdr(pvals)
    assert out["p_bh_fdr"].between(0, 1).all()
    pvals_sig = np.array([1e-6, 1e-5, 0.001, 0.2, 0.9])
    out = bh_fdr(pvals_sig)
    assert out["reject_bh_fdr"].sum() >= 3


def test_kruskal_wallis_detects_shift():
    rng = np.random.default_rng(1)
    g = ["A"] * 50 + ["B"] * 50 + ["C"] * 50
    v = np.concatenate([rng.normal(0, 1, 50), rng.normal(1, 1, 50), rng.normal(2, 1, 50)])
    stat, p, k, n = kruskal_wallis_across_groups(v, g)
    assert p < 0.01
    assert k == 3 and n == 150


def test_pairwise_wilcoxon_fdr():
    rng = np.random.default_rng(2)
    g = np.array(["A"] * 30 + ["B"] * 30 + ["C"] * 30)
    v = np.concatenate([rng.normal(0, 1, 30), rng.normal(1.5, 1, 30), rng.normal(0.2, 1, 30)])
    out = pairwise_wilcoxon(v, g, groups_order=["A", "B", "C"])
    assert {"A", "B", "C"} == set(out[["group_a", "group_b"]].melt()["value"].unique())
    # A vs B should be strongest
    ab = out[(out["group_a"] == "A") & (out["group_b"] == "B")].iloc[0]
    assert ab["p_raw"] < 0.01
    assert "p_bh_fdr" in out.columns


def test_cohens_kappa_perfect_and_mid():
    lab_a = ["X"] * 20 + ["Y"] * 20
    lab_b = ["X"] * 20 + ["Y"] * 20
    res = cohens_kappa_table(lab_a, lab_b)
    assert float(res["cohens_kappa"].iloc[0]) == pytest.approx(1.0, abs=1e-6)
    rng = np.random.default_rng(3)
    lab_c = rng.choice(["X", "Y"], size=40)
    lab_d = rng.choice(["X", "Y"], size=40)
    res2 = cohens_kappa_table(lab_c, lab_d)
    assert -0.3 <= float(res2["cohens_kappa"].iloc[0]) <= 0.3
    assert 0.0 <= float(res2["cramers_v"].iloc[0]) <= 1.0


def test_direction_concordance_same_and_opposite():
    base = pd.Series({"P1": 0.5, "P2": -0.3, "P3": 1.2, "P4": -0.8})
    same = pd.Series({"P1": 0.4, "P2": -0.1, "P3": 0.9, "P4": -1.0})
    conc = direction_concordance(base, same)
    assert int(conc["n_same_direction"].iloc[0]) == 4
    opp = pd.Series({"P1": -0.4, "P2": 0.1, "P3": -0.9, "P4": 1.0})
    conc2 = direction_concordance(base, opp)
    assert int(conc2["n_same_direction"].iloc[0]) == 0


def test_ssgsea_detects_enrichment():
    rng = np.random.default_rng(42)
    samples = [f"s{i}" for i in range(30)]
    genes = [f"g{i}" for i in range(100)]
    expr = pd.DataFrame(rng.normal(0, 1, (30, 100)), index=samples, columns=genes)
    # Upregulate signature in first 15 samples
    sig = ["g0", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"]
    expr.loc[samples[:15], sig] += 1.5
    scores = ssgsea_score(expr, sig)
    assert scores.loc[samples[:15]].mean() > scores.loc[samples[15:]].mean()


def test_score_all_programs_shape_and_z():
    rng = np.random.default_rng(5)
    samples = [f"s{i}" for i in range(50)]
    genes = sorted({g for genes in ALL_PROGRAMS.values() for g in genes}) + [f"zz{i}" for i in range(50)]
    expr = pd.DataFrame(rng.normal(5, 1, (50, len(genes))), index=samples, columns=genes)
    scores = score_all_programs(expr)
    assert scores.shape == (50, len(ALL_PROGRAMS))
    z = zscore_per_program(scores)
    assert z.mean().abs().max() < 1e-10
    assert np.isclose(z.std(ddof=0).mean(), 1.0, atol=0.1)


def test_synthetic_cohort_biome_ground_truth_exists():
    cohort = build_synthetic_test_cohort(n=220, seed=7)
    meta, expr = cohort["metadata"], cohort["expression"]
    assert meta.shape[0] == 220
    assert expr.shape[0] == 220
    assert "HDAC1" in expr.columns and "PSMB5" in expr.columns
    # Subgroup differences: check directionally ordered means (ground truth seeded)
    # Hypermetabolic >> NF2-inactivated canonical > Immunogenic for HDAC1/2, PSMB5
    order = ["Hypermetabolic", "NF2-inactivated canonical", "Immunogenic"]
    means = []
    for g in order:
        vals = expr.loc[meta["nassiri_group"] == g, "HDAC1"].values
        means.append(float(vals.mean()))
    # Strictly decreasing (Hypermetabolic highest, Immunogenic lowest)
    assert means[0] > means[1] > means[2], f"HDAC1 means not ordered as seeded: {means}"
    proteo_means = []
    for g in order:
        vals = expr.loc[meta["nassiri_group"] == g, "PSMB5"].values
        proteo_means.append(float(vals.mean()))
    assert proteo_means[0] > proteo_means[1] > proteo_means[2], f"PSMB5 means not ordered: {proteo_means}"
    # Immunogenic >> others for the HDAC8/TGFb/EMT resistance axis
    res = []
    for g in ["Immunogenic", "Hypermetabolic", "NF2-inactivated canonical"]:
        vals = expr.loc[meta["nassiri_group"] == g, "HDAC8"].values
        res.append(float(vals.mean()))
    assert res[0] > res[1] and res[0] > res[2], f"HDAC8 Immunogenic should be highest: {res}"
    # NF2 mutant/loss >> intact for HDAC1, PSMB5 (seeded delta +0.6)
    nf2_mut = expr.loc[meta["nf2_status"] == "Mutant/Loss", "HDAC1"].mean()
    nf2_int = expr.loc[meta["nf2_status"] == "Intact", "HDAC1"].mean()
    assert nf2_mut > nf2_int + 0.2, f"NF2 effect on HDAC1 not in seeded direction: {nf2_mut} vs {nf2_int}"


def test_aim2_end_to_end_runs_and_returns_tables():
    from meningeal_extension.target_association import run_aim2
    cohort = build_synthetic_test_cohort(n=100, seed=9)
    out = run_aim2(cohort["expression"], cohort["metadata"], cohort_name="test")
    assert out["h2_subgroup_kw"].shape[0] == len(ALL_PROGRAMS)
    assert "kruskal_wallis_p_bh_fdr" in out["h2_subgroup_kw"].columns
    assert not out["h3_nf2"].empty
    assert "p_bh_fdr_family" in out["h3_nf2"].columns
    assert "grade_adjusted" in out
    for col in ["subgroup_p_given_grade", "subgroup_p_bh_fdr_given_grade"]:
        assert col in out["grade_adjusted"].columns


def test_replication_end_to_end_runs():
    from meningeal_extension.target_association import run_aim2
    from meningeal_extension.replication import run_replication
    d = build_synthetic_test_cohort(n=120, seed=11)
    r = build_synthetic_test_cohort(n=140, seed=12)
    od = run_aim2(d["expression"], d["metadata"], cohort_name="disc")
    orr = run_aim2(r["expression"], r["metadata"], cohort_name="rep")
    res = run_replication(od, orr)
    assert "concordance_summary" in res
    assert "per_effect_forest" in res
    assert len(res["per_effect_forest"]) > 0
