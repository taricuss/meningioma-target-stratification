from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np

from .config import NF2_EXPORT_COL, NF2_INTERNAL_COL, PATHS, SEED
from .datasets import (
    CohortRecord,
    build_synthetic_test_cohort,
    cbioportal_list_meningioma_studies,
    geo_screen_meningioma_expressions,
    load_cohort,
    save_aim0_lock,
    save_cohort,
)
from .concordance import run_concordance, save_concordance_tables
from .target_association import run_aim2, save_aim2_tables
from .replication import run_replication, save_replication_tables
from .survival import run_survival_analysis, save_survival_tables
from .visualization import all_figures


def _resolve_nf2_col(meta: pd.DataFrame) -> Optional[str]:
    if NF2_INTERNAL_COL in meta.columns:
        return NF2_INTERNAL_COL
    if NF2_EXPORT_COL in meta.columns:
        return NF2_EXPORT_COL
    return None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-20s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stderr),
        logging.FileHandler(PATHS.logs / f"run_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.log", mode="w"),
    ],
)
logger = logging.getLogger("meningeal_extension.cli")

_REQUIRED_COHORTS = {
    "discovery_nassiri": {
        "expected_accession": "GSE180061",
        "fetch_script": "scripts/fetch_geo_gse180061.py",
        "description": "Nassiri 2021 discovery cohort (GSE180061 methylation + cBioPortal mRNA)",
    },
    "discovery_bi": {
        "expected_accession": "GSE212666",
        "fetch_script": "scripts/fetch_geo_gse212666.py",
        "description": "Bi-lab 2023 discovery cohort (GSE212666 three-group classification)",
    },
    "replication_gse136661": {
        "expected_accession": "GSE136661",
        "fetch_script": "scripts/fetch_geo_replication.py",
        "description": "Replication cohort (GSE136661 or other suitable GEO dataset)",
    },
}


def _setup_standard_parser(prog: str, desc: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog, description=desc)
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed (default: 20260317)")
    parser.add_argument(
        "--unit-test-synthetic-only",
        action="store_true",
        help=(
            "[DEVELOPMENT ONLY — NEVER FOR RESULTS] Generate synthetic seeded cohorts "
            "for unit-testing the pipeline. Outputs must NOT be cited, copied to "
            "manuscript/, or presented as findings. Synthetic data has ground-truth "
            "effects baked in — tests the STATISTICAL MACHINERY only, never biology."
        ),
    )
    return parser


def _check_cohort_presence(cohort_name: str) -> tuple[bool, str]:
    cohort_dir = PATHS.data_processed / cohort_name
    meta_path = cohort_dir / "metadata.csv"
    expr_path = cohort_dir / "expression.csv"
    has_meta = meta_path.is_file()
    has_expr = expr_path.is_file()

    if not has_meta and not has_expr:
        return False, f"Cohort directory {cohort_dir} has no metadata.csv or expression.csv"
    if not has_meta:
        return False, f"Missing {meta_path}"
    if not has_expr:
        return False, f"Missing {expr_path}"

    meta = pd.read_csv(meta_path, index_col=0)
    sample_ids = meta.index.tolist()
    if len(sample_ids) == 0:
        return False, f"{meta_path} has 0 rows"

    any_synthetic = any(
        "synthetic" in str(sid).lower() or "synth" in str(sid).lower()
        for sid in sample_ids
    )
    if any_synthetic:
        hint = (
            f"[INTEGRITY] {cohort_name}: sample_id column contains synthetic-labeled "
            f"rows. If you intended to run real data, delete {cohort_dir} and run "
            f"{_REQUIRED_COHORTS[cohort_name]['fetch_script']}."
        )
        return False, hint

    return True, f"OK: {len(meta)} samples, {len(pd.read_csv(expr_path, index_col=0).columns)} genes"


def _require_real_cohorts() -> None:
    failures = []
    for cohort_name, info in _REQUIRED_COHORTS.items():
        ok, msg = _check_cohort_presence(cohort_name)
        if not ok:
            failures.append(f"  [{cohort_name}] {info['description']}")
            failures.append(f"       -> {msg}")
            failures.append(f"       -> Run: python {info['fetch_script']}")
    if failures:
        error_block = "\n".join(
            [
                "",
                "=" * 78,
                "ERROR: Real data has not been acquired.",
                "=" * 78,
                "The pipeline will NOT silently fall back to synthetic data.",
                "Synthetic seeded outputs were quarantined to dryrun_synthetic_DO_NOT_CITE/.",
                "",
                "Missing/failed cohorts:",
                *failures,
                "",
                "If you genuinely want to run unit-tests against synthetic data",
                "(e.g. to debug stats.py or ssgsea.py logic against known ground truth),",
                "re-run with the explicit opt-in flag:",
                "    --unit-test-synthetic-only",
                "and accept that all outputs are engineering fixtures, not results.",
                "=" * 78,
                "",
            ]
        )
        logger.error(error_block)
        sys.exit(3)


def _write_synthetic_for_unit_tests(seed: int) -> None:
    logger.warning("=" * 72)
    logger.warning("EXPLICIT --unit-test-synthetic-only FLAG SET.")
    logger.warning("All output from this run is SEEDED SYNTHETIC DATA.")
    logger.warning("Nothing generated below is a finding or result.")
    logger.warning("Copying to results/ or manuscript/ is a fabrication risk.")
    logger.warning("=" * 72)
    disc_nassiri = build_synthetic_test_cohort(n=185, cohort_name="nassiri_synthetic_utest", seed=seed)
    save_cohort("discovery_nassiri", disc_nassiri["metadata"], disc_nassiri["expression"])
    disc_bi = build_synthetic_test_cohort(n=565, cohort_name="bi_synthetic_utest", seed=seed + 1)
    save_cohort("discovery_bi", disc_bi["metadata"], disc_bi["expression"])
    rep = build_synthetic_test_cohort(n=242, cohort_name="replication_synthetic_utest", seed=seed + 7)
    save_cohort("replication_gse136661", rep["metadata"], rep["expression"])


def aim0_main(argv: Optional[list] = None) -> int:
    parser = _setup_standard_parser("aim0-cohort-lock", "Aim 0 — Cohort lock-in (real data only unless --unit-test-synthetic-only).")
    args = parser.parse_args(argv)

    logger.info("Aim 0 started: screening GEO candidates and cBioPortal meningioma studies")

    geo_screen = geo_screen_meningioma_expressions()
    cbio = cbioportal_list_meningioma_studies()
    geo_screen.to_csv(PATHS.results_tables / "aim0_geo_screening.csv", index=False)
    cbio.to_csv(PATHS.results_tables / "aim0_cbio_studies.csv", index=False)
    logger.info("GEO screening: %d candidates", len(geo_screen))
    logger.info("cBioPortal studies: %d meningioma entries", len(cbio))

    records = [
        CohortRecord(
            accession="GSE180061",
            title="Nassiri et al. 2021 — DNA methylation profiling and multi-omics of 185 meningiomas",
            source="Nassiri 2021 Nature; cBioPortal: mng_utoronto_2021",
            role="discovery",
            platform="Illumina Infinium HumanMethylation450 + bulk RNA-seq + WES + snRNA-seq",
            n_samples=185,
            has_nf2=True,
            has_grade=True,
            has_recurrence=True,
            notes="EGAS00001004982 (WES/mRNA/snRNA) controlled access; processed data via cBioPortal or GEO idat (GSE180061).",
        ),
        CohortRecord(
            accession="GSE212666",
            title="Choudhury/Bi lab Cancer Cell 2023 — 565 meningiomas (Merlin-intact / immune-enriched / hypermitotic)",
            source="Bi-lab Cancer Cell 2023",
            role="discovery",
            platform="Illumina EPIC methylation + RNA-seq + proteomics + scRNA-seq",
            n_samples=565,
            has_nf2=True,
            has_grade=True,
            has_recurrence=False,
            notes="Open GEO accession. Subgroup labels are the primary 3-tier classification.",
        ),
        CohortRecord(
            accession="GSE136661",
            title="Meningioma expression GEO candidate 1 — replication",
            source="GEO",
            role="replication_candidate",
            platform="GPL570 [HG-U133_Plus_2]",
            n_samples=geo_screen.set_index("accession").loc["GSE136661", "n_samples_expected"]
            if "GSE136661" in geo_screen.index else 0,
            has_nf2=False,
            has_grade=False,
            has_recurrence=False,
            notes="Requires real GEO download via scripts/fetch_geo_replication.py. No synthetic fallback.",
        ),
    ]
    save_aim0_lock(records)
    with open(PATHS.results_tables / "aim0_cohort_lock.json", "w") as fh:
        json.dump([r.__dict__ for r in records], fh, indent=2, default=str)

    if args.unit_test_synthetic_only:
        _write_synthetic_for_unit_tests(args.seed)
    else:
        _require_real_cohorts()

    logger.info("Aim 0 complete — cohort lock written. Cohort presence verified (real).")
    return 0


def aim1_main(argv: Optional[list] = None) -> int:
    import os
    # Force UTF-8 stdout on Windows to avoid UnicodeEncodeError for any
    # remaining Unicode characters in logs. If PYTHONIOENCODING isn't set
    # at shell startup, reconfiguring sys.stdout here is the next-best
    # fallback. (cp1252 on Windows cmd.exe is the primary failure mode.)
    import sys as _sys
    if hasattr(_sys.stdout, "reconfigure"):
        try:
            _sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            _sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    parser = _setup_standard_parser("aim1-concordance", "Aim 1 — Classifier concordance (H1).")
    args = parser.parse_args(argv)

    if not args.unit_test_synthetic_only:
        # Aim 1 ONLY requires Nassiri discovery cohort (with both nassiri_group
        # and bi_group columns). The Bi-lab cohorts are NOT needed for H1:
        #   - discovery_bi was scoped out (no phenotype labels on GSE212666).
        #   - replication cohorts are for Aim 3, not Aim 1.
        # So we require ONLY discovery_nassiri to pass the real-data gate,
        # not the full set. This fixes a false-positive failure of
        # _require_real_cohorts() on the official H1 entrypoint.
        ok, msg = _check_cohort_presence("discovery_nassiri")
        if not ok:
            info = _REQUIRED_COHORTS["discovery_nassiri"]
            logger.error(
                "\n================================================================\n"
                "Aim 1 blocked: real Nassiri cohort not acquired.\n"
                "  %s\n"
                "  -> %s\n"
                "  -> Run: python %s\n"
                "Synthetic fallback requires --unit-test-synthetic-only (explicit opt-in).\n"
                "================================================================",
                info["description"], msg, info["fetch_script"],
            )
            _sys.exit(3)

    nassiri = load_cohort("discovery_nassiri")
    meta = nassiri["metadata"]
    if "bi_group" not in meta.columns or "nassiri_group" not in meta.columns:
        logger.error(
            "Aim 1 requires both nassiri_group and bi_group in metadata. "
            "Current columns: %s. Run the bridging classifier first "
            "(python -m meningeal_extension.bridging_classifier or "
            "scripts/00_* to populate bi_group from the bridging predictions).",
            list(meta.columns),
        )
        return 2
    cohort_label = "nassiri_synthetic_utest" if args.unit_test_synthetic_only else "nassiri_real"

    # --- H1A: Concordance -----------------------------------------------------
    # Pass the pre-registered Nassiri->Bi mapping explicitly so that
    # run_concordance can compute a legitimate 3x3 Cohen's kappa on the
    # same-label-space table (otherwise label sets are disjoint, kappa was
    # silently returning 0.0 before the guard in stats.py was added).
    from .concordance import test_pairing_specificity, save_pairing_specificity_tables
    from .config import NASSIRI_BI_MATCH_HYPOTHESIS, NASSIRI_GROUP_ORDER, BI_GROUP_ORDER

    result = run_concordance(
        meta, "nassiri_group", "bi_group",
        cohort_name=cohort_label,
        label_a_to_b_mapping=dict(NASSIRI_BI_MATCH_HYPOTHESIS),
    )
    save_concordance_tables(result)
    logger.info(
        "Aim1 H1A concordance: n=%d, 3x3 mapped Cohen's kappa=%.4f (p=%.3g), "
        "4x3 table Cramer's V=%.3f (p=%.3g)",
        result.n_paired, result.cohens_kappa, result.p_cohens_kappa,
        result.cramers_v, result.p_chi2,
    )
    if result.cohens_kappa is None or not np.isfinite(result.cohens_kappa):
        logger.error(
            "H1A kappa is NaN. This means label_a_to_b_mapping projection "
            "failed. Check NASSIRI_BI_MATCH_HYPOTHESIS keys match "
            "nassiri_group levels in metadata: %s",
            sorted(meta["nassiri_group"].dropna().unique().tolist()),
        )
        return 4

    # --- H1B: Pairing-specificity falsifiability ----------------------------
    # Scores the literature-registered true mapping against the full
    # enumeration of 12 structurally equivalent 4-to-3 pairings (Scenario C:
    # Hypermitotic is the merge target, 6 merge choices x 2 singleton perms =
    # 12, minus the true one = 11 wrong alternatives).
    psr = test_pairing_specificity(
        meta,
        nassiri_col="nassiri_group",
        true_bi_col="bi_group",
        true_mapping=dict(NASSIRI_BI_MATCH_HYPOTHESIS),
        merge_target="Hypermitotic",
        cohort_name=cohort_label,
    )
    save_pairing_specificity_tables(psr)
    logger.info(
        "Aim1 H1B pairing-specificity | tested=%d pairings | true rank(kappa)=#%d "
        "rank(agr)=#%d | delta_kappa_vs_best_alt=%+.4f (best_alt kappa=%.4f, true kappa=%.4f) | "
        "singleton-swap placebo kappa=%.4f | worst-alt-merge kappa=%.4f | "
        "pass_strict=%s pass_any_alt=%s",
        psr.n_pairings_tested,
        psr.true_mapping_rank_kappa, psr.true_mapping_rank_agreement,
        psr.delta_kappa_vs_best_alt, psr.best_alt_kappa, psr.true_mapping_kappa,
        psr.singleton_swap_kappa, psr.worst_merge_kappa,
        psr.pass_strict, psr.pass_any_alt,
    )
    # Fail loudly if the strict falsifiability criterion is not met — H1 must
    # rank #1 among all 12 pairings by BOTH kappa and raw agreement.
    if not psr.pass_strict:
        logger.error(
            "Aim1 H1 FAILS strict pairing-specificity gate. True mapping rank "
            "(kappa)=#%d, rank(agreement)=#%d (expected both = #1). "
            "Full ranked pairings written to results/tables/"
            "table_s1b_aim1_pairing_specificity_all_pairings_ranked.csv.",
            psr.true_mapping_rank_kappa, psr.true_mapping_rank_agreement,
        )
        return 5
    logger.info("Aim1 H1B: STRICT PASS (true pairing #1 by kappa AND agreement).")

    if args.unit_test_synthetic_only:
        logger.warning("[UNIT-TEST ONLY] Aim1 numbers above come from synthetic seeded data, not GSE180061.")
    return 0


def aim2_main(argv: Optional[list] = None) -> int:
    parser = _setup_standard_parser("aim2-target-assoc", "Aim 2 — Target-gene association (H2, H3).")
    args = parser.parse_args(argv)

    import sys as _sys

    if not args.unit_test_synthetic_only:
        # Aim 2 is fully testable on the Nassiri discovery cohort alone:
        #   H2 (subgroup-differential programs) runs on Nassiri native 4-group labels AND on
        #     the bridging-predicted Bi 3-group labels (both live in Nassiri metadata.csv).
        #   H3 (NF2 CNA-proxy association with merlin-intact / NF2-loss programs) uses the
        #     Nassiri CNA-derived NF2 proxy column — also independent of any Bi expression data.
        # The Bi-lab discovery cohort is NOT required for the manuscript H2/H3 claims.
        # If discovery_bi is later acquired, we'll re-run and compare the within-Bi
        # subgroup enrichments as a validation layer.
        ok_n, msg_n = _check_cohort_presence("discovery_nassiri")
        if not ok_n:
            info_n = _REQUIRED_COHORTS["discovery_nassiri"]
            logger.error(
                "\n================================================================\n"
                "Aim 2 blocked: real Nassiri cohort not acquired.\n"
                "  %s\n"
                "  -> %s\n"
                "  -> Run: python %s\n"
                "Synthetic fallback requires --unit-test-synthetic-only (explicit opt-in).\n"
                "================================================================",
                info_n["description"], msg_n, info_n["fetch_script"],
            )
            _sys.exit(3)
        ok_b, _ = _check_cohort_presence("discovery_bi")
        bi_ok = bool(ok_b)
        if not bi_ok:
            logger.warning(
                "Aim 2: Bi-lab discovery cohort (discovery_bi) not yet acquired. "
                "Nassiri-native H2+ are computed on both native 4-group labels and the "
                "bridging-predicted 3-group labels. Bi-lab within-cohort H2 is omitted "
                "from tables (declared not-testable in manuscript Limitations)."
            )
    else:
        # Synthetic mode: both cohorts are generated via build_synthetic_test_cohort
        # inside load_cohort, so we can safely attempt both.
        bi_ok = True

    nassiri = load_cohort("discovery_nassiri")
    meta_n = nassiri["metadata"]
    expr_n = nassiri["expression"]
    label_n = "nassiri_synthetic_utest" if args.unit_test_synthetic_only else "nassiri_real"
    label_b = "bi_synthetic_utest" if args.unit_test_synthetic_only else "bi_real"

    res_n_native = run_aim2(expr_n, meta_n, subgroup_col="nassiri_group", cohort_name=f"{label_n}_native4grp")
    save_aim2_tables(res_n_native, prefix="table_s2a_aim2_nassiri_native4grp")

    if "bi_group" in meta_n.columns:
        res_n_bi = run_aim2(expr_n, meta_n, subgroup_col="bi_group", cohort_name=f"{label_n}_predicted_bi3grp")
        save_aim2_tables(res_n_bi, prefix="table_s2b_aim2_nassiri_predicted_bi3grp")
    else:
        logger.warning("No bi_group column in Nassiri metadata — skipping 3-group-predicted Aim 2 analysis.")
        res_n_bi = None

    logger.info(
        "Aim2 Nassiri H2 (native 4-group) KW BH-FDR rejects at alpha=0.05: %d/%d programs",
        int(res_n_native["h2_subgroup_kw"]["kruskal_wallis_reject_bh_fdr"].sum()),
        len(res_n_native["h2_subgroup_kw"]),
    )
    if res_n_bi is not None:
        logger.info(
            "Aim2 Nassiri H2 (predicted Bi 3-group) KW BH-FDR rejects: %d/%d programs",
            int(res_n_bi["h2_subgroup_kw"]["kruskal_wallis_reject_bh_fdr"].sum()),
            len(res_n_bi["h2_subgroup_kw"]),
        )
        h3_df = res_n_bi.get("h3_nf2", pd.DataFrame())
        if len(h3_df) > 0 and "reject_bh_fdr_family" in h3_df.columns:
            n_h3_rej = int(h3_df["reject_bh_fdr_family"].sum())
            logger.info("Aim2 Nassiri H3 (NF2 CNA proxy) MW BH-FDR rejects: %d/%d", n_h3_rej, len(h3_df))
        else:
            logger.warning(
                "Aim2 Nassiri H3: no valid NF2 testable comparisons in h3_nf2 output. "
                "Declared not-testable on Nassiri metadata. h3_nf2 keys present: %s",
                list(res_n_bi.keys()),
            )

    if bi_ok:
        try:
            bi = load_cohort("discovery_bi")
            res_b = run_aim2(bi["expression"], bi["metadata"], subgroup_col="bi_group", cohort_name=label_b)
            save_aim2_tables(res_b, prefix="table_s2_aim2_bi")
            logger.info(
                "Aim2 Bi-lab H2 KW BH-FDR rejects: %d/%d",
                int(res_b["h2_subgroup_kw"]["kruskal_wallis_reject_bh_fdr"].sum()),
                len(res_b["h2_subgroup_kw"]),
            )
        except Exception as exc:
            logger.warning("Aim2 Bi-lab analysis failed (graceful fallback; declared not-testable): %s", exc)

    if args.unit_test_synthetic_only:
        logger.warning("[UNIT-TEST ONLY] Aim2 numbers above come from synthetic seeded data.")
    return 0


def aim3_main(argv: Optional[list] = None) -> int:
    parser = _setup_standard_parser("aim3-replication", "Aim 3 — External replication (H4).")
    args = parser.parse_args(argv)

    if not args.unit_test_synthetic_only:
        _require_real_cohorts()

    disc = load_cohort("discovery_nassiri")
    rep = load_cohort("replication_gse136661")
    label_disc = "nassiri_synthetic_utest" if args.unit_test_synthetic_only else "nassiri_real"
    label_rep = "replication_synthetic_utest" if args.unit_test_synthetic_only else "replication_real"
    res_disc = run_aim2(disc["expression"], disc["metadata"],
                        subgroup_col="nassiri_group", cohort_name=label_disc)
    res_rep = run_aim2(rep["expression"], rep["metadata"],
                       subgroup_col="nassiri_group", cohort_name=label_rep)
    replication = run_replication(res_disc, res_rep)
    save_replication_tables(replication)
    logger.info("Aim3 replication concordance summary:\n%s",
                replication["concordance_summary"].to_string(index=False))
    if args.unit_test_synthetic_only:
        logger.warning("[UNIT-TEST ONLY] Aim3 replication above is synthetic vs synthetic.")
    return 0


def aim4_main(argv: Optional[list] = None) -> int:
    parser = _setup_standard_parser("aim4-survival", "Aim 4 — Optional exploratory survival (H5).")
    args = parser.parse_args(argv)

    if not args.unit_test_synthetic_only:
        _require_real_cohorts()

    from .ssgsea import score_all_programs, zscore_per_program

    nassiri = load_cohort("discovery_nassiri")
    meta = nassiri["metadata"]
    scores_z = zscore_per_program(score_all_programs(nassiri["expression"]))
    surv = run_survival_analysis(meta, scores_z)
    save_survival_tables(surv)
    logger.info("Aim4 H5 survival: powered=%s, n_events=%d. Declaration: %s",
                surv["powered"], surv["n_events"], surv["message"])
    if args.unit_test_synthetic_only:
        logger.warning("[UNIT-TEST ONLY] Aim4 survival above uses synthetic event times.")
    return 0


def run_all(argv: Optional[list] = None) -> int:
    parser = _setup_standard_parser("run-all", "End-to-end execution (Aims 0-4 + figures + tables). REAL DATA ONLY unless --unit-test-synthetic-only.")
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args(argv)

    extra_flags = ["--seed", str(args.seed)]
    if args.unit_test_synthetic_only:
        extra_flags.append("--unit-test-synthetic-only")

    ret = aim0_main(extra_flags)
    if ret != 0:
        return ret
    aim1_main(extra_flags)
    aim2_main(extra_flags)
    aim3_main(extra_flags)
    aim4_main(extra_flags)

    from .concordance import run_concordance as _rc
    from .config import NASSIRI_BI_MATCH_HYPOTHESIS
    from .ssgsea import score_all_programs, zscore_per_program as _z

    ctx = {}
    nassiri = load_cohort("discovery_nassiri")
    bi = load_cohort("discovery_bi")
    rep = load_cohort("replication_gse136661")
    ctx["metadata_discovery"] = nassiri["metadata"]
    ctx["metadata_discovery_bi"] = bi["metadata"]
    ctx["aim2"] = {
        "nassiri": run_aim2(nassiri["expression"], nassiri["metadata"], "nassiri_group", cohort_name="nassiri"),
        "bi": run_aim2(bi["expression"], bi["metadata"], "bi_group", cohort_name="bi"),
    }
    ctx["concordance"] = _rc(
        nassiri["metadata"], "nassiri_group", "bi_group", "nassiri",
        label_a_to_b_mapping=dict(NASSIRI_BI_MATCH_HYPOTHESIS),
    )
    res_disc = ctx["aim2"]["nassiri"]
    res_rep = run_aim2(rep["expression"], rep["metadata"], "nassiri_group", cohort_name="replication")
    ctx["replication"] = run_replication(res_disc, res_rep)

    label_cohort_0 = "Nassiri 2021 (n=185, discovery)"
    label_cohort_1 = "Bi-lab 2023 (n=565, discovery)"
    label_cohort_2 = "Replication cohort (GSE136661)"
    if args.unit_test_synthetic_only:
        label_cohort_0 += " [SYNTHETIC UTEST — NOT REAL]"
        label_cohort_1 += " [SYNTHETIC UTEST — NOT REAL]"
        label_cohort_2 += " [SYNTHETIC UTEST — NOT REAL]"
    _write_cohort_demographics_table(
        nassiri["metadata"], bi["metadata"], rep["metadata"],
        labels=[label_cohort_0, label_cohort_1, label_cohort_2],
    )

    if not args.skip_figures:
        figs = all_figures(ctx)
        total = sum(len(v) for v in figs.values())
        logger.info("Generated %d figure assets across %d panels", total, len(figs))
        if args.unit_test_synthetic_only:
            logger.warning("[UNIT-TEST ONLY] Figures depict synthetic seeded distributions. Do NOT cite.")

    manifest = {
        "generated_at_utc": pd.Timestamp.utcnow().isoformat(),
        "seed": args.seed,
        "is_synthetic_unit_test_only": args.unit_test_synthetic_only,
        "data_provenance": "SYNTHETIC SEEDED — NOT REAL DATA" if args.unit_test_synthetic_only else (
            "Real data — verified via sample_id integrity checks (no '*synthetic*' IDs). "
            "Source: GEO GSE180061 / GSE212666 / cBioPortal mng_utoronto_2021."
        ),
        "tables": sorted([p.name for p in PATHS.results_tables.glob("*.csv")]),
        "figures": sorted([p.name for p in PATHS.results_figures.glob("*")]),
    }
    with open(PATHS.project_root / "manuscript" / "artifact_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info("Pipeline complete. is_synthetic_unit_test_only=%s", args.unit_test_synthetic_only)
    return 0


def _write_cohort_demographics_table(m_nas, m_bi, m_rep, labels=None) -> None:
    if labels is None:
        labels = [
            "Nassiri 2021 (n=185, discovery)",
            "Bi-lab 2023 (n=565, discovery)",
            "Replication cohort",
        ]
    blocks = []
    for (name, meta) in zip(labels, [m_nas, m_bi, m_rep]):
        n = len(meta)
        row = {"Cohort": name, "N": n}
        groups_series = meta.get("nassiri_group", pd.Series([], dtype="category"))
        if hasattr(groups_series, "cat"):
            for grp in groups_series.cat.categories.tolist():
                count = int((groups_series == grp).sum())
                row[f"Nassiri: {grp} n(%)"] = f"{count} ({100 * count / n:.1f})"
        if "bi_group" in meta.columns:
            bg = meta["bi_group"]
            cat_list = bg.cat.categories.tolist() if hasattr(bg, "cat") else sorted(bg.dropna().unique().tolist())
            for grp in cat_list:
                count = int((bg == grp).sum())
                row[f"Bi-lab: {grp} n(%)"] = f"{count} ({100 * count / n:.1f})"
        nf2_col = _resolve_nf2_col(meta)
        if nf2_col is not None:
            for lv in ["Intact", "Mutant/Loss"]:
                c = int((meta[nf2_col] == lv).sum())
                row[f"NF2 (CNA proxy) {lv} n(%)"] = f"{c} ({100 * c / n:.1f})"
        if "who_grade" in meta.columns:
            for lv in ["I", "II", "III"]:
                c = int((meta["who_grade"] == lv).sum()) if lv in meta["who_grade"].values else 0
                row[f"WHO {lv} n(%)"] = f"{c} ({100 * c / n:.1f})"
        if "recurrence_event" in meta.columns:
            ev = int(meta["recurrence_event"].astype(bool).sum())
            row[f"Recurrence events n(%)"] = f"{ev} ({100 * ev / n:.1f})"
        blocks.append(row)
    pd.DataFrame(blocks).to_csv(PATHS.results_tables / "table_s0_cohort_demographics.csv", index=False)


if __name__ == "__main__":
    sys.exit(run_all())
