#!/usr/bin/env python
"""Aggregate Aim 0 report: honest state of ALL real-data cohorts.

Run this AFTER running the individual fetch scripts:
    python scripts/fetch_geo_gse180061.py
    python scripts/fetch_geo_gse212666.py
    python scripts/fetch_cbio_mng_utoronto_2021.py
    python scripts/fetch_geo_replication.py   [one or more accessions]

Output:
    data/processed/aim0_cohort_status_report.txt   (human readable)
    data/processed/aim0_cohort_status_report.json  (machine readable)
    Exit code 0 if at least 1 hypothesis is testable across cohorts;
    exit code 4 if nothing is testable anywhere (still a valid finding — just a
    strong signal that H1-H5 need to be pre-registered as "not evaluable with
    currently available open data").

This script NEVER writes to results/tables or results/figures. It is purely a
diagnostic for Aim 0.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from meningeal_extension.config import PATHS  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-22s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger("aim0_report")

COHORTS_TO_SCAN: List[Dict] = [
    {
        "cohort_key": "discovery_nassiri",
        "expected_label": "Nassiri 2021 (GSE180061 + mng_utoronto_2021)",
        "processed_dir": PATHS.data_processed / "discovery_nassiri",
        "fetch_cmds": [
            "python scripts/fetch_geo_gse180061.py",
            "python scripts/fetch_cbio_mng_utoronto_2021.py",
        ],
    },
    {
        "cohort_key": "discovery_bi",
        "expected_label": "Bi lab 2023 (GSE212666)",
        "processed_dir": PATHS.data_processed / "discovery_bi",
        "fetch_cmds": ["python scripts/fetch_geo_gse212666.py"],
    },
    {
        "cohort_key": "replication_gse136661",
        "expected_label": "Replication (GSE136661)",
        "processed_dir": PATHS.data_processed / "replication_gse136661",
        "fetch_cmds": ["python scripts/fetch_geo_replication.py --accession GSE136661"],
    },
]


def scan_cohort(info: Dict) -> Dict:
    d = Path(info["processed_dir"])
    meta_path = d / "metadata.csv"
    expr_path = d / "expression.csv"
    prov_path = d / "metadata_PROVENANCE.csv"
    report: Dict = {
        "cohort_key": info["cohort_key"],
        "label": info["expected_label"],
        "processed_dir": str(d),
        "exists_dir": d.is_dir(),
        "metadata_csv_exists": meta_path.is_file(),
        "expression_csv_exists": expr_path.is_file(),
        "provenance_csv_exists": prov_path.is_file(),
        "fetch_commands": info["fetch_cmds"],
    }
    if meta_path.is_file():
        meta = pd.read_csv(meta_path, index_col=0)
        report["n_samples_metadata"] = len(meta)
        report["fields_present"] = list(meta.columns)
        report["field_non_nulls"] = {
            str(c): int(meta[c].notna().sum()) for c in meta.columns
        }
        report["field_missing_pct"] = {
            str(c): round(100 * (1 - float(meta[c].notna().sum()) / len(meta)), 2) if len(meta) else 100.0
            for c in meta.columns
        }
        # Per-hypothesis quick assessment
        report["hypothesis_quick"] = _quick_hypotheses(meta)
    else:
        report["note"] = (
            "metadata.csv not present. Run fetch commands listed above BEFORE "
            "running the rest of the pipeline. Synthetic fallback is disabled."
        )
    if expr_path.is_file():
        try:
            expr = pd.read_csv(expr_path, index_col=0, nrows=0)
            report["n_genes_in_expression"] = len(expr.columns)
        except Exception as e:  # noqa: BLE001
            report["n_genes_in_expression"] = f"ERROR reading: {e}"
    if prov_path.is_file():
        try:
            prov = pd.read_csv(prov_path)
            report["provenance"] = prov.to_dict(orient="records")
        except Exception as e:  # noqa: BLE001
            report["provenance"] = f"ERROR reading provenance CSV: {e}"
    return report


def _quick_hypotheses(meta: pd.DataFrame) -> Dict[str, Dict]:
    n = len(meta)

    def _field(name: str, min_non_null: int = 20):
        if name not in meta.columns:
            return {"present": False, "n_non_null": 0, "meets_min": False}
        nn = int(meta[name].notna().sum())
        return {"present": True, "n_non_null": nn, "meets_min": nn >= min_non_null}

    return {
        "H1_subgroup_concordance": {
            "requires": ["nassiri_group", "bi_group"],
            "nassiri_group": _field("nassiri_group"),
            "bi_group": _field("bi_group"),
            "testable_minimum": (
                _field("nassiri_group")["meets_min"] and _field("bi_group")["meets_min"]
            ),
            "note_if_missing": (
                "bi_group on Nassiri samples is the OUTPUT of the Bi-lab classifier "
                "run on mRNA expression. That classifier application step has NOT "
                "happened yet in this repo (it's step between Aim 0 and Aim 1). "
                "If only nassiri_group is present but bi_group is not, H2 is still "
                "testable on the Nassiri cohort using nassiri_group alone — H1 is "
                "simply deferred until both classifier outputs exist."
            ),
        },
        "H2_subgroup_target_enrichment": {
            "requires": ["either nassiri_group or bi_group"],
            "nassiri_group": _field("nassiri_group"),
            "bi_group": _field("bi_group"),
            "testable_minimum": (
                _field("nassiri_group")["meets_min"] or _field("bi_group")["meets_min"]
            ),
        },
        "H3_nf2_target_association": {
            "requires": ["nf2_status"],
            "nf2_status": _field("nf2_status", min_non_null=30),
            "testable_minimum": _field("nf2_status", min_non_null=30)["meets_min"],
        },
        "H5_survival": {
            "requires": ["recurrence_event AND recurrence_months"],
            "recurrence_event": _field("recurrence_event", min_non_null=20),
            "recurrence_months": _field("recurrence_months", min_non_null=20),
            "testable_minimum": (
                _field("recurrence_event", 20)["meets_min"]
                and _field("recurrence_months", 20)["meets_min"]
            ),
        },
    }


def format_report_txt(all_reports: List[Dict]) -> str:
    lines: List[str] = []
    lines.append("=" * 78)
    lines.append("AIM 0 COHORT STATUS REPORT  —  REAL DATA ONLY  —  NO SYNTHETIC FALLBACK")
    lines.append(f"Generated: {pd.Timestamp.now().isoformat()}")
    lines.append("=" * 78)
    overall_any_testable = False
    for r in all_reports:
        lines.append("")
        lines.append(f"─── Cohort: {r['label']} ───")
        lines.append(f"    processed_dir: {r['processed_dir']}")
        lines.append(f"    metadata.csv present: {r['metadata_csv_exists']}")
        lines.append(f"    expression.csv present: {r['expression_csv_exists']}"
                     + (f" ({r.get('n_genes_in_expression','?')} genes)" if r.get("n_genes_in_expression") else ""))
        lines.append(f"    provenance.csv present: {r['provenance_csv_exists']}")
        if not r["metadata_csv_exists"]:
            lines.append(f"    ⚠  MISSING — Run fetch commands:")
            for cmd in r["fetch_commands"]:
                lines.append(f"        $ {cmd}")
            continue
        lines.append(f"    N samples (metadata): {r.get('n_samples_metadata','?')}")
        lines.append(f"    Fields present: {', '.join(str(f) for f in r.get('fields_present',[]))}")
        lines.append("    Field completeness (% non-null):")
        for field, pct in (r.get("field_missing_pct") or {}).items():
            nn = (r.get("field_non_nulls") or {}).get(field, 0)
            comp = 100 - float(pct)
            bar = "█" * int(comp // 5) + "░" * (20 - int(comp // 5))
            lines.append(f"      {field:<24s} {bar} {comp:5.1f}%  ({nn}/{r.get('n_samples_metadata','?')})")
        hq = r.get("hypothesis_quick") or {}
        lines.append("    Hypothesis quick assessment:")
        for h_name, h_info in hq.items():
            ok = h_info.get("testable_minimum", False)
            overall_any_testable = overall_any_testable or bool(ok)
            tag = "✅ TESTABLE" if ok else "❌ NOT TESTABLE"
            lines.append(f"      [{tag}] {h_name}")
            for req_field in h_info.get("requires", []):
                lines.append(f"          requires: {req_field}")
            for k, v in h_info.items():
                if k in {"requires", "testable_minimum", "note_if_missing"}:
                    continue
                if isinstance(v, dict):
                    status = "OK" if v.get("meets_min") else (f"{v.get('n_non_null',0)} non-null" if v.get("present") else "MISSING")
                    lines.append(f"          {k}: {status}")
            if "note_if_missing" in h_info and not ok:
                for ln in str(h_info["note_if_missing"]).splitlines():
                    lines.append(f"          ℹ {ln}")
    lines.append("")
    lines.append("=" * 78)
    lines.append("Overall assessment:")
    lines.append(f"  Any hypothesis testable across any cohort: {overall_any_testable}")
    lines.append("=" * 78)
    lines.append("")
    lines.append("If a field says NOT TESTABLE that is a REAL, REPORTABLE finding.")
    lines.append("Do NOT plug synthetic data into the gap. Pre-register:")
    lines.append('  "H[X] was not evaluable in open data; deferring to authors\' '
                 "Supplementary Tables or controlled-access repositories.\"")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    reports = [scan_cohort(info) for info in COHORTS_TO_SCAN]
    out_txt = format_report_txt(reports)
    PATHS.data_processed.mkdir(parents=True, exist_ok=True)
    txt_path = PATHS.data_processed / "aim0_cohort_status_report.txt"
    txt_path.write_text(out_txt, encoding="utf-8")
    logger.info("\n%s", out_txt)
    json_path = PATHS.data_processed / "aim0_cohort_status_report.json"
    json_path.write_text(json.dumps(reports, indent=2, default=str), encoding="utf-8")
    logger.info("Reports written:\n  %s\n  %s", txt_path, json_path)
    any_testable = any(
        any(v.get("testable_minimum", False)
            for v in (rep.get("hypothesis_quick") or {}).values())
        for rep in reports
        if rep.get("metadata_csv_exists")
    )
    return 0 if any_testable else 4


if __name__ == "__main__":
    sys.exit(main())
