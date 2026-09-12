#!/usr/bin/env python
"""Fetch and audit replication GEO candidates (PLAN.md Aim 0 screening shortlist).

AIM 0 REAL DATA ACQUISITION — NO SYNTHETIC FALLBACK.

Default target: GSE136661 (listed first in REPLICATION_CANDIDATE_STUDIES).
The user can override with --accession GSEXXXXX to try any other shortlisted
candidate (GSE77259, GSE115966, GSE94474, etc).

This script:
  1. Downloads one GEO Series SOFT record via GEOparse.
  2. Extracts sample-level GSM metadata (characteristics_ch1 + source_name_ch1).
  3. Runs the SAME canonical-column mapping as the discovery scripts, and
     the SAME hypothesis-feasibility audit.
  4. If ANY of H2/H3/H4 are testable, writes to data/processed/replication_gseXXXXXX/.
     Otherwise, the raw audit is saved to data/raw/ only and the processed
     directory stays empty. This is a real finding ("this cohort can't replicate").

Honesty rule: a replication cohort can't invent subgroup labels. If a GEO dataset
was generated BEFORE the Nassiri/Bi subgroup classifiers existed (most of them),
the labels must be produced by RUNNING the published classifiers on the expression
matrix. That's a separate step — NOT this script's job. This script only records
what the GEO record actually says, and reports honestly what fields are missing.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from meningeal_extension.config import PATHS, REPLICATION_CANDIDATE_STUDIES  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-22s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger("fetch_geo_replication")


@dataclass
class FieldAudit:
    field: str
    n_non_null: int
    n_total: int
    missing_pct: float
    n_unique: int
    examples: List[str]


HYPOTHESIS_MINIMUM_FIELDS = {
    "H4 replication (same-sign direction)": ["nassiri_group OR bi_group"],
    "H2 replication (subgroup enrichment)": ["nassiri_group OR bi_group"],
    "H3 replication (NF2 association)": ["nf2_status"],
    "H5 replication (survival)": ["recurrence_event", "recurrence_months"],
}


def _import_geoparse():
    try:
        import GEOparse  # type: ignore
        return GEOparse
    except ImportError as e:
        logger.critical("GEOparse missing — `pip install -r requirements.txt`")
        raise SystemExit(2) from e


def fetch_gse(accession: str, raw_dir: Path, skip_download: bool):
    GEOparse = _import_geoparse()
    raw_dir.mkdir(parents=True, exist_ok=True)
    if skip_download and any(raw_dir.glob(f"{accession}_family.soft*")):
        logger.info("--skip-download: reusing cached SOFT in %s", raw_dir)
        return GEOparse.get_GEO(filepath=str(next(raw_dir.glob(f"{accession}_family.soft*"))), silent=False)
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(raw_dir), silent=False)
    except Exception as e:
        logger.critical(
            "Fetch failed for %s: %s\n"
            "Manual URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n"
            "Save SOFT to %s and rerun with --skip-download.\n"
            "IMPORTANT: If this candidate cohort is unreachable, TRY THE NEXT ONE in "
            "REPLICATION_CANDIDATE_STUDIES (config.py). Shortlist: %s\n"
            "DO NOT fall back to a synthetic substitute.",
            accession, e, accession, raw_dir, REPLICATION_CANDIDATE_STUDIES,
        )
        raise SystemExit(3) from e
    return gse


def extract_metadata(gse) -> pd.DataFrame:
    rows = []
    for gsm_id, gsm in gse.gsms.items():
        row = {"sample_id": gsm_id, "gsm_title": gsm.metadata.get("title", [""])[0]}
        for c in gsm.metadata.get("characteristics_ch1", []):
            if ":" in c:
                k, v = c.split(":", 1)
                row[k.strip()] = v.strip()
            else:
                row.setdefault("_unparsed", []).append(c)
        for i, s in enumerate(gsm.metadata.get("source_name_ch1", [])):
            row[f"source_name_{i}"] = s
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample_id")


def audit(meta: pd.DataFrame) -> pd.DataFrame:
    n = len(meta)
    audits = []
    for col in meta.columns:
        s = meta[col]
        nn = s.notna() & (s.astype(str).str.strip() != "") & (s.astype(str).str.lower() != "nan")
        audits.append(FieldAudit(
            field=col,
            n_non_null=int(nn.sum()), n_total=n,
            missing_pct=round(100 * (1 - int(nn.sum()) / n), 2) if n else 100,
            n_unique=int(s[nn].nunique()),
            examples=s[nn].astype(str).unique().tolist()[:5],
        ))
    return pd.DataFrame([asdict(a) for a in audits]).sort_values("missing_pct").reset_index(drop=True)


def map_canonical(meta: pd.DataFrame) -> tuple[pd.DataFrame, Dict[str, str]]:
    col_lower = {c.lower().replace(" ", "_").replace("-", "_"): c for c in meta.columns}
    can = pd.DataFrame(index=meta.index)
    maplist = [
        ("nassiri_group", ["subgroup", "group", "classification", "methylation_subgroup",
                           "nassiri_group", "nassiri_subgroup", "cluster", "class",
                           "dna_methylation_class", "molecular_subtype", "molecular_group"]),
        ("bi_group", ["bi_group", "bi_subgroup", "bi_classification", "merlin_group",
                      "immune_group", "subgroup_bi", "epigenomic_subgroup"]),
        ("who_grade", ["who_grade", "grade", "histologic_grade", "histological_grade",
                       "pathology", "tumor_grade", "pathological_grade", "histopathology",
                       "who_2016_grade", "who_2021_grade", "histology_grade"]),
        ("nf2_status", ["nf2_status", "nf2", "nf2_mutation", "merlin_status",
                        "nf2_mutation_status", "nf2_alteration", "merlin"]),
        ("recurrence_event", ["recurrence", "recurred", "progression", "progressed",
                              "rfs_event", "dfs_event", "outcome", "event"]),
        ("recurrence_months", ["rfs_months", "recurrence_months", "follow_up_months",
                               "survival_months", "time_to_recurrence", "dfs_months",
                               "months_to_recurrence", "time_to_event_months", "followup_months"]),
        ("age_at_surgery", ["age", "age_at_surgery", "diagnosis_age", "patient_age",
                            "age_at_diagnosis", "years", "age_years"]),
        ("sex", ["sex", "gender"]),
    ]
    mapper_log: Dict[str, str] = {}
    for canon, cands in maplist:
        found = next((col_lower[c] for c in cands if c in col_lower), None)
        if found:
            can[canon] = meta[found].values
            mapper_log[canon] = f"mapped from raw column '{found}'"
        else:
            can[canon] = np.nan
            mapper_log[canon] = f"NOT MAPPED — no alias match among raw cols: {sorted(col_lower.keys())}"
    if "who_grade" in can.columns:
        can["who_grade"] = can["who_grade"].astype(str).str.upper().str.extract(r"\b([I]{1,3})\b", expand=False)
    if "nf2_status" in can.columns:
        norm = can["nf2_status"].astype(str).str.strip().str.lower().replace({
            "wt": "Intact", "wild": "Intact", "wild_type": "Intact", "intact": "Intact", "normal": "Intact",
            "mutant": "Mutant/Loss", "loss": "Mutant/Loss", "del": "Mutant/Loss",
            "deleted": "Mutant/Loss", "mutated": "Mutant/Loss", "altered": "Mutant/Loss", "loh": "Mutant/Loss",
        })
        norm[~norm.isin(["Intact", "Mutant/Loss"])] = np.nan
        can["nf2_status"] = norm
    if "recurrence_event" in can.columns:
        can["recurrence_event"] = can["recurrence_event"].astype(str).str.strip().str.lower().map({
            "yes": True, "true": True, "1": True, "recurred": True, "progressed": True, "event": True,
            "no": False, "false": False, "0": False, "none": False,
        })
    for nc in ["age_at_surgery", "recurrence_months"]:
        if nc in can.columns:
            can[nc] = pd.to_numeric(can[nc], errors="coerce")
    return can, mapper_log


def feasibility(can: pd.DataFrame) -> Dict[str, Dict]:
    n = len(can)
    out = {}
    for h_name, req in HYPOTHESIS_MINIMUM_FIELDS.items():
        status = {}
        ok = True
        for f in req:
            alt_fields = f.split(" OR ")
            any_present = False
            for af in alt_fields:
                af = af.strip()
                if af in can.columns and can[af].notna().sum() >= 10:
                    any_present = True
                    status[af] = f"{int(can[af].notna().sum())}/{n} non-null"
            if not any_present:
                ok = False
                status[f] = f"insufficient (<10 non-null in any of {alt_fields})"
        out[h_name] = {"testable": ok, "fields": status}
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default="GSE136661",
                        help=f"GEO accession to fetch. Default GSE136661. Shortlist: {REPLICATION_CANDIDATE_STUDIES}")
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()

    accession = args.accession.upper()
    if not accession.startswith("GSE"):
        logger.critical("--accession must start with GSE (e.g. GSE136661). Got: %s", accession)
        return 2

    raw_dir = PATHS.data_raw / accession
    gse = fetch_gse(accession, raw_dir, args.skip_download)

    meta_raw = extract_metadata(gse)
    meta_raw.to_csv(raw_dir / "gsm_all_metadata.csv")
    logger.info("Metadata: %d samples × %d columns.", *meta_raw.shape)

    aud = audit(meta_raw)
    aud.to_csv(raw_dir / "field_completeness_audit.csv", index=False)
    logger.info(
        "Raw GSM metadata: %d/%d characteristic fields have ≥90%% non-missing "
        "(this is a descriptor of GEO's GSM schema, NOT of hypothesis-critical fields).",
        int((aud["missing_pct"] <= 10).sum()), len(aud),
    )

    can, mapper_diag = map_canonical(meta_raw)
    logger.info("Canonical-column mapper diagnostic (raw → mapped OR NOT MAPPED):")
    for canon, diag in mapper_diag.items():
        short_diag = diag[:200] + ("…" if len(diag) > 200 else "")
        nn = int(can[canon].notna().sum()) if canon in can.columns else 0
        logger.info("  %s <- %s [non-null after mapping: %d/%d]",
                    canon, short_diag, nn, len(can))

    n_canon_ok = sum(1 for c in can.columns if can[c].notna().sum() >= 20)
    logger.info(
        "Hypothesis-critical canonical fields with ≥20 non-null: %d/%d.",
        n_canon_ok, len(can.columns),
    )
    can.to_csv(raw_dir / "metadata_canonical_unfiltered.csv")

    feas = feasibility(can)
    with open(raw_dir / "hypothesis_feasibility.json", "w") as fh:
        json.dump(feas, fh, indent=2, default=str)
    for h_name, info in feas.items():
        logger.info("  %s: testable=%s — %s", h_name, info["testable"], info["fields"])

    # A replication cohort needs at least one testable outcome AND a classifier
    # label column (nassiri_group or bi_group). Classifier labels will almost
    # always be missing because these datasets predate the classifiers — that's
    # OK and expected. The report should say so explicitly. We still write
    # the processed output with a clear provenance note stating the labels are
    # absent and need to be generated by running the classifier.
    proc_dir = PATHS.data_processed / f"replication_{accession.lower()}"
    proc_dir.mkdir(parents=True, exist_ok=True)
    can.to_csv(proc_dir / "metadata.csv")
    pd.DataFrame([{
        "dataset": accession,
        "n_samples": len(can),
        "h4_testable": feas.get("H4 replication (same-sign direction)", {}).get("testable", False),
        "h2_rep_testable": feas.get("H2 replication (subgroup enrichment)", {}).get("testable", False),
        "h3_rep_testable": feas.get("H3 replication (NF2 association)", {}).get("testable", False),
        "note": (
            "Subgroup labels (nassiri_group / bi_group) will almost always be "
            "absent in raw GEO metadata for replication cohorts published before "
            "2021-2023. They must be generated by running the Nassiri / Bi lab "
            "classifiers on the expression matrix. DO NOT invent labels to "
            "make H4 'testable'."
        ),
    }]).to_csv(proc_dir / "metadata_PROVENANCE.csv", index=False)
    logger.info("Processed metadata written to %s. expression.csv must be populated separately.", proc_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
