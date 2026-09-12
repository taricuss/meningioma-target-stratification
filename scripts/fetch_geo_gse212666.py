#!/usr/bin/env python
"""Fetch and audit GSE212666 (Choudhury/Raleigh/Bi lab Nature Genetics 2022 — 565 meningiomas).

AIM 0 REAL DATA ACQUISITION — NO SYNTHETIC FALLBACK.

Source:
  - GEO accession: GSE212666
  - Paper: "Meningioma DNA methylation groups identify biological drivers
    and therapeutic vulnerabilities" — Choudhury et al. / Raleigh & Bi labs,
    Nature Genetics 54(5):649–659, 2022 (PMID 35534562).
    Follow-up subgroup refinement: Choudhury et al., Neuro-Oncology 25(3):
    520–530, 2023 — Hypermitotic subdivided into FOXM1+/Proliferative vs
    FOXM1–/Hypermetabolic, concordant with Nassiri MG3/MG4.
  - Subgroup labels expected: 3-tier classification originally derived from
      Illumina EPIC DNA methylation (n=565 total; n=200 discovery + n=365
      validation). Subgroup labels were also propagated to a subset of
      samples with RNA-seq (n=200 discovery + n=302 validation), proteomics,
      and scRNA-seq through the integrated multi-omic model (see
      BI_CLASSIFIER_METHODS_NOTE in config.py).
      Merlin-intact / Immune-enriched / Hypermitotic
  - Platforms listed: Illumina EPIC methylation, bulk RNA-seq, proteomics, scRNA-seq
  - Open GEO accession.

This script:
  1. Downloads GEO Series SOFT via GEOparse.
  2. Extracts sample-level metadata from every GSM record (characteristics_ch1,
     source_name_ch1).
  3. Exhaustive field-completeness audit with per-column missingness rates.
  4. Maps to canonical columns (bi_group, nassiri_group, nf2_status, who_grade,
     recurrence). NO IMPUTATION — only exact-case-insensitive matches are mapped.
  5. Honest hypothesis-feasibility report for H1-H5, per field, with minimum
     cell-size checks.
  6. Writes to data/processed/discovery_bi/ ONLY if at least one hypothesis is
     testable given fields actually present.

Errors loudly; does not fabricate a single row if the fetch or metadata is
insufficient.
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
from meningeal_extension.config import PATHS  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-22s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger("fetch_GSE212666")

GSE_ACCESSION = "GSE212666"
EXPECTED_N = 565
COHORT_NAME = "discovery_bi"


@dataclass
class FieldAudit:
    field: str
    n_non_null: int
    n_total: int
    missing_pct: float
    n_unique: int
    examples: List[str]


# Minimum fields per hypothesis (same structure as Nassiri script — used in report).
HYPOTHESIS_MINIMUM_FIELDS = {
    "H1 (subgroup concordance)": ["bi_group", "nassiri_group"],
    "H2 (subgroup × target-program)": ["bi_group"],
    "H3 (NF2 × target-program)": ["nf2_status"],
    "H4 (replication)": ["external cohort, not this one"],
    "H5 (survival)": ["recurrence_event", "recurrence_months"],
}


def _import_geoparse():
    try:
        import GEOparse  # type: ignore
        return GEOparse
    except ImportError as e:
        logger.critical(
            "GEOparse not installed. Run: pip install -r requirements.txt"
        )
        raise SystemExit(2) from e


def fetch_gse_soft(accession: str, dest_dir: Path):
    GEOparse = _import_geoparse()
    dest_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Fetching %s via GEOparse → %s", accession, dest_dir)
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(dest_dir), silent=False)
    except Exception as e:
        logger.critical(
            "GSE%s fetch failed: %s\n"
            "Manual URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n"
            "Save *_family.soft.gz to %s and rerun with --skip-download.\n"
            "DO NOT substitute synthetic data — a failed fetch is a real finding.",
            accession, e, accession, dest_dir,
        )
        raise SystemExit(3) from e
    logger.info("GSE%s downloaded. GSM count in SOFT: %d", accession, len(gse.gsms))
    return gse


def extract_gsm_metadata(gse) -> pd.DataFrame:
    rows = []
    for gsm_id, gsm in gse.gsms.items():
        row = {"sample_id": gsm_id, "gsm_title": gsm.metadata.get("title", [""])[0]}
        for c in gsm.metadata.get("characteristics_ch1", []):
            if ":" in c:
                k, v = c.split(":", 1)
                row[k.strip()] = v.strip()
            else:
                row.setdefault("characteristics_unparsed", []).append(c)
        for i, s in enumerate(gsm.metadata.get("source_name_ch1", [])):
            row[f"source_name_{i}"] = s
        for i, s in enumerate(gsm.metadata.get("supplementary_file", [])):
            row[f"supp_file_{i}"] = s
        for k, v in gsm.metadata.items():
            if k not in {"characteristics_ch1", "source_name_ch1", "supplementary_file", "title"}:
                if isinstance(v, list) and len(v) == 1:
                    row[k] = v[0]
                elif isinstance(v, list):
                    row[k] = "|".join(str(x) for x in v)
                else:
                    row[k] = v
        rows.append(row)
    meta = pd.DataFrame(rows).set_index("sample_id")
    logger.info("Extracted %d cols × %d GSM samples.", meta.shape[1], meta.shape[0])
    return meta


def audit_field_completeness(meta: pd.DataFrame) -> pd.DataFrame:
    n_total = len(meta)
    audits = []
    for col in meta.columns:
        s = meta[col]
        non_null = s.notna() & (s.astype(str).str.strip() != "") & (s.astype(str).str.lower() != "nan")
        n_non = int(non_null.sum())
        examples = s[non_null].astype(str).unique().tolist()[:5]
        audits.append(FieldAudit(
            field=col,
            n_non_null=n_non,
            n_total=n_total,
            missing_pct=round(100.0 * (1 - n_non / n_total), 2) if n_total else 100.0,
            n_unique=int(s[non_null].nunique()),
            examples=examples,
        ))
    return pd.DataFrame([asdict(a) for a in audits]).sort_values("missing_pct").reset_index(drop=True)


BI_GROUP_CANONICAL = {"merlin-intact", "immune-enriched", "hypermitotic"}
# Case-insensitive + hyphen/underscore variants accepted for exact label matching:
BI_GROUP_LABEL_MAP = {
    "merlin-intact": "Merlin-intact",
    "merlin_intact": "Merlin-intact",
    "merlin intact": "Merlin-intact",
    "immune-enriched": "Immune-enriched",
    "immune_enriched": "Immune-enriched",
    "immune enriched": "Immune-enriched",
    "hypermitotic": "Hypermitotic",
    "hyper_mitotic": "Hypermitotic",
    "proliferative": "Hypermitotic",  # Sometimes labelled this in preprints
}


def map_to_canonical_columns(meta: pd.DataFrame) -> pd.DataFrame:
    """Map GEO metadata columns → canonical names. NO IMPUTATION.

    Exact case-insensitive / underscore-insensitive matches only.
    Bi-group labels are verified against the 3 published labels and anything
    unrecognized is set to NaN (never coerced).
    """
    col_lower = {c.lower().replace(" ", "_").replace("-", "_"): c for c in meta.columns}
    canonical = pd.DataFrame(index=meta.index)

    mappings = [
        ("bi_group", ["subgroup", "group", "classification", "epigenomic_subgroup",
                      "methylation_subgroup", "bi_subgroup", "cluster", "class",
                      "meningioma_subgroup", "dna_methylation_class", "molecular_subtype",
                      "molecular_group"]),
        ("nassiri_group", ["nassiri_subgroup", "nassiri_group", "nassiri_classification"]),
        ("who_grade", ["who_grade", "grade", "histologic_grade", "histological_grade",
                       "who_2016_grade", "who_2021_grade", "pathology", "tumor_grade",
                       "pathological_grade", "histopathology", "histology_grade"]),
        ("nf2_status", ["nf2", "nf2_status", "nf2_mutation", "merlin_status",
                        "nf2_mutation_status", "nf2_alteration", "merlin"]),
        ("recurrence_event", ["recurrence", "recurred", "progression", "progressed",
                              "rfs_event", "dfs_event", "outcome", "event"]),
        ("recurrence_months", ["rfs_months", "recurrence_months", "dfs_months",
                               "follow_up_months", "time_to_recurrence", "survival_months",
                               "months_to_recurrence", "time_to_event_months", "followup_months"]),
        ("age_at_surgery", ["age", "age_at_surgery", "diagnosis_age", "patient_age",
                            "age_at_diagnosis", "years", "age_years"]),
        ("sex", ["sex", "gender"]),
    ]

    applied = []
    skipped = []
    mapper_diag = {}
    for canon_name, candidates in mappings:
        found = None
        for cand in candidates:
            if cand in col_lower:
                found = col_lower[cand]
                break
        if found is not None:
            canonical[canon_name] = meta[found].values
            applied.append((canon_name, found))
            mapper_diag[canon_name] = f"mapped from raw column '{found}'"
        else:
            skipped.append(canon_name)
            mapper_diag[canon_name] = f"NOT MAPPED — no alias match in raw cols: {sorted(col_lower.keys())[:30]}{'…' if len(col_lower) > 30 else ''}"
    logger.info("Mapped columns: %s\nMissing canonical columns: %s", applied, skipped)
    logger.info("Canonical mapper diagnostic (with post-map non-null counts):")
    for cn, diag in mapper_diag.items():
        nn = int(canonical[cn].notna().sum()) if cn in canonical.columns else 0
        logger.info("  %s <- %s [non-null: %d/%d]",
                    cn, diag[:200] + ("…" if len(diag) > 200 else ""), nn, len(canonical))

    # --- Bi-group label sanitization (only on actually-present values):
    if "bi_group" in canonical.columns:
        raw_labels = canonical["bi_group"].astype(str).str.strip().str.lower().str.replace(" ", "_").str.replace("-", "_")
        mapped_labels = raw_labels.map(lambda x: BI_GROUP_LABEL_MAP.get(x, np.nan))
        unknowns = canonical.loc[mapped_labels.isna() & canonical["bi_group"].notna(), "bi_group"].astype(str).unique().tolist()
        if unknowns:
            logger.warning(
                "bi_group column has %d unrecognized label(s): %s. "
                "These are set to NaN — NOT coerced to a known group. "
                "If these are valid subgroup labels from the paper, add them to "
                "BI_GROUP_LABEL_MAP in this script AFTER verifying against the "
                "paper's supplementary tables.",
                len(unknowns), unknowns[:10],
            )
        canonical["bi_group"] = mapped_labels

    # WHO grade: pull roman numeral (I/II/III) from whatever string is there; fail to NaN otherwise.
    if "who_grade" in canonical.columns:
        extracted = (
            canonical["who_grade"]
            .astype(str)
            .str.upper()
            .str.extract(r"\b([I]{1,3})\b", expand=False)
        )
        canonical["who_grade"] = extracted

    # NF2: only accept unambiguous values → NaN otherwise.
    if "nf2_status" in canonical.columns:
        normalized = (
            canonical["nf2_status"]
            .astype(str)
            .str.strip()
            .str.lower()
            .replace({
                "mutant": "Mutant/Loss", "loss": "Mutant/Loss", "del": "Mutant/Loss",
                "deleted": "Mutant/Loss", "mutated": "Mutant/Loss",
                "abnormal": "Mutant/Loss", "altered": "Mutant/Loss",
                "loh": "Mutant/Loss",
                "wt": "Intact", "wild": "Intact", "wild_type": "Intact",
                "wild-type": "Intact", "intact": "Intact", "normal": "Intact",
                "no_alteration": "Intact",
            })
        )
        unknown_mask = ~normalized.isin(["Mutant/Loss", "Intact"])
        normalized[unknown_mask] = np.nan
        canonical["nf2_status"] = normalized

    # Recurrence event: strict yes/no booleanization.
    if "recurrence_event" in canonical.columns:
        ce = canonical["recurrence_event"].astype(str).str.strip().str.lower()
        canonical["recurrence_event"] = ce.map({
            "yes": True, "true": True, "1": True, "recurred": True,
            "progressed": True, "event": True, "recurrence": True,
            "no": False, "false": False, "0": False, "none": False,
        })

    for numeric_col in ["age_at_surgery", "recurrence_months"]:
        if numeric_col in canonical.columns:
            canonical[numeric_col] = pd.to_numeric(canonical[numeric_col], errors="coerce")

    return canonical


def assess_hypothesis_feasibility(canonical: pd.DataFrame) -> Dict[str, Dict]:
    n = len(canonical)
    out: Dict[str, Dict] = {}
    for h_name, required_fields in HYPOTHESIS_MINIMUM_FIELDS.items():
        if h_name == "H4 (replication)":
            out[h_name] = {
                "testable": False,
                "reason": "Requires separate replication cohort. Not assessed here.",
            }
            continue
        field_status = {}
        all_ok = True
        for f in required_fields:
            if f not in canonical.columns:
                field_status[f] = "NOT PRESENT"
                all_ok = False
                continue
            non_null = int(canonical[f].notna().sum())
            pct = 100 * non_null / n if n else 0
            field_status[f] = f"{non_null}/{n} ({pct:.1f}% non-null)"
            if non_null < 20:
                all_ok = False
        out[h_name] = {
            "testable": all_ok,
            "n_total": n,
            "fields": field_status,
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true",
                        help="Use existing soft file in data/raw/GSE212666/ (no re-fetch).")
    args = parser.parse_args()

    raw_dir = PATHS.data_raw / "GSE212666"
    raw_dir.mkdir(parents=True, exist_ok=True)

    if args.skip_download and any(raw_dir.glob(f"{GSE_ACCESSION}_family.soft*")):
        logger.info("--skip-download; reusing existing SOFT in %s", raw_dir)
        GEOparse = _import_geoparse()
        gse = GEOparse.get_GEO(filepath=str(next(raw_dir.glob(f"{GSE_ACCESSION}_family.soft*"))), silent=False)
    else:
        gse = fetch_gse_soft(GSE_ACCESSION, raw_dir)

    meta_raw = extract_gsm_metadata(gse)
    meta_raw.to_csv(raw_dir / "gsm_all_metadata.csv")
    logger.info("Raw GSM metadata: %s", str(meta_raw.shape))

    audit = audit_field_completeness(meta_raw)
    audit.to_csv(raw_dir / "field_completeness_audit.csv", index=False)
    high_coverage = int((audit["missing_pct"] <= 10).sum())
    logger.info(
        "Raw GSM metadata: %d/%d total fields; %d/%d with ≥90%% non-missing "
        "(GEO schema completeness only — hypothesis-critical field completeness reported separately after mapping).",
        len(audit), len(audit), high_coverage, len(audit),
    )

    if len(meta_raw) < EXPECTED_N:
        logger.warning(
            "GSE212666 expected N=%d (paper) but SOFT has N=%d GSM entries. "
            "This may reflect GEO hosting a subset; check if multi-platform runs "
            "have separate GSM pools. DO NOT pad the shortfall with synthetic rows.",
            EXPECTED_N, len(meta_raw),
        )

    canonical = map_to_canonical_columns(meta_raw)
    canonical.to_csv(raw_dir / "metadata_canonical_unfiltered.csv")

    # Report bi_group label distribution IF PRESENT
    if "bi_group" in canonical.columns and canonical["bi_group"].notna().any():
        dist = canonical["bi_group"].value_counts(dropna=False)
        logger.info("Bi-group label distribution (real values only; NaN=unmapped):\n%s", dist.to_string())

    feasibility = assess_hypothesis_feasibility(canonical)
    with open(raw_dir / "hypothesis_feasibility.json", "w") as fh:
        json.dump(feasibility, fh, indent=2, default=str)
    for h_name, info in feasibility.items():
        logger.info("  %s: testable=%s — %s",
                    h_name, info.get("testable"),
                    info.get("reason") or str(info.get("fields")))

    any_testable = any(v.get("testable") for v in feasibility.values())
    if not any_testable:
        logger.critical(
            "=" * 78 + "\n"
            "GSE212666: NO HYPOTHESIS MEETS MINIMUM FIELD REQUIREMENTS.\n"
            "Nothing written to data/processed/discovery_bi/.\n"
            "This is a real, reportable finding: the open GEO metadata alone\n"
            "is insufficient. Next step: check paper Supplementary Tables\n"
            "for subgroup labels / grade / NF2 and add them as a sidecar CSV,\n"
            "then re-merge. DO NOT fabricate labels.\n"
            + "=" * 78
        )
        return 0

    out_dir = PATHS.data_processed / COHORT_NAME
    out_dir.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(out_dir / "metadata.csv")

    provenance = pd.DataFrame([{
        "dataset": GSE_ACCESSION,
        "n_metadata_rows": len(canonical),
        "n_expected_paper": EXPECTED_N,
        "metadata_source": "GEO SOFT GSM characteristics_ch1 + source_name_ch1",
        "expression_csv_generated_here": False,
        "note": (
            "expression.csv is NOT yet populated. If GSE212666 has a GEO Series "
            "Matrix file with RNA-seq counts / normalized expression, it should "
            "be parsed separately. DO NOT fabricate an expression matrix."
        ),
        "bi_group_testable": feasibility.get("H2 (subgroup × target-program)", {}).get("testable", False),
        "nf2_testable": feasibility.get("H3 (NF2 × target-program)", {}).get("testable", False),
        "h1_testable": feasibility.get("H1 (subgroup concordance)", {}).get("testable", False),
    }])
    provenance.to_csv(out_dir / "metadata_PROVENANCE.csv", index=False)
    logger.info("Processed metadata written to %s. Run an expression-matrix fetcher next.", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
