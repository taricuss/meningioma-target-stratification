#!/usr/bin/env python
"""Fetch and audit GSE180061 (Nassiri et al. 2021 Nature — 185 meningiomas).

AIM 0 REAL DATA ACQUISITION — NO SYNTHETIC FALLBACK.

Source:
  - GEO accession: GSE180061 (Illumina Infinium HumanMethylation450, idat files + SOFT)
  - Controlled multi-omics: EGAS00001004982 (WES, bulk mRNA, snRNA-seq)
  - Open processed mRNA / clinical / mutations via cBioPortal: mng_utoronto_2021
    (see companion script: fetch_cbio_mng_utoronto_2021.py)

This script:
  1. Downloads the GEO Series SOFT record via GEOparse (metadata, NOT raw idats —
     those are ~20GB and require separate processing with minfi).
  2. Extracts every sample-level metadata field available in the GSM records.
  3. Writes an HONEST completeness report: per-field non-missing counts,
     missingness rates, and which hypotheses (H1-H4) are actually testable given
     the fields present.
  4. Saves extracted metadata to data/raw/GSE180061/ and a cleaned version
     to data/processed/discovery_nassiri/ ONLY if the minimum required fields
     for at least one hypothesis exist.

If the GEO fetch fails, this script ERRORS OUT with exit code != 0 and
a message telling the user exactly what to do next. It does NOT invent data.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from meningeal_extension.config import PATHS  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-22s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger("fetch_GSE180061")

GSE_ACCESSION = "GSE180061"
EXPECTED_N = 185
COHORT_NAME = "discovery_nassiri"


@dataclass
class FieldAudit:
    field: str
    n_non_null: int
    n_total: int
    missing_pct: float
    n_unique: int
    examples: List[str]
    pass_minimum_for_any_hypothesis: bool


HYPOTHESIS_MINIMUM_FIELDS = {
    "H1 (subgroup concordance)": ["nassiri_group", "bi_group"],
    "H2 (subgroup × target-program)": ["nassiri_group"],
    "H3 (NF2 × target-program)": ["nf2_status"],
    "H4 (replication, external cohort)": ["requires separate cohort, not this one"],
    "H5 (survival)": ["recurrence_event", "recurrence_months"],
}


def _import_geoparse():
    try:
        import GEOparse  # type: ignore
        return GEOparse
    except ImportError as e:
        logger.critical(
            "GEOparse not installed. Run: pip install GEOparse>=2.0.3\n"
            "(Already listed in requirements.txt — `pip install -r requirements.txt` should work.)"
        )
        raise SystemExit(2) from e


def fetch_gse_soft(accession: str, dest_dir: Path) -> object:
    """Download GSE SOFT file to dest_dir. Raise on failure; no silent fallback."""
    GEOparse = _import_geoparse()
    dest_dir.mkdir(parents=True, exist_ok=True)
    soft_path = dest_dir / f"{accession}_family.soft.gz"
    logger.info("Fetching %s via GEOparse → %s", accession, soft_path)
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(dest_dir), silent=False)
    except Exception as e:
        logger.critical(
            "Failed to fetch GSE%s from NCBI GEO.\n"
            "Error: %s\n"
            "Next steps:\n"
            "  1. Check network connectivity to ncbi.nlm.nih.gov\n"
            "  2. Manually download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n"
            "     and save the *_family.soft.gz file to %s\n"
            "  3. Retry this script.\n"
            "DO NOT substitute synthetic data — a failed fetch is a real result.",
            accession, e, accession, dest_dir,
        )
        raise SystemExit(3) from e
    logger.info("GSE%s fetched. %d GSM samples in SOFT.", accession, len(gse.gsms))
    return gse


def extract_gsm_metadata(gse) -> pd.DataFrame:
    """Flatten every characteristics_ch1 + source_name_ch1 + title field from GSM records."""
    rows: List[Dict] = []
    for gsm_id, gsm in gse.gsms.items():
        row: Dict = {"sample_id": gsm_id, "gsm_title": gsm.metadata.get("title", [""])[0]}
        # characteristics_ch1 is the primary structured metadata list
        chars = gsm.metadata.get("characteristics_ch1", [])
        for c in chars:
            if ":" in c:
                k, v = c.split(":", 1)
                row[k.strip()] = v.strip()
            else:
                row.setdefault("characteristics_unparsed", []).append(c)
        # source_name_ch1 usually contains subgroup or tissue label
        sources = gsm.metadata.get("source_name_ch1", [])
        for i, s in enumerate(sources):
            row[f"source_name_{i}"] = s
        # supplemental_file — we record this so we know if idats are linked
        supp = gsm.metadata.get("supplementary_file", [])
        for i, s in enumerate(supp):
            row[f"supp_file_{i}"] = s
        rows.append(row)
    meta = pd.DataFrame(rows).set_index("sample_id")
    logger.info("Extracted %d metadata columns across %d GSM samples.", meta.shape[1] - 0, meta.shape[0])
    return meta


def audit_field_completeness(meta: pd.DataFrame) -> pd.DataFrame:
    n_total = len(meta)
    audits: List[FieldAudit] = []
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
            pass_minimum_for_any_hypothesis=False,
        ))
    return pd.DataFrame([asdict(a) for a in audits]).sort_values("missing_pct").reset_index(drop=True)


def map_to_canonical_columns(meta: pd.DataFrame) -> pd.DataFrame:
    """Try to map GEO metadata column names to our canonical column names.

    CRITICAL: This function never IMPUTES values. It only RENAMES columns when
    a near-exact case-insensitive match exists. Unknown values remain NaN
    with a full audit trail of what was NOT found.
    """
    # Lowercase + underscore lookup of all existing columns
    col_lower = {c.lower().replace(" ", "_").replace("-", "_"): c for c in meta.columns}
    canonical = pd.DataFrame(index=meta.index)

    mappings = [
        # (canonical_name, list_of_possible_lowercase_source_names)
        ("nassiri_group", ["dna_methylation_subgroup", "methylation_subgroup",
                           "subgroup", "nassiri_subgroup", "group", "classification",
                           "cluster", "class", "molecular_subtype", "molecular_group",
                           "molecule_subtype"]),
        ("bi_group", ["bi_group", "bi_subgroup", "bi_classification", "merlin_group",
                      "immune_group", "subgroup_bi", "epigenomic_subgroup"]),
        ("who_grade", ["who_grade", "grade", "who_2016_grade", "histological_grade",
                       "pathology", "tumor_grade", "pathological_grade", "histopathology",
                       "who_2021_grade", "histology_grade", "histologic_grade"]),
        ("nf2_status", ["nf2_status", "nf2_mutation", "nf2", "merlin_status",
                        "nf2_mutation_status", "nf2_alteration", "merlin"]),
        ("recurrence_event", ["recurrence", "recurred", "tumor_recurrence",
                              "progression", "progressed", "rfs_event",
                              "dfs_event", "outcome", "event"]),
        ("recurrence_months", ["recurrence_free_months", "rfs_months", "follow_up_months",
                               "time_to_recurrence_months", "survival_months",
                               "dfs_months", "months_to_recurrence",
                               "time_to_event_months", "followup_months"]),
        ("age_at_surgery", ["age", "age_at_surgery", "patient_age", "diagnosis_age",
                            "age_at_diagnosis", "years", "age_years"]),
        ("sex", ["sex", "gender"]),
        ("sample_name", ["sample_name", "specimen_id", "patient_id"]),
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
            mapper_diag[canon_name] = f"NOT MAPPED — no alias match among raw cols: {sorted(col_lower.keys())[:25]}{'…' if len(col_lower) > 25 else ''}"

    logger.info(
        "Column mapping applied: %s\nColumns NOT found in GEO metadata (will be NaN): %s",
        applied, skipped,
    )
    logger.info("Canonical mapper diagnostic (per field — with post-map non-null count):")
    for cn, diag in mapper_diag.items():
        nn = int(canonical[cn].notna().sum()) if cn in canonical.columns else 0
        logger.info("  %s <- %s [non-null: %d/%d]",
                    cn, diag[:200] + ("…" if len(diag) > 200 else ""), nn, len(canonical))

    # Type coercion — only on columns that actually exist, with errors → NaN
    if "who_grade" in canonical.columns:
        canonical["who_grade"] = (
            canonical["who_grade"]
            .astype(str)
            .str.strip()
            .str.upper()
            .str.extract(r"([I]+)", expand=False)  # Pull roman numeral I/II/III
        )
    if "nf2_status" in canonical.columns:
        canonical["nf2_status"] = (
            canonical["nf2_status"]
            .astype(str)
            .str.strip()
            .str.lower()
            .replace(
                {
                    "mutant": "Mutant/Loss", "loss": "Mutant/Loss",
                    "deleted": "Mutant/Loss", "mutated": "Mutant/Loss",
                    "abnormal": "Mutant/Loss", "altered": "Mutant/Loss",
                    "wt": "Intact", "wild": "Intact", "wild-type": "Intact",
                    "intact": "Intact", "normal": "Intact",
                }
            )
        )
        # Unmapped values → NaN, NOT silently coerced to Intact
        mask_unknown = ~canonical["nf2_status"].isin(["Mutant/Loss", "Intact"])
        canonical.loc[mask_unknown, "nf2_status"] = np.nan
    if "recurrence_event" in canonical.columns:
        ce = canonical["recurrence_event"].astype(str).str.strip().str.lower()
        canonical["recurrence_event"] = ce.map(
            {
                "yes": True, "true": True, "1": True, "recurred": True,
                "progressed": True, "event": True, "recurrence": True,
                "no": False, "false": False, "0": False, "none": False,
            }
        )
    if "age_at_surgery" in canonical.columns:
        canonical["age_at_surgery"] = pd.to_numeric(canonical["age_at_surgery"], errors="coerce")
    if "recurrence_months" in canonical.columns:
        canonical["recurrence_months"] = pd.to_numeric(canonical["recurrence_months"], errors="coerce")

    return canonical


def assess_hypothesis_feasibility(canonical: pd.DataFrame) -> Dict[str, Dict]:
    n = len(canonical)
    results: Dict[str, Dict] = {}
    for h_name, required_fields in HYPOTHESIS_MINIMUM_FIELDS.items():
        if h_name == "H4 (replication, external cohort)":
            results[h_name] = {
                "testable": False,
                "reason": "Requires separate replication cohort (see fetch_geo_replication.py), not GSE180061.",
            }
            continue
        field_status = {}
        all_present = True
        for f in required_fields:
            if f not in canonical.columns:
                field_status[f] = "NOT PRESENT in any GSM metadata"
                all_present = False
                continue
            non_null = canonical[f].notna().sum()
            pct = 100 * non_null / n if n else 0
            field_status[f] = f"{non_null}/{n} non-null ({pct:.1f}%)"
            if non_null < 20:  # minimum cell size for any sensible test
                all_present = False
        results[h_name] = {
            "testable": all_present,
            "n_samples_total": n,
            "fields": field_status,
        }
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true",
                        help="Use already-downloaded SOFT file in data/raw/GSE180061/ without re-fetching.")
    args = parser.parse_args()

    raw_dir = PATHS.data_raw / "GSE180061"
    raw_dir.mkdir(parents=True, exist_ok=True)

    # 1) Fetch
    if args.skip_download and any(raw_dir.glob(f"{GSE_ACCESSION}_family.soft*")):
        logger.info("--skip-download set; loading existing SOFT from %s", raw_dir)
        GEOparse = _import_geoparse()
        soft_file = next(raw_dir.glob(f"{GSE_ACCESSION}_family.soft*"))
        gse = GEOparse.get_GEO(filepath=str(soft_file), silent=False)
    else:
        gse = fetch_gse_soft(GSE_ACCESSION, raw_dir)

    # 2) Extract full GSM metadata
    meta_raw = extract_gsm_metadata(gse)
    meta_raw.to_csv(raw_dir / "gsm_all_metadata.csv")
    logger.info("Raw GSM metadata saved to %s (shape %s)", raw_dir / "gsm_all_metadata.csv", meta_raw.shape)

    # 3) Field-level audit
    audit = audit_field_completeness(meta_raw)
    audit.to_csv(raw_dir / "field_completeness_audit.csv", index=False)
    logger.info(
        "Raw GSM metadata audit: %d/%d total characteristic fields; %d/%d have ≥90%% non-missing "
        "(GSM schema completeness ≠ hypothesis-critical field completeness.",
        len(audit), len(audit),
        int((audit["missing_pct"] <= 10).sum()), len(audit),
    )
    if len(meta_raw) < EXPECTED_N:
        logger.warning(
            "GSE180061 expected N=%d but found N=%d in SOFT. This may or may not be a "
            "problem — the paper reports N=185; GEO may include technical replicates "
            "or a subset. Check audit above. DO NOT top up with synthetic samples.",
            EXPECTED_N, len(meta_raw),
        )

    # 4) Map to canonical columns — NO IMPUTATION
    canonical = map_to_canonical_columns(meta_raw)
    canonical.to_csv(raw_dir / "metadata_canonical_unfiltered.csv")

    # 5) Hypothesis feasibility
    feasibility = assess_hypothesis_feasibility(canonical)
    feasibility_path = raw_dir / "hypothesis_feasibility.json"
    with open(feasibility_path, "w") as fh:
        json.dump(feasibility, fh, indent=2, default=str)
    logger.info("Hypothesis feasibility report written to %s", feasibility_path)
    for h_name, info in feasibility.items():
        logger.info("  %s: testable=%s — %s", h_name, info.get("testable"), info.get("reason") or list(info.get("fields", {}).values()))

    # 6) Only write to data/processed if ≥ 1 hypothesis is actually testable.
    testable_any = any(v.get("testable") for v in feasibility.values())
    if not testable_any:
        logger.critical(
            "=" * 78 + "\n"
            "GSE180061 METADATA INSUFFICIENT FOR ANY PRE-SPECIFIED HYPOTHESIS.\n"
            "Canonical columns saved to raw_dir only. NOTHING written to\n"
            "data/processed/discovery_nassiri/. This is a real result:\n"
            "  \"GSE180061 open-access GEO metadata alone cannot support H1-H5.\"\n"
            "  Next step: run scripts/fetch_cbio_mng_utoronto_2021.py to pull\n"
            "  processed clinical + mutation + mRNA data from the cBioPortal entry\n"
            "  mng_utoronto_2021 (the Nassiri paper's own processed multi-omics).\n"
            "  DO NOT invent subgroup/NF2 labels to make hypotheses testable.\n"
            + "=" * 78
        )
        # Still exit 0 because the audit itself succeeded — but the processed dir stays empty.
        return 0

    # Need at least subgroup labels to populate nassiri_group (the classifier-output column
    # that Bi lab groups will be compared against for H1). If only methylation subgroup is
    # available (which it should be for Nassiri), that IS the nassiri_group — no relabeling.
    if "nassiri_group" not in canonical.columns or canonical["nassiri_group"].notna().sum() < 20:
        logger.warning(
            "nassiri_group not populated from GEO metadata. Subgroup labels for the "
            "Nassiri classification may only be in the paper Supplementary Tables or "
            "cBioPortal. H1 cannot be run until nassiri_group exists. NOT imputing."
        )

    # H1 requires BI-LAB GROUP labels on THESE samples, which we do NOT have from
    # GSE180061 alone. The Bi-lab classifier must be run on the expression data to
    # produce bi_group. That's a separate step — NOT this script's job, and NEVER
    # an excuse to plug in synthetic labels.
    if "bi_group" not in canonical.columns:
        logger.info(
            "bi_group not in GSE180061 metadata. Expected: the Bi-lab classifier must "
            "be applied to the cBioPortal mRNA matrix to generate bi_group for H1. "
            "H2 (subgroup program enrichment using Nassiri's own nassiri_group) is still "
            "testable if nassiri_group is present."
        )

    out_dir = PATHS.data_processed / COHORT_NAME
    out_dir.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(out_dir / "metadata.csv")

    # Expression: this script does NOT process the 20GB idat methylation files.
    # It records the fact, and tells the user to pull mRNA from cBioPortal.
    provenance_note = pd.DataFrame([{
        "dataset": "GSE180061",
        "n_metadata_rows": len(canonical),
        "n_expected": EXPECTED_N,
        "metadata_source": "GEO SOFT GSM characteristics_ch1 + source_name_ch1",
        "methylation_idat_available": "YES in GEO (use minfi/ChAMP to process idats)",
        "mrna_processed_available": "NO on GEO — run fetch_cbio_mng_utoronto_2021.py for z-scores",
        "wes_snRNA_available": "Controlled access EGAS00001004982",
        "expression_csv_generated_here": False,
        "note": (
            "expression.csv in this directory will be populated by "
            "fetch_cbio_mng_utoronto_2021.py. DO NOT fabricate a matrix."
        ),
    }])
    provenance_note.to_csv(out_dir / "metadata_PROVENANCE.csv", index=False)
    logger.info(
        "Metadata written to %s. Expression matrix NOT yet populated — "
        "run fetch_cbio_mng_utoronto_2021.py next.", out_dir,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
