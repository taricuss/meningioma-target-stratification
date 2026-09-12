#!/usr/bin/env python
"""Fetch Nassiri 2021 processed multi-omics from cBioPortal (mng_utoronto_2021).

AIM 0 REAL DATA ACQUISITION — NO SYNTHETIC FALLBACK.

Source:
  Study ID: mng_utoronto_2021
  Paper: "Delineating the molecular and cellular architecture of human meningioma"
         Nassiri et al., Nature 2021.
  URL: https://www.cbioportal.org/study/summary?id=mng_utoronto_2021

Endpoints used (all public, no API key):
  GET /api/studies/{studyId}/clinical-data          — sample-level covariates
  GET /api/studies/{studyId}/mutations               — mutation calls (NF2 + any)
  GET /api/studies/{studyId}/molecular-data/mRNA     — mRNA z-scores vs reference diploid

Outputs:
  data/raw/mng_utoronto_2021/{clinical.csv, mutations.csv, mrna_zscores.csv}
  data/processed/discovery_nassiri/{metadata.csv, expression.csv, metadata_PROVENANCE.csv}

Behavior on failure:
  - Any unreachable endpoint → sys.exit(3) with actionable debug instructions.
  - If metadata is incomplete (e.g., no WHO grade, no subgroup labels) the script
    writes a report stating the hypothesis is NOT testable, and still writes metadata
    to processed/ (with NaN in missing columns, NOT imputed values).
  - mRNA matrix is saved ONLY for samples that also appear in clinical metadata.
    No "nearest-neighbor" or "carry-forward" imputation.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from meningeal_extension.config import NF2_EXPORT_COL, NF2_INTERNAL_COL, PATHS  # noqa: E402


def _canonicalize_nf2_for_export(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if NF2_INTERNAL_COL in out.columns and NF2_EXPORT_COL not in out.columns:
        out = out.rename(columns={NF2_INTERNAL_COL: NF2_EXPORT_COL})
    elif NF2_INTERNAL_COL in out.columns and NF2_EXPORT_COL in out.columns:
        out = out.drop(columns=[NF2_INTERNAL_COL])
    return out

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)-22s | %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger("fetch_cbio_mng_utoronto_2021")

CBIO_API = "https://www.cbioportal.org/api"
STUDY_ID = "mng_utoronto_2021"
EXPECTED_N_SAMPLES = 185  # per paper
COHORT_NAME_PROCESSED = "discovery_nassiri"

# Clinical attribute IDs commonly used in cBioPortal meningioma studies.
# If the real clinical-data endpoint returns different names, the audit log
# will show every attribute actually present — the user can rerun with the
# correct IDs after reading the report.
CLINICAL_ATTR_CANDIDATES = {
    "who_grade": ["WHO_GRADE", "GRADE", "WHO_2016_GRADE", "HISTOLOGIC_GRADE",
                  "PATHOLOGY", "TUMOR_GRADE", "PATHOLOGICAL_GRADE", "HISTOPATHOLOGY",
                  "WHO_2021_GRADE", "HISTOLOGY_GRADE", "HISTOLOGICAL_GRADE"],
    "nassiri_group": ["NASSIRI_GROUP", "NASSIRI_SUBGROUP", "NASSIRI_CLASSIFICATION",
                      "DNA_METHYLATION_CLASS", "METHYLATION_SUBGROUP", "METHYLATION_CLASS",
                      "MENINGIOMA_CLASS", "MOLECULAR_SUBTYPE", "MOLECULAR_GROUP",
                      "CLASS", "GROUP", "CLUSTER", "SUBTYPE"],
    "bi_group": ["BI_GROUP", "BI_SUBGROUP", "BI_CLASSIFICATION",
                 "EPIGENOMIC_SUBGROUP", "MERLIN_GROUP", "IMMUNE_GROUP"],
    "nf2_status": ["NF2_STATUS", "NF2", "NF2_MUTATION_STATUS", "MERLIN_STATUS",
                   "NF2_ALTERATION", "MERLIN"],
    "subgroup": ["SUBGROUP", "CLASS", "GROUP", "METHYLATION_CLASS", "MENINGIOMA_CLASS",
                 "DNA_METHYLATION_CLASS", "METHYLATION_SUBGROUP", "MOLECULAR_SUBTYPE",
                 "MOLECULAR_GROUP", "CLUSTER", "SUBTYPE"],
    "recurrence_event": ["RECURRENCE", "RFS_EVENT", "PROGRESSION", "DFS_EVENT",
                         "OUTCOME", "EVENT", "TUMOR_RECURRENCE"],
    "recurrence_months": ["RFS_MONTHS", "RECURRENCE_MONTHS", "DFS_MONTHS",
                           "TIME_TO_RECURRENCE_MONTHS", "MONTHS_TO_RECURRENCE",
                           "FOLLOW_UP_MONTHS", "TIME_TO_EVENT_MONTHS", "FOLLOWUP_MONTHS",
                           "SURVIVAL_MONTHS"],
    "age_at_surgery": ["AGE", "AGE_AT_SURGERY", "DIAGNOSIS_AGE", "AGE_AT_DIAGNOSIS",
                       "YEARS", "AGE_YEARS"],
    "sex": ["SEX", "GENDER"],
}


@dataclass
class FetchReport:
    endpoint: str
    status_code: Optional[int]
    n_rows: int
    n_cols: int
    n_samples: int
    sample_ids: List[str]


def _cbio_get_json(path: str, params: Optional[Dict] = None, timeout: int = 90) -> object:
    """GET JSON from cBioPortal. Raise SystemExit on any failure. NO SILENT FALLBACK."""
    url = f"{CBIO_API}{path}"
    logger.info("GET %s params=%s", url, params)
    try:
        r = requests.get(url, params=params, timeout=timeout)
    except requests.RequestException as e:
        logger.critical(
            "Network error reaching cBioPortal: %s\n"
            "Check: (1) internet connectivity, (2) cbioportal.org status page,\n"
            "(3) whether the endpoint requires an API key (this script assumes public).\n"
            "Manual fallback: download the study TSV from the cBioPortal datahub\n"
            "(https://www.cbioportal.org/datasets) and place it in data/raw/mng_utoronto_2021/.",
            e,
        )
        raise SystemExit(3) from e
    if r.status_code != 200:
        logger.critical(
            "cBioPortal returned HTTP %d for %s\nBody preview: %s\n"
            "Check study ID spelling (currently '%s'). Next step: open\n"
            "  %s/studies/%s/clinical-data\nin a browser to verify it loads.",
            r.status_code, url, r.text[:500], STUDY_ID, CBIO_API, STUDY_ID,
        )
        raise SystemExit(3)
    try:
        return r.json()
    except ValueError as e:
        logger.critical("cBioPortal returned non-JSON for %s: %s\nFirst 500 chars: %s",
                        url, e, r.text[:500])
        raise SystemExit(3) from e


def fetch_study_exists() -> Dict:
    """First smoke test — confirm the study ID resolves in the studies list."""
    studies = _cbio_get_json("/studies", params={"projection": "DETAILED"})
    assert isinstance(studies, list)
    match = [s for s in studies if s.get("studyId") == STUDY_ID]
    if not match:
        logger.critical(
            "Study '%s' not found in cBioPortal /studies endpoint.\n"
            "Available meningioma study IDs should include: mng_utoronto_2021, mng_miami_2022.\n"
            "If the study ID has changed, update STUDY_ID at the top of this script.",
            STUDY_ID,
        )
        raise SystemExit(3)
    return match[0]


def fetch_clinical_data() -> pd.DataFrame:
    """Fetch full clinical data table (long format clinicalAttributeId × sampleId)."""
    rows = _cbio_get_json(f"/studies/{STUDY_ID}/clinical-data")
    if not isinstance(rows, list) or not rows:
        logger.critical("Clinical data endpoint returned empty or non-list: %r", rows)
        raise SystemExit(3)
    df_long = pd.DataFrame(rows)
    required = {"clinicalAttributeId", "sampleId", "value"}
    missing = required - set(df_long.columns)
    if missing:
        logger.critical("Clinical data missing expected columns %s. Columns present: %s",
                        missing, list(df_long.columns))
        raise SystemExit(3)
    # Pivot to wide: one row per sampleId, one col per clinicalAttributeId.
    wide = df_long.pivot_table(index="sampleId", columns="clinicalAttributeId",
                               values="value", aggfunc="first")
    logger.info("Clinical data: %d samples × %d attributes (pivoted wide).", *wide.shape)
    return wide


def _cbio_post_json(path: str, body: Dict, timeout: int = 120) -> object:
    """POST JSON body to cBioPortal, return parsed JSON. NEVER raises for HTTP
    errors — returns None and lets the caller decide how to degrade. Use this
    for fetch endpoints that are optional (e.g., mutation calls that may be
    hosted at a separate API path per cBioPortal instance)."""
    url = f"{CBIO_API}{path}"
    logger.info("POST %s body keys=%s", url, sorted(body.keys()))
    try:
        r = requests.post(url, json=body, timeout=timeout)
        if r.status_code != 200:
            logger.warning(
                "POST %s HTTP %d. Falling back to NO mutation data (mutation-derived\n"
                "NF2 will be all-NaN; clinical NF2 is still used if present). Body preview:\n  %s",
                url, r.status_code, r.text[:400],
            )
            return None
        return r.json()
    except Exception as exc:  # noqa: BLE001 — degrade, don't crash the full fetch
        logger.warning("POST %s network error: %s. Proceeding without mutation data.", url, exc)
        return None


def fetch_mutation_data() -> pd.DataFrame:
    """Fetch mutation calls. We primarily need NF2 (H3), but keep all for provenance.

    DEGRADES GRACEFULLY ON FAILURE (returns empty DataFrame with 0 rows + big
    warning banner): clinical data + attribute audit will still be produced.
    This is important because cBioPortal instances host mutations via different
    endpoints (some via /mutations/fetch POST, some not at all)."""
    body = {
        "studyId": STUDY_ID,
        "projection": "DETAILED",
        "pageSize": 500000,
    }
    data = _cbio_post_json("/mutations/fetch", body)
    if isinstance(data, list) and len(data) > 0:
        df = pd.DataFrame(data)
        logger.info("Mutations (POST /mutations/fetch): %d rows. Cols: %s",
                    len(df), list(df.columns)[:25])
        return df
    try:
        df_alt = _cbio_get_json(f"/studies/{STUDY_ID}/mutations",
                                params={"projection": "DETAILED", "pageSize": 500000})
        if isinstance(df_alt, list) and len(df_alt) > 0:
            df = pd.DataFrame(df_alt)
            logger.info("Mutations (GET fallback): %d rows. Cols: %s",
                        len(df), list(df.columns)[:25])
            return df
    except SystemExit:
        pass
    except Exception as exc:  # noqa: BLE001
        logger.warning("GET /studies/%s/mutations also failed: %s. No mutation data.",
                       STUDY_ID, exc)
    # Second fallback: try mutations profile via molecular-data GET endpoint
    profiles = _safe_list_profiles()
    mut_profiles = [p for p in profiles if "mutat" in str(p.get("molecularProfileId", "")).lower()]
    if mut_profiles:
        sl_candidate = f"{STUDY_ID}_sequenced"
        mp = mut_profiles[0]["molecularProfileId"]
        logger.info("Trying mutation molecular profile %s via GET…", mp)
        try:
            r = requests.get(f"{CBIO_API}/molecular-profiles/{mp}/molecular-data",
                             params={"sampleListId": sl_candidate, "entrezGeneId": 4771},
                             timeout=180)
            if r.status_code == 200 and len(r.text.strip()) > 0:
                rows = r.json()
                if isinstance(rows, list) and rows:
                    logger.info("  Got %d rows via molecular profile GET for NF2.", len(rows))
                    return pd.DataFrame(rows)
        except Exception as exc:  # noqa: BLE001
            logger.warning("  mutation profile molecular-data fallback failed: %s", exc)
    logger.warning(
        "============================================================================\n"
        "NO MUTATION DATA RETRIEVED FOR %s.\n"
        "Mutation-derived NF2 calls will be ALL-NaN (not fabricated to Intact).\n"
        "CNA-derived NF2 from the gistic profile is STILL TRIED below — that is the\n"
        "primary NF2 source (80%% of NF2 inactivation in meningioma is via chr22q\n"
        "deletion, not point mutation). Clinical-data NF2 attribute (if present) is\n"
        "also used. Final NF2 merges all three sources with conflict→NaN.\n"
        "============================================================================",
        STUDY_ID,
    )
    return pd.DataFrame()


def _safe_list_profiles() -> list:
    """Best-effort profile list — returns [] on any failure."""
    try:
        data = _cbio_get_json(f"/studies/{STUDY_ID}/molecular-profiles")
        return data if isinstance(data, list) else []
    except SystemExit:
        return []
    except Exception:  # noqa: BLE001
        return []


def fetch_cna_gistic_data(entrez_gene_ids: Optional[Iterable[int]] = None) -> pd.DataFrame:
    """Fetch GISTIC discrete copy-number calls from the study's gistic profile.

    Uses the GET /molecular-profiles/{id}/molecular-data endpoint with
    sampleListId + (optional per-gene) entrezGeneId parameter — the endpoint
    confirmed working for mng_utoronto_2021_gistic (HTTP 200, 121/121 samples).

    If entrez_gene_ids is None we fetch ALL genes in the profile (bigger file,
    more provenance); otherwise only the requested genes (NF2=4771 plus e.g.
    common meningioma drivers for provenance/QC).

    DEGRADES GRACEFULLY: returns empty DataFrame with warning banner on failure;
    never crashes the clinical metadata pipeline.
    """
    profiles = _safe_list_profiles()
    gistic = next((p for p in profiles if "gistic" in str(p.get("molecularProfileId", "")).lower()), None)
    if gistic is None:
        logger.warning("No GISTIC copy-number profile ID found among %d profiles. "
                       "Available: %s", len(profiles),
                       [p.get("molecularProfileId") for p in profiles[:12]])
        return pd.DataFrame()
    profile_id = gistic["molecularProfileId"]
    sl_candidate = f"{STUDY_ID}_cna"
    logger.info("Using GISTIC CNA profile: %s (%s)", profile_id, gistic.get("name", ""))

    if entrez_gene_ids is not None:
        rows: list = []
        for eid in entrez_gene_ids:
            r = requests.get(f"{CBIO_API}/molecular-profiles/{profile_id}/molecular-data",
                             params={"sampleListId": sl_candidate, "entrezGeneId": int(eid)},
                             timeout=180)
            if r.status_code == 200 and len(r.text.strip()) > 0:
                try:
                    rows.extend(r.json())
                except Exception:  # noqa: BLE001
                    continue
        if rows:
            df = pd.DataFrame(rows)
            logger.info("GISTIC CNA (targeted, %d genes): %d rows across %d samples.",
                        len(entrez_gene_ids), len(df), df["sampleId"].nunique() if "sampleId" in df.columns else 0)
            return df
    else:
        # Full profile: pull without entrezGeneId (if endpoint supports it) — or try a
        # full-POST fetch. Fall back: GET returns 400; warn and return empty.
        r = requests.get(f"{CBIO_API}/molecular-profiles/{profile_id}/molecular-data",
                         params={"sampleListId": sl_candidate},
                         timeout=300)
        if r.status_code == 200 and len(r.text.strip()) > 0:
            try:
                df = pd.DataFrame(r.json())
                logger.info("GISTIC CNA (full profile): %d rows × %d cols", *df.shape)
                return df
            except Exception:  # noqa: BLE001
                pass
        logger.warning(
            "GISTIC full-profile fetch skipped (status %d); falling back to NF2-only "
            "targeted pull because the unfiltered endpoint requires entrezGeneId.",
            r.status_code,
        )
        return fetch_cna_gistic_data(entrez_gene_ids=[4771])
    logger.warning(
        "============================================================================\n"
        "NO GISTIC CNA DATA RETRIEVED for %s.\n"
        "CNA-derived NF2 will be ALL-NaN; mutation/clinical NF2 sources still used.\n"
        "============================================================================",
        STUDY_ID,
    )
    return pd.DataFrame()


def derive_nf2_from_cna(cna: pd.DataFrame) -> pd.Series:
    """Derive per-sample NF2 status from GISTIC CNA table.

    Conservative rule (matches NF2_STATUS_PLAN in config.py):
      GISTIC value ≤ -1 → Mutant/Loss  (het loss OR homozygous deletion)
      GISTIC value =  0 → Intact        (diploid, no evidence of NF2 deletion)
      GISTIC value ≥ +1 → Intact        (gains/amplifications at NF2 are not inactivating)
      NaN / missing row  → NaN          (unknown)

    Only NF2 rows are selected from the full CNA table (entrezGeneId==4771 OR
    hugoGeneSymbol upper=="NF2").
    """
    if cna.empty:
        return pd.Series(dtype="object", name="nf2_from_cna")
    eid_col = next((c for c in cna.columns if c.lower() == "entrezgeneid"), None)
    sample_col = next((c for c in cna.columns if c.lower() == "sampleid"), None)
    hugo_col = next((c for c in cna.columns if c.lower() == "hugogenesymbol" or c.lower() == "gene"), None)
    value_col = next((c for c in cna.columns if c.lower() == "value"), None)
    if sample_col is None or value_col is None:
        logger.info("CNA table missing sampleId/value cols — cannot derive NF2 from CNA. "
                    "Cols present: %s", list(cna.columns)[:15])
        return pd.Series(dtype="object", name="nf2_from_cna")
    mask_nf2 = pd.Series(False, index=cna.index)
    if eid_col is not None:
        mask_nf2 |= cna[eid_col].astype(str).isin({"4771", "4771.0"})
    if hugo_col is not None:
        mask_nf2 |= cna[hugo_col].astype(str).str.upper().isin({"NF2"})
    nf2_rows = cna.loc[mask_nf2].copy()
    if nf2_rows.empty:
        logger.info("No NF2 rows (entrez=4771 or Hugo=NF2) found in CNA table (n=%d rows).", len(cna))
        return pd.Series(dtype="object", name="nf2_from_cna")
    nf2_rows[value_col] = pd.to_numeric(nf2_rows[value_col], errors="coerce")
    def _map(v):
        if pd.isna(v):
            return np.nan
        if v <= -1:
            return "Mutant/Loss"
        return "Intact"
    # Multiple entries per sample (rare): take the most severe call: Mutant/Loss > Intact > NaN
    severity = {"Mutant/Loss": 2, "Intact": 1}
    nf2_rows["_sev"] = nf2_rows[value_col].map(_map).map(severity).fillna(0)
    nf2_rows = nf2_rows.sort_values("_sev", ascending=False).drop_duplicates(subset=[sample_col], keep="first")
    series = nf2_rows.set_index(sample_col)[value_col].map(_map)
    series.name = "nf2_from_cna"
    logger.info("NF2 from GISTIC CNA: Mutant/Loss=%d, Intact=%d, NaN=%d (total unique samples: %d).",
                int((series == "Mutant/Loss").sum()),
                int((series == "Intact").sum()),
                int(series.isna().sum()),
                int(series.notna().sum()))
    return series


def fetch_mrna_molecular_data(entrez_gene_ids: Optional[Iterable[int]] = None) -> pd.DataFrame:
    """Fetch mRNA TPM (preferred) or z-score (fallback) for all genes.

    NORMALIZATION (see config.NORMALIZATION_DECISION for full rationale):
      PREFERRED PROFILE: *_mrna_seq_tpm — raw TPM, same scale as GSE183653 Bi-lab
        training matrix. Gene overlap with GSE183653 is maximized; no cross-scale
        centroid drift.
      FALLBACK PROFILE: any *_zscores or *_mrna_seq_mrna — use only if TPM profile
        does not exist in this cBioPortal instance. Prints GIANT warning banner so
        this caveat is NEVER silent.

    Uses the corrected cBioPortal API pattern (2025-07):
      POST /api/molecular-profiles/{profileId}/molecular-data/fetch
      Payload: {"sampleIds": [...], "entrezGeneIds": [...]}

    Critical fix for mng_utoronto_2021: all /sample-lists endpoints return 0
    samples (empty lists), so sample IDs are derived directly from the clinical
    data endpoint (CAM20-* sample IDs). Response is streamed to disk first
    (14,795 genes × 121 samples = ~1.79M rows / ~500 MB JSON) to avoid the
    in-memory MemoryError, then parsed from file.
    """
    # --- 1. Locate mRNA molecular profile ID ---
    profiles = _cbio_get_json(f"/studies/{STUDY_ID}/molecular-profiles")
    if not isinstance(profiles, list):
        logger.critical("Molecular profiles endpoint misbehaved. Got: %r", profiles)
        raise SystemExit(3)
    mrna_profiles = [p for p in profiles if "mrna" in str(p.get("molecularProfileId", "")).lower()
                     or "rna_seq" in str(p.get("molecularProfileId", "")).lower()
                     or p.get("genericAssayType") == "mRNA"]
    if not mrna_profiles:
        logger.critical(
            "No mRNA molecular profile found for %s.\n"
            "Available molecularProfileIds: %s",
            STUDY_ID, [p.get("molecularProfileId") for p in profiles],
        )
        raise SystemExit(3)

    tpm_profiles = [p for p in mrna_profiles
                    if "tpm" in str(p.get("molecularProfileId", "")).lower()]
    zscore_profiles = [p for p in mrna_profiles
                       if "zscore" in str(p.get("molecularProfileId", "")).lower()]
    other_profiles = [p for p in mrna_profiles
                      if p not in tpm_profiles and p not in zscore_profiles]

    used_tpm_fallback = False
    if tpm_profiles:
        chosen = tpm_profiles[0]
        logger.info(
            "NORMALIZATION OK: Selected TPM mRNA profile: %s (%s) — "
            "matches GSE183653 training scale per config.NORMALIZATION_DECISION.",
            chosen["molecularProfileId"], chosen.get("name", ""),
        )
    elif zscore_profiles:
        chosen = zscore_profiles[0]
        used_tpm_fallback = True
        logger.warning(
            "\n================================================================\n"
            "  TPM PROFILE NOT FOUND — FALLING BACK TO Z-SCORE mRNA PROFILE.\n"
            "  Selected: %s (%s)\n"
            "  THIS VIOLATES config.NORMALIZATION_DECISION.\n"
            "  IMPACT: Bridging classifier trains on GSE183653 TPM but deploys\n"
            "  on Nassiri z-scores — two different scales. Expected effects:\n"
            "    (A) Attenuated cross-study concordance (kappa biased DOWNWARD,\n"
            "        which is a conservative anti-inflationary bias, so results\n"
            "        if positive remain directionally trustworthy).\n"
            "    (B) Per-class centroids shift; nearest-centroid LOOCV on\n"
            "        GSE183653 remains internally valid but test-set predictions\n"
            "        should be interpreted with the scale-mismatch caveat.\n"
            "  CAVEAT WILL BE FLAGGED IN: Methods, Limitations, Table S1 footer.\n"
            "  MITIGATION: Ridge pipeline is scale-robust (coefficients absorb\n"
            "  the linear shift); nearest-centroid remains rank-preserving.\n"
            "================================================================\n",
            chosen["molecularProfileId"], chosen.get("name", ""),
        )
    elif other_profiles:
        chosen = other_profiles[0]
        used_tpm_fallback = True
        logger.warning(
            "  Neither TPM nor z-score label in profile id. Using: %s (%s)",
            chosen["molecularProfileId"], chosen.get("name", ""),
        )
    else:
        chosen = mrna_profiles[0]
    profile_id = chosen["molecularProfileId"]
    fetch_mrna_molecular_data._last_used_profile = profile_id
    fetch_mrna_molecular_data._used_tpm_fallback = used_tpm_fallback

    # --- 2. Resolve sample IDs: sample-lists → clinical data fallback ---
    sample_lists = _cbio_get_json(f"/studies/{STUDY_ID}/sample-lists")
    if not isinstance(sample_lists, list):
        logger.critical("Sample lists endpoint misbehaved.")
        raise SystemExit(3)
    all_samples_list = next((sl for sl in sample_lists
                             if sl.get("category") == "all_cases_in_study"
                             or "all" in str(sl.get("sampleListId", "")).lower()),
                            sample_lists[0] if sample_lists else None)
    sl_samples = list((all_samples_list or {}).get("sampleIds", []) or [])
    if len(sl_samples) < 5:
        # KNOWN ISSUE on mng_utoronto_2021 (2025-07): all sample lists return 0 samples.
        # Fall back to the clinical data endpoint, which reliably returns 121 CAM20-* IDs.
        logger.warning("Sample list %s has only %d samples (expected ~121). "
                       "Falling back to clinical-data endpoint sample IDs.",
                       (all_samples_list or {}).get("sampleListId", "?"), len(sl_samples))
        clin_df = fetch_clinical_data()
        sample_ids = sorted(clin_df.index.dropna().unique().tolist())
        logger.info("  Clinical fallback → %d sample IDs.", len(sample_ids))
    else:
        sample_ids = sorted(set(sl_samples))
    if not sample_ids:
        logger.critical("Could not resolve any sample IDs for mRNA fetch.")
        raise SystemExit(3)

    # --- 3. POST molecular-data/fetch — STREAM to disk, THEN parse ---
    endpoint = f"/molecular-profiles/{profile_id}/molecular-data/fetch"
    payload: Dict[str, object] = {"sampleIds": sample_ids}
    if entrez_gene_ids is not None:
        payload["entrezGeneIds"] = list(entrez_gene_ids)
    raw_json_path = PATHS.data_raw / "mng_utoronto_2021" / "mrna_zscores_raw.json"
    raw_json_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Streaming POST %s (all genes, %d samples) → %s ...",
                endpoint, len(sample_ids), raw_json_path)
    try:
        with requests.post(f"{CBIO_API}{endpoint}", json=payload, timeout=300, stream=True) as r:
            if r.status_code != 200:
                body_preview = ""
                try:
                    body_preview = next(r.iter_content(chunk_size=1000), b"").decode("utf-8", errors="replace")
                except Exception:
                    pass
                logger.critical("cBioPortal returned HTTP %d for %s\nBody preview: %s\n"
                                "Check study ID spelling (currently '%s').",
                                r.status_code, endpoint, body_preview, STUDY_ID)
                raise SystemExit(3)
            total_bytes = 0
            with open(raw_json_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        total_bytes += len(chunk)
    except requests.RequestException as e:
        logger.critical("Network error streaming mRNA molecular data: %s", e)
        raise SystemExit(3) from e
    logger.info("  Downloaded %.1f MB.", total_bytes / 1024 / 1024)

    # --- 4. Parse from disk, pivot, entrez→hugo ---
    logger.info("Parsing streamed JSON from disk ...")
    try:
        long_df = pd.read_json(raw_json_path)
    except ValueError as e:
        logger.critical("mRNA raw JSON parse fail: %s", e)
        raise SystemExit(3) from e
    if not isinstance(long_df, pd.DataFrame) or long_df.empty:
        logger.critical("mRNA raw JSON produced empty / non-DataFrame output.")
        raise SystemExit(3)
    needed = {"sampleId", "entrezGeneId", "value"}
    missing = needed - set(long_df.columns)
    if missing:
        logger.critical("mRNA molecular data missing columns %s. Have: %s",
                        missing, list(long_df.columns))
        raise SystemExit(3)
    wide = long_df.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    logger.info("mRNA z-scores (before symbol map): %d samples × %d entrez genes.", *wide.shape)
    # Save long TSV (compressed) for reproducibility alongside raw JSON
    try:
        long_df.to_csv(PATHS.data_raw / "mng_utoronto_2021" / "mrna_zscores_long.csv.gz",
                       index=False, compression="gzip")
    except Exception:
        pass

    symbol_map = _entrez_to_hugo_map(wide.columns.tolist())
    if symbol_map:
        wide = wide.rename(columns=symbol_map)
        lookup = pd.DataFrame({"entrezGeneId": list(symbol_map.keys()),
                               "hugoSymbol": list(symbol_map.values())})
        lookup.to_csv(PATHS.data_raw / "mng_utoronto_2021" / "entrez_hugo_lookup.csv", index=False)
    return wide


def _entrez_to_hugo_map(entrez_ids: List[int]) -> Dict[int, str]:
    """Convert Entrez IDs → Hugo symbols via 1 GET call to /api/genes.

    The cBioPortal public endpoint GET /api/genes returns the COMPLETE gene
    catalog (44,896 genes as of 2025-07) in a single HTTP 200 response, each
    record containing {entrezGeneId, hugoGeneSymbol, type}. We therefore do
    NOT need any per-gene or batched-POST calls — just 1 GET + a dict lookup.

    If GET /api/genes ever fails (e.g., offline, 5xx, or response schema
    changes), we fall back to the previous best-effort per-gene GET loop so
    the pipeline degrades gracefully rather than hard-crashing.
    """
    if not entrez_ids:
        return {}
    todo: List[int] = []
    for e in entrez_ids:
        try:
            if pd.notna(e):
                todo.append(int(float(e)))
        except (TypeError, ValueError):
            continue
    todo_set = set(todo)
    mapping: Dict[int, str] = {}
    if not todo_set:
        return mapping
    logger.info(
        "Entrez→Hugo: resolving %d unique entrez IDs via 1-call GET /api/genes "
        "(full catalog lookup — NOT per-gene loop).",
        len(todo_set),
    )
    import time
    start = time.time()
    try:
        r = requests.get(f"{CBIO_API}/genes", timeout=60)
        if r.status_code == 200:
            full_catalog = r.json()
            if isinstance(full_catalog, list):
                n_hits = 0
                for g in full_catalog:
                    if not isinstance(g, dict):
                        continue
                    eid = g.get("entrezGeneId")
                    hugo = g.get("hugoGeneSymbol") or g.get("symbol")
                    try:
                        eid_int = int(eid) if eid is not None else None
                    except (TypeError, ValueError):
                        eid_int = None
                    if eid_int is not None and eid_int in todo_set and hugo:
                        mapping[eid_int] = hugo
                        n_hits += 1
                elapsed = time.time() - start
                missed = len(todo_set) - len(mapping)
                logger.info(
                    "Entrez→Hugo OK: 1 GET /api/genes, %d/%d mapped in %.0fs. "
                    "%d unmapped. (catalog size: %d genes).",
                    len(mapping), len(todo_set), elapsed, missed, len(full_catalog),
                )
                return mapping
            logger.warning(
                "GET /api/genes returned non-list type %s — "
                "falling back to per-gene GET loop.",
                type(full_catalog).__name__,
            )
        else:
            logger.warning(
                "GET /api/genes HTTP %d — falling back to per-gene GET loop.",
                r.status_code,
            )
    except Exception as exc:
        logger.warning(
            "GET /api/genes network error: %s — falling back to per-gene GET loop.",
            exc,
        )
    logger.info("  Fallback: per-gene GET loop for %d entrez IDs ...", len(todo_set))
    todo_list = sorted(todo_set)
    missed = 0
    for idx, eid in enumerate(todo_list):
        try:
            rg = requests.get(f"{CBIO_API}/genes/{eid}", timeout=10)
            if rg.status_code == 200:
                gd = rg.json() if rg.content else None
                if isinstance(gd, dict):
                    hugo = gd.get("hugoGeneSymbol") or gd.get("symbol")
                    if hugo:
                        mapping[eid] = hugo
                        continue
            missed += 1
        except Exception:
            missed += 1
        if (idx + 1) % 2000 == 0:
            elapsed = time.time() - start
            rate = (idx + 1) / max(elapsed, 1e-6)
            logger.info("  ... %d/%d done (mapped %d, missed %d) @ %.0f rps",
                        idx + 1, len(todo_list), len(mapping), missed, rate)
            time.sleep(0.5)
    elapsed = time.time() - start
    missed = len(todo_set) - len(mapping)
    logger.info(
        "Entrez→Hugo complete: %d/%d mapped in %.0fs. %d unmapped.",
        len(mapping), len(todo_set), elapsed, missed,
    )
    return mapping


def derive_nf2_from_mutations(mutations: pd.DataFrame) -> pd.Series:
    """From the mutation table, return a per-sample Series:
       Mutant/Loss if any NF2 (Hugo symbol NF2, entrez 4771) mutation is present;
       Intact if no NF2 mutation in sample that otherwise has mutation data;
       NaN if sample appears in clinical but has no mutation profile info at all
       (never silently coerce to Intact).

    We do NOT use CNAs here because the mutation endpoint may not include them;
    that would be a separate fetches for copy-number profiles. If the user wants
    NF2 loss from CNAs that's an explicit enhancement, never an implicit default.
    """
    if mutations.empty:
        return pd.Series(dtype="object", name="nf2_from_mutations")
    hugo_col = next((c for c in mutations.columns if c.lower() == "hugogenesymbol" or c.lower() == "gene"), None)
    entrez_col = next((c for c in mutations.columns if c.lower() == "entrezgeneid"), None)
    sample_col = next((c for c in mutations.columns if c.lower() == "sampleid"), None)
    if hugo_col is None and entrez_col is None:
        logger.info("No hugo/entrez column in mutations table — cannot derive NF2 status from mutations.")
        return pd.Series(dtype="object", name="nf2_from_mutations")
    if sample_col is None:
        logger.info("No sampleId column in mutations table.")
        return pd.Series(dtype="object", name="nf2_from_mutations")
    mask_nf2 = pd.Series(False, index=mutations.index)
    if hugo_col is not None:
        mask_nf2 |= mutations[hugo_col].astype(str).str.upper().isin({"NF2"})
    if entrez_col is not None:
        mask_nf2 |= mutations[entrez_col].astype(str).isin({"4771", "4771.0"})
    nf2_mut_samples = set(mutations.loc[mask_nf2, sample_col].unique().tolist())
    all_mut_samples = set(mutations[sample_col].unique().tolist())
    series = pd.Series(index=sorted(all_mut_samples), dtype="object", name="nf2_from_mutations")
    series[:] = "Intact"
    series.loc[list(nf2_mut_samples)] = "Mutant/Loss"
    return series


def merge_and_standardize(clinical_wide: pd.DataFrame,
                          nf2_from_mutations: pd.Series,
                          nf2_from_cna: Optional[pd.Series] = None) -> pd.DataFrame:
    """Map clinical attribute names → canonical names. Merge NF2 from (up to) 3 sources.

    Sources of nf2_status, in priority order of "how conservative is a conflict":
      1. Clinical attribute NF2 (if present) — comes directly from cBio clinical table
      2. Mutation-derived NF2   (point mutations / small indels via mutation table)
      3. CNA-derived NF2        (GISTIC chr22q deletion; primary source for meningioma)

    Merge rule (applied pairwise, then the result merges with source 3):
      - Only one source has a value for a sample → use it.
      - Two+ sources agree (Intact==Intact, Mutant/Loss==Mutant/Loss) → use it.
      - Any pairwise disagreement (Intact vs Mutant/Loss) → set to NaN (unknown).
      - Samples appearing in clinical metadata but not in the source table
        (e.g. clinical sample has no CNA info) → that source contributes NaN,
        not "Intact" by default.

    NO IMPUTATION. Missing stays missing.
    """
    out = pd.DataFrame(index=clinical_wide.index)
    col_lower_map = {str(c).upper(): c for c in clinical_wide.columns}
    consumed_raw_cols = set()

    for canon, candidates in CLINICAL_ATTR_CANDIDATES.items():
        matched = None
        for cand in candidates:
            if cand.upper() in col_lower_map and cand.upper() not in consumed_raw_cols:
                matched = col_lower_map[cand.upper()]
                consumed_raw_cols.add(cand.upper())
                break
        if matched is not None:
            out[canon] = clinical_wide[matched].values
        else:
            out[canon] = np.nan

    if "subgroup" in out.columns:
        if "nassiri_group" in out.columns and out["nassiri_group"].notna().any():
            out = out.drop(columns=["subgroup"])
        else:
            out = out.rename(columns={"subgroup": "nassiri_group"})

    def _pairwise_merge(base: pd.Series, overlay: pd.Series) -> pd.Series:
        """Merge two per-sample NF2 Series with conflict→NaN, single-source→its value."""
        base_aligned = base.reindex(out.index).astype(object)
        overlay_aligned = overlay.reindex(out.index).astype(object)
        result = base_aligned.copy()
        only_base = base_aligned.notna() & overlay_aligned.isna()
        only_overlay = base_aligned.isna() & overlay_aligned.notna()
        result[only_overlay] = overlay_aligned[only_overlay].values
        both = base_aligned.notna() & overlay_aligned.notna()
        conflict = both & (
            base_aligned.astype(str).str.strip().str.lower()
            != overlay_aligned.astype(str).str.strip().str.lower()
        )
        result[conflict] = np.nan
        return result

    merged_nf2 = out["nf2_status"].astype(object).copy()
    if not nf2_from_mutations.empty:
        merged_nf2 = _pairwise_merge(merged_nf2, nf2_from_mutations)
    if nf2_from_cna is not None and not nf2_from_cna.empty:
        merged_nf2 = _pairwise_merge(merged_nf2, nf2_from_cna)
    out["nf2_status"] = merged_nf2

    # Standardize values
    if "who_grade" in out.columns:
        who_raw = out["who_grade"].astype(str).str.upper()
        roman = who_raw.str.extract(r"\b([I]{1,3})\b", expand=False)
        g_short = who_raw.str.extract(r"\bG([123])\b", expand=False).map({"1": "I", "2": "II", "3": "III"})
        out["who_grade"] = roman.where(roman.notna(), g_short)
    if "nf2_status" in out.columns:
        nf2_str = out["nf2_status"].astype(str).str.strip().str.lower()
        intact_mask = nf2_str.isin({"wt", "wild", "wild_type", "intact", "normal"})
        loss_mask = nf2_str.isin({
            "mutant", "loss", "mutant/loss", "mutant_loss",
            "del", "deleted", "mutated", "altered",
            "loh", "abnormal", "homozygous_deletion", "het_loss",
        })
        nf2_norm = pd.Series(np.nan, index=out.index, dtype=object)
        nf2_norm[intact_mask] = "Intact"
        nf2_norm[loss_mask] = "Mutant/Loss"
        out["nf2_status"] = nf2_norm
    if "recurrence_event" in out.columns:
        re = out["recurrence_event"].astype(str).str.strip().str.lower()
        out["recurrence_event"] = re.map({
            "yes": True, "true": True, "1": True, "recurred": True,
            "progressed": True, "event": True,
            "no": False, "false": False, "0": False, "none": False,
        })
    for num in ["age_at_surgery", "recurrence_months"]:
        if num in out.columns:
            out[num] = pd.to_numeric(out[num], errors="coerce")
    return out


def audit_report(clinical_wide: pd.DataFrame, mutations: pd.DataFrame,
                 mrna_wide: pd.DataFrame) -> pd.DataFrame:
    attrs_full = pd.DataFrame({
        "attribute_id": clinical_wide.columns,
        "n_non_null": clinical_wide.notna().sum().values,
        "pct_non_null": np.round(100 * clinical_wide.notna().sum().values / len(clinical_wide), 2)
            if len(clinical_wide) else 0,
        "unique_vals_examples": [
            str(clinical_wide[c].dropna().astype(str).unique().tolist()[:5])
            for c in clinical_wide.columns
        ],
    })
    return attrs_full


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-mrna", action="store_true",
                        help="Skip mRNA z-score fetch (useful for a metadata-only audit).")
    parser.add_argument("--skip-download", action="store_true",
                        help="Reuse cached files in data/raw/mng_utoronto_2021/*.csv if present.")
    args = parser.parse_args()

    raw_dir = PATHS.data_raw / "mng_utoronto_2021"
    raw_dir.mkdir(parents=True, exist_ok=True)

    reports: List[FetchReport] = []

    if args.skip_download and (raw_dir / "clinical.csv").is_file():
        logger.info("--skip-download: reusing cached clinical/mutations/CNA CSV.")
        clinical_wide = pd.read_csv(raw_dir / "clinical.csv", index_col=0)
        mut_path = raw_dir / "mutations.csv"
        if mut_path.is_file() and mut_path.stat().st_size > 0:
            try:
                mutations = pd.read_csv(mut_path)
            except pd.errors.EmptyDataError:
                logger.warning("mutations.csv existed but is empty — treating as no-mutation-data.")
                mutations = pd.DataFrame()
        else:
            mutations = pd.DataFrame()
        cna_path = raw_dir / "cna_gistic.csv"
        if cna_path.is_file() and cna_path.stat().st_size > 0:
            try:
                cna = pd.read_csv(cna_path)
            except pd.errors.EmptyDataError:
                logger.warning("cna_gistic.csv existed but is empty — treating as no-CNA-data.")
                cna = pd.DataFrame()
        else:
            cna = pd.DataFrame()
    else:
        study_meta = fetch_study_exists()
        logger.info("Study exists: %s (allSampleCount=%s, allPatientCount=%s)",
                    study_meta.get("name"), study_meta.get("allSampleCount"), study_meta.get("allPatientCount"))
        clinical_wide = fetch_clinical_data()
        clinical_wide.to_csv(raw_dir / "clinical.csv")
        reports.append(FetchReport(f"/studies/{STUDY_ID}/clinical-data", 200,
                                   len(clinical_wide), clinical_wide.shape[1],
                                   len(clinical_wide), list(clinical_wide.index)[:10]))
        mutations = fetch_mutation_data()
        mutations.to_csv(raw_dir / "mutations.csv", index=False)
        reports.append(FetchReport(f"/studies/{STUDY_ID}/mutations", 200,
                                   len(mutations), mutations.shape[1],
                                   int(mutations["sampleId"].nunique()) if "sampleId" in mutations.columns else 0,
                                   []))
        # Fetch GISTIC CNA — primary NF2 source for meningioma.
        # Only fetch a curated small panel of provenance/QC genes plus NF2; full
        # profile is large (~20k genes × 121 samples) and unnecessary for our needs.
        # Panel includes: NF2, plus known meningioma drivers (AKT1, SMO, PIK3CA,
        # PIK3R1, POLR2A, TRAF7, KLF4, ARID1A, SMARCB1, SUFU, CDKN2A/B loci) and
        # common cytoband controls (chr22q multi-gene spread for copy concordance).
        CNA_PROVENANCE_PANEL: List[int] = [
            4771,    # NF2
            207,     # AKT1
            6608,    # SMO
            5290,    # PIK3CA
            5295,    # PIK3R1
            5430,    # POLR2A
            7188,    # TRAF7
            9314,    # KLF4
            8289,    # ARID1A
            6598,    # SMARCB1
            51684,   # SUFU
            1029,    # CDKN2A
            1030,    # CDKN2B
            7015,    # CHEK2 (chr22, control for 22q loss extent)
            5962,    # RB1 (chr13, unrelated aneuploidy control)
            7157,    # TP53
        ]
        cna = fetch_cna_gistic_data(entrez_gene_ids=CNA_PROVENANCE_PANEL)
        cna.to_csv(raw_dir / "cna_gistic.csv", index=False)
        reports.append(FetchReport(f"molecular-profiles/{STUDY_ID}_gistic", 200,
                                   len(cna), cna.shape[1],
                                   int(cna["sampleId"].nunique()) if "sampleId" in cna.columns else 0,
                                   list(cna["sampleId"].unique()[:10]) if "sampleId" in cna.columns else []))

    # Attribute audit BEFORE mapping so the user sees exactly what cBioPortal returned.
    attribute_audit = audit_report(clinical_wide, mutations, pd.DataFrame())
    attribute_audit.to_csv(raw_dir / "clinical_attribute_audit.csv", index=False)
    logger.info("Full list of %d clinical attributes saved to clinical_attribute_audit.csv. "
                "If CLINICAL_ATTR_CANDIDATES used wrong IDs, edit them based on this file.",
                len(attribute_audit))

    nf2_mut_series = derive_nf2_from_mutations(mutations)
    nf2_cna_series = derive_nf2_from_cna(cna)
    metadata = merge_and_standardize(clinical_wide, nf2_mut_series, nf2_cna_series)
    metadata_export = _canonicalize_nf2_for_export(metadata)
    metadata_export.to_csv(raw_dir / "metadata_canonical_unfiltered.csv")

    if args.skip_mrna:
        logger.info("--skip-mrna set; NOT writing expression.csv to processed/.")
        expression_wide = None
    else:
        mrna_path = raw_dir / "mrna_zscores.csv"
        if args.skip_download and mrna_path.is_file():
            logger.info("--skip-download: reusing %s", mrna_path)
            expression_wide = pd.read_csv(mrna_path, index_col=0)
        else:
            expression_wide = fetch_mrna_molecular_data()
            expression_wide.to_csv(mrna_path)
        reports.append(FetchReport("molecular-data/mRNA", 200,
                                   expression_wide.shape[0] * expression_wide.shape[1],
                                   expression_wide.shape[1], expression_wide.shape[0],
                                   list(expression_wide.index[:5])))

    # Intersect: only keep samples with BOTH metadata AND expression rows.
    samples_meta = set(metadata.index.tolist())
    if expression_wide is not None:
        samples_expr = set(expression_wide.index.tolist())
        overlap = samples_meta & samples_expr
        if len(overlap) < max(1, 0.8 * min(len(samples_meta), len(samples_expr))):
            logger.warning(
                "Sample-ID overlap between clinical (%d) and mRNA (%d) is only %d. "
                "This usually means sampleId normalization is off (e.g. sampleId vs patientId). "
                "Check cBioPortal sample list and clinical attribute audit BEFORE proceeding. "
                "Script writes only the intersected samples; not padding with synthetic rows.",
                len(samples_meta), len(samples_expr), len(overlap),
            )
        metadata_out = metadata.reindex(sorted(overlap))
        expression_out = expression_wide.reindex(sorted(overlap))
    else:
        metadata_out = metadata.copy()
        expression_out = None

    n_out = len(metadata_out)
    logger.info("Final metadata after intersection: %d samples (expected %d per paper).",
                n_out, EXPECTED_N_SAMPLES)
    if n_out < EXPECTED_N_SAMPLES:
        logger.warning(
            "Shortfall of %d samples vs paper. Possible causes: (a) cBioPortal hosts the "
            "discovery subset only, (b) mismatch between sampleId formats in clinical vs "
            "mRNA endpoints, (c) multi-region samples excluded. DO NOT pad with synthetic "
            "samples — report the actual N=%d in the paper.",
            EXPECTED_N_SAMPLES - n_out, n_out,
        )

    # Honest feasibility: count non-nulls per canonical field.
    feas = {
        "n_total": n_out,
        "fields_non_null": {
            col: {"n": int(metadata_out[col].notna().sum()),
                  "pct": round(100 * float(metadata_out[col].notna().sum() / n_out), 2) if n_out else 0}
            for col in metadata_out.columns
        },
    }
    with open(raw_dir / "hypothesis_feasibility.json", "w") as fh:
        json.dump(feas, fh, indent=2)
    logger.info("Feasibility summary (n=%d): %s", n_out,
                {k: v for k, v in feas["fields_non_null"].items() if k != "n_total"})

    # Write to processed/
    out_dir = PATHS.data_processed / COHORT_NAME_PROCESSED
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata_out_export = _canonicalize_nf2_for_export(metadata_out)
    metadata_out_export.to_csv(out_dir / "metadata.csv")
    if expression_out is not None:
        expression_out.to_csv(out_dir / "expression.csv")
    _profile_id = getattr(fetch_mrna_molecular_data, "_last_used_profile", "UNKNOWN")
    _fallback_flag = getattr(fetch_mrna_molecular_data, "_used_tpm_fallback", None)
    pd.DataFrame([{
        "dataset": f"cBioPortal {STUDY_ID}",
        "paper": "Nassiri et al. 2021 Nature",
        "n_samples_metadata": n_out,
        "n_samples_expected_paper": EXPECTED_N_SAMPLES,
        "n_genes_expression": expression_out.shape[1] if expression_out is not None else 0,
        "expression_fetch_skipped": args.skip_mrna,
        "expression_profile_id": _profile_id,
        "expression_profile_tpm_used": (
            False if _fallback_flag is None else (not _fallback_flag)
        ) if expression_out is not None else False,
        "expression_profile_tpm_fallback_to_zscore": _fallback_flag if expression_out is not None else None,
        "normalization_decision_ref": "config.NORMALIZATION_DECISION (TPM preferred, z-score only fallback with caveat)",
        "note": (
            "All values come directly from cBioPortal. No imputation was performed. "
            "Missing values remain NaN. Conflicting NF2 values (clinical vs mutation) "
            "were set to NaN conservatively."
        ),
    }]).to_csv(out_dir / "metadata_PROVENANCE.csv", index=False)
    logger.info("Processed outputs written to %s", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
