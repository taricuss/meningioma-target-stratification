from __future__ import annotations

import logging
import pickle
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import requests

from .config import (
    CBIO_MENINGIOMA_STUDY_IDS,
    NF2_EXPORT_COL,
    NF2_INTERNAL_COL,
    PATHS,
    REPLICATION_CANDIDATE_STUDIES,
)

logger = logging.getLogger(__name__)

GEO_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CBIO_API = "https://www.cbioportal.org/api"


@dataclass
class CohortRecord:
    accession: str
    title: str
    source: str
    role: str
    platform: str
    n_samples: int
    has_nf2: bool
    has_grade: bool
    has_recurrence: bool
    notes: str = ""
    download_url: Optional[str] = None


def _geo_esearch(term: str, retmax: int = 50) -> List[str]:
    params = {"db": "gds", "term": term, "retmax": retmax, "retmode": "json", "sort": "relevance"}
    r = requests.get(f"{GEO_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    data = r.json().get("esearchresult", {})
    return data.get("idlist", []) or []


def _geo_esummary(gds_ids: List[str]) -> List[dict]:
    if not gds_ids:
        return []
    params = {"db": "gds", "id": ",".join(gds_ids), "retmode": "json"}
    r = requests.get(f"{GEO_BASE}/esummary.fcgi", params=params, timeout=60)
    r.raise_for_status()
    result = r.json().get("result", {})
    return [result[gid] for gid in gds_ids if gid in result]


def geo_screen_meningioma_expressions() -> pd.DataFrame:
    rows = []
    for accession in REPLICATION_CANDIDATE_STUDIES:
        try:
            summary = _geo_esummary([accession.replace("GSE", "")])
            info = summary[0] if summary else {}
            title = info.get("title", "")
            n_samples = info.get("n_samples", np.nan)
            platform = info.get("gpl", "")
            rows.append(
                {
                    "accession": accession,
                    "title": title,
                    "n_samples_expected": int(n_samples) if pd.notna(n_samples) else np.nan,
                    "platform": platform,
                    "role": "replication_candidate",
                    "source": "GEO",
                    "open_access": True,
                }
            )
            time.sleep(0.3)
        except Exception as e:
            logger.warning("Failed fetching %s: %s", accession, e)
            rows.append(
                {
                    "accession": accession,
                    "title": "",
                    "n_samples_expected": np.nan,
                    "platform": "",
                    "role": "replication_candidate",
                    "source": "GEO",
                    "open_access": np.nan,
                }
            )
    return pd.DataFrame(rows)


def cbioportal_list_meningioma_studies() -> pd.DataFrame:
    try:
        r = requests.get(f"{CBIO_API}/studies?projection=DETAILED", timeout=30)
        r.raise_for_status()
        studies = r.json()
        rows = []
        for s in studies:
            name = (s.get("name") or "").lower()
            if "mening" in name or s.get("studyId") in CBIO_MENINGIOMA_STUDY_IDS:
                rows.append(
                    {
                        "study_id": s.get("studyId"),
                        "name": s.get("name"),
                        "description": s.get("description"),
                        "citation": s.get("citation"),
                        "n_samples": s.get("allSampleCount"),
                        "n_patients": s.get("allPatientCount"),
                        "has_mutation": s.get("mutation"),
                        "has_cna": s.get("cna"),
                        "has_mrna": s.get("mrna"),
                        "has_methylation": s.get("methylation_hm27") or s.get("methylation_hm450"),
                    }
                )
        return pd.DataFrame(rows)
    except Exception as e:
        logger.warning("cBioPortal fetch failed: %s", e)
        return pd.DataFrame(
            [
                {
                    "study_id": sid,
                    "n_samples": np.nan,
                    "note": "cBioPortal API unreachable; ID hard-coded from plan.",
                }
                for sid in CBIO_MENINGIOMA_STUDY_IDS
            ]
        )


def build_synthetic_test_cohort(
    n: int = 185,
    cohort_name: str = "unit_test_synthetic",
    seed: int = 20260317,
) -> Dict[str, pd.DataFrame]:
    """[UNIT-TEST FIXTURE — NEVER FOR RESULTS] Produce a synthetic but biologically
    plausible meningioma cohort for testing that stats.py, ssgsea.py, FDR logic, etc.
    all execute correctly against KNOWN SEEDED GROUND TRUTH.

    Ground-truth subgroup differences are DELIBERATELY SEEDED to match PLAN.md
    hypotheses — statistical "significance" is guaranteed by construction. This
    function validates PIPELINE MACHINERY only, never biology.

    CRITICAL: Outputs from this function must NEVER be copied to results/tables,
    results/figures, manuscript/, or presented to anyone as findings.
    """
    rng = np.random.default_rng(seed)
    sample_ids = [f"{cohort_name}_{i:04d}" for i in range(n)]

    nassiri_probs = np.array([0.30, 0.45, 0.25])
    bi_probs = np.array([0.30, 0.35, 0.35])
    nassiri_groups = rng.choice(
        ["Immunogenic", "NF2-inactivated canonical", "Hypermetabolic"], size=n, p=nassiri_probs
    )
    # Make Bi-lab groups concordant but not identical (kappa ~ 0.55)
    transition = {
        "Immunogenic": [0.10, 0.75, 0.15],
        "NF2-inactivated canonical": [0.70, 0.10, 0.20],
        "Hypermetabolic": [0.15, 0.15, 0.70],
    }
    bi_labels = ["Merlin-intact", "Immune-enriched", "Hypermitotic"]
    bi_groups = np.array(
        [rng.choice(bi_labels, p=transition[g]) for g in nassiri_groups]
    )

    # NF2 status: high in NF2-inactivated canonical, low elsewhere
    nf2 = []
    for g in nassiri_groups:
        if g == "NF2-inactivated canonical":
            nf2.append(rng.choice(["Mutant/Loss", "Intact"], p=[0.82, 0.18]))
        else:
            nf2.append(rng.choice(["Mutant/Loss", "Intact"], p=[0.22, 0.78]))
    nf2 = np.array(nf2)

    # WHO grade correlated with hypermetabolic/hypermitotic
    grade_probs = {
        "Hypermetabolic": [0.20, 0.45, 0.35],
        "Immunogenic": [0.55, 0.35, 0.10],
        "NF2-inactivated canonical": [0.45, 0.40, 0.15],
    }
    grade_order = ["I", "II", "III"]
    who_grade = np.array(
        [rng.choice(grade_order, p=grade_probs[g]) for g in nassiri_groups]
    )

    # Recurrence: 25% events
    has_event = rng.binomial(1, 0.25, size=n).astype(bool)
    time_to_event = rng.uniform(6, 96, size=n)
    time_to_event = np.where(has_event, time_to_event, time_to_event * rng.uniform(0.8, 1.5, size=n))

    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "nassiri_group": pd.Categorical(
                nassiri_groups,
                categories=["Immunogenic", "NF2-inactivated canonical", "Hypermetabolic"],
                ordered=False,
            ),
            "bi_group": pd.Categorical(
                bi_groups, categories=bi_labels, ordered=False
            ),
            NF2_INTERNAL_COL: pd.Categorical(nf2, categories=["Intact", "Mutant/Loss"]),
            "who_grade": pd.Categorical(who_grade, categories=grade_order, ordered=True),
            "recurrence_event": has_event,
            "recurrence_months": time_to_event,
            "age_at_surgery": rng.integers(28, 84, size=n).astype(float),
            "sex": rng.choice(["F", "M"], size=n, p=[0.62, 0.38]),
        }
    ).set_index("sample_id")

    # Expression matrix: ~1500 genes covering the target programs
    from .gene_programs import ALL_PROGRAMS

    all_program_genes = sorted({g for genes in ALL_PROGRAMS.values() for g in genes})
    # Add ~1400 random context genes
    np.random.seed(seed)
    context_genes = [f"GENE_{i:04d}" for i in range(1400)]
    all_genes = all_program_genes + context_genes
    rng_expr = np.random.default_rng(seed + 1)
    expr = pd.DataFrame(
        rng_expr.normal(loc=6.0, scale=1.2, size=(n, len(all_genes))),
        index=sample_ids,
        columns=all_genes,
    )

    # Seed true subgroup effects per program (H2 ground truth)
    program_effect_sizes = {
        "HDAC_Panobinostat_Romidepsin": {
            "Hypermetabolic": +0.9,
            "NF2-inactivated canonical": +0.2,
            "Immunogenic": -0.4,
        },
        "Proteasome_Carfilzomib": {
            "Hypermetabolic": +1.2,
            "NF2-inactivated canonical": +0.3,
            "Immunogenic": -0.5,
        },
        "Tubulin_Ixabepilone": {
            "Hypermetabolic": +1.1,
            "NF2-inactivated canonical": +0.1,
            "Immunogenic": -0.3,
        },
        "Translation_Omacetaxine": {
            "Hypermetabolic": +0.8,
            "NF2-inactivated canonical": +0.5,
            "Immunogenic": -0.3,
        },
        "HDAC8_TGFb_EMT_Resistance": {
            "Immunogenic": +1.0,
            "Hypermetabolic": +0.2,
            "NF2-inactivated canonical": -0.3,
        },
    }

    for program, effect in program_effect_sizes.items():
        genes = ALL_PROGRAMS[program]
        for group, delta in effect.items():
            mask = metadata["nassiri_group"] == group
            for g in genes:
                if g in expr.columns and mask.any():
                    expr.loc[mask, g] += delta + rng_expr.normal(0, 0.2, size=int(mask.sum()))

    # NF2 effect on HDAC1/2 and proteasome
    nf2_mask = metadata[NF2_INTERNAL_COL] == "Mutant/Loss"
    for g in ["HDAC1", "HDAC2", "PSMB5", "EEF2"]:
        if g in expr.columns and nf2_mask.any():
            expr.loc[nf2_mask, g] += 0.6 + rng_expr.normal(0, 0.25, size=int(nf2_mask.sum()))

    return {
        "metadata": metadata,
        "expression": expr,
    }


def _export_nf2_column(meta: pd.DataFrame) -> pd.DataFrame:
    out = meta.copy()
    if NF2_INTERNAL_COL in out.columns and NF2_EXPORT_COL not in out.columns:
        out = out.rename(columns={NF2_INTERNAL_COL: NF2_EXPORT_COL})
    elif NF2_EXPORT_COL in out.columns and NF2_INTERNAL_COL in out.columns:
        out = out.drop(columns=[NF2_INTERNAL_COL])
    return out


def _load_nf2_columns_for_internal(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if NF2_EXPORT_COL in out.columns and NF2_INTERNAL_COL not in out.columns:
        out[NF2_INTERNAL_COL] = out[NF2_EXPORT_COL]
    elif NF2_INTERNAL_COL in out.columns and NF2_EXPORT_COL not in out.columns:
        out[NF2_EXPORT_COL] = out[NF2_INTERNAL_COL]
    return out


def save_cohort(
    cohort_name: str,
    metadata: pd.DataFrame,
    expression: Optional[pd.DataFrame] = None,
    extra: Optional[Dict[str, pd.DataFrame]] = None,
) -> Path:
    out_dir = PATHS.data_processed / cohort_name
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_export = _export_nf2_column(metadata)
    meta_export.to_csv(out_dir / "metadata.csv")
    if expression is not None:
        expression.to_csv(out_dir / "expression.csv")
    if extra:
        for name, df in extra.items():
            if name == "metadata":
                df = _export_nf2_column(df)
            df.to_csv(out_dir / f"{name}.csv")
    return out_dir


def load_cohort(cohort_name: str) -> Dict[str, pd.DataFrame]:
    d = PATHS.data_processed / cohort_name
    result = {}
    for f in d.glob("*.csv"):
        df = pd.read_csv(f, index_col=0)
        if f.stem == "metadata":
            df = _load_nf2_columns_for_internal(df)
        result[f.stem] = df
    return result


def save_aim0_lock(records: List[CohortRecord], path: Optional[Path] = None) -> Path:
    if path is None:
        path = PATHS.results_tables / "aim0_cohort_lock.csv"
    df = pd.DataFrame([asdict(r) for r in records])
    df.to_csv(path, index=False)
    (PATHS.results_tables / "aim0_cohort_lock.pkl").write_bytes(pickle.dumps(records))
    return path
