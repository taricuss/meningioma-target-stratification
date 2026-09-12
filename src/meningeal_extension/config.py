from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np

SEED: int = 20260317

PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Paths:
    project_root: Path = PROJECT_ROOT
    data_raw: Path = PROJECT_ROOT / "data" / "raw"
    data_interim: Path = PROJECT_ROOT / "data" / "interim"
    data_processed: Path = PROJECT_ROOT / "data" / "processed"
    results_tables: Path = PROJECT_ROOT / "results" / "tables"
    results_figures: Path = PROJECT_ROOT / "results" / "figures"
    manuscript: Path = PROJECT_ROOT / "manuscript"
    logs: Path = PROJECT_ROOT / "logs"

    def __post_init__(self) -> None:
        for p in (
            self.data_raw,
            self.data_interim,
            self.data_processed,
            self.results_tables,
            self.results_figures,
            self.manuscript,
            self.logs,
        ):
            p.mkdir(parents=True, exist_ok=True)


PATHS = Paths()


@dataclass(frozen=True)
class CohortSpec:
    accession: str
    role: str
    platform: str
    source: str
    n_expected: int
    notes: str = ""


COHORTS: Dict[str, CohortSpec] = {
    "nassiri_2021": CohortSpec(
        accession="GSE180061",
        role="discovery",
        platform="Illumina Infinium HumanMethylation450 + RNA-seq + WES",
        source="Nassiri et al. 2021 Nature (GEO methylation idat; cBioPortal mng_utoronto_2021 for mRNA)",
        n_expected=185,
        notes="EGAS00001004982 (WES/mRNA/snRNA) is controlled-access; cBioPortal entry used as default entry point.",
    ),
    "bi_lab_2023": CohortSpec(
        accession="GSE212666",
        role="validation_unlabeled_scoped_out",
        platform="Illumina EPIC methylation + RNA-seq (Kallisto) + proteomics + scRNA-seq",
        source="Choudhury/Bi lab Cancer Cell 2023 three-group classification (Merlin-intact / immune-enriched / hypermitotic)",
        n_expected=302,
        notes="SCOPED OUT 2026-07-31 (see GSE212666_SCOPE_OUT rationale). GSM-level SOFT: 0% phenotype metadata "
              "(302/302 = tissue: Meningioma only; no subgroup/grade/NF2/age/sex/recurrence). Data format: 302 per-sample "
              "Output_QM Kallisto tarballs (abundance.h5/tsv with ENST transcript IDs + TPM/counts). Engineering lift "
              "is moderate (tximport aggregation tractable) but produces ZERO testable info without ground-truth labels. "
              "Retained for provenance; re-activate if Choudhury/Bi/Raleigh lab shares per-sample subgroup labels via "
              "correspondence. Labeled training data for bridging classifier comes from sibling accession GSE183653 "
              "(n=185 processed TPM from 2022 Nature Genetics discovery cohort).",
    ),
}

REPLICATION_CANDIDATE_GEO_TERMS: List[str] = [
    "meningioma[tiab] AND expression profiling by array[pt]",
    "meningioma[tiab] AND RNA-seq[tiab]",
]

REPLICATION_CANDIDATE_STUDIES = [
    "GSE136661",
    "GSE77259",
    "GSE115966",
    "GSE94474",
    "GSE16581",
    "GSE32062",
    "GSE40516",
    "GSE43290",
    "GSE48841",
    "GSE55018",
]

CBIO_MENINGIOMA_STUDY_IDS = [
    "mng_utoronto_2021",
    "mng_miami_2022",
]

PLOT_STYLE: Dict[str, object] = {
    "context": "talk",
    "style": "whitegrid",
    "palette": "viridis",
    "font_scale": 1.1,
    "figure_dpi": 300,
    "figure_format": ["png", "pdf", "svg"],
    "cmap_heatmap": "RdBu_r",
    "cmap_categorical": "Set2",
}

SIGNIFICANCE_THRESHOLD: float = 0.05
FDR_ALPHA: float = 0.10

GRADE_ORDER: List[str] = ["I", "II", "III"]
NF2_ORDER: List[str] = ["Intact", "Mutant/Loss"]
NF2_INTERNAL_COL: str = "nf2_status"
NF2_EXPORT_COL: str = "nf2_cna_loss_proxy"
NF2_CATEGORIES: List[str] = ["Intact", "Mutant/Loss"]

NASSIRI_GROUP_ORDER: List[str] = [
    "Immunogenic",
    "MG2",
    "hypermetabolic",
    "proliferative",
]
NASSIRI_GROUP_PAPER_ALIASES: Dict[str, str] = {
    "Immunogenic": (
        "MG1 Immunogenic (Nassiri Nature 2021 Fig 1; immune-infiltrated, favorable). "
        "LITERATURE CLAIM (independent validation, replicated in our data): 'MG1 tumors invariably "
        "featured NF2 biallelic inactivation from co-occurring NF2 point mutation AND "
        "22q (NF2 copy-number loss' (Oxford Neuro-Oncology / Frontiers Oncol 2022 reviews). 26789377; "
        "Frontiers Oncol 12:962702; npj Precis Oncol; Lynes et al 2022). OUR DATA "
        "REPLICATES THIS FALSIFIABLE CLAIM: in cBioPortal mng_utoronto_2021 GISTIC CNA "
        "nf2_cna_loss_proxy shows 0/17 (0.0%) Immunogenic samples are NF2-Intact, i.e. "
        "17/17 carry evidence of 22q/NF2 deletion, directionally consistent with 'invariable "
        "biallelic loss' (point-mutation data from cBioPortal mutation endpoint is broken; the "
        "CNA arm captures the 22q arm-deletion half of biallelic inactivation; the 0% intact rate "
        "serves as a necessary, not sufficient, validation of the derivation method's specificity "
        "for MG1-level biallelic calls)."
    ),
    "MG2": (
        "MG2 Benign NF2-wildtype (Nassiri Nature 2021 Fig 1; low recurrence, clinically benign). "
        "LITERATURE CLAIM (5-source consensus): explicitly NF2-wildtype with benign clinical course, "
        "low recurrence. OUR DATA REPLICATES: 27/32 (84.4%) nf2_cna_loss_proxy = Intact (CNA-diploid "
        "at NF2 locus), consistent with wildtype. 5/32 residual CNA-loss calls in MG2 may represent "
        "isolated 22q loss without the concurrent NF2 point mutation required for true biallelic "
        "inactivation — our proxy cannot distinguish isolated shallow 22q loss from biallelic MG1-type "
        "inactivation (mutation endpoint broken in this cBioPortal instance; Limitations section of "
        "manuscript states this explicitly)."
    ),
    "hypermetabolic": (
        "MG3 Hypermetabolic (Nassiri Nature 2021 Fig 1; metabolic pathway upregulation, "
        "intermediate prognosis. Maps to Choudhury/Raleigh 2022 Nature Genetics 'Hypermitotic' "
        "subgroup combined with MG4 per the 2023 Oxford Neuro-Oncology and Frontiers reviews "
        "(see NASSIRI_BI_MATCH_HYPOTHESIS). The Hypermitotic bin pools Nassiri MG3+MG4 per the "
        "2023 Neuro-Oncology subgroup-refinement paper (Choudhury et al.): FOXM1− Hypermetabolic = MG3, "
        "FOXM1+ Proliferative = MG4."
    ),
    "proliferative": (
        "MG4 Proliferative (Nassiri Nature 2021 Fig 1; high proliferation / Ki-67, least favorable "
        "prognosis. Maps to Choudhury/Raleigh 2022 Nature Genetics 'Hypermitotic' subgroup combined "
        "with MG3 per literature crosswalk (see NASSIRI_BI_MATCH_HYPOTHESIS)."
    ),
}

BI_CLASSIFIER_METHODS_NOTE: str = (
    "METHODS (for Bi mRNA cohort classification, reproducibility pre-registration):\n"
    "  The 3-tier Bi-lab subgroup labels {Merlin-intact, Immune-enriched, Hypermitotic} were\n"
    "  originally derived in Choudhury et al. 2022 (Nature Genetics 54:649-659, PMID 35534562)\n"
    "  from Illumina EPIC DNA methylation profiling (N=565 total; N=200 discovery + N=365\n"
    "  validation) using an unsupervised consensus clustering + nearest-centroid classifier\n"
    "  trained on methylation beta-values at the most variably methylated CpG probes.\n"
    "  The associated GitHub/Zenodo release (abrarc/meningioma-svm v1, doi:10.5281/zenodo.6353877)\n"
    "  contains the R SVM training script skeleton (svm_classifier_ng.R) but does NOT include\n"
    "  the per-sample subgroups.csv or beta_values.csv needed to reproduce the classifier;\n"
    "  those files are referenced by the script but are not open-access via GitHub.\n"
    "\n"
    "  The original Nature Genetics 2022 study also performed paired bulk RNA-seq on the\n"
    "  N=200-sample DISCOVERY COHORT only (the validation cohort of N=365 was methylation-only\n"
    "  in the 2022 publication). A subsequent 2023 subgroup-refinement paper from the same\n"
    "  group (Choudhury et al., Neuro-Oncology 25:520-530) added RNA-seq for N=302 validation\n"
    "  samples, giving N=502 total with paired RNA-seq across both manuscripts.\n"
    "  FOR THIS REPO'S DATA AVAILABILITY REALITY CHECK (verified by direct GEO FTP audit, 2026):\n"
    "    - GSE183653 = original 2022 discovery RNA-seq: N=185 samples with a PROCESSED TPM\n"
    "      expression matrix available via GEO supplementary (GSE183653_meningioma_tpm.csv.gz).\n"
    "      185 < 200 because 15 discovery samples failed RNA-seq QC. This is the training set.\n"
    "    - GSE212666 = 2023 validation RNA-seq: N=302 samples on GPL24676. GSM-level\n"
    "      characteristics_ch1 are EMPTY except for a single 'tissue: Meningioma' key on all\n"
    "      302 records (0/302 = 0% phenotype metadata coverage, confirmed by direct SOFT parse).\n"
    "      Only a RAW.tar with 302 per-sample Output_QM quant tarballs is on GEO; no series-level\n"
    "      processed expression matrix exists and no GSM-level subgroup labels are present.\n"
    "  Neither study deposited mRNA subgroup calls on GEO. The 2022 paper did NOT publish a\n"
    "  standalone mRNA-only 3-group classifier (no 'molecular-grade mRNA centroid'); the 3-group\n"
    "  labels were propagated from the methylation classifier onto the paired RNA-seq samples via\n"
    "  the paper's integrated multi-omic model.\n"
    "\n"
    "  For THIS study's Aim 1 concordance analysis (H1: Nassiri 4-group <-> Bi 3-group mapping),\n"
    "  the classification strategy is as follows, with A and B explicitly recorded as NON-VIABLE\n"
    "  or SCOPED-OUT a priori so downstream readers are not surprised by their absence:\n"
    "\n"
    "    (A) [PRE-REGISTERED AS NON-VIABLE, CONFIRMED 0/302 BY GEO SOFT AUDIT]\n"
    "        Plan A: use ground-truth methylation-based Bi subgroup labels directly from GEO GSM\n"
    "        characteristics_ch1. OUTCOME: GSM-level characteristics on GSE212666 contain only\n"
    "        'tissue = Meningioma' (N=302/302). NO subgroup, grade, age, NF2-status, recurrence,\n"
    "        or sex metadata is present. Plan A therefore provides 0% coverage and is abandoned.\n"
    "        This is a known, documented contingency, not a post-hoc failure mode.\n"
    "\n"
    "    (B) [SCOPED OUT 2026-07-31, SEE GSE212666_SCOPE_OUT RATIONALE IN THIS FILE]\n"
    "        Plan B (originally proposed): apply trained bridging classifier to GSE212666\n"
    "        N=302 to derive predicted bi_group labels. OUTCOME: Scoped out on cost-benefit\n"
    "        grounds. GSE212666 has 0% phenotype metadata (confirmed above), so predicted\n"
    "        labels cannot be validated against any ground truth. Engineering lift to convert\n"
    "        302 Output_QM Kallisto tarballs (abundance.h5/tsv, confirmed per-tarball content\n"
    "        audit) into a gene-level expression matrix is moderate (tximport-style ENST->gene\n"
    "        aggregation, 302 tarball extractions) but produces ZERO testable new information\n"
    "        for H1 (Nassiri<->Bi concordance only needs Nassiri cohort labeled by bridging\n"
    "        classifier), H2 (subgroup/target enrichment in Bi needs labeled Bi set from\n"
    "        GSE183653 only), or H3 (NF2 association needs labeled Bi set from GSE183653 only).\n"
    "        The n=302 unlabeled set would add only weak 'transcriptomic cluster-consistency'\n"
    "        evidence that cannot confirm label prediction accuracy. Cost-benefit ratio is\n"
    "        unambiguously negative. GSE212666 retained in raw data folder and listed as\n"
    "        future direction if Choudhury/Bi lab provides per-sample subgroup labels via\n"
    "        correspondence.\n"
    "\n"
    "    (C) [PRIMARY, ONLY VIABLE PATH — mRNA-to-methylation-label bridging classifier]\n"
    "        Train a 3-class nearest-centroid (or robust regularized alternative, e.g.\n"
    "        Ridge-penalized multinomial) RNA expression BRIDGING CLASSIFIER on the GSE183653\n"
    "        N=185 discovery samples, using top-N most variably expressed protein-coding genes\n"
    "        that overlap between the Nassiri expression matrix (mRNA-seq from cBioPortal\n"
    "        mng_utoronto_2021) and the Bi TPM matrix (GSE183653_meningioma_tpm.csv.gz from\n"
    "        GEO series supplementary). The per-sample training labels for the bridging\n"
    "        classifier are obtained from the Choudhury 2022 Nature Genetics paper's\n"
    "        Supplementary Tables (Supp Table 1 or Supp Table 2 carry per-sample clinical +\n"
    "        DNA-methylation-subgroup metadata for the N=200 discovery cohort, of which\n"
    "        N=185 passed RNA-seq QC and are present in the TPM matrix). If the printed\n"
    "        Supp Tables are figure-rendered and not machine-readable, labels are extracted\n"
    "        by manual double-entry with cross-validation against the paper's reported\n"
    "        subgroup fractions for the N=200 discovery set (Merlin-intact 34%,\n"
    "        Immune-enriched 38%, Hypermitotic 28%, expected N=68/76/56 in N=200,\n"
    "        approximately N=63/70/52 in N=185 after QC attrition). Expected label fractions\n"
    "        are enforced as a sanity-check assertion; any label-extraction CSV that deviates\n"
    "        by >5 absolute count per bin from the paper's reported fractions is REJECTED\n"
    "        with a loud pipeline failure (exit code 3) to force manual re-verification.\n"
    "        The bridging classifier's leave-one-out (LOOCV) OR 5-fold-repeated x 10-fold\n"
    "        cross-validation accuracy on the GSE183653 training set is reported ALONGSIDE\n"
    "        the Aim 1 concordance results so readers can judge label-recovery fidelity\n"
    "        independently of the cross-study Nassiri/Bi mapping claim.\n"
    "        CRITICAL: This bridging classifier is an INTERNAL ADAPTATION STEP, NOT a\n"
    "        claim that the 2022/2023 Bi-lab papers validated a standalone mRNA-input\n"
    "        3-group classifier. This caveat is flagged explicitly in the manuscript\n"
    "        Methods under 'Cross-study subgroup concordance (Aim 1)'.\n"
    "\n"
    "  Post-training deployment: the validated bridging classifier is applied ONLY to the\n"
    "  Nassiri mng_utoronto_2021 discovery mRNA-seq cohort (N=121, from cBioPortal) to\n"
    "  derive predicted bi_group labels per sample, which are then used as the input to\n"
    "  Aim 1's concordance test (Nassiri 4-group <-> predicted Bi 3-group cross-tab,\n"
    "  evaluated with Cohen's kappa + permutation-based pairing-specificity check in\n"
    "  concordance.py::test_pairing_specificity()).\n"
    "\n"
    "  Open-data caveat: if Choudhury et al.'s 2022 per-sample Supp Table cannot be\n"
    "  machine-extracted and the authors cannot provide a machine-readable subgroups.csv\n"
    "  on request, the bridging classifier's training labels are the best-effort manual\n"
    "  double-entry extraction, and this limitation is stated explicitly in the Aim 1\n"
    "  results paragraph BEFORE any concordance metric is reported (never buried in\n"
    "  Limitations). PIPELINE FAILS LOUDLY (exit 3) if label sidecar is missing."
)

GSE212666_SCOPE_OUT: str = (
    "GSE212666 (N=302 Bi-lab 2023 validation RNA-seq, GPL24676) — SCOPED-OUT RATIONALE:\n"
    "  Audited 2026-07-31. Key facts:\n"
    "    (1) 0% phenotype metadata coverage in GEO GSM-level characteristics_ch1:\n"
    "        302/302 GSMs have ONLY 'tissue: Meningioma'; NO subgroup, grade, age,\n"
    "        NF2-status, recurrence, or sex metadata. Confirmed by direct SOFT parse.\n"
    "    (2) Data format is tractable but not immediately usable: GEO supplies only\n"
    "        GSE212666_RAW.tar + 302 per-sample Output_QM*.tar.gz files. Each tarball\n"
    "        contains standard Kallisto output (abundance.h5, abundance.tsv with ENST\n"
    "        transcript IDs + TPM/counts, run_info.json). Already quantified (NOT raw\n"
    "        sequencing reads), so alignment/index-building is unnecessary. Engineering\n"
    "        lift = moderate: 302 tarball extractions + tximport-style ENST->gene map\n"
    "        (Ensembl BioMart tx2gene) + matrix construction. No blockers; all tractable.\n"
    "    (3) Zero testable information gain, given (1). H1 (Nassiri<->Bi concordance)\n"
    "        needs only Nassiri N=121 labeled via bridging classifier, which does not\n"
    "        involve GSE212666 at all. H2 (Bi subgroup target enrichment) and H3\n"
    "        (Bi subgroup NF2 association) need labeled Bi samples, which GSE183653\n"
    "        N=185 provides after label extraction from 2022 paper Supp Tables.\n"
    "        GSE212666 N=302 has NO labels, so predicted labels from bridging classifier\n"
    "        can only be checked for 3-cluster transcriptomic consistency (silhouette\n"
    "        score, PCA separation) — NOT against ground-truth subgroup assignments.\n"
    "        Cluster consistency is a necessary but EXTREMELY non-specific check that\n"
    "        cannot validate that the CORRECT samples are assigned to the CORRECT groups.\n"
    "    (4) Cost-benefit is unambiguously negative: moderate engineering work yields\n"
    "        evidence too weak to move any hypothesis needle. Scoped out to keep the\n"
    "        pipeline focused on H1/H2/H3 tests that have actual falsifiability.\n"
    "  Retention policy: raw SOFT file retained in data/raw/GSE212666/ for provenance\n"
    "  and future re-evaluation. GSE212666 accession is NOT removed from COHORTS dict\n"
    "  (preserves audit trail) but its role field is treated as 'validation_unlabeled'\n"
    "  by all pipeline scripts, which skip it for hypothesis-testing analyses.\n"
    "  Re-entry criterion: if Choudhury/Bi/Raleigh lab provides per-sample subgroup\n"
    "  labels for the 302 GSMs via correspondence, re-scope in. Re-activation checklist:\n"
    "  (a) write tximport aggregation step, (b) merge labels into metadata,\n"
    "  (c) add GSE212666 as independent labeled replication for H2/H3 only\n"
    "  (still not useful for H1 label-validation, since bridging classifier is trained\n"
    "  on GSE183653 from the SAME lab — would be same-lab validation not independent)."
)

BRIDGING_CLASSIFIER_LABELS_REQUIREMENT: str = (
    "BRIDGING CLASSIFIER TRAINING LABELS — PREREQUISITE FOR Aim 1, FAIL LOUDLY IF MISSING:\n"
    "  File location: data/raw/GSE183653/GSE183653_sample_bi_labels.csv\n"
    "  Required columns (CSV with header, comma-delimited):\n"
    "    - sample_id: str, matches TPM matrix column names EXACTLY (e.g. 'M1', 'M10',\n"
    "      'M100', case-sensitive, no GSM prefix). Must be a 1:1 match to the non-gene\n"
    "      column headers in GSE183653_meningioma_tpm.csv.gz (expected N=185).\n"
    "    - bi_group: str, one of {Merlin-intact, Immune-enriched, Hypermitotic}\n"
    "      exactly matching BI_GROUP_ORDER capitalization and hyphenation.\n"
    "  Optional additional columns (retained but not required): age, grade, sex,\n"
    "    nf2_status, recurrence — any such covariates will be used to confirm GSM<->M##\n"
    "    covariate match against GEO GSM characteristics if probabilistic matching is\n"
    "    needed to fill partial label gaps.\n"
    "  Sanity-check assertions (pipeline REJECTS CSV with exit 3 if violated):\n"
    "    (A) set(sample_id) must equal set(TPM matrix columns) exactly — no missing,\n"
    "        no extra samples, no typos.\n"
    "    (B) bi_group values must be subset of BI_GROUP_ORDER — no typos, no 4th group.\n"
    "    (C) Per-bin count sanity: |n(Merlin-intact)  - 63| <= 8,\n"
    "                               |n(Immune-enriched) - 70| <= 8,\n"
    "                               |n(Hypermitotic)    - 52| <= 8\n"
    "        Derived from paper's N=200 discovery fractions (34%/38%/28%) × QC attrition\n"
    "        (N=185/200 = 92.5%). Tolerance of ±8 allows for ~4% QC-differential loss\n"
    "        per group; tighter than any plausible manual-extraction error without\n"
    "        being so tight as to reject a correct extraction. If labels are rejected,\n"
    "        the error message prints actual counts vs expected and suggests verifying\n"
    "        group labels against Fig 1 of Choudhury 2022 Nature Genetics.\n"
    "  Label source (priority order, documented verbatim in Methods):\n"
    "    1. Choudhury et al. 2022 Nature Genetics Supplementary Table 1 or\n"
    "       Supplementary Table 2 (per-sample clinical metadata with DNA methylation\n"
    "       subgroup assignments), machine-extracted if PDF-native vector table,\n"
    "       or manual double-entry by two independent readers if figure-rendered.\n"
    "    2. If #1 fails: Abrarc/meningioma-svm GitHub repo (Zenodo 10.5281/zenodo.6353877)\n"
    "       if the subgroups.csv file referenced in svm_classifier_ng.R becomes available\n"
    "       via author correspondence.\n"
    "    3. If #1 and #2 both fail: Aim 1 (and any Bi-lab-based H2/H3 analyses) are\n"
    "       declared NOT TESTABLE and documented as such in Aim 0 status report.\n"
    "       Pipeline continues with Nassiri-only H2/H3 analyses; Bi-lab hypotheses\n"
    "       are explicitly listed as deferred to future work with author data-sharing.\n"
    "  Provenance: the provenance field in metadata_PROVENANCE.csv for any processed\n"
    "  bi_lab_2023-derived cohort MUST cite the exact label source used (priority 1/2/3),\n"
    "  including date of extraction, name of extractor, and (for manual entry) confirmation\n"
    "  of inter-reader agreement (Cohen's kappa > 0.95 required for manual double-entry)."
)

BI_GROUP_ORDER: List[str] = [
    "Merlin-intact",
    "Immune-enriched",
    "Hypermitotic",
]

NASSIRI_BI_MATCH_HYPOTHESIS: Dict[str, str] = {
    "Immunogenic": "Immune-enriched",
    "MG2": "Merlin-intact",
    "hypermetabolic": "Hypermitotic",
    "proliferative": "Hypermitotic",
}
NASSIRI_BI_MATCH_RATIONALE: str = (
    "PRECISE DIRECTIONAL 1:1(:1) PRE-REGISTERED HYPOTHESIS — literature crosswalk confirmed "
    "across 5 independent sources (2023 Oxford Neuro-Oncology Lynes et al. review, 2 Oxford Neuro-Oncology "
    "primary papers, npj Precision Oncology, Frontiers Oncol 12:962702). Exact directional pairing: "
    "(1) Bi Merlin-intact ←→ Nassiri MG2 (benign NF2-wildtype) — both explicitly named as the "
    "clinically benign, NF2-diploid/wildtype subgroup; (2) Bi Immune-enriched ←→ Nassiri MG1 "
    "(Immunogenic) — both named as the immune-infiltrated subgroup with INVARIABLE biallelic NF2 "
    "inactivation (22q loss + concurrent NF2 point mutation); (3) Bi Hypermitotic ←→ Nassiri "
    "{MG3 Hypermetabolic ∪ MG4 Proliferative} — the two higher-grade, higher-proliferation "
    "Nassiri subgroups pooled into the single Bi Hypermitotic malignant bucket. "
    "Pooling of MG3∪MG4 into Bi Hypermitotic is registered in config BEFORE any classifier "
    "application. H1 is falsifiable via test_pairing_specificity() in concordance.py: "
    "the true mapping is scored against the full enumeration of 11 alternative pairings that "
    "also route a merged Nassiri pair to Bi Hypermitotic (C(4,2)=6 choices of which 2 Nassiri "
    "groups merge × 2! = 2 permutations of the remaining 2 singleton assignments = 12 total "
    "pairings, 1 correct + 11 wrong). The true 1:1:1 pairing must rank #1 in Cohen's κ and/or "
    "raw agreement among all 12; this is a stronger check than the strict 2-pairing Scenario A "
    "(which only swaps the two singletons for 1 wrong alternative) because it additionally "
    "verifies that the specific MG3∪MG4 merge is empirically better than any of the 5 "
    "alternative merge choices (MG1∪MG2, MG1∪MG3, MG1∪MG4, MG2∪MG3, MG2∪MG4) — i.e., the "
    "literature's choice of which 2 groups pool into Hypermitotic is recoverable from data alone. "
    "The wrong singleton-swap (MG1↔Merlin-intact, MG2↔Immune-enriched) is the most adversarial "
    "single placebo and always included; the 10 additional alternatives span alternative merge "
    "choices to rule out trivial 'any 2-group merge works' explanations."
)

H2_MIN_CELL_SIZE: int = 10
H2_MIN_CELL_SIZE_RATIONALE: str = (
    "PRE-REGISTRATION DECISION: minimum 10 observations per subgroup cell for the "
    "H2 Kruskal-Wallis omnibus test. Any level of nassiri_group with N < H2_MIN_CELL_SIZE "
    "is pooled into a single 'pooled_small_groups' cell; if the pooled cell is also <10, "
    "the subgroup is excluded from the test (degrees of freedom adjusted accordingly). "
    "Rationale: K-W rank-sum distribution asymptotics are unreliable below 5-8 per cell; "
    "N=10 is a conservative floor that preserves the K-W type-I error rate at alpha=0.05 "
    "in small-cell Monte-Carlo simulations. Current Nassiri discovery Ns: Immunogenic=17, "
    "MG2=32, hypermetabolic=43, proliferative=29 — all 4 groups exceed the floor, so no "
    "pooling is triggered by current data. Rule retained to constrain future sub-analyses."
)

H5_STATUS: str = "not_testable"
H5_STATUS_RATIONALE: str = (
    "PRE-REGISTRATION DECLARATION (Aim 0). The originally specified H5 endpoint "
    "(Cox PH / Kaplan-Meier time-to-recurrence analysis) is NOT TESTABLE in the "
    "public cBioPortal mng_utoronto_2021 release. recurrence_event is present at 100% "
    "non-null (64 events / 57 censored) but recurrence_months / time-to-event / follow-up "
    "duration is not present in the clinical attributes (verified 2026-07-31 against "
    "clinical_attribute_audit.csv for mng_utoronto_2021, which contains only 7 attributes "
    "total, none a duration-type field). H5 is therefore declared not testable BEFORE any "
    "analysis of the event variable itself, to eliminate any subsequent appearance of "
    "result-driven endpoint selection. If the paper's Supplementary survival table is "
    "obtained as a sidecar CSV with per-sample months-to-recurrence, H5 may be re-enabled "
    "as a pre-registered amendment citing that table as the source."
)

NF2_STATUS_PLAN: Dict[str, str] = {
    "internal_colname": "nf2_status (legacy/internal; never used in manuscript tables or exported CSV headers).",
    "exported_colname": "nf2_cna_loss_proxy (canonical name for all manuscript-facing CSVs and tables).",
    "limitation_statement": "COPY-NUMBER-LOSS PROXY ONLY — not a true biallelic NF2 inactivation call. "
        "GISTIC ≤ -1 captures only one hit (het loss or hom del at 22q/NF2); true functional NF2 inactivation "
        "in meningioma specifically requires BOTH 22q copy-number loss AND a concurrent NF2 point mutation "
        "(MG1-type biallelic inactivation per Nassiri 2021). The mutation endpoint returns 404 in this "
        "cBioPortal instance, so the second hit is invisible. Interpretation: this proxy correctly captures "
        "MG1-type biallelic cases (which uniformly carry 22q loss) but may overcall isolated shallow 22q loss "
        "without a point mutation as 'Mutant/Loss' even though Merlin may remain functionally expressed. "
        "This limitation is stated explicitly in the manuscript Limitations section. Results using this "
        "variable are labeled as directional-association evidence, not definitive inactivation calls.",
    "primary_source": "cBioPortal mng_utoronto_2021 GISTIC copy-number profile (mng_utoronto_2021_gistic, NF2 entrez 4771).",
    "derivation": "Discrete GISTIC values ≤ -1 (het loss or hom del) → Mutant/Loss; 0 → Intact; ≥+1 → Intact. "
        "Per-sample severity tie-break: if multiple CNA entries per sample, most severe call wins (Mutant/Loss > Intact > NaN).",
    "coverage_expected": "121/121 samples in mng_utoronto_2021 (verified via GET molecular-profiles/gistic/molecular-data).",
    "mutation_overlay": "Mutations endpoint 404 in this cBioPortal instance. If a MAF supplementary file is obtained from the paper, point mutations will be OR-merged with CNA calls (any mutation OR any CNA loss → Mutant/Loss) under the same conflict-conservative rule used elsewhere — clinical/mutation disagreement → NaN. IF mutations are ever merged in, the variable suffix will be updated to nf2_biallelic_combined rather than silently reusing the _cna_loss_proxy label.",
    "fallback_if_cna_missing": "If GISTIC profile is ever removed: (a) try the log2CNA profile for NF2 log ratio thresholds; (b) fall back to GEO GSE180061 with a paper-supplementary NF2 sidecar; (c) if none work, H3 is declared not testable rather than imputed.",
}

MULTIPLE_TESTING_FAMILY_DESCRIPTION: str = (
    "BH-FDR applied across the full family of "
    "program × subgroup × cohort tests."
)

NORMALIZATION_DECISION: str = (
    "EXPRESSION NORMALIZATION — PRE-REGISTERED, EXPLICIT DECISION (Phase 0, 2026-08-01):\n"
    "  Both cohorts enter the bridging classifier on the SAME SCALE: raw (non-z-scored)\n"
    "  TPM (transcripts-per-million) with log2(TPM+1) transform, followed by per-gene\n"
    "  z-scoring on the GSE183653 TRAINING SET ONLY. The StandardScaler inside the\n"
    "  sklearn Pipeline then re-standarizes internally; this double-z is statistically\n"
    "  harmless for nearest-centroid/Ridge and preserves relative rank ordering.\n"
    "\n"
    "  RATIONALE (why TPM, not cBioPortal z-scores):\n"
    "    (A) GSE183653 (Bi-lab training set) is only available as processed TPM from\n"
    "        GEO series supplementary — there is NO corresponding Bi-lab z-score profile.\n"
    "        Training on TPM → deploying on cBioPortal z-scores would mix scales and\n"
    "        systematically shift centroids, introducing avoidable distribution shift.\n"
    "    (B) cBioPortal's default mRNA profile (mng_utoronto_2021_mrna_seq_mrna) returns\n"
    "        z-scores vs a reference diploid pool. This study explicitly tests for\n"
    "        subgroup-LEVEL mean differences between 3 groups; if z-scoring already\n"
    "        centers every gene at 0, between-subgroup signal is attenuated by\n"
    "        construction — a real but silent bias that we eliminate by using TPM.\n"
    "    (C) Both TPM matrices use HUGO gene symbols (Entrez→Hugo resolved via\n"
    "        cBioPortal POST /genes/fetch batch API, ≤30 calls). Gene-symbol overlap\n"
    "        is maximized on the same coordinate system.\n"
    "    (D) If the TPM molecular profile (mng_utoronto_2021_mrna_seq_tpm) is not\n"
    "        available in this cBioPortal instance, the script FALLS BACK to the\n"
    "        z-score profile and prints a GIANT warning banner rather than failing\n"
    "        silently — this caveat is then flagged in Limitations and Methods.\n"
    "\n"
    "  Implementation: scripts/fetch_cbio_mng_utoronto_2021.py selects the molecular\n"
    "  profile whose molecularProfileId matches '*mrna_seq_tpm*' first; if absent,\n"
    "  falls back to '*mrna*' with a banner. The selection is logged and recorded\n"
    "  in metadata_PROVENANCE.csv under 'expression_profile_id' so the exact choice\n"
    "  is auditable.\n"
)
