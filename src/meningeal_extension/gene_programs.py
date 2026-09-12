from __future__ import annotations

from typing import Dict, List

GENE_PROGRAMS: Dict[str, List[str]] = {
    "HDAC_Panobinostat_Romidepsin": [
        "HDAC1",
        "HDAC2",
    ],
    "Proteasome_Carfilzomib": [
        "PSMB5",
        "PSMB1",
        "PSMB2",
        "PSMA1",
        "PSMA2",
        "PSMA3",
        "PSMA4",
        "PSMA5",
        "PSMA6",
        "PSMA7",
        "PSMC1",
        "PSMC2",
        "PSMC3",
        "PSMC4",
        "PSMC5",
        "PSMC6",
        "PSMD1",
        "PSMD2",
        "PSMD11",
        "PSMD14",
    ],
    "Tubulin_Ixabepilone": [
        "TUBB",
        "TUBB3",
        "TUBB2A",
        "TUBB2B",
        "TUBB4A",
        "TUBB4B",
        "TUBB6",
        "TUBA1A",
        "TUBA1B",
        "TUBA1C",
        "TUBA3C",
        "TUBA3D",
        "TUBA4A",
    ],
    "Translation_Omacetaxine": [
        "EEF2",
        "EEF1A1",
        "EEF1A2",
        "EEF1B2",
        "EEF1D",
        "EEF1G",
        "EIF4A1",
        "EIF4A2",
        "EIF4A3",
        "EIF4E",
        "EIF4G1",
        "EIF4G2",
        "RPL3",
        "RPL5",
        "RPL10A",
        "RPL13",
        "RPL13A",
        "RPL19",
        "RPL23",
        "RPL27A",
        "RPL28",
        "RPL29",
        "RPL30",
        "RPL31",
        "RPL34",
        "RPL35",
        "RPL36A",
        "RPL37A",
        "RPS3",
        "RPS5",
        "RPS6",
        "RPS7",
        "RPS8",
        "RPS10",
        "RPS11",
        "RPS13",
        "RPS15A",
        "RPS16",
        "RPS18",
        "RPS19",
        "RPS20",
        "RPS23",
        "RPS24",
        "RPS25",
        "RPS27",
        "RPS27A",
        "RPS29",
    ],
}

HDAC_RESISTANCE_PROGRAM: Dict[str, List[str]] = {
    "HDAC8_TGFb_EMT_Resistance": [
        "HDAC8",
        "TGFB1",
        "TGFB2",
        "TGFBR1",
        "TGFBR2",
        "SMAD2",
        "SMAD3",
        "SMAD4",
        "SNAI1",
        "SNAI2",
        "TWIST1",
        "TWIST2",
        "ZEB1",
        "ZEB2",
        "FN1",
        "VIM",
        "CDH2",
        "MMP2",
        "MMP9",
        "COL1A1",
        "COL3A1",
    ],
}

ALL_PROGRAMS: Dict[str, List[str]] = {**GENE_PROGRAMS, **HDAC_RESISTANCE_PROGRAM}

DRUG_PROGRAM_MAP: Dict[str, str] = {
    "Panobinostat": "HDAC_Panobinostat_Romidepsin",
    "Romidepsin": "HDAC_Panobinostat_Romidepsin",
    "Carfilzomib": "Proteasome_Carfilzomib",
    "Ixabepilone": "Tubulin_Ixabepilone",
    "Omacetaxine": "Translation_Omacetaxine",
}

PROGRAM_DESCRIPTIONS: Dict[str, str] = {
    "HDAC_Panobinostat_Romidepsin": (
        "Class-I HDACs HDAC1/2 — direct shared molecular target of panobinostat "
        "(pan-HDAC inhibitor, Jungwirth Sci Transl Med 2026) and romidepsin "
        "(class-I selective HDAC inhibitor, Jungwirth Clin Cancer Res 2023)."
    ),
    "Proteasome_Carfilzomib": (
        "26S proteasome core (PSM-family alpha, beta, regulatory subunits) — "
        "target of carfilzomib, a selective irreversible proteasome inhibitor "
        "(IC50 0.12 nmol/L in CCR 2023 screen)."
    ),
    "Tubulin_Ixabepilone": (
        "Alpha/beta tubulin gene family including TUBB3 (class-III beta-tubulin, "
        "established epothilone/taxane resistance marker) — target of ixabepilone "
        "(epothilone-B microtubule stabilizer)."
    ),
    "Translation_Omacetaxine": (
        "Eukaryotic elongation factors (EEF2, EEF1-family) and 40S/60S ribosomal "
        "protein genes — targets of omacetaxine/homoharringtonine (translation "
        "elongation inhibitor)."
    ),
    "HDAC8_TGFb_EMT_Resistance": (
        "HDAC8 + TGFβ/SMAD axis + canonical EMT transcription factors and "
        "mesenchymal markers — the panobinostat-resistance axis identified in "
        "Jungwirth Sci Transl Med 2026 (HDAC8 depletion restores sensitivity)."
    ),
}

PROGRAM_TO_DRUGS: Dict[str, List[str]] = {
    "HDAC_Panobinostat_Romidepsin": ["Panobinostat (StM 2026)", "Romidepsin (CCR 2023)"],
    "Proteasome_Carfilzomib": ["Carfilzomib (CCR 2023)"],
    "Tubulin_Ixabepilone": ["Ixabepilone (CCR 2023)"],
    "Translation_Omacetaxine": ["Omacetaxine (CCR 2023)"],
    "HDAC8_TGFb_EMT_Resistance": ["Panobinostat resistance (StM 2026 mechanism)"],
}
