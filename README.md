# Meningeal Extension — convergent classifier stratification of meningioma drug programs

Locked computational pipeline extending Jungwirth *Clin Cancer Res* 2023 and Jungwirth *Sci Transl Med* 2026 Heidelberg meningioma drug screens with two independent molecular classification systems (Nassiri 2021 Nature; Bi-lab 2023 Cancer Cell), external replication, and a pre-specified exploratory survival module.

## Quick start

```bash
pip install -r requirements.txt
PYTHONPATH=src python scripts/99_run_all.py
```

This runs Aims 0–4, writes 29 output tables to `results/tables/`, 21 figure assets (PNG/PDF/SVG) to `results/figures/`, and produces the compliance manifest at `manuscript/artifact_manifest.json`.

## Tests

```bash
PYTHONPATH=src pytest tests/test_all_aims.py -v
```

## Manuscript, 1-page summary, pre-submission checklist

All under `manuscript/`:
- `paper_manuscript.md` — full manuscript
- `one_page_summary_for_warta.md` — Heidelberg-facing 1-page summary (attach Figures 2A + 5)
- `pre_submission_checklist.md` — 43-point numbered checklist across 3 blocks

## Protocol lock

Plan, hypotheses, and frozen gene-program lists are detailed in `PLAN.md` and enforced literally in `src/meningeal_extension/`.
