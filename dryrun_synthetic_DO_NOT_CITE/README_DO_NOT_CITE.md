###############################################################################
#                          DO NOT CITE — SYNTHETIC DRY-RUN                    #
###############################################################################
#
# EVERY FILE IN THIS TREE WAS GENERATED FROM SYNTHETIC SEEDED DATA, NOT REAL
# GEO/cBioPortal COHORTS.
#
# Source: `build_synthetic_cohort()` in datasets.py
# - Seed: 20260317
# - Sample IDs: nassiri_synthetic_XXXX / bi_synthetic_XXXX / replication_synthetic_XXXX
# - Ground-truth subgroup/NF2 effects were DELIBERATELY SEEDED into the data
#   to match PLAN.md hypotheses, so statistical "significance" here is
#   guaranteed by construction — it tests the pipeline machinery, not biology.
#
# Contents:
#   results/           — Tables 1-S5 and Figures 1-5 + S1 (all synthetic)
#   manuscript/        — paper_manuscript.md, one_page_summary_for_warta.md,
#                        pre_submission_checklist.md (all describe synthetic
#                        numbers as if they were real findings)
#   data_processed_discovery_nassiri/   — Synthetic Nassiri-format cohort (n=185)
#   data_processed_discovery_bi/        — Synthetic Bi-lab-format cohort (n=565)
#   data_processed_replication_gse136661/ — Synthetic replication cohort (n=242)
#
# Allowed use:
#   ✅ Comparing pipeline output format (CSV shapes, column names)
#   ✅ Testing stats.py / ssgsea.py logic against known ground truth
#   ❌ Reading any number as a "result" or "finding"
#   ❌ Sending manuscript/, one_page_summary, or any figure to Warta / journals / collaborators
#   ❌ Copying any file back into results/, manuscript/, or data/processed/
#
# Regeneration of results/manuscript is permitted ONLY from real data acquired
# via:
#   scripts/fetch_geo_gse180061.py   (Nassiri — GSE180061)
#   scripts/fetch_geo_gse212666.py   (Bi lab — GSE212666)
#   scripts/fetch_cbio_mng_utoronto_2021.py   (cBioPortal mRNA)
# followed by:
#   scripts/01_aim0_cohort_lock.py --require-real-data
#
###############################################################################
