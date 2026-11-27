# Consolidated RWDID analysis code: source provenance

The public repository contains exactly one analysis source file:

`analysis_code/RWDID_complete_analysis.R`

It contains 11 ordered source blocks. The split-generator block has one documented path normalization to the shared downstream `results/Step_03/ukb_70_30_split` directory; the split algorithm and seed are unchanged.

1. **Step_02_1_1** — `Step_02_1_1_data_preparation_and_wrangling_food_group(5).R` (733 lines; local current source)
2. **Step_03_1_0_SPLIT** — `Step_03_1_0_ukb_70_30_split_manifest_and_cohort_generation(1).R` (869 lines; Library source dated 2026-11-23; SPLIT_DIR normalized to the shared path used by all downstream Step_03/Step_04 scripts)
3. **Step_03_1_0_MODEL2_REFERENCE** — `Step_03_1_0_association_between_food_group_and_depression_70_30_MAIN_ROTATING_LOWEST_RISK_REFERENCE_V1.R` (3150 lines; Library source dated 2026-11-05)
4. **Step_03_2_7_FINAL15_AUDIT** — `Step_03_2_7_final15_original_table_item_level_portion.R` (1667 lines; local current source)
5. **Step_04_1_1_1_RCS_BASE** — `Step_04_1_1_1_rcs_legacy_anchored.R` (3986 lines; project source archive data_analysis_v2.zip)
6. **Step_04_1_1_2_INTERVAL_RECONCILIATION** — `Step_04_1_1_2_data_derived_interval_reconciliation.R` (1016 lines; project source archive data_analysis_v2.zip)
7. **Step_04_1_1_3_RCS20_OVERLAP** — `Step_04_1_1_3__rcs_plateau_overlap.R` (835 lines; project source archive data_analysis_v2.zip)
8. **Step_04_1_1_4_FINAL_RCS_BINARY** — `Step_04_1_1_4_RCS_analysis_between_food_and_depression_FINAL_B500_PARALLEL_BOOTSTRAP_5_10_20_FONT_REVISION.R` (5476 lines; Library source dated 2026-11-23)
9. **Step_04_1_2_RWDID_C** — `Step_04_1_2_RCS_continuous_proportional_scoring_colizzi(2).R` (1428 lines; Library source dated 2026-11-23)
10. **Step_04_2_1_BETA_WEIGHTED_TRAINING** — `Step_04_2_1_beta_weighted_diet_score_training_V1_3_RCS_DATADIST_FIX(20261118-041346).R` (3239 lines; Library source dated 2026-11-18; current Step_04_2_1 V1.3 RCS/DATADIST fix)
11. **Step_04_2_2_HELDOUT_VALIDATION** — `Step_04_2_2_beta_weighted_diet_score_test_validation_V1_4_RCS_DATADIST_FIX(6).R` (2959 lines; Library source dated 2026-11-16; current Step_04_2_2 V1.4 RCS/DATADIST fix)

## Upstream dependencies

The analysis code requires authorised UK Biobank inputs and preprocessed project data objects. The final RCS stage also uses the preceding RCS anchor outputs. These inputs are listed in the consolidated script header; no individual-level participant data are included in the public repository.
