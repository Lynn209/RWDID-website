# RWDID Website V6 — QA report

## Scope

This revision corrects the linked UK Biobank bread-roll type selector for the active refined-grains pseudo-item. It does not change the 15 components, 82 scored intake fields, score intervals, beta coefficients, continuous targets, P1–P99 support, or score equations.

## Source rule reproduced

- Amount field: UKB Field 101160, Data-Coding 100006.
- Type field: UKB Field 20094, categorical multiple, Data-Coding 81.
- Type values: 1=White; 2=Mixed; 3=Wholemeal; 4=Seeded; 5=Other.
- The original Step_02_1_1 transformation uses `p20094_iN == "[1]"` before passing the Field 101160 amount into the retained white-bread pseudo-item.
- The calculator therefore includes the amount only when White is the sole selected type. Other single types and multi-type selections (including White plus another type) contribute zero to this derived item.
- Field 20094 is a linked selector, not an additional scored food item or a new dietary component.

## Automated validation

`node scripts/validate_questionnaire_mapping.js` completed successfully and confirmed:

1. Exactly 82 active intake fields: 76 categorical Oxford WebQ fields and 6 continuous grams/day fields.
2. Exactly 15 frozen dietary components.
3. All categorical item-specific options and the six continuous field IDs remain valid.
4. Field 101160 links to Field 20094 with Data-Coding 81 and all five correct type options.
5. Exact White-only selector gating is enforced; non-White and multiple selections do not contribute.
6. The cols4all `friendly13` asset contains 13 unique HEX colors and the scoring-engine smoke test passes.

## Chromium interaction tests

The production `app.js`, `scoring.js`, CSS, palette asset and JSON parameter/mapping files were loaded into headless Chromium. Only the three static JSON fetches were stubbed to avoid network access. The following cases passed:

| Scenario | Refined-grains intake | Binary component score | Result |
|---|---:|---:|---|
| White only, 1 roll | 1.000 | 0 | PASS |
| White only, 3 rolls | 3.000 | 1 | PASS |
| Wholemeal only, 1 roll | 0.000 | 0 | PASS |
| White + Mixed, 1 roll | 0.000 | 0 | PASS |
| No bread roll + 1 naan | 1.000 | 0 | PASS |

Additional UI checks passed: amount is required when a type is selected; an unconfirmed empty selector cannot complete the group; group-level no-intake completes the group; Clear resets the linked type and amount; Reset returns progress to 0/15. No JavaScript page errors or console errors were detected.

A second synthetic end-to-end profile was constructed with every component intake placed inside its frozen binary target interval. All 15 component indicators were 1, yielding RWDID **1.0000** and RWDID-C **15.00**. This is a deliberately target-aligned software test profile, not a real participant or a clinical low-risk classification. An all-groups-no-intake test yielded RWDID **0.2425** and RWDID-C **4.00**, because zero lies within the frozen target intervals for four components; this is an expected property of the scoring rules, not a calculation error.

## Frozen file integrity

- `analysis_code/RWDID_complete_analysis.R` unchanged.
- `data/rwdid_parameters.json` unchanged (SHA-256: `048359e910291a9ec6abd4c52c5008412efc8678387e07c4c7cbd1d06c9421d5`).
- `data/friendly13_palette.json` unchanged.
- Binary and continuous score parameter files unchanged.

## Publication note

These are local Chromium and static-file tests. They do not assert that GitHub Pages has deployed this revision. After uploading the package, confirm that the GitHub Actions Pages workflow succeeds and hard-refresh the live site before treating the public deployment as updated.
