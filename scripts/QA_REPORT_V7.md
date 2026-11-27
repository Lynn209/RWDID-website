# RWDID Website V7.0.0 — Bread-roll selector deployment fix

## Reason for revision

The live acceptance check reported that the Refined grains item showed only “Bread roll intake (white)” and did not visibly expose the Wholemeal/Mixed choices. The V6 local source contained the selector, but the public page had not demonstrated that the same build was live. This revision makes the selector self-contained, explicit and cache-versioned.

## Visible UI change

The Refined grains item is now labelled **Bread roll intake** (rather than “Bread roll intake (white)”). Its linked Field 20094 block explicitly displays five checkbox choices before the amount control:

1. White
2. Mixed
3. Wholemeal
4. Seeded
5. Other

The block identifies UKB Field 20094, Data-Coding 81, and says “Select all that apply”. The Field 101160 amount selector becomes visible when at least one type is selected.

## Frozen scoring rule

The original Step_02_1_1 transformation includes the Field 101160 amount only when the Field 20094 response is exactly the single code `[1]`. Therefore:

- White only: amount included.
- Mixed only, Wholemeal only, Seeded only, Other only: amount contributes zero to the derived white-bread item.
- White plus any other type: amount contributes zero.

Field 20094 is a linked selector and is not an 83rd scored intake field.

## Deployment and cache protection

- The application normalizes the Field 101160 selector metadata to the frozen Field 20094 / Data-Coding 81 definition at runtime, so an older cached mapping JSON cannot silently remove the five visible options.
- CSS, JavaScript, palette JSON, field mapping JSON and score parameter JSON requests use the `v=7.0.0` cache key.
- The calculator header visibly displays `V7.0.0`, allowing the live build to be identified without relying on a cached screenshot.
- The validation script checks the five selector labels, runtime normalization and build/cache markers.

## Tests run in this revision

- `node --check assets/app.js`: PASS.
- `node --check assets/palette.js`: PASS.
- `node --check assets/methods.js`: PASS.
- `node scripts/validate_questionnaire_mapping.js`: PASS.
- `node scripts/test_v7_bread_selector.js`: PASS.
- Scoring engine: White only (1 roll) => derived intake 1.000; Mixed only, Wholemeal only and White + Mixed => derived intake 0.000.
- With all other components set to zero, White only (1 roll) => RWDID 0.2425; RWDID-C 4.33. Mixed only, Wholemeal only and White + Mixed => RWDID 0.2425; RWDID-C 4.00.
- White only (3 rolls) => refined-grains binary score 1; RWDID 0.2711; RWDID-C 5.00.
- All-zero synthetic profile => RWDID 0.2425; RWDID-C 4.00.
- All-target synthetic profile => RWDID 1.0000; RWDID-C 15.00.

These are source, mapping and scoring-engine tests. They do not claim that GitHub Pages has deployed V7.0.0 or that a live browser session has been operated. Live acceptance requires a successful latest Pages workflow and visual confirmation that the page header says V7.0.0 and all five Field 20094 options are visible.

## Frozen files preserved

- 15 components and 82 active scored fields.
- Training-derived beta coefficients and beta denominator.
- Binary component intervals.
- RWDID-C continuous parameters and training P1–P99 support.
- Both score equations.
- `analysis_code/RWDID_complete_analysis.R`.
- cols4all `friendly13` palette asset and exact ordered HEX vector.
