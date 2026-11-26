# RWDID GitHub Pages — UKB questionnaire-aligned calculator

A static, browser-based calculator for the frozen 15-component RWDID and RWDID-C scoring framework.

## What changed in this revision
- The questionnaire now follows the Oxford WebQ interaction pattern: tick an item as consumed, then select its amount category.
- Unchecked categorical items contribute zero; the amount selector does not contain a synthetic zero option.
- Each field shows its item-specific UKB response labels rather than a generic numeric range.
- The 82 active inputs are correctly described as 76 categorical Oxford WebQ fields plus 6 continuous UKB Showcase total-weight fields.
- The six continuous fields (26102, 26131, 26133, 26151, 26152, 26153) accept grams/day; they are not categorical WebQ response fields.
- The discrepancy for Data-Coding 100001 is documented: lists 2 and 3+ separately, but the released Showcase coding has 200=2+ and the frozen R pipeline maps 200 to 2.0. The calculator follows the released coding and the frozen analysis.
- UKB questionnaire cautions are included for mixed vegetables/salads and sushi to reduce double counting.
- The exact 13-color `cols4all::friendly13` asset is included.

## Pages
- `index.html`: calculator-first homepage and 15-component questionnaire.
- `methods.html`: frozen scoring formulas and parameter summary.
- `resources.html`: consolidated analysis code, field mapping and score parameters.
- `input-guide.html`: concise input instructions and source notes.
- `docs/ukb_questionnaire_response_mapping.md`: complete 82-field response option mapping and coding provenance.

## Exact Friendly13 palette
`data/friendly13_palette.json` contains the ordered HEX vector returned by the user's verified `cols4all::c4a("friendly13")` output using `cols4all` 0.11. The GitHub Actions workflow regenerates and validates the same palette before deployment. The site does not substitute another palette.

## Scoring
- RWDID: 15 binary components weighted by frozen training-derived β coefficients; range 0–1.
- RWDID-C: 15 continuous proportional component scores summed with equal weights; range 0–15.
- The website applies frozen parameters; it does not refit Cox or RCS models.
- No imputation is performed.

## Privacy and use
Questionnaire inputs are processed locally in the browser and are not uploaded. RWDID and RWDID-C are research indices, not clinical risk predictions or diagnostic tools. Public files do not include individual-level UK Biobank data.


## Questionnaire linkage update
The Field 101160 bread-roll amount is now linked to the UKB Field 20094 multiple-response type selector (Data-Coding 81). The amount contributes to the frozen refined-grains white-bread pseudo-item only when White (code 1) is the sole selected type, exactly matching the original Step_02_1_1 `p20094_iN == "[1]"` rule. This linked selector is not an additional scored input or component. The frozen score parameters and formulas are unchanged. See `docs/ukb_questionnaire_response_mapping.md` and `CHANGELOG_UKB_QUESTIONNAIRE_FIX.md`.
