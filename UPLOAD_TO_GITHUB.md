# Upload to GitHub Pages — V7.0.0

1. Unzip this package. The ZIP is arranged with `index.html`, `assets/`, `data/`, `docs/` and `scripts/` at its top level.
2. In `Lynn209/RWDID-website`, upload the **contents** of the extracted package directly to the repository root. Do not upload the enclosing folder as a nested subfolder. Replace same-named files.
3. Preserve the repository's existing `.github/workflows/pages.yml`; this package intentionally does not replace the existing Pages workflow. Do not delete `.github/`.
4. Commit the files to `main`.
5. Open Actions and wait for the latest Pages build/deployment run to finish successfully.
6. Open `https://lynn209.github.io/RWDID-website/` and hard-refresh with Ctrl+F5.
7. Confirm that the header shows **V7.0.0**. Then open Refined grains and confirm that all five Field 20094 choices (White, Mixed, Wholemeal, Seeded, Other) are visible.

## What V7 fixes

- The bread-roll card is displayed as “Bread roll intake”.
- UKB Field 20094 is shown explicitly as a multiple-response selector with five visible choices.
- Field 101160 amount selection appears after at least one type is selected.
- The frozen analysis rule is unchanged: Field 101160 contributes only when Field 20094 is exactly the single selection `[1]` (White only). Mixed, Wholemeal, Seeded, Other, and White plus any other type contribute zero to the derived white-bread item.
- Field 20094 remains a linked questionnaire selector, not an 83rd scored input.
- CSS, JavaScript and JSON requests use `v=7.0.0` cache-busting; the header build marker lets you identify the deployed version.

## Validation

Run these commands in the extracted folder if you want to repeat the local checks:

```bash
node --check assets/app.js
node --check assets/palette.js
node --check assets/methods.js
node scripts/validate_questionnaire_mapping.js
node scripts/test_v7_bread_selector.js
```

See `docs/QA_REPORT_V7.md` for scope and results. These checks do not by themselves prove that GitHub Pages has deployed V7; the header marker and visible five-choice selector must be checked after deployment.
