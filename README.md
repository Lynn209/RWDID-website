# RWDID GitHub Pages — Minimal Calculator Edition

A static, browser-based calculator for the frozen 15-component RWDID and RWDID-C scoring framework.

## Pages
- `index.html`: concise calculator-first homepage with the 15-group questionnaire, progress, results, and local-processing notice.
- `methods.html`: scoring formulas and the frozen 15-component parameter table.
- `resources.html`: one consolidated research analysis script and essential mapping/parameter downloads.
- `input-guide.html`: short questionnaire input instructions.

The homepage intentionally excludes the three-step cards, full parameter table, palette swatches, audit documents, and long research explanations. Those details are separated into dedicated pages.

## Exact Friendly13 palette
Categorical accents are generated from `cols4all::friendly13` using `cols4all` 0.11. The GitHub Actions workflow runs `scripts/export_friendly13_palette.R` and validates the generated 13 unique HEX values before deployment. The site does not substitute another palette.

## Deploy
1. Upload all package contents to the root of `Lynn209/RWDID-website`.
2. In GitHub **Settings → Pages**, set the source to **GitHub Actions**.
3. Push to `main`.
4. Confirm both build and deploy jobs are successful in **Actions**.
5. Hard-refresh the published site.

## Scoring
- RWDID: 15 binary components weighted by the frozen training-derived β coefficients; range 0–1.
- RWDID-C: 15 continuous proportional component scores summed with equal weights; range 0–15.
- The website applies frozen parameters; it does not refit Cox or RCS models.
- Blank responses remain missing; the calculator does not impute.

## Privacy and use
Questionnaire inputs are processed locally in the browser and are not uploaded. RWDID and RWDID-C are research indices, not clinical risk predictions or diagnostic tools. Public files do not include individual-level UK Biobank data.
