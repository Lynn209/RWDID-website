# Upload to GitHub Pages

1. Unzip this package.
2. Upload **all files and folders** to the root of `Lynn209/RWDID-website`, replacing files with the same names.
3. Check that `index.html` and `.github/workflows/pages.yml` are both at the repository root / expected path.
4. Open **Settings → Pages** and set the source to **GitHub Actions**.
5. Push a commit to `main`.
6. Open **Actions**. Both `build` and `deploy` must finish successfully.
7. Open `https://lynn209.github.io/RWDID-website/` and hard-refresh with Ctrl+F5.

The homepage is now calculator-first. Methods, analysis code and input instructions are on separate pages. The exact `cols4all::friendly13` palette is generated and validated during the Pages build; no alternative palette is substituted.


The questionnaire now mirrors the UKB Oxford WebQ checkbox → amount workflow for 76 categorical items. Six Showcase total-weight fields are entered in g/day. The complete response-category crosswalk and coding exceptions are documented in `docs/ukb_questionnaire_response_mapping.md`.


## v5 UKB response-alignment checks
- The bread-roll amount field 101160 is linked to the multiple-response bread-type selector 20094.
- The website includes the amount in the refined-grains component only when Field 20094 is exactly White alone (`[1]`), matching the frozen R pipeline.
- The five Data-Coding 100001 grain fields retain the Showcase `2+` option mapped to 2.0.
- The build checks all 82 mappings, 15 scoring parameters, exact Friendly13 HEX order, and a synthetic maximum-score scenario.
