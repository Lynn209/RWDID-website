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

## V6 bread-roll selector update

This package also corrects the linked UKB questionnaire fields:

- Field 101160: bread-roll amount (Data-Coding 100006).
- Field 20094: bread-roll type (categorical multiple; Data-Coding 81).
- The amount is included in the derived white-bread item only when White (code 1) is the sole selected type, reproducing the frozen R preprocessing condition exactly.

The linked type selector is not an additional scored intake item. The 82 active scoring inputs and all frozen score parameters are unchanged.

**Important:** this ZIP does not replace the repository's existing `.github/workflows/pages.yml`. Keep that workflow in the repository. Upload the ZIP contents over the existing site files; do not delete the `.github` directory. After commit, verify the Actions run and then hard-refresh the deployed page.

For the detailed audit, see `docs/QA_REPORT_V6.md`.
