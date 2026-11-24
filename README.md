# RWDID GitHub Pages calculator

Static, browser-only research calculator for the 15-component Risk-Weighted Dietary Index (RWDID) and continuous equal-weight RWDID-C.

## Repository layout

```text
index.html
input-guide.html
assets/
  app.js
  scoring.js
  styles.css
data/
  field_mapping_active82.json
  field_mapping_active82.csv
  field_mapping_final15_85.csv
  rwdid_parameters.json
  friendly13_palette.json
  ...
analysis_code/
  RWDID_complete_analysis.R
scripts/
  export_friendly13_palette.R
docs/
  ...
.nojekyll
README.md
README_CN.md
```

`index.html` is at the repository root. This is a static site: the browser performs the score calculation locally. No Shiny server, Python server, Node build, API key, database, or GitHub Actions workflow is required to serve the published website.

## Publish on GitHub Pages

1. Create a public GitHub repository, for example `RWDID-website`.
2. Upload **the contents of this folder** to the repository root. `index.html` must be in the root.
3. Commit/push to the `main` branch.
4. Open **Settings → Pages**.
5. Under **Build and deployment**, select **Deploy from a branch**.
6. Choose branch **main** and folder **/(root)**, then click **Save**.
7. GitHub will provide a URL such as `https://YOUR-USERNAME.github.io/RWDID-website/`.

This follows the same branch-root publishing approach used by the public betaSTEER_web repository. No `.github/workflows` directory is needed.

## One-time exact friendly13 palette export

The source ZIP includes an empty palette placeholder rather than guessed HEX values. To ensure the page uses the exact ordered output of cols4all, generate the palette once before uploading:

```bash
Rscript scripts/export_friendly13_palette.R
```

The script calls:

```r
cols4all::c4a(
  palette = "friendly13",
  n = 13,
  type = "cat",
  format = "hex",
  verbose = FALSE
)
```

It writes the exact 13 returned HEX colors and the installed package version to `data/friendly13_palette.json`. Commit this generated JSON with the website. The browser validates the palette identity, count, HEX format and uniqueness. It does not silently substitute another palette. If the palette file has not been generated, the page displays a visible warning and uses neutral presentation colors; it does not claim to be using friendly13.

## Local preview

Because the page uses browser `fetch()`, preview through a local static server rather than opening `index.html` as `file://`:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000/`.

## Calculator inputs and scoring

- 15 frozen dietary components.
- 82 active questionnaire inputs, mapped from the 85-row Table S1 item audit.
- RWDID: binary component scores weighted by the 15 training-set Cox β coefficients and divided by their algebraic sum; range 0–1.
- RWDID-C: sum of 15 continuous proportional component scores, each bounded 0–1; range 0–15.
- No imputation. Missing item responses are not silently converted to zero.
- The browser applies frozen parameters; it does not refit the cohort Cox or RCS models.

## Complete underlying analysis code

`analysis_code/RWDID_complete_analysis.R` consolidates the supplied research scripts into one downloadable source file. It contains dietary preparation, final-15 item/portion QC, RCS and binary interval scoring, RWDID-C, β-weighted training analysis, and the supplied held-out test validation source.

**Important reproducibility boundary:** individual-level UK Biobank data are not included. The frozen 70/30 split and the upstream Model 2 / RCS20 anchor artifacts must already exist in the authorised project environment. The supplied Step_04_2_1 V5.3.1 and Step_04_2_2 V1.4 scripts also contain different versioned output-path conventions; the consolidated file documents this rather than silently altering model definitions or claiming an unverified one-command rerun. See the header of the consolidated analysis code and `docs/analysis_code_provenance.json`.

## Data protection and intended use

No participant-level UK Biobank data, EIDs, or participant-level training scores are included. The calculator processes entered responses in the browser and does not transmit them to a server. For research use only; not a clinical risk prediction or diagnostic tool.
