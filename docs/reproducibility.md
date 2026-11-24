# Reproducibility and public release

This repository separates the user-facing calculator from the underlying research analysis code. The calculator applies frozen scoring parameters. Core R scripts document food-group construction, field mapping audits, RCS interval derivation, continuous proportional scoring and beta estimation.

No individual-level UK Biobank records, EIDs or participant-level score files are included. The participant-level `beta_weighted_diet_score_training.csv` was used locally to validate the score formula and is deliberately excluded from the public repository. The analysis scripts require authorised access to the corresponding UK Biobank data and project inputs.

## Local preview

```bash
python -m http.server 8000
```

Open `http://localhost:8000`. Opening `index.html` directly may block JSON `fetch()` calls.

## GitHub Pages

1. Create a public repository and upload the contents of this folder.
2. Set Settings → Pages → Source to **GitHub Actions**.
3. Push to `main`; the included workflow exports the exact `cols4all::c4a("friendly13", 13)` palette and deploys the static site.
4. Run `Rscript scripts/export_friendly13_palette.R` locally if you want to regenerate the palette file before deployment.
