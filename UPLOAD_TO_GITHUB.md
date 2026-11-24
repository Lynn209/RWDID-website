# Upload to GitHub Pages

1. Unzip this package.
2. Upload **all files and folders** to the root of `Lynn209/RWDID-website`, replacing files with the same names.
3. Check that `index.html` and `.github/workflows/pages.yml` are both at the repository root / expected path.
4. Open **Settings → Pages** and set the source to **GitHub Actions**.
5. Push a commit to `main`.
6. Open **Actions**. Both `build` and `deploy` must finish successfully.
7. Open `https://lynn209.github.io/RWDID-website/` and hard-refresh with Ctrl+F5.

The homepage is now calculator-first. Methods, analysis code and input instructions are on separate pages. The exact `cols4all::friendly13` palette is generated and validated during the Pages build; no alternative palette is substituted.
