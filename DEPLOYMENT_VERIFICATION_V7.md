# RWDID Website V7.0.0 — Deployment verification

## What changed

- The Refined grains card now displays the item as **Bread roll intake** (without the misleading suffix “(white)”).
- The linked UKB Field 20094 selector is displayed as a clearly labelled, always-visible five-option multiple-response checklist: White, Mixed, Wholemeal, Seeded and Other.
- The Field 101160 amount selector appears after one or more types are selected.
- The amount contributes to the frozen refined-grains white-bread pseudo-item only when the selected type is exactly White alone (code 1).
- Non-White-only responses and all multi-type combinations contribute zero to this derived item.
- The selector is normalized from the fixed UKB Field 20094 / Data-Coding 81 definition at runtime, so a stale mapping JSON cannot silently remove Mixed/Wholemeal options.
- CSS, JavaScript, palette and JSON requests use the V7.0.0 cache-busting version. The header displays V7.0.0.

## Frozen items not changed

No change to the 15 dietary components, 82 scored intake fields, binary cut-points, training beta coefficients, RWDID-C continuous parameters, P1–P99 support, score equations, or consolidated R analysis code. Field 20094 remains a linked questionnaire selector and is not counted as an 83rd scored item.

## Required post-upload verification

1. Upload the *contents* of this ZIP to the repository root, not the enclosing folder.
2. Preserve the repository's existing Pages workflow. Commit to `main`.
3. Wait for the latest Pages workflow run to finish successfully.
4. Hard-refresh the site. Confirm the header shows `V7.0.0`.
5. Open Refined grains. Confirm all five Field 20094 options are visible before selecting an amount.
6. Test White only, Mixed only, Wholemeal only, and White + Mixed. The amount should be included only for White alone.

If the header still does not show V7.0.0, the V7 index page has not been deployed or the wrong repository folder was published.
