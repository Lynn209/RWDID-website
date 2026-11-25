# Synthetic scoring-engine test profile

This is a deterministic software test, not a real participant record, dietary recommendation, or clinical risk assessment.

The test chooses component intakes that fall inside the current frozen RWDID binary target intervals and the RWDID-C full-score targets. It is designed to verify that the website aggregates the 82 field inputs correctly and that the two scoring functions return their expected upper bounds.

| Component | Synthetic intake (portions/day) | Rule satisfied |
|---|---:|---|
| Citrus fruits | 1.00 | 0.034338–1.682583 |
| Red/orange vegetables | 1.00 | ≥0.010 |
| Other vegetables | 2.00 | 0.010–4.620 |
| Potatoes | 1.25 | 1.00–1.50 |
| Refined grains | 3.00 | ≥3.00 |
| Eggs | 1.00 | 0.030151–1.869347 |
| Low-fat dairy | 1.00 | 0.135135–2.196911 |
| Red meats | 1.00 | 1.00–1.75 |
| Processed meats | 3.00 | ≥2.70 |
| Fish and seafood | 1.00 | 0–1.98 |
| Coffee | 2.00 | 0.01–4.50 |
| Fruit juices | 1.00 | 0–1.37 |
| Sweetened beverages | 0.00 | 0–2.08 |
| Water | 2.00 | 0.01–3.08 |
| Wine | 0.00 | 0–2.00 |

Field-level test selections are defined in `scripts/test_scoring_engine.js`. Unselected categorical items contribute zero. Low-fat yogurt is entered as 125 g/day (one Table S1 portion); the other five continuous Showcase fields are entered as 0 g/day.

Expected deterministic output using the bundled frozen parameters:

- RWDID = 1.0000
- RWDID-C = 15.00
- All 15 binary component indicators = 1
- All 15 continuous component scores = 1

**Important interpretation:** the test profile is constructed to maximize the index under the study's frozen scoring rules. In particular, some component directions reflect the fitted study-specific score definitions and should not be interpreted as general healthy-eating guidance. Do not use this profile as dietary advice.
