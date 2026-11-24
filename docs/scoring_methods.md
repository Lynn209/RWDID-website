# Scoring methods

## RWDID: binary, beta-weighted score

For component i, the binary indicator S_i is 1 when the food-group intake meets its frozen direction-specific interval and 0 otherwise. The composite score is:

\[ RWDID = \frac{\sum_{i=1}^{15} S_i \beta_i}{\sum_{i=1}^{15} \beta_i} \]

The current uploaded 15-component beta file has an algebraic sum of **-1.621792082127**. All 15 beta coefficients are negative in this frozen file, so normalised weights are positive and RWDID is bounded from 0 to 1.

## RWDID-C: continuous, equal-weight score

Each component intake is clamped to its training-derived P1–P99 support interval, then mapped to a 0–1 component score using the frozen function:

- `high`: linear adequacy ramp from training P1 to `value1`, then 1.
- `low`: 1 up to `value1`, then a linear decline to training P99.
- `middle`: 1 within [`value1`, `value2`], with linear ramps to training P1 and P99.
- `point`: triangular score centred on `value1`.

\[ RWDID-C = \sum_{i=1}^{15} C_i \]

The theoretical range is 0–15. RWDID-C is not beta-weighted.

## Missing values

No component is imputed. The website requires every item to be explicitly answered or the user to confirm no intake for the food group. Blank answers are not silently converted to zero.

## Frozen parameter files

- `data/RWDID_binary_scoring_rules_15_components.csv`
- `data/beta_coefficients_training.csv`
- `data/Table_Sx_continuous_scoring_parameters_15.csv`
- `data/training_only_P1_P99_support.csv`
- `data/food_rcs_range_all_final.csv`

The website does not refit Cox models, estimate new beta coefficients, select new intervals, or recalculate P1/P99 limits.
