# Change log — UKB questionnaire alignment

## Revision: UKB questionnaire-aligned input flow

### User-facing changes
- Categorical Oxford WebQ fields now reproduce the original item-checkbox → amount-response sequence.
- Amount selectors are hidden until the corresponding categorical item is ticked as consumed.
- Unchecked categorical items contribute zero. A synthetic `0 / Not consumed` option has been removed.
- Item-specific response wording is shown for every categorical Field ID, including units such as servings, whole fruit, eggs, slices, sticks, cups, glasses, cans, bottles, inches, cloves and rasher counts.
- The six continuous Showcase fields display direct grams/day inputs and are not presented as categorical WebQ items.
- Group-level “No intake” remains as a calculator convenience and sets all items in that component to zero.
- Mixed-vegetable/side-salad and sushi instructions are displayed to reduce duplicate recording.

### Source and coding reconciliation
- 76 categorical fields use the Oxford WebQ v1.5 response wording and their UKB Showcase Data-Coding IDs.
- 6 continuous fields (26102, 26131, 26133, 26151, 26152, 26153) use Showcase Continuous, g values.
- Data-Coding 100001 retains the released `200 = 2+` category and maps it to 2.0, matching the frozen R pipeline. The PDF v1.5 display discrepancy is documented.
- Other egg dishes uses Showcase Field ID 102980. The PDF v1.5 appendix appears to print 120980; the released Showcase field, Table S1 mapping and frozen analysis use 102980.

### Research calculations
- No binary scoring intervals, β coefficients, continuous scoring targets, P1–P99 support values, or score formulas were changed.
- The consolidated R analysis code was retained unchanged.
- GitHub Actions now validates the 82-field mapping, option counts, six continuous fields, the Data-Coding 100001 exception, palette, and score-engine smoke test before deployment.


## Revision: V6 — linked bread-roll type selector (Field 20094)

### User-facing changes
- Field 101160 now displays the linked Field 20094 type selector with all five UKB Data-Coding 81 choices: White, Mixed, Wholemeal, Seeded and Other.
- The amount response is required when any bread-roll type is selected; users can explicitly confirm that no bread roll was consumed.
- The group-level no-intake and clear/reset controls also reset and disable the linked selector correctly.

### Frozen analysis compatibility
- The calculator applies the original Step_02_1_1 condition exactly: Field 101160 contributes only when Field 20094 is the exact single selection White (`[1]`).
- Other single types and multi-type selections (including White plus another type) contribute zero to the derived white bread-roll item, as in the frozen source code.
- The 82 scored input fields, 15 components, RWDID intervals, β coefficients, RWDID-C targets, P1–P99 supports and score equations are unchanged. Field 20094 is a linked selector, not an extra scored item.
- Added machine-readable selector metadata to both mapping CSVs and the active 82-field JSON.

### Validation
- Mapping, selector coding, selector gating, palette and score-engine tests are run before packaging.
- Browser interaction tests cover White only, non-White only, White plus another type, no bread-roll intake, group-level no intake, clear and reset.


## V7.0.0 — Explicit Field 20094 selector and deployment versioning

- Made the bread-roll type block explicitly visible with all five options: White, Mixed, Wholemeal, Seeded and Other.
- Renamed the displayed item to “Bread roll intake” to avoid implying the type selector is already fixed to White.
- Kept the exact frozen White-only scoring gate.
- Added runtime selector normalization for Field 101160 → Field 20094 (Data-Coding 81).
- Added V7.0.0 cache-busting to CSS, JavaScript and JSON requests, plus a visible build marker.
- Added V7 mapping/deployment validation and a scoring regression test.
- No score definitions or research analysis code changed.
