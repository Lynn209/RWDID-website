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
