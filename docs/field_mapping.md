# Questionnaire field mapping and input logic

The website mapping starts from the final 15 food components in Table S1, selects the associated Field IDs, and validates them against `field.csv`.

- Table S1 contains 206 food-item rows across 34 groups.
- The final 15 components contribute 85 Table S1 rows.
- All 85/85 numeric Field IDs matched the source dictionary.
- Three of the 85 rows are not active calculator inputs; the frozen Step_02_1_1 implementation retains 82 active inputs.
- The 82 active inputs comprise **76 categorical Oxford WebQ response fields** and **6 continuous UKB Showcase total-weight fields**. One categorical input (`w-101160`) is a derived white-bread row and requires the dependent multi-response type selector, Field 20094.

## Correct questionnaire interaction

The Oxford WebQ presents categorical food items as a checkbox and, when selected, a single-choice amount response. The website mirrors that sequence:

1. Tick `Consumed yesterday` for each item consumed.
2. Choose the corresponding amount category.
3. Leave an unconsumed item unticked; it contributes zero.
4. Alternatively, select `No intake from this food group` to set all fields in that component to zero.

The amount selector does not include a synthetic `0 / Not consumed` category. This is important because zero intake is represented by an unchecked item in the original questionnaire, not by one of the categorical amount codes.

For the derived white-bread item `w-101160`, the UI also requests Field 20094 (bread-roll type). The amount is included in refined grains only when the recorded type response is exactly White alone (`[1]`), reproducing the original Step_02_1_1 code. Multiple selected bread types are not treated as White alone.

## Coded response values

Categorical response options are item-specific and use the exact wording from the UK Biobank Oxford WebQ v1.5 Appendix 1. Stored code-to-value conversion follows the UKB Showcase Data-Coding ID associated with each field and the frozen R pipeline.

- 100017: ¼, ½, 1, 2, 3+.
- 100004: ½, 1, 2, 3, 4+.
- 100006: ½, 1, 2, 3, 4, 5, 6+.
- 100003: ½, 1, 2, 3+.
- 100016: ½, 1, 2, 3, 4, 5+.
- 100001: Showcase codebook has ½, 1, 2+.

### Data-Coding 100001 discrepancy

The Oxford WebQ v1.5 PDF appendix shows separate 2 and 3+ labels for the five selected grain fields, but Showcase Data-Coding 100001 lists only three released values (`555=half`, `1=1`, `200=2+`). The frozen RWDID code maps 200 to 2.0. To reproduce the released UKB fields and the frozen analysis, the website retains one `2+` option mapped to 2.0. It does not fabricate a separate 3+ stored category.

## Continuous grams/day fields

Fields 26102, 26131, 26133, 26151, 26152 and 26153 are UKB Showcase `Total weight by food group yesterday` fields with value type `Continuous, g`. They are not categorical Oxford WebQ responses. The website accepts grams/day and divides by the corresponding Table S1 standard portion size.

## Repeated recalls

The original questionnaire asks about the previous 24 hours. The research pipeline averages available repeated recall instances at the item level before constructing component intakes. The browser calculator accepts one representative-day response per item; it does not reproduce repeated-instance averaging or impute missing recalls.

See [`ukb_questionnaire_response_mapping.md`](ukb_questionnaire_response_mapping.md) for the complete 82-field option table and source links.
