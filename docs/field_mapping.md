# Questionnaire field mapping

The mapping was built by selecting the final 15 food groups from the uploaded Table S1, extracting all associated FieldIDs, and matching the numeric FieldID to `field.csv`.

- Table S1 contains 206 food-item rows across 34 groups.
- The final 15 components contribute **85 Table S1 item rows**.
- All **85/85** numeric FieldIDs matched `field.csv`.
- 81 are direct fields; four are white-bread pseudo-fields.
- The frozen Step_02_1_1 script writes all four white-bread pseudo-fields sequentially to the same `pw100950_i0–i4` variable. The final assignment is white bread roll (`w-101160`), selected through Field ID 20094. Thus, the frozen score retains 81 direct fields plus this one pseudo-field: **82 active inputs**. The other three pseudo-fields remain in the audit CSV but are not scored.

## Unit conversion

- Categorical fields use the encoded UK Biobank response-to-portion mappings.
- `encoding_id = 0` fields are entered in grams/day and divided by the Table S1 portion size in grams.
- Item-level standard portions/day are summed within each food group.
- The original research pipeline averaged available repeated recall instances for each item before summing the food group. The website therefore asks for an average daily amount across typical days.

## Legacy Shiny mapping

The old `field_data_final.csv` contains 98 rows. It includes non-final groups (Green leafy vegetables, High-fat dairy products) and omits final groups (Fruit juices, Water). It is superseded by this Table S1 + field.csv mapping.
