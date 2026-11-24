# UK Biobank questionnaire response mapping

## Source and scope

This file documents the 82 active inputs used by the frozen 15-component RWDID calculator. The 82 inputs comprise **76 categorical Oxford WebQ item fields** and **6 continuous UKB Showcase total-weight fields**. The latter are not categorical questionnaire-response items.

Primary questionnaire source: [UK Biobank 24-hour dietary recall questionnaire (Oxford WebQ), v1.5, Appendix 1](https://biobank.ctsu.ox.ac.uk/showcase/showcase/docs/DietWebQ.pdf).

Each categorical WebQ food item is presented as a checkbox followed by a single-choice amount response. The website therefore uses an item checkbox (`Consumed yesterday`) and reveals the amount selector only after the checkbox is selected. An unchecked item contributes zero. Zero is not inserted as a fake UKB coding category.

UKB Showcase coding sources: [100017](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100017), [100004](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100004), [100006](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100006), [100003](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100003), [100016](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100016), and [100001](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100001).

## Response coding and numeric conversion

| Showcase Data-Coding | Original stored code → portions/day |
|---|---|
| 100017 | 444 → 0.25; 555 → 0.5; 1 → 1; 2 → 2; 300 → 3 |
| 100004 | 555 → 0.5; 1 → 1; 2 → 2; 3 → 3; 400 → 4 |
| 100006 | 555 → 0.5; 1 → 1; 2 → 2; 3 → 3; 4 → 4; 5 → 5; 600 → 6 |
| 100003 | 555 → 0.5; 1 → 1; 2 → 2; 300 → 3 |
| 100016 | 555 → 0.5; 1 → 1; 2 → 2; 3 → 3; 4 → 4; 500 → 5 |
| 100001 | 555 → 0.5; 1 → 1; 200 → 2 |

These numeric values reproduce the recoding in the frozen RWDID preparation code (`444 → 0.25`, `555 → 0.5`, `200 → 2`, `300 → 3`, `400 → 4`, `500 → 5`, `600 → 6`). For open-ended categories, the frozen analysis uses the lower-bound numeric value shown above.

## Field-by-field response options

| Component | Field ID | Food item | Coding / input type | Portion (g) | Questionnaire response options | Source page |
|---|---:|---|---|---:|---|---:|
| Citrus fruits | 104490 | Grapefruit intake | Data-Coding 100004 | 148 | Half a grapefruit or half a serving; 1 grapefruit or serving; 2 grapefruits or servings; 3 grapefruits or servings; 4+ grapefruits or servings | 89 |
| Citrus fruits | 104530 | Orange intake | Data-Coding 100004 | 120 | Half an orange; 1 orange; 2 oranges; 3 oranges; 4+ oranges | 90 |
| Citrus fruits | 104540 | Satsuma intake | Data-Coding 100004 | 70 | Half a small citrus fruit or half a serving; 1 small citrus fruit or serving; 2 small citrus fruits or servings; 3 small citrus fruits or servings; 4+ small citrus fruits or servings | 90 |
| Red/Orange vegetables | 104150 | Butternut squash intake | Data-Coding 100017 | 130 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 82 |
| Red/Orange vegetables | 104170 | Carrot intake | Data-Coding 100017 | 60 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 83 |
| Red/Orange vegetables | 104290 | Sweet pepper intake | Data-Coding 100017 | 133 | Quarter of a whole pepper; Half a whole pepper; 1 whole pepper; 2 whole peppers; 3+ whole peppers | 85 |
| Red/Orange vegetables | 104330 | Sweet potato intake | Data-Coding 100017 | 130 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 85 |
| Red/Orange vegetables | 104340 | Fresh tomato intake | Data-Coding 100017 | 85 | Quarter of a medium tomato; Half a medium tomato; 1 medium tomato; 2 medium tomatoes; 3+ medium tomatoes | 85 |
| Red/Orange vegetables | 104350 | Tinned tomato intake | Data-Coding 100017 | 135 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 86 |
| Other vegetables | 104060 | Mixed vegetable intake | Data-Coding 100017 | 90 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 81 |
| Other vegetables | 104070 | Vegetable pieces intake | Data-Coding 100017 | 90 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 81 |
| Other vegetables | 104080 | Coleslaw intake | Data-Coding 100017 | 90 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 81 |
| Other vegetables | 104090 | Side salad intake | Data-Coding 100017 | 66 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 81 |
| Other vegetables | 104100 | Avocado intake | Data-Coding 100017 | 136 | Quarter of a medium avocado; Half a medium avocado; 1 medium avocado; 2 medium avocados; 3+ medium avocados | 81 |
| Other vegetables | 104130 | Beetroot intake | Data-Coding 100017 | 48 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 82 |
| Other vegetables | 104180 | Cauliflower intake | Data-Coding 100017 | 90 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 83 |
| Other vegetables | 104190 | Celery intake | Data-Coding 100017 | 30 | Quarter of a stick; Half a stick; 1 stick; 2 sticks; 3+ sticks | 83 |
| Other vegetables | 104200 | Courgette intake | Data-Coding 100017 | 90 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 83 |
| Other vegetables | 104210 | Cucumber intake | Data-Coding 100017 | 60 | Quarter of an inch; Half an inch; 1 inch (2.5 cm); 2 inches; 3+ inches | 83 |
| Other vegetables | 104220 | Garlic intake | Data-Coding 100017 | 5 | Quarter of a clove; Half a clove; 1 clove; 2 cloves; 3+ cloves | 83 |
| Other vegetables | 104230 | Leek intake | Data-Coding 100017 | 80 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 84 |
| Other vegetables | 104260 | Onion intake | Data-Coding 100017 | 60 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 84 |
| Other vegetables | 104270 | Parsnip intake | Data-Coding 100017 | 65 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 84 |
| Other vegetables | 104360 | Turnip/swede intake | Data-Coding 100017 | 60 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 86 |
| Other vegetables | 104380 | Other vegetables intake | Data-Coding 100017 | 60 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 86 |
| Potatoes | 104030 | Boiled/baked potatoes intake | Data-Coding 100017 | 175 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 80 |
| Potatoes | 104050 | Mashed potato intake | Data-Coding 100017 | 180 | Quarter of a serving; Half a serving; 1 serving; 2 servings; 3+ servings | 81 |
| Refined grains | 101230 | Naan bread intake | Data-Coding 100006 | 160 | Half a naan; 1 naan; 2 naan; 3 naan; 4 naan; 5 naan; 6+ naan | 48 |
| Refined grains | 101240 | Garlic bread intake | Data-Coding 100006 | 20 | Half a slice; 1 slice; 2 slices; 3 slices; 4 slices; 5 slices; 6+ slices | 48 |
| Refined grains | 101250 | Crispbread intake | Data-Coding 100006 | 10 | Half a cracker/crispbread/cake; 1 cracker/crispbread/cake; 2 crackers/crispbreads/cakes; 3 crackers/crispbreads/cakes; 4 crackers/crispbreads/cakes; 5 crackers/crispbreads/cakes; 6+ crackers/crispbreads/cakes | 48 |
| Refined grains | 101260 | Oatcakes intake | Data-Coding 100006 | 13 | Half an oatcake; 1 oatcake; 2 oatcakes; 3 oatcakes; 4 oatcakes; 5 oatcakes; 6+ oatcakes | 48 |
| Refined grains | 101270 | Other bread intake | Data-Coding 100006 | 45 | Half a slice/wrap/stick/item; 1 slice/wrap/stick/item; 2 slices/wraps/sticks/items; 3 slices/wraps/sticks/items; 4 slices/wraps/sticks/items; 5 slices/wraps/sticks/items; 6+ slices/wraps/sticks/items | 49 |
| Refined grains | 102710 | White pasta intake | Data-Coding 100001 | 230 | Half a serving; 1 serving; 2+ servings | 65 |
| Refined grains | 102730 | White rice intake | Data-Coding 100001 | 150 | Half a serving; 1 serving; 2+ servings | 66 |
| Refined grains | 102750 | Sushi intake | Data-Coding 100001 | 278 | Half a sushi; 1 sushi; 2+ sushi | 66 |
| Refined grains | 102760 | Snackpot intake | Data-Coding 100001 | 280 | Half a pot; 1 pot; 2+ pots | 66 |
| Refined grains | 102770 | Couscous intake | Data-Coding 100001 | 150 | Half a serving; 1 serving; 2+ servings | 66 |
| Refined grains | 101160 | Bread roll intake (white) | Data-Coding 100006 | 60 | Half a roll/bap/bun/bagel; 1 roll/bap/bun/bagel; 2 rolls/baps/buns/bagels; 3 rolls/baps/buns/bagels; 4 rolls/baps/buns/bagels; 5 rolls/baps/buns/bagels; 6+ rolls/baps/buns/bagels | 47 |
| Eggs | 102940 | Whole egg intake | Data-Coding 100003 | 50 | Half an egg; 1 egg; 2 eggs; 3+ eggs | 69 |
| Eggs | 102950 | Omelette intake | Data-Coding 100003 | 120 | Half a serving; 1 serving; 2 servings; 3+ servings | 69 |
| Eggs | 102960 | Eggs in sandwiches intake | Data-Coding 100003 | 120 | Half a serving; 1 serving; 2 servings; 3+ servings | 69 |
| Eggs | 102970 | Scotch egg intake | Data-Coding 100003 | 120 | Half a Scotch egg; 1 Scotch egg; 2 Scotch eggs; 3+ Scotch eggs | 69 |
| Eggs | 102980 | Other egg intake | Data-Coding 100003 | 75 | Half a serving; 1 serving; 2 servings; 3+ servings | 69 |
| Low-fat dairy products | 26102 | Low fat yogurt | Showcase continuous g/day | 125 | Enter grams/day; no categorical response options | Showcase field metadata |
| Low-fat dairy products | 26131 | Semi skimmed milk | Showcase continuous g/day | 259 | Enter grams/day; no categorical response options | Showcase field metadata |
| Low-fat dairy products | 26133 | Skimmed milk and cholesterol-lowering milk | Showcase continuous g/day | 259 | Enter grams/day; no categorical response options | Showcase field metadata |
| Low-fat dairy products | 102810 | Low fat hard cheese intake | Data-Coding 100003 | 40 | Half a serving; 1 serving; 2 servings; 3+ servings | 68 |
| Low-fat dairy products | 102850 | Low fat cheese spread intake | Data-Coding 100003 | 15 | Half a serving; 1 serving; 2 servings; 3+ servings | 68 |
| Low-fat dairy products | 102870 | Cottage cheese intake | Data-Coding 100003 | 60 | Half a serving; 1 serving; 2 servings; 3+ servings | 68 |
| Red meats | 103020 | Beef intake | Data-Coding 100016 | 120 | Half a serving; 1 serving; 2 servings; 3 servings; 4 servings; 5+ servings | 70 |
| Red meats | 103030 | Pork intake | Data-Coding 100016 | 120 | Half a serving; 1 serving; 2 servings; 3 servings; 4 servings; 5+ servings | 70 |
| Red meats | 103040 | Lamb intake | Data-Coding 100016 | 120 | Half a serving; 1 serving; 2 servings; 3 servings; 4 servings; 5+ servings | 71 |
| Processed meats | 103010 | Sausage intake | Data-Coding 100016 | 30 | Half a sausage; 1 sausage; 2 sausages; 3 sausages; 4 sausages; 5+ sausages | 70 |
| Processed meats | 103070 | Bacon intake | Data-Coding 100016 | 46 | Half a rasher; 1 rasher; 2 rashers; 3 rashers; 4 rashers; 5+ rashers | 72 |
| Processed meats | 103080 | Ham intake | Data-Coding 100016 | 23 | Half a slice; 1 slice; 2 slices; 3 slices; 4 slices; 5+ slices | 72 |
| Fish and other seafood | 103150 | Tinned tuna intake | Data-Coding 100004 | 92 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 74 |
| Fish and other seafood | 103160 | Oily fish intake | Data-Coding 100004 | 100 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 74 |
| Fish and other seafood | 103190 | White fish intake | Data-Coding 100004 | 120 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 74 |
| Fish and other seafood | 103230 | Other fish intake | Data-Coding 100004 | 100 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 75 |
| Fish and other seafood | 103200 | Prawns intake | Data-Coding 100004 | 60 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 74 |
| Fish and other seafood | 103210 | Lobster/crab intake | Data-Coding 100004 | 85 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 75 |
| Fish and other seafood | 103220 | Shellfish intake | Data-Coding 100004 | 40 | Half a serving; 1 serving; 2 servings; 3 servings; 4+ servings | 75 |
| Coffee | 100250 | Instant coffee intake | Data-Coding 100006 | 190 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 23 |
| Coffee | 100270 | Filtered coffee intake | Data-Coding 100006 | 190 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 24 |
| Coffee | 100290 | Cappuccino intake | Data-Coding 100006 | 30 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 25 |
| Coffee | 100300 | Latte intake | Data-Coding 100006 | 30 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 25 |
| Coffee | 100310 | Espresso intake | Data-Coding 100006 | 30 | Half a cup; 1 cup; 2 cups; 3 cups; 4 cups; 5 cups; 6+ cups | 26 |
| Coffee | 100330 | Other coffee type | Data-Coding 100006 | 190 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 27 |
| Fruit juices | 100190 | Orange juice intake | Data-Coding 100006 | 250 | Half a glass/carton; 1 glass/carton (250 ml); 2 glasses/cartons; 3 glasses/cartons; 4 glasses/cartons; 5 glasses/cartons; 6+ glasses/cartons | 33 |
| Fruit juices | 100200 | Grapefruit juice intake | Data-Coding 100006 | 250 | Half a glass/carton; 1 glass/carton (250 ml); 2 glasses/cartons; 3 glasses/cartons; 4 glasses/cartons; 5 glasses/cartons; 6+ glasses/cartons | 33 |
| Fruit juices | 100210 | Pure fruit/vegetable juice intake | Data-Coding 100006 | 250 | Half a glass/carton; 1 glass/carton (250 ml); 2 glasses/cartons; 3 glasses/cartons; 4 glasses/cartons; 5 glasses/cartons; 6+ glasses/cartons | 34 |
| Sweetened beverages | 100170 | Fizzy drink intake | Data-Coding 100006 | 330 | Half a glass/can; 1 glass/can; 2 glasses/cans; 3 glasses/cans; 4 glasses/cans; 5 glasses/cans; 6+ glasses/cans | 33 |
| Sweetened beverages | 100180 | Squash intake | Data-Coding 100006 | 250 | Half a glass/carton; 1 glass/carton (250 ml); 2 glasses/cartons; 3 glasses/cartons; 4 glasses/cartons; 5 glasses/cartons; 6+ glasses/cartons | 33 |
| Sweetened beverages | 100220 | Fruit smoothie intake | Data-Coding 100006 | 260 | Half a glass/bottle; 1 glass/bottle (250 ml); 2 glasses/bottles; 3 glasses/bottles; 4 glasses/bottles; 5 glasses/bottles; 6+ glasses/bottles | 34 |
| Sweetened beverages | 100230 | Dairy smoothie intake | Data-Coding 100006 | 260 | Half a glass/bottle; 1 glass/bottle (250 ml); 2 glasses/bottles; 3 glasses/bottles; 4 glasses/bottles; 5 glasses/bottles; 6+ glasses/bottles | 34 |
| Sweetened beverages | 100550 | Hot chocolate intake | Data-Coding 100006 | 143 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 35 |
| Sweetened beverages | 100160 | Low calorie drink intake | Data-Coding 100006 | 330 | Half a glass/can; 1 glass/can; 2 glasses/cans; 3 glasses/cans; 4 glasses/cans; 5 glasses/cans; 6+ glasses/cans | 32 |
| Sweetened beverages | 100540 | Low calorie hot chocolate intake | Data-Coding 100006 | 143 | Half a mug/cup; 1 mug/cup; 2 mugs/cups; 3 mugs/cups; 4 mugs/cups; 5 mugs/cups; 6+ mugs/cups | 35 |
| Water | 100150 | Drinking water intake | Data-Coding 100006 | 250 | Half a glass/beaker; 1 glass/beaker (250 ml); 2 glasses/beakers; 3 glasses/beakers; 4 glasses/beakers; 5 glasses/beakers; 6+ glasses/beakers | 32 |
| Wine | 26151 | Fortified wine intake | Showcase continuous g/day | 50 | Enter grams/day; no categorical response options | Showcase field metadata |
| Wine | 26152 | Red wine intake | Showcase continuous g/day | 175 | Enter grams/day; no categorical response options | Showcase field metadata |
| Wine | 26153 | White wine intake | Showcase continuous g/day | 175 | Enter grams/day; no categorical response options | Showcase field metadata |

## Important discrepancy: Data-Coding 100001

The Oxford WebQ v1.5 PDF appendix displays four response labels for the five grain fields in this component set: half, 1, 2 and 3+. However, the UKB Showcase Data-Coding 100001 attached to the released fields lists only three stored categories: `555 = half`, `1 = 1`, and `200 = 2+`. The frozen RWDID R code maps `200` to `2.0`. The website therefore presents `2+` as one category and maps it to `2.0`; it does not create a separate 3+ category that cannot be recovered from the released coded field.

Relevant Showcase fields: [102710 White pasta](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102710), [102730 White rice](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102730), [102750 Sushi](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102750), [102760 Snackpot](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102760), [102770 Couscous](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102770).

## Field ID discrepancy: Other egg dishes

The Oxford WebQ v1.5 PDF appendix prints **120980** for “Other egg dishes”. The UKB Showcase released data field is **102980** (“Other egg intake”), and it is included in Data-Coding 100003. The Table S1 mapping and the frozen research data pipeline also use 102980. The calculator therefore uses **102980**, matching the released Showcase field rather than the apparent PDF typographical discrepancy.

Showcase reference: [Data-Field 102980 — Other egg intake](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=102980). The PDF discrepancy is visible on page 69 of the [Oxford WebQ v1.5 questionnaire](https://biobank.ctsu.ox.ac.uk/showcase/showcase/docs/DietWebQ.pdf).

## Six continuous Showcase fields

The following six fields belong to UKB Showcase's **Total weight by food group yesterday** category and have value type **Continuous, g**. They are not checkbox-and-category WebQ item fields. The website asks for grams/day and converts to standard portions by dividing by the Table S1 portion size.

| Field ID | Showcase field | Type | Table S1 portion (g) |
|---:|---|---|---:|
| 26102 | Low fat yogurt | Continuous, g | 125 |
| 26131 | Semi skimmed milk | Continuous, g | 259 |
| 26133 | Skimmed milk and cholesterol-lowering milk | Continuous, g | 259 |
| 26151 | Fortified wine | Continuous, g | 50 |
| 26152 | Red wine | Continuous, g | 175 |
| 26153 | White wine | Continuous, g | 175 |

Showcase references: [26102](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26102), [26131](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26131), [26133](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26133), [26151](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26151), [26152](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26152), [26153](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26153).

## Recall-period note

The original Oxford WebQ asks about the previous 24 hours. The study pipeline used available repeated recall instances and averaged item-level values before summing each food group. The browser calculator accepts one representative-day response per item and does not reproduce the repeated-instance averaging or perform imputation.
