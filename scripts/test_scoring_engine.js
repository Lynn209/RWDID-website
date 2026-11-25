#!/usr/bin/env node
"use strict";

const fs = require("fs");
const path = require("path");
const root = path.resolve(__dirname, "..");
const mapping = JSON.parse(fs.readFileSync(path.join(root, "data/field_mapping_active82.json"), "utf8"));
const params = JSON.parse(fs.readFileSync(path.join(root, "data/rwdid_parameters.json"), "utf8"));
const scoring = require(path.join(root, "assets/scoring.js"));

const chosenCodes = {
  104530: "1",    // orange: 1 fruit
  104170: "1",    // carrot: 1 serving
  104060: "2",    // mixed vegetables: 2 servings
  104030: "1",    // boiled/baked potato: 1 serving
  104050: "444",  // mashed potato: 1/4 serving; total potatoes = 1.25
  101230: "3",    // naan: 3 items/servings
  102940: "1",    // whole egg: 1 egg
  103020: "1",    // beef: 1 serving
  103010: "3",    // sausage: 3 units; algorithmic target for processed meats >= 2.7
  103160: "1",    // oily fish: 1 serving
  100250: "2",    // instant coffee: 2 mugs/cups
  100190: "1",    // orange juice: 1 glass/carton
  100150: "2"     // water: 2 glasses/beakers
};

const gramsByField = {
  26102: 125,     // low-fat yogurt: 1 Table S1 portion
  26131: 0,       // no semi-skimmed milk
  26133: 0,       // no skimmed milk
  26151: 0,       // no fortified wine
  26152: 0,       // no red wine
  26153: 0        // no white wine
};

const groupIntakes = {};
for (const item of mapping) {
  const key = item.component_key;
  if (!(key in groupIntakes)) groupIntakes[key] = 0;

  let contribution = 0;
  if (Number(item.encoding_id) === 0) {
    const grams = gramsByField[Number(item.field_id)];
    if (grams === undefined) throw new Error(`Missing synthetic grams input for Field ${item.field_id}`);
    contribution = scoring.decodeItemValue(item, String(grams));
  } else if (Object.prototype.hasOwnProperty.call(chosenCodes, Number(item.field_id))) {
    const rawCode = chosenCodes[Number(item.field_id)];
    contribution = scoring.decodeItemValue(item, rawCode);
  }

  groupIntakes[key] += contribution;
}

const result = scoring.calculateScores(groupIntakes, params);
const expected = {
  citrus_fruits: 1.0,
  red_orange_veg: 1.0,
  other_veg: 2.0,
  potatoes: 1.25,
  refined_grains: 3.0,
  eggs: 1.0,
  low_fat_dairy: 1.0,
  red_meats: 1.0,
  processed_meats: 3.0,
  fish_seafood: 1.0,
  coffee: 2.0,
  fruit_juices: 1.0,
  sweetened_bev: 0.0,
  water: 2.0,
  wine: 0.0
};

for (const [food, value] of Object.entries(expected)) {
  if (Math.abs(groupIntakes[food] - value) > 1e-10) {
    throw new Error(`${food}: expected ${value}, got ${groupIntakes[food]}`);
  }
}
if (Math.abs(result.rwdid - 1) > 1e-10) throw new Error(`Expected RWDID=1.0, got ${result.rwdid}`);
if (Math.abs(result.rwdid_c - 15) > 1e-10) throw new Error(`Expected RWDID-C=15.0, got ${result.rwdid_c}`);

const bread = mapping.find(x => x.is_derived_bread_item === true);
const whiteOnly = scoring.isDerivedItemIncluded(bread, ["1"])
  ? scoring.decodeItemValue(bread, "1") : 0;
const whiteAndWholemeal = scoring.isDerivedItemIncluded(bread, ["1", "3"])
  ? scoring.decodeItemValue(bread, "1") : 0;
const wholemealOnly = scoring.isDerivedItemIncluded(bread, ["3"])
  ? scoring.decodeItemValue(bread, "1") : 0;
if (whiteOnly !== 1 || whiteAndWholemeal !== 0 || wholemealOnly !== 0) {
  throw new Error("Derived bread-roll selector does not reproduce the frozen White-only inclusion rule.");
}

console.log("PASS: synthetic 15-component target-profile input is aggregated from the 82-field map.");
console.log("PASS: component intakes match all 15 expected frozen target values.");
console.log(`PASS: RWDID=${result.rwdid.toFixed(4)}; RWDID-C=${result.rwdid_c.toFixed(2)}.`);
console.log("PASS: bread-roll amount is included only for a White-only Field 20094 response.");
console.log("NOTE: This is a deterministic algorithm test, not dietary advice or a clinical risk assessment.");
