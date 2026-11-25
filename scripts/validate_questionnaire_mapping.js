#!/usr/bin/env node
"use strict";

const fs = require("fs");
const path = require("path");
const root = path.resolve(__dirname, "..");
const mapping = JSON.parse(fs.readFileSync(path.join(root, "data/field_mapping_active82.json"), "utf8"));
const params = JSON.parse(fs.readFileSync(path.join(root, "data/rwdid_parameters.json"), "utf8"));
const palette = JSON.parse(fs.readFileSync(path.join(root, "data/friendly13_palette.json"), "utf8"));
const scoring = require(path.join(root, "assets/scoring.js"));

function fail(message) {
  console.error(`FAIL: ${message}`);
  process.exit(1);
}
function assert(condition, message) {
  if (!condition) fail(message);
}
function sameArray(a, b) {
  return JSON.stringify(a) === JSON.stringify(b);
}

assert(Array.isArray(mapping) && mapping.length === 82, "Expected exactly 82 active mapping rows.");
assert(new Set(mapping.map(x => Number(x.field_id))).size === 82, "Field IDs must be unique.");
assert(new Set(mapping.map(x => x.component_key)).size === 15, "Expected 15 dietary components.");
assert(Array.isArray(params.components) && params.components.length === 15, "Expected 15 frozen score components.");

const categorical = mapping.filter(x => Number(x.encoding_id) !== 0);
const continuous = mapping.filter(x => Number(x.encoding_id) === 0);
assert(categorical.length === 76, "Expected 76 categorical Oxford WebQ fields.");
assert(continuous.length === 6, "Expected 6 continuous UKB Showcase gram fields.");

const expectedCodes = {
  100017: [["444", 0.25], ["555", 0.5], ["1", 1], ["2", 2], ["300", 3]],
  100004: [["555", 0.5], ["1", 1], ["2", 2], ["3", 3], ["400", 4]],
  100006: [["555", 0.5], ["1", 1], ["2", 2], ["3", 3], ["4", 4], ["5", 5], ["600", 6]],
  100003: [["555", 0.5], ["1", 1], ["2", 2], ["300", 3]],
  100016: [["555", 0.5], ["1", 1], ["2", 2], ["3", 3], ["4", 4], ["500", 5]],
  100001: [["555", 0.5], ["1", 1], ["200", 2]]
};

for (const item of categorical) {
  const coding = Number(item.encoding_id);
  const expected = expectedCodes[coding];
  const options = item.encoding_options;
  assert(expected, `Unexpected Data-Coding ${coding} on Field ID ${item.field_id}.`);
  assert(Array.isArray(options) && options.length === expected.length,
    `Incorrect response-option count for Field ID ${item.field_id}.`);
  assert(sameArray(options.map(o => [String(o.code), Number(o.value)]), expected),
    `Incorrect stored-code to portion conversion for Field ID ${item.field_id}.`);
  assert(options.every(o => String(o.code) !== "0"), `Synthetic zero category found for Field ID ${item.field_id}.`);
  assert(options.every(o => typeof o.label === "string" && o.label.length > 0),
    `Missing item-specific UKB response wording for Field ID ${item.field_id}.`);
}

const continuousIds = continuous.map(x => Number(x.field_id)).sort((a,b) => a-b);
const expectedContinuousIds = [26102,26131,26133,26151,26152,26153];
assert(sameArray(continuousIds, expectedContinuousIds),
  "Continuous Showcase field IDs do not match the frozen six-field mapping.");

const grains = mapping.filter(x => Number(x.encoding_id) === 100001);
assert(grains.length === 5, "Expected five fields using Data-Coding 100001.");
assert(grains.every(x => x.coding_discrepancy_note &&
  x.coding_discrepancy_note.includes("200=2+") &&
  x.encoding_options.some(o => String(o.code) === "200" && Number(o.value) === 2)),
  "Data-Coding 100001 discrepancy must be documented and mapped as 2+ -> 2.0.");

const derivedBread = mapping.filter(x => x.is_derived_bread_item === true);
assert(derivedBread.length === 1, "Expected exactly one active derived white-bread item.");
const bread = derivedBread[0];
assert(Number(bread.field_id) === 101160, "Derived bread amount must be Field 101160.");
assert(Number(bread.selector_field_id) === 20094, "Derived bread type selector must be Field 20094.");
assert(Number(bread.selector_field_encoding_id) === 81, "Bread type selector must use Data-Coding 81.");
assert(bread.selector_selection_mode === "multiple", "Field 20094 must be presented as multiple response.");
assert(bread.response_control === "checkbox_then_select_plus_dependent_multiselect",
  "Derived bread roll must expose the dependent Field 20094 selector in the UI.");
assert(sameArray(
  bread.selector_response_options.map(o => [String(o.code), o.label]),
  [
    ["1", "White"],
    ["2", "Granary, brown, mixed flours/grains, white with extras (e.g. Best of Both, 50/50)"],
    ["3", "Wholemeal"],
    ["5", "Other flour - please describe"],
    ["4", "Seeded (inside and/or on top)"]
  ]
), "Field 20094 selector options do not match the Oxford WebQ / Showcase mapping.");
assert(scoring.isDerivedItemIncluded(bread, ["1"]) === true, "White-only bread selector should be included.");
assert(scoring.isDerivedItemIncluded(bread, ["1", "3"]) === false, "Multiple bread types must not be treated as White-only.");
assert(scoring.isDerivedItemIncluded(bread, ["3"]) === false, "Non-white bread must not be included in the white-bread row.");

const paletteColors = [
  "#E73F74", "#F1CE63", "#77AADD", "#009988", "#9467BD",
  "#FF9D9A", "#99DDFF", "#AAAA00", "#225522", "#882255",
  "#997700", "#4B4B8F", "#8C564B"
];
assert(palette.palette === "cols4all::friendly13", "Incorrect palette identity.");
assert(String(palette.package_version) === "0.11", "Expected cols4all 0.11.");
assert(palette.generated === true && Array.isArray(palette.colors) && palette.colors.length === 13,
  "Friendly13 asset must contain 13 generated colors.");
assert(sameArray(palette.colors.map(x => String(x).toUpperCase()), paletteColors),
  "Friendly13 HEX values or their order differ from the verified cols4all output.");
assert(new Set(palette.colors.map(x => String(x).toUpperCase())).size === 13,
  "Friendly13 colors must be unique.");
assert(palette.colors.every(x => /^#[0-9A-F]{6}$/i.test(x)), "Friendly13 contains an invalid HEX color.");

const betaSum = params.components.reduce((sum, c) => sum + Number(c.beta), 0);
assert(Math.abs(betaSum - Number(params.beta_sum)) < 1e-10,
  "The beta_sum denominator does not equal the sum of the 15 published beta coefficients.");
assert(params.components.every(c => Number(c.beta) < 0),
  "The current frozen RWDID beta parameter set is expected to contain 15 negative coefficients.");
assert(sameArray(
  params.components.map(c => c.food),
  [...new Set(mapping.slice().sort((a,b) => Number(a.component_order)-Number(b.component_order))
    .map(x => x.component_key))]
), "Mapping component order does not match the frozen score parameter order.");

const zeroIntakes = Object.fromEntries(params.components.map(x => [x.food, 0]));
const zeroSmoke = scoring.calculateScores(zeroIntakes, params);
assert(zeroSmoke.rwdid >= 0 && zeroSmoke.rwdid <= 1, "RWDID zero-intake smoke score outside 0–1.");
assert(zeroSmoke.rwdid_c >= 0 && zeroSmoke.rwdid_c <= 15, "RWDID-C zero-intake smoke score outside 0–15.");

console.log("PASS: 82 active fields = 76 categorical Oxford WebQ + 6 continuous Showcase g/day.");
console.log("PASS: all six Data-Coding option sets and stored-code-to-portion conversions.");
console.log("PASS: derived Field 101160 correctly depends on multiple-response Field 20094.");
console.log("PASS: 100001 grain discrepancy is preserved as Showcase 2+ -> 2.0.");
console.log("PASS: exact ordered cols4all::friendly13 HEX vector.");
console.log(`PASS: 15 frozen beta coefficients reconcile to beta_sum ${params.beta_sum}.`);
console.log("PASS: scoring-engine smoke test.");
