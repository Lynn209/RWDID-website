#!/usr/bin/env node
"use strict";

const fs = require("fs");
const path = require("path");
const root = path.resolve(__dirname, "..");
const mapping = JSON.parse(fs.readFileSync(path.join(root, "data/field_mapping_active82.json"), "utf8"));
const params = JSON.parse(fs.readFileSync(path.join(root, "data/rwdid_parameters.json"), "utf8"));
const palette = JSON.parse(fs.readFileSync(path.join(root, "data/friendly13_palette.json"), "utf8"));

function fail(message) {
  console.error(`FAIL: ${message}`);
  process.exit(1);
}
function assert(condition, message) {
  if (!condition) fail(message);
}

assert(Array.isArray(mapping) && mapping.length === 82, "Expected exactly 82 active mapping rows.");
assert(new Set(mapping.map(x => Number(x.field_id))).size === 82, "Field IDs must be unique.");
assert(new Set(mapping.map(x => x.component_key)).size === 15, "Expected 15 dietary components.");
assert(Array.isArray(params.components) && params.components.length === 15, "Expected 15 frozen score components.");

const categorical = mapping.filter(x => Number(x.encoding_id) !== 0);
const continuous = mapping.filter(x => Number(x.encoding_id) === 0);
assert(categorical.length === 76, "Expected 76 categorical Oxford WebQ fields.");
assert(continuous.length === 6, "Expected 6 continuous UKB Showcase gram fields.");

const expectedOptionCounts = {
  100017: 5,
  100004: 5,
  100006: 7,
  100003: 4,
  100016: 6,
  100001: 3
};

for (const item of categorical) {
  const coding = Number(item.encoding_id);
  const options = item.encoding_options;
  assert(expectedOptionCounts[coding] !== undefined, `Unexpected Data-Coding ${coding} on Field ID ${item.field_id}.`);
  assert(Array.isArray(options) && options.length === expectedOptionCounts[coding],
    `Incorrect response-option count for Field ID ${item.field_id}.`);
  assert(options.every(o => String(o.code) !== "0"), `Synthetic zero category found for Field ID ${item.field_id}.`);
  assert(options.every(o => typeof o.label === "string" && o.label.length > 0),
    `Missing item-specific UKB response wording for Field ID ${item.field_id}.`);
  assert(options.every(o => Number.isFinite(Number(o.value))),
    `Invalid numeric portion mapping for Field ID ${item.field_id}.`);
}

const continuousIds = continuous.map(x => Number(x.field_id)).sort((a,b) => a-b);
const expectedContinuousIds = [26102,26131,26133,26151,26152,26153];
assert(JSON.stringify(continuousIds) === JSON.stringify(expectedContinuousIds),
  "Continuous Showcase field IDs do not match the frozen six-field mapping.");

const grains = mapping.filter(x => Number(x.encoding_id) === 100001);
assert(grains.length === 5, "Expected five fields using Data-Coding 100001.");
assert(grains.every(x => x.coding_discrepancy_note && x.encoding_options.some(o => String(o.code) === "200" && Number(o.value) === 2)),
  "Data-Coding 100001 2+ -> 2.0 discrepancy must be documented and retained.");

const breadRoll = mapping.find(x => Number(x.field_id) === 101160);
assert(breadRoll, "The active white bread-roll amount field 101160 is missing.");
assert(Number(breadRoll.selector_field_id) === 20094, "Field 101160 must be linked to selector Field 20094.");
assert(Number(breadRoll.selector_encoding_id) === 81, "Bread-roll selector must use Data-Coding 81.");
assert(breadRoll.selector_response_type === "categorical_multiple", "Field 20094 must be treated as multiple response.");
assert(breadRoll.selector_match_rule === "exactly_one_selected_code", "Frozen exact-single White selector rule is missing.");
assert(JSON.stringify((breadRoll.selector_options || []).map(x => [String(x.code), x.label])) === JSON.stringify([
  ["1", "White"], ["2", "Mixed"], ["3", "Wholemeal"], ["4", "Seeded"], ["5", "Other"]
]), "Field 20094 options do not match UKB Data-Coding 81.");

assert(palette.palette === "cols4all::friendly13", "Incorrect palette identity.");
assert(String(palette.package_version) === "0.11", "Expected cols4all 0.11.");
assert(palette.generated === true && Array.isArray(palette.colors) && palette.colors.length === 13,
  "Friendly13 asset must contain 13 generated colors.");
assert(new Set(palette.colors.map(x => String(x).toUpperCase())).size === 13,
  "Friendly13 colors must be unique.");
assert(palette.colors.every(x => /^#[0-9A-F]{6}$/i.test(x)), "Friendly13 contains an invalid HEX color.");

const scoring = require(path.join(root, "assets/scoring.js"));
const zeroIntakes = Object.fromEntries(params.components.map(x => [x.food, 0]));
const smoke = scoring.calculateScores(zeroIntakes, params);
assert(smoke.rwdid >= 0 && smoke.rwdid <= 1, "RWDID smoke score outside 0–1.");
assert(smoke.rwdid_c >= 0 && smoke.rwdid_c <= 15, "RWDID-C smoke score outside 0–15.");

assert(scoring.selectorGatedIntake(breadRoll, "1", ["1"]) === 1,
  "Exact White-only selection must include one roll amount.");
assert(scoring.selectorGatedIntake(breadRoll, "555", ["1"]) === 0.5,
  "Exact White-only selection must preserve the half-roll coding.");
assert(scoring.selectorGatedIntake(breadRoll, "1", ["2"]) === 0,
  "A non-White-only selection must contribute zero to the derived white-bread item.");
assert(scoring.selectorGatedIntake(breadRoll, "1", ["1", "2"]) === 0,
  "White plus another type must contribute zero, matching the frozen exact '[1]' rule.");
assert(scoring.selectorGatedIntake(breadRoll, "1", []) === 0,
  "No selected type must contribute zero to the derived white-bread item.");
assert(scoring.selectorGatedIntake(breadRoll, "", [], true) === 0,
  "Explicit no bread-roll intake must contribute zero.");
let missingWhiteAmountRejected = false;
try { scoring.selectorGatedIntake(breadRoll, "", ["1"]); }
catch (e) { missingWhiteAmountRejected = /Missing amount/.test(e.message); }
assert(missingWhiteAmountRejected, "White-only selection must require an amount response.");

console.log("PASS: 82 fields = 76 categorical WebQ + 6 continuous Showcase g/day.");
console.log("PASS: 15 components, item-specific UKB options, no synthetic zero response.");
console.log("PASS: Data-Coding 100001 discrepancy documented and mapped as frozen R code.");
console.log("PASS: Field 101160 is linked to Field 20094 (Data-Coding 81) with all five multiple-response types.");
console.log("PASS: White-only exact selector gate; non-White and multi-type responses are not counted.");
console.log("PASS: exact cols4all friendly13 asset and scoring-engine smoke test.");
