(() => {
  "use strict";

  const $ = (s, root = document) => root.querySelector(s);
  const $$ = (s, root = document) => Array.from(root.querySelectorAll(s));

  let items = [];
  let params = null;
  let latestResult = null;

  const escapeHTML = v => String(v ?? "").replace(/[&<>"']/g, c => ({
    "&":"&amp;", "<":"&lt;", ">":"&gt;", '"':"&quot;", "'":"&#39;"
  }[c]));

  const fmt = (x, d = 3) => Number(x).toLocaleString(undefined, {
    minimumFractionDigits: d, maximumFractionDigits: d
  });

  function setStatus(message, kind = "") {
    const el = $("#form-status");
    if (!el) return;
    el.textContent = message;
    el.className = `status-text ${kind}`.trim();
  }

  async function loadJSON(path) {
    const response = await fetch(path);
    if (!response.ok) throw new Error(`Could not load ${path} (${response.status}).`);
    return response.json();
  }

  function foodGroups() {
    const map = new Map();
    items.forEach(item => {
      if (!map.has(item.component_key)) {
        map.set(item.component_key, {
          key: item.component_key,
          order: Number(item.component_order),
          name: item.component_name,
          items: []
        });
      }
      map.get(item.component_key).items.push(item);
    });
    return [...map.values()].sort((a, b) => a.order - b.order);
  }

  function cleanNote(item) {
    let note = String(item.condition_note || item.field_notes || "");
    const prompt = String(item.question_prompt || "");
    if (prompt) note = note.replace(`Question asked: "${prompt}"`, "").trim();
    const help = String(item.questionnaire_help_note || "").trim();
    return [note, help].filter(Boolean).join(" ");
  }

  function renderField(item) {
    const id = `field-${item.active_item_id}`;
    const itemId = escapeHTML(item.active_item_id);
    const isContinuous = Number(item.encoding_id) === 0;
    const isDerivedBread = item.is_derived_bread_item === true;

    if (isContinuous) {
      return `<div class="item-response continuous-response">
        <div class="amount-control">
          <label class="amount-label" for="${id}">Amount consumed (g/day)</label>
          <input class="answer response-input" id="${id}"
            data-item-id="${itemId}" type="number" min="0" step="any"
            inputmode="decimal" placeholder="Enter 0 if none">
        </div>
      </div>`;
    }

    const options = (item.encoding_options || []).map(o =>
      `<option value="${escapeHTML(o.code)}">${escapeHTML(o.label)}</option>`
    ).join("");

    const amountControl = `<div class="amount-control" hidden>
      <label class="amount-label" for="${id}">Amount consumed</label>
      <select class="answer response-input" id="${id}"
        data-item-id="${itemId}" disabled>
        <option value="">Select UKB response</option>
        ${options}
      </select>
    </div>`;

    if (isDerivedBread) {
      const selectorOptions = (item.selector_response_options || []).map(option => `
        <label class="selector-option">
          <input type="checkbox" class="bread-type-check"
            data-item-id="${itemId}" value="${escapeHTML(option.code)}" disabled>
          <span>${escapeHTML(option.label)}</span>
        </label>
      `).join("");

      return `<div class="item-response derived-bread-response">
        <label class="item-consumed">
          <input type="checkbox" class="item-consumed-check" data-item-id="${itemId}">
          <span>Consumed yesterday</span>
        </label>
        ${amountControl}
        <div class="selector-control" hidden>
          <span class="amount-label">Bread type · Field ID ${escapeHTML(item.selector_field_id)} · Select all that apply</span>
          <div class="selector-option-list">${selectorOptions}</div>
          <p class="field-caution">The frozen research pipeline includes this bread-roll amount in the refined-grains component only when the recorded bread type is White alone (code 1). If multiple types are selected, this row contributes zero, matching the original analysis code.</p>
        </div>
      </div>`;
    }

    return `<div class="item-response">
      <label class="item-consumed">
        <input type="checkbox" class="item-consumed-check" data-item-id="${itemId}">
        <span>Consumed yesterday</span>
      </label>
      ${amountControl}
    </div>`;
  }

  function renderQuestionnaire() {
    const host = $("#questionnaire-container");
    host.innerHTML = "";

    foodGroups().forEach((group, index) => {
      const details = document.createElement("details");
      details.className = "food-group";
      details.dataset.group = group.key;
      if (index === 0) details.open = true;
      details.addEventListener("toggle", () => {
        if (details.open) {
          $$(".food-group").forEach(other => {
            if (other !== details) other.open = false;
          });
        }
      });

      const paletteIndex = (index % 13) + 1;
      details.style.setProperty("--group-color", `var(--friendly-${paletteIndex})`);
      details.style.setProperty("--group-contrast", `var(--friendly-${paletteIndex}-contrast)`);

      const fields = group.items.map(item => {
        const isContinuous = Number(item.encoding_id) === 0;
        const prompt = item.question_prompt || (
          isContinuous
            ? "Enter the grams consumed per day for this UK Biobank continuous total-weight field."
            : `How much ${item.food_item} did you consume yesterday?`
        );
        const note = cleanNote(item);
        const portion = item.portion_size_g == null ? "" : ` · Table S1 portion ${fmt(item.portion_size_g, 1)} g`;
        const meta = isContinuous
          ? `UKB Showcase Field ID ${escapeHTML(item.table_s1_field_id)} · Continuous g/day${portion}`
          : item.is_derived_bread_item
            ? `UKB Field ID ${escapeHTML(item.field_id)} · type selector ${escapeHTML(item.selector_field_id)} · Data-Coding ${escapeHTML(item.encoding_id)}${portion}`
            : `Field ID ${escapeHTML(item.table_s1_field_id)} · Data-Coding ${escapeHTML(item.encoding_id)}${portion}`;
        const codingNote = Number(item.encoding_id) === 100001
          ? `<p class="field-source-note">Showcase Data-Coding 100001 stores the upper category as 2+ (code 200); the frozen RWDID pipeline maps it to 2.0. This calculator follows that released coding for all five fields using this encoding.</p>`
          : "";
        const questionnaireHelp = item.questionnaire_help_note
          ? `<p class="field-caution">${escapeHTML(item.questionnaire_help_note)}</p>`
          : "";

        return `<div class="field-card">
          <div class="field-label">
            ${escapeHTML(item.food_item)}
            ${note ? `<span class="field-info" tabindex="0" title="${escapeHTML(note)}" aria-label="Additional question information">i</span>` : ""}
          </div>
          <p class="field-prompt">${escapeHTML(prompt)}</p>
          ${renderField(item)}
          <small class="field-meta">${meta}</small>
          ${questionnaireHelp}
          ${codingNote}
        </div>`;
      }).join("");

      details.innerHTML = `
        <summary>
          <span class="group-index">${String(group.order).padStart(2, "0")}</span>
          <span class="group-title">${escapeHTML(group.name)}
            <small>${group.items.length} dietary fields</small>
          </span>
          <span class="group-progress" data-group-status>Incomplete</span>
          <span class="chevron" aria-hidden="true">⌄</span>
        </summary>
        <div class="group-content">
          <div class="group-tools">
            <label class="no-intake">
              <input type="checkbox" class="no-intake-check">
              No intake from this food group
            </label>
            <span class="group-hint">Tick categorical items; enter grams/day for continuous fields</span>
            <button type="button" class="clear-group">Clear</button>
          </div>
          <div class="field-grid">${fields}</div>
        </div>`;
      host.appendChild(details);
    });

    $$(".response-input").forEach(input => {
      input.addEventListener("input", updateProgress);
      input.addEventListener("change", updateProgress);
    });

    $$(".item-consumed-check").forEach(check => check.addEventListener("change", e => {
      const card = e.target.closest(".field-card");
      const group = e.target.closest(".food-group");
      const input = $(".response-input", card);
      const amount = $(".amount-control", card);
      const selector = $(".selector-control", card);

      if (e.target.checked) {
        $(".no-intake-check", group).checked = false;
        input.disabled = false;
        amount.hidden = false;
        if (selector) {
          selector.hidden = false;
          $$(".bread-type-check", card).forEach(typeCheck => {
            typeCheck.disabled = false;
          });
        }
      } else {
        input.value = "";
        input.disabled = true;
        amount.hidden = true;
        if (selector) {
          selector.hidden = true;
          $$(".bread-type-check", card).forEach(typeCheck => {
            typeCheck.checked = false;
            typeCheck.disabled = true;
          });
        }
      }
      updateProgress();
    }));

    $$(".bread-type-check").forEach(check => {
      check.addEventListener("change", updateProgress);
    });

    $$(".no-intake-check").forEach(check => check.addEventListener("change", e => {
      const group = e.target.closest(".food-group");
      const noIntake = e.target.checked;

      $$(".item-consumed-check", group).forEach(itemCheck => {
        itemCheck.checked = false;
        itemCheck.disabled = noIntake;
      });

      $$(".bread-type-check", group).forEach(typeCheck => {
        typeCheck.checked = false;
        typeCheck.disabled = true;
      });

      $$(".selector-control", group).forEach(selector => {
        selector.hidden = true;
      });

      $$(".response-input", group).forEach(input => {
        const item = items.find(x => String(x.active_item_id) === String(input.dataset.itemId));
        const amount = $(".amount-control", input.closest(".field-card"));
        if (noIntake) {
          input.value = "";
          input.disabled = true;
          amount.hidden = true;
        } else if (item && Number(item.encoding_id) === 0) {
          input.disabled = false;
          amount.hidden = false;
        } else {
          input.value = "";
          input.disabled = true;
          amount.hidden = true;
        }
      });
      updateProgress();
    }));

    $$(".clear-group").forEach(button => button.addEventListener("click", e => {
      const group = e.target.closest(".food-group");
      $(".no-intake-check", group).checked = false;

      $$(".item-consumed-check", group).forEach(check => {
        check.checked = false;
        check.disabled = false;
      });

      $$(".bread-type-check", group).forEach(check => {
        check.checked = false;
        check.disabled = true;
      });
      $$(".selector-control", group).forEach(selector => { selector.hidden = true; });

      $$(".response-input", group).forEach(input => {
        const item = items.find(x => String(x.active_item_id) === String(input.dataset.itemId));
        const amount = $(".amount-control", input.closest(".field-card"));
        input.value = "";
        if (item && Number(item.encoding_id) === 0) {
          input.disabled = false;
          amount.hidden = false;
        } else {
          input.disabled = true;
          amount.hidden = true;
        }
      });
      updateProgress();
    }));
  }

  function isGroupComplete(el) {
    if ($(".no-intake-check", el).checked) return true;

    let anyIntake = false;
    for (const item of items.filter(x => x.component_key === el.dataset.group)) {
      const input = $(`.response-input[data-item-id="${item.active_item_id}"]`, el);
      if (Number(item.encoding_id) === 0) {
        if (!input || input.value === "") return false;
        const grams = Number(input.value);
        if (!Number.isFinite(grams) || grams < 0) return false;
        if (grams > 0) anyIntake = true;
      } else {
        const check = $(`.item-consumed-check[data-item-id="${item.active_item_id}"]`, el);
        if (check && check.checked) {
          if (!input || input.value === "") return false;
          if (item.is_derived_bread_item) {
            const card = input.closest(".field-card");
            const selectedTypes = $$(".bread-type-check:checked", card).map(x => x.value);
            if (selectedTypes.length === 0) return false;
          }
          anyIntake = true;
        }
      }
    }
    return anyIntake;
  }

  function updateProgress() {
    const groups = $$(".food-group");
    const done = groups.filter(isGroupComplete).length;
    const pct = groups.length ? Math.round(done / groups.length * 100) : 0;
    $("#progress-count").textContent = `${done} / ${groups.length} completed`;
    $("#progress-fill").style.width = `${pct}%`;
    $("#progress-percent").textContent = `${pct}%`;
    $("#calculate-button").disabled = done !== 15;
    $("#progress-hint").textContent = done === 15
      ? "All 15 components are complete. You can calculate both scores."
      : `${groups.length - done} groups remaining`;
    groups.forEach(group => {
      const complete = isGroupComplete(group);
      const status = $("[data-group-status]", group);
      status.textContent = complete ? "Complete" : "Incomplete";
      status.classList.toggle("complete", complete);
    });
    setStatus(done === 15
      ? "All required components are complete."
      : "For each component, confirm no intake or complete all gram fields and select amounts for consumed categorical items.");
  }

  function collectIntakes() {
    const intakes = {};
    for (const group of foodGroups()) {
      const el = $(`.food-group[data-group="${group.key}"]`);
      const noIntake = $(".no-intake-check", el).checked;
      let total = 0;

      for (const item of group.items) {
        const input = $(`.response-input[data-item-id="${item.active_item_id}"]`, el);
        let raw = "0";
        if (!noIntake) {
          if (Number(item.encoding_id) === 0) {
            raw = input.value;
          } else {
            const check = $(`.item-consumed-check[data-item-id="${item.active_item_id}"]`, el);
            if (check && check.checked) {
              if (item.is_derived_bread_item) {
                const card = input.closest(".field-card");
                const selectedTypes = $$(".bread-type-check:checked", card).map(x => x.value);
                const includeAmount = window.RWDIDScoring.isDerivedItemIncluded(item, selectedTypes);
                raw = includeAmount ? input.value : "0";
              } else {
                raw = input.value;
              }
            } else {
              raw = "0";
            }
          }
        }
        const portion = window.RWDIDScoring.decodeItemValue(item, raw);
        if (portion === null) return null;
        total += portion;
      }
      intakes[group.key] = total;
    }
    return intakes;
  }

  function renderResults(result) {
    latestResult = result;
    $("#results-panel").hidden = false;
    $("#rwdid-value").textContent = fmt(result.rwdid, 4);
    $("#rwdidc-value").textContent = fmt(result.rwdid_c, 2);
    $("#rwdid-bar").style.width = `${result.rwdid * 100}%`;
    $("#rwdidc-bar").style.width = `${result.rwdid_c / 15 * 100}%`;

    const tbody = $("#component-results-body");
    tbody.innerHTML = "";
    result.components.forEach((c, i) => {
      const row = document.createElement("tr");
      row.innerHTML = `
        <td><span class="table-dot" style="background:var(--friendly-${(i % 13) + 1})"></span>${escapeHTML(c.display_name)}</td>
        <td>${fmt(c.intake, 3)}</td>
        <td>${c.binary_score}</td>
        <td>${fmt(c.beta, 4)}</td>
        <td>${fmt(c.normalised_weight * 100, 2)}%</td>
        <td>${fmt(c.continuous_score, 3)}</td>`;
      tbody.appendChild(row);
    });
    $("#results-panel").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function downloadResults() {
    if (!latestResult) return;
    const rows = [
      ["component_order", "component", "intake_portions_per_day", "RWDID_binary_score", "beta", "normalised_weight", "RWDID_C_component_score"],
      ...latestResult.components.map(c => [
        c.order, c.food, c.intake, c.binary_score, c.beta, c.normalised_weight, c.continuous_score
      ]),
      [],
      ["RWDID", latestResult.rwdid],
      ["RWDID-C", latestResult.rwdid_c],
      ["beta_sum", latestResult.beta_sum],
      ["parameter_version", params.version]
    ];
    const csv = rows.map(row => row.map(v => `"${String(v ?? "").replace(/"/g, '""')}"`).join(",")).join("\r\n");
    const blob = new Blob(["\ufeff" + csv], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "RWDID_and_RWDID-C_result.csv";
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  async function init() {
    await window.RWDIDPaletteReady;
    const [mapping, parameters] = await Promise.all([
      loadJSON("data/field_mapping_active82.json"),
      loadJSON("data/rwdid_parameters.json")
    ]);
    items = mapping;
    params = parameters;

    if (!Array.isArray(items) || items.length !== 82) {
      throw new Error(`Expected 82 active questionnaire fields; found ${items?.length ?? "none"}.`);
    }
    if (!Array.isArray(params.components) || params.components.length !== 15) {
      throw new Error("Frozen parameter file must contain 15 components.");
    }

    renderQuestionnaire();
    updateProgress();

    $("#calculate-button").addEventListener("click", () => {
      try {
        const intakes = collectIntakes();
        if (!intakes) throw new Error("Complete all questionnaire items before calculating.");
        const result = window.RWDIDScoring.calculateScores(intakes, params);
        setStatus("Both dietary scores have been calculated.", "success");
        renderResults(result);
      } catch (error) {
        setStatus(error.message, "error");
      }
    });

    $("#reset-button").addEventListener("click", () => {
      $$(".no-intake-check").forEach(check => { check.checked = false; });
      $$(".item-consumed-check").forEach(check => {
        check.checked = false;
        check.disabled = false;
      });
      $$(".bread-type-check").forEach(check => {
        check.checked = false;
        check.disabled = true;
      });
      $$(".selector-control").forEach(selector => { selector.hidden = true; });
      $$(".response-input").forEach(input => {
        const item = items.find(x => String(x.active_item_id) === String(input.dataset.itemId));
        const amount = $(".amount-control", input.closest(".field-card"));
        input.value = "";
        if (item && Number(item.encoding_id) === 0) {
          input.disabled = false;
          amount.hidden = false;
        } else {
          input.disabled = true;
          amount.hidden = true;
        }
      });
      $("#results-panel").hidden = true;
      latestResult = null;
      updateProgress();
      $("#calculator").scrollIntoView({ behavior: "smooth", block: "start" });
    });

    $("#download-results").addEventListener("click", downloadResults);
    $("#print-results").addEventListener("click", () => window.print());

    $("#boot-state").hidden = true;
    $("#app-shell").hidden = false;
  }

  document.addEventListener("DOMContentLoaded", () => {
    init().catch(error => {
      const boot = $("#boot-state");
      if (boot) {
        boot.classList.add("boot-error");
        boot.innerHTML = `<strong>Calculator could not be initialized.</strong><p>${escapeHTML(error.message)}</p>`;
      }
      $("#app-shell").hidden = true;
    });
  });
})();
