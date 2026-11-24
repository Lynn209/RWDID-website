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
    return note;
  }

  function renderField(item) {
    const id = `field-${item.active_item_id}`;
    if (Number(item.encoding_id) === 0) {
      return `<input class="answer response-input" id="${id}"
        data-item-id="${escapeHTML(item.active_item_id)}" type="number"
        min="0" step="any" inputmode="decimal" placeholder="Grams per day">`;
    }
    const options = (item.encoding_options || []).map(o =>
      `<option value="${escapeHTML(o.code)}">${escapeHTML(o.label)}</option>`
    ).join("");
    return `<select class="answer response-input" id="${id}"
      data-item-id="${escapeHTML(item.active_item_id)}">
      <option value="">Select response</option>
      <option value="0">0 · Not consumed</option>
      ${options}
    </select>`;
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
        const prompt = item.question_prompt || `Enter the amount consumed for ${item.food_item}.`;
        const note = cleanNote(item);
        const portion = item.portion_size_g == null ? "" : ` · Portion ${fmt(item.portion_size_g, 1)} g`;
        const units = item.units ? ` · ${escapeHTML(item.units)}` : "";
        const meta = `Field ID ${escapeHTML(item.table_s1_field_id)}${units}${portion}`;

        return `<div class="field-card">
          <label class="field-label" for="${escapeHTML(`field-${item.active_item_id}`)}">
            ${escapeHTML(item.food_item)}
            ${note ? `<span class="field-info" tabindex="0" title="${escapeHTML(note)}" aria-label="Additional question information">i</span>` : ""}
          </label>
          <p class="field-prompt">${escapeHTML(prompt)}</p>
          ${renderField(item)}
          <small class="field-meta">${meta}</small>
        </div>`;
      }).join("");

      details.innerHTML = `
        <summary>
          <span class="group-index">${String(group.order).padStart(2, "0")}</span>
          <span class="group-title">${escapeHTML(group.name)}
            <small>${group.items.length} questionnaire items</small>
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

    $$(".no-intake-check").forEach(check => check.addEventListener("change", e => {
      const group = e.target.closest(".food-group");
      $$(".response-input", group).forEach(input => {
        input.value = e.target.checked ? "0" : "";
        input.disabled = e.target.checked;
      });
      updateProgress();
    }));

    $$(".clear-group").forEach(button => button.addEventListener("click", e => {
      const group = e.target.closest(".food-group");
      $(".no-intake-check", group).checked = false;
      $$(".response-input", group).forEach(input => {
        input.value = "";
        input.disabled = false;
      });
      updateProgress();
    }));
  }

  function isGroupComplete(el) {
    if ($(".no-intake-check", el).checked) return true;
    return $$(".response-input", el).every(input => input.value !== "");
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
      ? "All required responses are complete."
      : "Complete all 15 groups, or explicitly confirm no intake.");
  }

  function collectIntakes() {
    const intakes = {};
    for (const group of foodGroups()) {
      const el = $(`.food-group[data-group="${group.key}"]`);
      const noIntake = $(".no-intake-check", el).checked;
      let total = 0;
      for (const item of group.items) {
        const input = $(`[data-item-id="${item.active_item_id}"]`, el);
        const raw = noIntake ? "0" : input.value;
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
      $$(".response-input").forEach(input => { input.value = ""; input.disabled = false; });
      $$(".no-intake-check").forEach(check => { check.checked = false; });
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
