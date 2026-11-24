(() => {
  "use strict";
  const esc = v => String(v ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const n = (x, d=3) => Number(x).toFixed(d);
  function binary(rule) {
    const a=n(rule.value1);
    if (rule.range_direction==="middle") return `${a}–${n(rule.value2)}`;
    if (rule.range_direction==="high") return `≥ ${a}`;
    if (rule.range_direction==="low") return `≤ ${a}`;
    return `Point ${a}`;
  }
  function continuous(rule) {
    if (rule.range_direction==="middle") return `${n(rule.value1)}–${n(rule.value2)}`;
    if (rule.range_direction==="high") return `Adequacy ≥ ${n(rule.value1)}`;
    if (rule.range_direction==="low") return `Moderation ≤ ${n(rule.value1)}`;
    return `Point target ${n(rule.value1)}`;
  }
  document.addEventListener("DOMContentLoaded", async () => {
    try {
      await window.RWDIDPaletteReady;
      const response = await fetch("data/rwdid_parameters.json");
      if (!response.ok) throw new Error("Parameter file could not be loaded.");
      const p = await response.json();
      const body = document.getElementById("parameter-table-body");
      body.innerHTML = p.components.map((c,i) => {
        const w=Number(c.beta)/Number(p.beta_sum)*100;
        return `<tr>
          <td><span class="table-dot" style="background:var(--friendly-${(i%13)+1})"></span>${esc(c.display_name)}</td>
          <td>${binary(c.binary)}</td>
          <td>${n(c.beta,4)}</td>
          <td>${n(w,2)}%</td>
          <td>${continuous(c.continuous)}</td>
          <td>${n(c.continuous.support_p1,2)}–${n(c.continuous.support_p99,2)}</td>
        </tr>`;
      }).join("");
      const status=document.getElementById("parameter-version");
      if(status) status.textContent=`Frozen parameter set · ${p.version}`;
      document.getElementById("boot-state").hidden=true;
      document.getElementById("methods-content").hidden=false;
    } catch(e) {
      const boot=document.getElementById("boot-state");
      boot.classList.add("boot-error");
      boot.innerHTML=`<strong>Methods page could not be initialized.</strong><p>${esc(e.message)}</p>`;
    }
  });
})();
