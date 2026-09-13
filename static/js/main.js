/* ==========================================================================
   main.js
   Shared front-end logic for every method page:
     - editable Time/Concentration table (add / remove / edit rows)
     - optional third column (C_R) for autocatalytic / irreversible pages
     - CSV upload
     - generic results renderer (equation, transformed table, Plotly graph,
       step-by-step explanation, result cards)
     - prediction calculator (concentration / time / rate)
   ========================================================================== */

function addRow(tableBodyId, withThirdCol) {
  const tbody = document.getElementById(tableBodyId);
  const row = document.createElement("tr");
  row.innerHTML = `
    <td><input type="number" step="any" class="time-input" placeholder="t"></td>
    <td><input type="number" step="any" class="conc-input" placeholder="C"></td>
    ${withThirdCol ? '<td><input type="number" step="any" class="conc-r-input" placeholder="C_R (optional)"></td>' : ""}
    <td><button type="button" class="remove-row-btn" onclick="this.closest('tr').remove()">&times;</button></td>
  `;
  tbody.appendChild(row);
}

function collectTableData(tableBodyId, withThirdCol) {
  const rows = document.querySelectorAll(`#${tableBodyId} tr`);
  const time = [], concentration = [], concentration_r = [];
  rows.forEach(r => {
    const t = r.querySelector(".time-input")?.value;
    const c = r.querySelector(".conc-input")?.value;
    const cr = r.querySelector(".conc-r-input")?.value;
    if (t !== "" && c !== "" && t !== undefined && c !== undefined) {
      time.push(parseFloat(t));
      concentration.push(parseFloat(c));
      if (withThirdCol && cr !== "" && cr !== undefined) {
        concentration_r.push(parseFloat(cr));
      }
    }
  });
  const out = { time, concentration };
  if (withThirdCol && concentration_r.length === time.length && concentration_r.length > 0) {
    out.concentration_r = concentration_r;
    out.time_r = time;
  }
  return out;
}

function fillTableFromCSV(tableBodyId, data, withThirdCol) {
  const tbody = document.getElementById(tableBodyId);
  tbody.innerHTML = "";
  const n = data.time.length;
  for (let i = 0; i < n; i++) {
    addRow(tableBodyId, withThirdCol);
    const rows = tbody.querySelectorAll("tr");
    const row = rows[rows.length - 1];
    row.querySelector(".time-input").value = data.time[i];
    row.querySelector(".conc-input").value = data.concentration[i];
    if (withThirdCol && data.concentration_r) {
      row.querySelector(".conc-r-input").value = data.concentration_r[i] ?? "";
    }
  }
}

function bindCSVUpload(inputId, tableBodyId, withThirdCol, errorBoxId) {
  const input = document.getElementById(inputId);
  if (!input) return;
  input.addEventListener("change", async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const formData = new FormData();
    formData.append("file", file);
    try {
      const resp = await fetch("/upload_csv", { method: "POST", body: formData });
      const data = await resp.json();
      if (data.error) {
        showError(errorBoxId, data.error);
        return;
      }
      fillTableFromCSV(tableBodyId, data, withThirdCol);
      hideError(errorBoxId);
    } catch (err) {
      showError(errorBoxId, "Failed to upload/parse CSV: " + err);
    }
  });
}

function showError(boxId, msg) {
  const box = document.getElementById(boxId);
  if (!box) return;
  box.style.display = "block";
  box.innerHTML = `<div class="alert-quiet-error">${msg}</div>`;
}
function hideError(boxId) {
  const box = document.getElementById(boxId);
  if (box) box.style.display = "none";
}

function toggleLoading(spinnerId, show) {
  const el = document.getElementById(spinnerId);
  if (el) el.style.display = show ? "block" : "none";
}

/* ---- Math / units formatting --------------------------------------------
   The backend returns plain-ASCII equation strings (e.g. "C_A0", "C_A^(1-n)",
   "k1", "log10(C_A)") since it has no notion of markup. formatMath() turns
   those into properly typeset HTML (real subscripts/superscripts) so the
   whole site reads like textbook notation instead of raw variable names.
   ========================================================================== */
function formatMath(str) {
  if (str === null || str === undefined) return "";
  let s = String(str);
  s = s
    .replace(/C_A0/g, "C<sub>A0</sub>")
    .replace(/C_R0/g, "C<sub>R0</sub>")
    .replace(/C_Ae/g, "C<sub>Ae</sub>")
    .replace(/C_Re/g, "C<sub>Re</sub>")
    .replace(/C_A(?![0e])/g, "C<sub>A</sub>")
    .replace(/C_R(?![0e])/g, "C<sub>R</sub>")
    .replace(/X_Ae/g, "X<sub>Ae</sub>")
    .replace(/X_A(?!e)/g, "X<sub>A</sub>")
    .replace(/r_A/g, "r<sub>A</sub>")
    .replace(/\bk1\b/g, "k<sub>1</sub>")
    .replace(/\bk2\b/g, "k<sub>2</sub>")
    .replace(/log10/g, "log<sub>10</sub>")
    .replace(/\^\(([^)]+)\)/g, "<sup>$1</sup>")
    .replace(/\^(-?[0-9A-Za-z.]+)/g, "<sup>$1</sup>");
  return s;
}

/* ---- Units --------------------------------------------------------------
   Units don't change the math (the analysis works on whatever numbers were
   entered) - they only change how k, rates, and axes are LABELED, so this
   is handled entirely client-side. injectUnitSelectors() drops a Time-unit /
   Concentration-unit picker above the data table; getUnits() reads the
   current selection; kUnitsForOrder()/kUnitsFirstOrder()/kUnitsSecondOrderLike()
   build the correct unit string for the rate constant given the reaction's
   functional form.
   ========================================================================== */
function injectUnitSelectors(beforeElId) {
  const target = document.getElementById(beforeElId);
  if (!target || document.getElementById('timeUnitSel')) return;
  const wrap = document.createElement('div');
  wrap.className = 'row g-3 mb-3';
  wrap.innerHTML = `
    <div class="col-md-4">
      <label class="form-label">Time unit</label>
      <select id="timeUnitSel" class="form-select">
        <option value="s">seconds (s)</option>
        <option value="min" selected>minutes (min)</option>
        <option value="hr">hours (hr)</option>
      </select>
    </div>
    <div class="col-md-4">
      <label class="form-label">Concentration unit</label>
      <select id="concUnitSel" class="form-select">
        <option value="mol&middot;L<sup>-1</sup>" selected>mol/L</option>
        <option value="mol&middot;dm<sup>-3</sup>">mol/dm&sup3;</option>
        <option value="mmol&middot;L<sup>-1</sup>">mmol/L</option>
      </select>
    </div>`;
  target.parentNode.insertBefore(wrap, target);
}

function getUnits() {
  const t = document.getElementById('timeUnitSel');
  const c = document.getElementById('concUnitSel');
  return { time: t ? t.value : 'time', conc: c ? c.value : 'conc' };
}

function kUnitsForOrder(n, concUnit, timeUnit) {
  const nn = Math.round(n * 100) / 100;
  if (Math.abs(nn - 1) < 1e-6) return `${timeUnit}<sup>-1</sup>`;
  const exp = 1 - nn;
  const expStr = Number.isInteger(exp) ? String(exp) : exp.toFixed(2);
  return `(${concUnit})<sup>${expStr}</sup>&middot;${timeUnit}<sup>-1</sup>`;
}
function kUnitsFirstOrder(timeUnit) {
  return `${timeUnit}<sup>-1</sup>`;
}
function kUnitsSecondOrderLike(concUnit, timeUnit) {
  return `${concUnit}<sup>-1</sup>&middot;${timeUnit}<sup>-1</sup>`;
}

/* ---- Order-by-order narrative ("check 0th, does it fit? no -> try 1st...")
   Builds a plain-language walk-through of the R² grid search so the reader
   sees WHY an order was chosen, not just the final winner. ------------------ */
function buildOrderNarrative(comparisonTable) {
  const sorted = [...comparisonTable].sort((a, b) => a.order - b.order);
  let best = null;
  sorted.forEach(c => {
    if (c.physically_valid && (!best || c.r2 > best.r2)) best = c;
  });
  let html = '<p>The app starts at zero order and works upward, testing how well a straight line fits the linearized data at each order:</p><ol>';
  sorted.forEach(c => {
    const isBest = best && c.order === best.order;
    let verdict;
    if (!c.physically_valid) verdict = 'rejected &mdash; gives a non-physical (negative) k';
    else if (c.r2 >= 0.98) verdict = 'excellent fit';
    else if (c.r2 >= 0.90) verdict = 'reasonable fit, but not the best';
    else verdict = 'poor fit &mdash; data is clearly not linear at this order';
    html += `<li${isBest ? ' style="font-weight:700;color:var(--success);"' : ''}>Order n = ${c.order}: R&sup2; = ${fmt(c.r2)} &rarr; ${verdict}${isBest ? ' ✅ <b>selected</b>' : ''}</li>`;
  });
  html += '</ol>';
  return { html, best };
}

/* ---- Generic results rendering ------------------------------------------ */

function renderInputTable(containerId, inputData) {
  let html = `<table class="kinetics-table"><thead><tr><th>Time</th><th>Concentration</th></tr></thead><tbody>`;
  inputData.forEach(r => {
    html += `<tr><td>${r.time}</td><td>${(+r.concentration).toPrecision(5)}</td></tr>`;
  });
  html += `</tbody></table>`;
  document.getElementById(containerId).innerHTML = html;
}

function renderTransformedTable(containerId, rows, columns) {
  // columns: [{key, label}]
  let html = `<div style="overflow-x:auto"><table class="kinetics-table"><thead><tr>`;
  columns.forEach(c => html += `<th>${formatMath(c.label)}</th>`);
  html += `</tr></thead><tbody>`;
  rows.forEach(r => {
    html += "<tr>";
    columns.forEach(c => {
      const v = r[c.key];
      html += `<td>${v === null || v === undefined ? "&mdash;" : (typeof v === "number" ? v.toPrecision(5) : v)}</td>`;
    });
    html += "</tr>";
  });
  html += `</tbody></table></div>`;
  document.getElementById(containerId).innerHTML = html;
}

function renderSteps(containerId, steps) {
  let html = "";
  steps.forEach(s => {
    html += `<div class="step-card"><h6>${formatMath(s.title)}</h6><div>${formatMath(s.text)}</div></div>`;
  });
  document.getElementById(containerId).innerHTML = html;
}

function renderGraph(divId, graph, title) {
  const points = {
    x: graph.points_x,
    y: graph.points_y,
    mode: "markers",
    type: "scatter",
    name: "Experimental data",
    marker: { color: "#AD9C8E", size: 9 }
  };
  const line = {
    x: graph.line_x,
    y: graph.line_y,
    mode: "lines",
    type: "scatter",
    name: "Best-fit line",
    line: { color: "#3a332c", width: 2 }
  };
  const layout = {
    title: title || "",
    xaxis: { title: graph.x_axis },
    yaxis: { title: graph.y_axis },
    plot_bgcolor: "#fffdf8",
    paper_bgcolor: "#fffdf8",
    margin: { t: 50 },
    legend: { orientation: "h", y: -0.2 }
  };
  Plotly.newPlot(divId, [points, line], layout, { responsive: true });
}

function renderResultCards(containerId, cards) {
  // cards: [{label, value}]
  let html = "";
  cards.forEach(c => {
    html += `
      <div class="col-6 col-md-3 mb-3">
        <div class="result-card">
          <div class="label">${formatMath(c.label)}</div>
          <div class="value">${formatMath(c.value)}</div>
        </div>
      </div>`;
  });
  document.getElementById(containerId).innerHTML = html;
}

function fmt(x, sig = 5) {
  if (x === null || x === undefined || isNaN(x)) return "&mdash;";
  return Number(x).toPrecision(sig);
}
