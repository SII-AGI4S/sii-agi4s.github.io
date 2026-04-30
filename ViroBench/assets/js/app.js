let RAW = { taxonomy: {}, host: {}, genome: {}, cds: {} };

const DATA_FILES = {
  taxonomy: "data/json/Task1_Taxonomy Classification.json",
  host: "data/json/Task2_Host Prediction.json",
  genome: "data/json/Tasks3_Genome Modeling.json",
  cds: "data/json/Tasks4_CDS Completion.json",
};

const TASK_DEFS = {
  taxonomy: { label: "Taxonomy", subtitle: "hierarchical Macro-F1", color: "#0f9f9a" },
  host: { label: "Host", subtitle: "host Macro-F1", color: "#2563eb" },
  genome: { label: "Genome", subtitle: "BPB / lower better", color: "#4f46e5" },
  cds: { label: "CDS", subtitle: "completion metrics", color: "#e76f51" },
};

const taxMetrics = ["kingdom_f1_macro", "phylum_f1_macro", "class_f1_macro", "order_f1_macro", "family_f1_macro"];
const metricLabel = {
  avg_taxon: "Avg. Macro-F1",
  kingdom_f1_macro: "Kingdom F1",
  phylum_f1_macro: "Phylum F1",
  class_f1_macro: "Class F1",
  order_f1_macro: "Order F1",
  family_f1_macro: "Family F1",
  f1_macro: "Macro-F1",
  mean: "Mean BPB",
  median: "Median BPB",
  min: "Min BPB",
  max: "Max BPB",
  edit_distance: "Edit distance",
  exact_match_acc: "Exact match",
  is_CDS: "CDS success",
  kmer_JSD: "K-mer JSD",
  kmer_KS: "K-mer KS",
};
const lowerBetter = new Set(["mean", "median", "min", "max", "edit_distance", "kmer_JSD", "kmer_KS"]);

let state = { task: "taxonomy", scenario: "ALL-taxon-genus", metric: "avg_taxon", search: "" };
let selectedModel = null;
let genBucket = "genome-short";
let cdsBucket = "cds-short";
const scatterState = {
  xTask: "taxonomy",
  yTask: "host",
  cdsScenario: "cds-short",
  cdsMetric: "kmer_JSD",
  hiddenFamilies: new Set(),
  hitPoints: [],
};

const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));

function parseVal(v) {
  if (v === null || v === undefined) return NaN;
  if (typeof v === "number") return v;
  return Number(String(v).split("(")[0]);
}

function normalizeModelKey(name = "") {
  return String(name).trim().toLowerCase();
}

function fmt(v, metric) {
  if (!Number.isFinite(v)) return "—";
  if (metric === "is_CDS") return (v * 100).toFixed(2) + "%";
  if (["edit_distance", "exact_match_acc", "kmer_JSD", "kmer_KS"].includes(metric)) return v.toFixed(4);
  return v.toFixed(2);
}

function scenarioLabel(k) {
  return k
    .replace("ALL", "All")
    .replace("DNA", "DNA")
    .replace("RNA", "RNA")
    .replace("taxon", "taxonomy")
    .replace("times", "temporal")
    .replace("genus", "genus-disjoint")
    .replace("host", "host")
    .replaceAll("-", " · ");
}

function familyOf(name = "") {
  const n = name.toLowerCase();
  if (n.includes("virohyena")) return ["ViroHyena", "#0f9f9a"];
  if (n.includes("evo2")) return ["Evo2", "#4f46e5"];
  if (n.includes("evo-1")) return ["Evo1", "#6366f1"];
  if (n.includes("aido.dna")) return ["AIDO.DNA", "#2563eb"];
  if (n.includes("aido.rna")) return ["AIDO.RNA", "#1d4ed8"];
  if (n.includes("aido")) return ["AIDO", "#2563eb"];
  if (n.includes("lucavirus")) return ["LucaVirus", "#10b981"];
  if (n.includes("lucaone")) return ["LucaOne", "#22c55e"];
  if (n.includes("genomeocean")) return ["GenomeOcean", "#06b6d4"];
  if (n.includes("genos")) return ["Genos", "#7c3aed"];
  if (n.includes("generator")) return ["GENERator", "#f97316"];
  if (n.includes("hyenadna")) return ["HyenaDNA", "#f59e0b"];
  if (n.includes("dnabert")) return ["DNABERT", "#ef4444"];
  if (n.startsWith("ntv3")) return ["Ntv3", "#b91c1c"];
  if (n.startsWith("ntv2")) return ["Ntv2", "#dc2626"];
  if (n.startsWith("nt-")) return ["Nt", "#ef4444"];
  if (n.includes("rna-fm") || n.includes("rinalmo") || n.includes("birna") || n.includes("rnabert") || n.includes("mp-rna")) return ["RNA Foundation", "#e76f51"];
  if (n.includes("gena-lm")) return ["Gena-lm", "#0284c7"];
  if (n.includes("caduceus") || n.includes("omnireg") || n.includes("grover")) return ["General Genomic", "#0ea5e9"];
  if (n.includes("blast") || n.includes("kraken")) return ["Classical", "#334155"];
  if (n.includes("cnn") || n.includes("bilstm")) return ["Neural baseline", "#64748b"];
  return ["General NFM", "#94a3b8"];
}

function avgTaxon(row) {
  return taxMetrics.reduce((a, m) => a + parseVal(row[m]), 0) / taxMetrics.length;
}

function scoreOf(row, metric = state.metric, task = state.task) {
  if (task === "taxonomy" && metric === "avg_taxon") return avgTaxon(row);
  return parseVal(row[metric]);
}

function currentRows() {
  let rows = [];
  if (state.task === "taxonomy") rows = RAW.taxonomy[state.scenario] || [];
  if (state.task === "host") rows = RAW.host[state.scenario] || [];
  if (state.task === "genome") rows = RAW.genome[state.scenario] || [];
  if (state.task === "cds") rows = RAW.cds[state.scenario] || [];
  const q = state.search.toLowerCase().trim();
  if (q) rows = rows.filter((r) => r.model.toLowerCase().includes(q));
  const lb = lowerBetter.has(state.metric);
  return rows
    .map((r) => ({ ...r, _score: scoreOf(r) }))
    .filter((r) => Number.isFinite(r._score))
    .sort((a, b) => (lb ? a._score - b._score : b._score - a._score));
}

function bestRowsFor(task, scenario, metric, n = 1) {
  const old = { ...state };
  state = { ...state, task, scenario, metric, search: "" };
  const rows = currentRows().slice(0, n);
  state = old;
  return rows;
}

function initStats() {
  const stats = [
    { v: "4", l: "Task types", s: "Taxonomy · Host · Genome · CDS", c: "rgba(15,159,154,.18)" },
    { v: "18", l: "Scenarios", s: "12 classification + 6 generation", c: "rgba(37,99,235,.16)" },
    { v: "74", l: "Models in JSON", s: "classification tables include classical tools", c: "rgba(79,70,229,.16)" },
    { v: "58K", l: "Viral samples", s: "paper-level curated corpus", c: "rgba(231,111,81,.16)" },
  ];
  $("#statGrid").innerHTML = stats
    .map(
      (x) =>
        `<div class="stat" style="--accent:${x.c}"><div class="val">${x.v}</div><div class="label">${x.l}</div><div class="sub">${x.s}</div></div>`
    )
    .join("");
}

function renderTaskCards() {
  const cards = [
    {
      id: "Task 1",
      t: "Taxonomy Classification",
      p: "Predict hierarchical viral taxonomy from Kingdom to Family under genus-disjoint and temporal shifts.",
      m: "6",
      l: "Scenarios",
      c1: "#0f9f9a",
      c2: "#2563eb",
      best: bestRowsFor("taxonomy", "ALL-taxon-genus", "avg_taxon", 1)[0]?.model || "—",
    },
    {
      id: "Task 2",
      t: "Host Prediction",
      p: "Classify standardized host categories and probe whether sequence representations encode ecological signals.",
      m: "6",
      l: "Scenarios",
      c1: "#2563eb",
      c2: "#4f46e5",
      best: bestRowsFor("host", "ALL-host-genus", "f1_macro", 1)[0]?.model || "—",
    },
    {
      id: "Task 3",
      t: "Genome Modeling",
      p: "Measure next-token sequence likelihood in short, medium, and long viral genome fragments using BPB.",
      m: "3",
      l: "Length buckets",
      c1: "#4f46e5",
      c2: "#7c3aed",
      best: bestRowsFor("genome", "genome-short", "mean", 1)[0]?.model || "—",
    },
    {
      id: "Task 4",
      t: "CDS Completion",
      p: "Evaluate protein-coding sequence continuation with edit distance, exact match, CDS success and k-mer statistics.",
      m: "3",
      l: "Length buckets",
      c1: "#e76f51",
      c2: "#f59e0b",
      best: bestRowsFor("cds", "cds-short", "is_CDS", 1)[0]?.model || "—",
    },
  ];
  const block = document.createElement("div");
  block.className = "cards4 fade-in delay-2";
  block.innerHTML = cards
    .map(
      (x) =>
        `<div class="task-card" style="--c1:${x.c1};--c2:${x.c2}"><div class="task-id">${x.id}</div><h3>${x.t}</h3><p>${x.p}</p><div class="task-metric"><div><strong>${x.m}</strong><span>${x.l}</span></div><div><strong style="font-size:16px;max-width:135px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis">${x.best}</strong><span>Current top</span></div></div></div>`
    )
    .join("");
  const ref = $(".section-title");
  ref.parentNode.insertBefore(block, ref);
}

function renderPills() {
  const tasks = Object.keys(TASK_DEFS);
  $("#taskPills").innerHTML = tasks
    .map((t) => `<button class="pill ${state.task === t ? "on" : ""}" data-task="${t}">${TASK_DEFS[t].label}</button>`)
    .join("");
  $$("#taskPills .pill").forEach((b) => (b.onclick = () => {
    state.task = b.dataset.task;
    state.search = "";
    if (state.task === "taxonomy") { state.scenario = "ALL-taxon-genus"; state.metric = "avg_taxon"; }
    if (state.task === "host") { state.scenario = "ALL-host-genus"; state.metric = "f1_macro"; }
    if (state.task === "genome") { state.scenario = "genome-short"; state.metric = "mean"; }
    if (state.task === "cds") { state.scenario = "cds-short"; state.metric = "kmer_JSD"; }
    selectedModel = null;
    renderAll();
  }));
}

function metricOptions() {
  if (state.task === "taxonomy") return ["avg_taxon", ...taxMetrics];
  if (state.task === "host") return ["f1_macro"];
  if (state.task === "genome") return ["mean", "median", "min", "max"];
  return ["kmer_JSD", "kmer_KS", "edit_distance", "exact_match_acc", "is_CDS"];
}

function sortScenarioKeys(keys) {
  const prefixRank = { ALL: 0, DNA: 1, RNA: 2 };
  return [...keys].sort((a, b) => {
    const pa = a.split("-")[0];
    const pb = b.split("-")[0];
    const ra = prefixRank[pa] ?? 99;
    const rb = prefixRank[pb] ?? 99;
    if (ra !== rb) return ra - rb;
    const sa = a.slice(pa.length + 1);
    const sb = b.slice(pb.length + 1);
    return sa.localeCompare(sb);
  });
}

function scenarioOptions() {
  if (state.task === "taxonomy") return sortScenarioKeys(Object.keys(RAW.taxonomy));
  if (state.task === "host") return sortScenarioKeys(Object.keys(RAW.host));
  if (state.task === "genome") return Object.keys(RAW.genome);
  return Object.keys(RAW.cds);
}

function scoreBarNorm(task, metric, value, min, max) {
  if ((task === "taxonomy" || task === "host") && !lowerBetter.has(metric)) {
    return Math.max(0, Math.min(1, value / 100));
  }
  if (lowerBetter.has(metric)) return (max - value) / (max - min || 1);
  return (value - min) / (max - min || 1);
}

function renderControls() {
  const scenarios = scenarioOptions();
  const metrics = metricOptions();
  $("#controls").innerHTML =
    `<div class="control"><label>Scenario</label><select id="scenarioSel">${scenarios.map((s) => `<option value="${s}" ${s === state.scenario ? "selected" : ""}>${scenarioLabel(s)}</option>`).join("")}</select></div>` +
    `<div class="control"><label>Metric</label><select id="metricSel">${metrics.map((m) => `<option value="${m}" ${m === state.metric ? "selected" : ""}>${metricLabel[m]}</option>`).join("")}</select></div>` +
    `<div class="control"><label>Search model</label><input id="searchInp" placeholder="e.g. Evo2, ViroHyena" value="${state.search}"></div>`;
  $("#scenarioSel").onchange = (e) => { state.scenario = e.target.value; selectedModel = null; renderAll(false); };
  $("#metricSel").onchange = (e) => { state.metric = e.target.value; selectedModel = null; renderAll(false); };
  $("#searchInp").oninput = (e) => { state.search = e.target.value; renderLeaderboard(); drawBarChart(); };
}

function renderLeaderboard() {
  const rows = currentRows();
  const max = Math.max(...rows.map((r) => r._score));
  const min = Math.min(...rows.map((r) => r._score));
  const task = state.task;
  $("#leaderHead").innerHTML = `<tr><th>#</th><th>Model</th><th>${metricLabel[state.metric]}</th><th>Notes</th></tr>`;
  $("#leaderBody").innerHTML = rows.slice(0, 74).map((r, i) => {
    const [fam, col] = familyOf(r.model);
    const norm = scoreBarNorm(task, state.metric, r._score, min, max);
    const rk = i + 1;
    return `<tr data-model="${r.model.replaceAll('"', "&quot;")}"><td class="rank ${rk === 1 ? "r1" : rk === 2 ? "r2" : rk === 3 ? "r3" : ""}">${rk}</td><td><div class="model-cell" style="--family-color:${col}"><span class="model-dot"></span><div><div class="model-name">${r.model}</div><div class="family">${fam}</div></div></div></td><td class="bar-cell"><div class="bar-line"><div class="bar-bg"><div class="bar-fill" style="width:${Math.max(3, norm * 100).toFixed(1)}%"></div></div><div class="score">${fmt(r._score, state.metric)}</div></div></td><td class="delta">${lowerBetter.has(state.metric) ? "lower is better" : "higher is better"}</td></tr>`;
  }).join("");
  $$("#leaderBody tr").forEach((tr) => (tr.onclick = () => {
    selectedModel = tr.dataset.model;
    renderDetail(rows.find((r) => r.model === selectedModel));
    drawBarChart();
  }));
  renderDetail(rows.find((r) => r.model === selectedModel) || rows[0]);
  const active = rows.find((r) => r.model === selectedModel) || rows[0];
  $("#chartHint").textContent = active ? `${active.model} · ${scenarioLabel(state.scenario)}` : `${TASK_DEFS[state.task].label} · ${scenarioLabel(state.scenario)}`;
}

function renderDetail(row) {
  if (!row) { $("#modelDetail").innerHTML = ""; return; }
  selectedModel = row.model;
  let extras = [];
  if (state.task === "taxonomy") extras = taxMetrics.map((m) => [metricLabel[m], fmt(parseVal(row[m]), m)]);
  else if (state.task === "host") extras = [["Macro-F1", fmt(parseVal(row.f1_macro), "f1_macro")]];
  else if (state.task === "genome") extras = ["mean", "median", "min", "max"].map((m) => [metricLabel[m], fmt(parseVal(row[m]), m)]);
  else extras = ["kmer_JSD", "kmer_KS", "edit_distance", "exact_match_acc", "is_CDS"].map((m) => [metricLabel[m], fmt(parseVal(row[m]), m)]);
  const [fam, col] = familyOf(row.model);
  $("#modelDetail").innerHTML = `<div class="detail-title" style="color:${col}">${row.model}</div><div class="family">${fam} · ${scenarioLabel(state.scenario)}</div><div class="detail-grid">${extras.slice(0, 6).map(([k, v]) => `<div class="kv"><small>${k}</small><strong>${v}</strong></div>`).join("")}</div>`;
}

function fitCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = rect.width * dpr;
  canvas.height = rect.height * dpr;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: rect.width, h: rect.height };
}

function roundRect(ctx, x, y, w, h, r) {
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.arcTo(x + w, y, x + w, y + h, r);
  ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r);
  ctx.arcTo(x, y, x + w, y, r);
  ctx.closePath();
}

function drawBars(canvas, rows, metric, opts = {}) {
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  const top = rows.slice(0, 10);
  if (!top.length) return;
  const task = opts.task || state.task;
  const vals = top.map((r) => r._score ?? scoreOf(r, metric, opts.task || state.task));
  const max = Math.max(...vals);
  const min = Math.min(...vals);
  const padL = 140;
  const padR = 28;
  const padT = 18;
  const rowH = (h - padT - 24) / top.length;
  const start = performance.now();
  function frame() {
    const p = Math.min(1, (performance.now() - start) / 700);
    ctx.clearRect(0, 0, w, h);
    ctx.font = "11px " + getComputedStyle(document.body).fontFamily;
    top.forEach((r, i) => {
      const val = vals[i];
      const best = scoreBarNorm(task, metric, val, min, max);
      const bw = (w - padL - padR) * (0.12 + 0.88 * best) * p;
      const y = padT + i * rowH + rowH * 0.22;
      const bh = Math.max(8, rowH * 0.43);
      ctx.fillStyle = "#64748b";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(r.model.length > 17 ? r.model.slice(0, 16) + "…" : r.model, padL - 12, y + bh / 2);
      ctx.fillStyle = "#eef3f8";
      roundRect(ctx, padL, y, w - padL - padR, bh, 999);
      ctx.fill();
      const grad = ctx.createLinearGradient(padL, 0, padL + bw, 0);
      grad.addColorStop(0, opts.c1 || TASK_DEFS[state.task].color);
      grad.addColorStop(1, opts.c2 || "#4f46e5");
      ctx.fillStyle = grad;
      roundRect(ctx, padL, y, bw, bh, 999);
      ctx.fill();
      ctx.fillStyle = "#0f172a";
      ctx.textAlign = "left";
      ctx.fillText(fmt(val, metric), padL + bw + 8, y + bh / 2);
    });
    if (p < 1) requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}

function selectedMetricItems(row) {
  if (!row) return [];
  if (state.task === "taxonomy") {
    return taxMetrics.map((m) => ({ label: metricLabel[m], value: parseVal(row[m]), metric: m }));
  }
  if (state.task === "host") {
    return [{ label: "Macro-F1", value: parseVal(row.f1_macro), metric: "f1_macro" }];
  }
  if (state.task === "genome") {
    return ["mean", "median", "min", "max"].map((m) => ({ label: metricLabel[m], value: parseVal(row[m]), metric: m }));
  }
  return ["kmer_JSD", "kmer_KS", "edit_distance", "exact_match_acc", "is_CDS"].map((m) => ({ label: metricLabel[m], value: parseVal(row[m]), metric: m }));
}

function drawSelectedMetrics(canvas, row) {
  const items = selectedMetricItems(row).filter((x) => Number.isFinite(x.value));
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  if (!items.length) return;

  const vals = items.map((x) => x.value);
  const max = Math.max(...vals);
  const pad = { l: 52, r: 18, t: 18, b: 56 };
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b;
  let yMax = max <= 0 ? 1 : max * 1.1;
  if (state.task === "taxonomy" || state.task === "host") yMax = 100;
  const start = performance.now();

  function frame() {
    const p = Math.min(1, (performance.now() - start) / 700);
    ctx.clearRect(0, 0, w, h);
    ctx.font = "11px " + getComputedStyle(document.body).fontFamily;
    ctx.textBaseline = "middle";

    // y-axis grid + labels
    const ticks = 4;
    for (let i = 0; i <= ticks; i++) {
      const ratio = i / ticks;
      const y = pad.t + plotH * (1 - ratio);
      const tv = yMax * ratio;
      ctx.strokeStyle = "#e5edf5";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(pad.l, y);
      ctx.lineTo(w - pad.r, y);
      ctx.stroke();

      ctx.fillStyle = "#94a3b8";
      ctx.textAlign = "right";
      const tlabel = yMax <= 1 ? tv.toFixed(2) : yMax <= 10 ? tv.toFixed(1) : Math.round(tv).toString();
      ctx.fillText(tlabel, pad.l - 8, y);
    }

    // axes
    ctx.strokeStyle = "#cfdbe8";
    ctx.beginPath();
    ctx.moveTo(pad.l, pad.t);
    ctx.lineTo(pad.l, h - pad.b);
    ctx.lineTo(w - pad.r, h - pad.b);
    ctx.stroke();

    const n = items.length;
    const gap = Math.min(18, plotW / (n * 3));
    const barW = Math.max(22, (plotW - gap * (n + 1)) / n);

    items.forEach((item, i) => {
      const x = pad.l + gap + i * (barW + gap);
      const vh = Math.max(2, (item.value / yMax) * plotH * p);
      const y = h - pad.b - vh;

      const grad = ctx.createLinearGradient(x, 0, x + barW, 0);
      grad.addColorStop(0, TASK_DEFS[state.task].color);
      grad.addColorStop(1, "#4f46e5");
      ctx.fillStyle = grad;
      roundRect(ctx, x, y, barW, vh, 8);
      ctx.fill();

      ctx.fillStyle = "#0f172a";
      ctx.textAlign = "center";
      ctx.fillText(fmt(item.value, item.metric), x + barW / 2, y - 10);

      ctx.fillStyle = "#64748b";
      // Do not strip "-F1" (e.g. Macro-F1); \s*F1$ would turn "Macro-F1" into "Macro-".
      const shortLabel = /-F1$/i.test(item.label) ? item.label : item.label.replace(/\s+F1$/i, "");
      ctx.fillText(shortLabel, x + barW / 2, h - pad.b + 16);
    });

    if (p < 1) requestAnimationFrame(frame);
  }

  requestAnimationFrame(frame);
}

function drawBarChart() {
  const rows = currentRows();
  const active = rows.find((r) => r.model === selectedModel) || rows[0];
  drawSelectedMetrics($("#barChart"), active);
}

function scatterTaskMeta(task) {
  const meta = {
    taxonomy: { label: "Task 1 · Taxonomy", scenario: "ALL-taxon-genus", metric: "avg_taxon", lowerBetter: false },
    host: { label: "Task 2 · Host", scenario: "ALL-host-genus", metric: "f1_macro", lowerBetter: false },
    genome: { label: "Task 3 · Genome", scenario: "genome-short", metric: "mean", lowerBetter: true },
    cds: { label: "Task 4 · CDS", scenario: scatterState.cdsScenario, metric: scatterState.cdsMetric, lowerBetter: lowerBetter.has(scatterState.cdsMetric) },
  };
  return meta[task];
}

function scatterScoreValue(task, row) {
  const m = scatterTaskMeta(task);
  if (task === "taxonomy") return avgTaxon(row);
  return parseVal(row[m.metric]);
}

function scatterSeries(task) {
  const m = scatterTaskMeta(task);
  const source = RAW[task]?.[m.scenario] || [];
  const entries = source
    .map((r) => ({ model: r.model, value: scatterScoreValue(task, r) }))
    .filter((x) => Number.isFinite(x.value));
  const map = new Map();
  entries.forEach((x) => map.set(normalizeModelKey(x.model), x));
  const metricName = metricLabel[m.metric] || m.metric;
  const suffix = m.lowerBetter ? " (lower better · axis reversed)" : " (higher better)";
  return { map, label: `${m.label} · ${metricName}${suffix}` };
}

function initScatterControls() {
  const tasks = ["taxonomy", "host", "genome", "cds"];
  const options = tasks.map((t) => `<option value="${t}">${scatterTaskMeta(t).label}</option>`).join("");
  $("#scatterXSel").innerHTML = options;
  $("#scatterYSel").innerHTML = options;
  $("#scatterXSel").value = scatterState.xTask;
  $("#scatterYSel").value = scatterState.yTask;
  $("#scatterXSel").onchange = (e) => { scatterState.xTask = e.target.value; updateScatterTask4ControlsVisibility(); drawScatter(); };
  $("#scatterYSel").onchange = (e) => { scatterState.yTask = e.target.value; updateScatterTask4ControlsVisibility(); drawScatter(); };

  const cdsScenarioOptions = Object.keys(RAW.cds || {}).map((s) => `<option value="${s}">${scenarioLabel(s)}</option>`).join("");
  const cdsMetricList = ["kmer_JSD", "kmer_KS", "edit_distance", "exact_match_acc", "is_CDS"];
  const cdsMetricOptions = cdsMetricList.map((m) => `<option value="${m}">${metricLabel[m]}</option>`).join("");
  $("#scatterCdsScenarioSel").innerHTML = cdsScenarioOptions;
  $("#scatterCdsMetricSel").innerHTML = cdsMetricOptions;
  if (!RAW.cds?.[scatterState.cdsScenario]) scatterState.cdsScenario = Object.keys(RAW.cds || {})[0] || "cds-short";
  $("#scatterCdsScenarioSel").value = scatterState.cdsScenario;
  $("#scatterCdsMetricSel").value = scatterState.cdsMetric;
  $("#scatterCdsScenarioSel").onchange = (e) => { scatterState.cdsScenario = e.target.value; drawScatter(); };
  $("#scatterCdsMetricSel").onchange = (e) => { scatterState.cdsMetric = e.target.value; drawScatter(); };
  updateScatterTask4ControlsVisibility();
}

function updateScatterTask4ControlsVisibility() {
  const show = scatterState.xTask === "cds" || scatterState.yTask === "cds";
  const scenarioCtl = $("#scatterCdsScenarioControl");
  const metricCtl = $("#scatterCdsMetricControl");
  if (scenarioCtl) scenarioCtl.style.display = show ? "" : "none";
  if (metricCtl) metricCtl.style.display = show ? "" : "none";
}

function renderScatterLegend(points) {
  const families = Array.from(new Map(points.map((p) => [p.family, p.color])).entries()).sort((a, b) => a[0].localeCompare(b[0]));
  const wrap = $("#scatterLegend");
  if (!wrap) return;
  wrap.innerHTML = families.map(([fam, col]) => {
    const off = scatterState.hiddenFamilies.has(fam) ? "off" : "";
    return `<button class="scatter-chip ${off}" data-family="${fam}"><span class="dot" style="background:${col}"></span><span>${fam}</span></button>`;
  }).join("");
  $$("#scatterLegend .scatter-chip").forEach((el) => {
    el.onclick = () => {
      const fam = el.dataset.family;
      if (scatterState.hiddenFamilies.has(fam)) scatterState.hiddenFamilies.delete(fam);
      else scatterState.hiddenFamilies.add(fam);
      drawScatter();
    };
  });
}

function hideScatterTooltip() {
  const tip = $("#scatterTooltip");
  if (!tip) return;
  tip.style.display = "none";
}

function initScatterHover() {
  const canvas = $("#scatterChart");
  if (!canvas || canvas.dataset.hoverBound === "1") return;
  canvas.dataset.hoverBound = "1";
  canvas.addEventListener("mousemove", (e) => {
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    let hit = null;
    let bestD2 = Infinity;
    scatterState.hitPoints.forEach((p) => {
      const dx = x - p.px;
      const dy = y - p.py;
      const d2 = dx * dx + dy * dy;
      const threshold = (p.r + 4) * (p.r + 4);
      if (d2 <= threshold && d2 < bestD2) {
        bestD2 = d2;
        hit = p;
      }
    });
    const tip = $("#scatterTooltip");
    if (!tip || !hit) {
      canvas.style.cursor = "default";
      hideScatterTooltip();
      return;
    }
    canvas.style.cursor = "pointer";
    const xMeta = scatterTaskMeta(scatterState.xTask);
    const yMeta = scatterTaskMeta(scatterState.yTask);
    const xMetricName = metricLabel[xMeta.metric] || xMeta.metric;
    const yMetricName = metricLabel[yMeta.metric] || yMeta.metric;
    tip.textContent = `${hit.model}\n${xMetricName}: ${fmt(hit.x, xMeta.metric)}\n${yMetricName}: ${fmt(hit.y, yMeta.metric)}`;
    tip.style.display = "block";
    const wrapRect = canvas.parentElement.getBoundingClientRect();
    const tx = e.clientX - wrapRect.left + 10;
    const ty = e.clientY - wrapRect.top - 10;
    tip.style.left = `${tx}px`;
    tip.style.top = `${ty}px`;
  });
  canvas.addEventListener("mouseleave", () => {
    canvas.style.cursor = "default";
    hideScatterTooltip();
  });
}

function drawScatter() {
  const canvas = $("#scatterChart");
  if (!canvas) return;
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  const xs = scatterSeries(scatterState.xTask);
  const ys = scatterSeries(scatterState.yTask);
  const pts = [];
  xs.map.forEach((xItem, modelKey) => {
    const yItem = ys.map.get(modelKey);
    if (!yItem || !Number.isFinite(yItem.value)) return;
    const displayModel = xItem.model || yItem.model;
    const [family, color] = familyOf(displayModel);
    pts.push({ model: displayModel, x: xItem.value, y: yItem.value, family, color });
  });
  renderScatterLegend(pts);
  const visible = pts.filter((p) => !scatterState.hiddenFamilies.has(p.family));
  if (!visible.length) {
    scatterState.hitPoints = [];
    hideScatterTooltip();
    $("#scatterHint").textContent = "No visible points. Click legend chips to re-enable families.";
    return;
  }
  const xv = visible.map((p) => p.x);
  const yv = visible.map((p) => p.y);
  let minX = Math.min(...xv);
  let maxX = Math.max(...xv);
  let minY = Math.min(...yv);
  let maxY = Math.max(...yv);
  if (minX === maxX) { minX -= 1; maxX += 1; }
  if (minY === maxY) { minY -= 1; maxY += 1; }
  const xPad = (maxX - minX) * 0.05;
  const yPad = (maxY - minY) * 0.05;
  minX -= xPad; maxX += xPad;
  minY -= yPad; maxY += yPad;

  const pad = { l: 58, r: 26, t: 24, b: 48 };
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b;
  const xLowerBetter = scatterTaskMeta(scatterState.xTask).lowerBetter;
  const yLowerBetter = scatterTaskMeta(scatterState.yTask).lowerBetter;
  function sx(x) {
    const t = (x - minX) / (maxX - minX || 1);
    const dir = xLowerBetter ? (1 - t) : t;
    return pad.l + dir * plotW;
  }
  function sy(y) {
    const t = (y - minY) / (maxY - minY || 1);
    const up = yLowerBetter ? (1 - t) : t;
    return h - pad.b - up * plotH;
  }

  // grid + ticks
  const tickCount = 5;
  const xTicks = Array.from({ length: tickCount + 1 }, (_, i) => minX + ((maxX - minX) * i) / tickCount);
  const yTicks = Array.from({ length: tickCount + 1 }, (_, i) => minY + ((maxY - minY) * i) / tickCount);
  ctx.strokeStyle = "#e5edf5";
  ctx.lineWidth = 1;
  xTicks.forEach((t) => {
    const x = sx(t);
    ctx.beginPath();
    ctx.moveTo(x, pad.t);
    ctx.lineTo(x, h - pad.b);
    ctx.stroke();
  });
  yTicks.forEach((t) => {
    const y = sy(t);
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(w - pad.r, y);
    ctx.stroke();
  });

  // axes
  ctx.strokeStyle = "#dbe5ef";
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, h - pad.b);
  ctx.lineTo(w - pad.r, h - pad.b);
  ctx.stroke();

  ctx.font = "10px " + getComputedStyle(document.body).fontFamily;
  ctx.fillStyle = "#94a3b8";
  ctx.textAlign = "center";
  xTicks.forEach((t) => ctx.fillText(t.toFixed(2).replace(/\.00$/, ""), sx(t), h - pad.b + 14));
  ctx.textAlign = "right";
  yTicks.forEach((t) => ctx.fillText(t.toFixed(2).replace(/\.00$/, ""), pad.l - 8, sy(t)));

  // axis titles (with metric info)
  ctx.font = "11px " + getComputedStyle(document.body).fontFamily;
  ctx.fillStyle = "#64748b";
  ctx.textAlign = "center";
  ctx.fillText(xs.label, w / 2, h - 12);
  ctx.save();
  ctx.translate(15, h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(ys.label, 0, 0);
  ctx.restore();

  scatterState.hitPoints = [];
  visible.forEach((p) => {
    const px = sx(p.x);
    const py = sy(p.y);
    const r = p.model === selectedModel ? 9 : 6;
    ctx.beginPath();
    ctx.fillStyle = p.color + "cc";
    ctx.strokeStyle = "#fff";
    ctx.lineWidth = 1.5;
    ctx.arc(px, py, r, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    scatterState.hitPoints.push({ ...p, px, py, r });
  });
  $("#scatterHint").textContent = `Visible models: ${visible.length}. Lower-better metrics use reversed axes so top-right remains best.`;
}

function renderAll(full = true) {
  if (full) { renderPills(); renderControls(); }
  else { renderControls(); }
  renderLeaderboard();
  drawBarChart();
  drawScatter();
}

function initNav() {
  const go = (p) => {
    $$(".page").forEach((x) => x.classList.toggle("on", x.id === "page-" + p));
    $$(".nav-btn").forEach((x) => x.classList.toggle("on", x.dataset.page === p));
    window.scrollTo({ top: 0, behavior: "smooth" });
    setTimeout(() => { drawBarChart(); drawScatter(); drawGenerationCharts(); }, 60);
  };
  $$("[data-page]").forEach((b) => (b.onclick = () => go(b.dataset.page)));
  $$("[data-page-go]").forEach((b) => (b.onclick = () => go(b.dataset.pageGo)));
  $("[data-scroll-leader]").onclick = () => $("#leaderboard").scrollIntoView({ behavior: "smooth", block: "start" });
  $(".brand").onclick = () => go("home");
}

function initTaskBrowser() {
  const tasks = [
    {
      id: "01",
      name: "Taxonomy Classification",
      cat: "Understanding axis",
      desc: "Predict hierarchical viral taxonomy labels from Kingdom to Family, across ALL/DNA/RNA cohorts and genus-disjoint/temporal splits.",
      badges: ["6 scenarios", "5 ranks", "Macro-F1"],
      metrics: [["6", "scenarios"], ["5", "taxonomy ranks"], ["74", "models in table"]],
      note: "This module connects taxonomy classification results with the leaderboard controls, enabling users to compare models across virus cohorts, split strategies, and taxonomic ranks through searchable tables and synchronized visualizations.",
    },
    {
      id: "02",
      name: "Host Prediction",
      cat: "Understanding axis",
      desc: "Predict standardized host categories from viral sequences to test whether models capture ecological virus–host signals.",
      badges: ["6 scenarios", "Macro-F1", "Genus/Temporal"],
      metrics: [["6", "scenarios"], ["1", "primary metric"], ["ALL/DNA/RNA", "cohorts"]],
      note: "This module presents host prediction performance across ALL, DNA, and RNA virus cohorts, allowing users to explore model robustness under genus-disjoint and temporal splits with interactive ranking, filtering, and chart updates.",
    },
    {
      id: "03",
      name: "Genome Modeling",
      cat: "Generation axis",
      desc: "Compute BPB for next-token modeling on full viral genome fragments stratified into short, medium, and long buckets.",
      badges: ["3 buckets", "BPB", "Lower better"],
      metrics: [["3", "length buckets"], ["28", "generative models"], ["BPB", "metric"]],
      note: "This module visualizes genome modeling results across short, medium, and long sequence regimes, helping users compare generative sequence modeling performance through BPB-centered rankings and dynamic metric views.",
    },
    {
      id: "04",
      name: "CDS Completion",
      cat: "Generation axis",
      desc: "Complete protein-coding sequences from a 129-bp prefix and evaluate exactness, edit distance, CDS validity, and k-mer distribution.",
      badges: ["3 buckets", "5 metrics", "CDS validity"],
      metrics: [["3", "length buckets"], ["5", "metrics"], ["129bp", "prompt"]],
      note: "This module summarizes CDS completion performance under length-stratified settings, linking sequence-level fidelity and coding-sequence validity metrics to interactive model comparison and diagnostic visualizations.",
    },
  ];
  let active = 0;
  function render() {
    const t = tasks[active];
    $("#taskBrowser").innerHTML =
      `<div class="task-list">${tasks.map((x, i) => `<div class="task-item ${i === active ? "on" : ""}" data-i="${i}"><div class="cat">${x.cat}</div><div class="name">${x.name}</div></div>`).join("")}</div>` +
      `<div class="task-detail"><div class="eyebrow">${t.cat}</div><h2>${t.name}</h2><p>${t.desc}</p><div class="badge-row">${t.badges.map((b) => `<span class="badge">${b}</span>`).join("")}</div><div class="metric-grid">${t.metrics.map((m) => `<div class="metric-box"><div class="big">${m[0]}</div><div class="txt">${m[1]}</div></div>`).join("")}</div><div class="section-title" style="margin-top:30px"><h2>INTERACTIVE TASK MODULE</h2></div><p style="color:var(--ink3)">${t.note}</p></div>`;
    $$(".task-item").forEach((el) => (el.onclick = () => { active = Number(el.dataset.i); render(); }));
  }
  render();
}

function drawGenerationCharts() {
  if (!$("#genomeChart")) return;
  const gRows = (RAW.genome[genBucket] || []).map((r) => ({ ...r, _score: parseVal(r.mean) })).sort((a, b) => a._score - b._score);
  drawBars($("#genomeChart"), gRows, "mean", { task: "genome", c1: "#4f46e5", c2: "#0f9f9a" });
  const cRows = (RAW.cds[cdsBucket] || []).map((r) => ({ ...r, _score: parseVal(r.kmer_JSD) })).sort((a, b) => a._score - b._score);
  drawBars($("#cdsChart"), cRows, "kmer_JSD", { task: "cds", c1: "#e76f51", c2: "#f59e0b" });
}

function initGenerationPills() {
  const gpWrap = $("#genomeBucketPills");
  const cpWrap = $("#cdsBucketPills");
  if (!gpWrap || !cpWrap) return;
  const gp = Object.keys(RAW.genome);
  const cp = Object.keys(RAW.cds);
  gpWrap.innerHTML = gp.map((b) => `<button class="pill ${b === genBucket ? "on" : ""}" data-b="${b}">${b.replace("genome-", "")}</button>`).join("");
  cpWrap.innerHTML = cp.map((b) => `<button class="pill ${b === cdsBucket ? "on" : ""}" data-b="${b}">${b.replace("cds-", "")}</button>`).join("");
  $$("#genomeBucketPills .pill").forEach((el) => (el.onclick = () => { genBucket = el.dataset.b; initGenerationPills(); drawGenerationCharts(); }));
  $$("#cdsBucketPills .pill").forEach((el) => (el.onclick = () => { cdsBucket = el.dataset.b; initGenerationPills(); drawGenerationCharts(); }));
}

function renderLoadError(err) {
  const msg = `Data loading failed: ${err?.message || String(err)}`;
  if ($("#statGrid")) {
    $("#statGrid").innerHTML = `<div class="stat"><div class="label">Error</div><div class="sub">${msg}</div></div>`;
  }
  console.error(msg, err);
}

async function fetchJson(path) {
  const res = await fetch(encodeURI(path));
  if (!res.ok) throw new Error(`${path} -> ${res.status}`);
  return res.json();
}

async function loadData() {
  const [taxonomy, host, genome, cds] = await Promise.all([
    fetchJson(DATA_FILES.taxonomy),
    fetchJson(DATA_FILES.host),
    fetchJson(DATA_FILES.genome),
    fetchJson(DATA_FILES.cds),
  ]);
  RAW = { taxonomy, host, genome, cds };
}

function initHeroDnaAnimation() {
  const container = $("#dnaLottie");
  if (!container || !window.lottie) return;
  const dnaAnim = window.lottie.loadAnimation({
    container,
    renderer: "svg",
    loop: true,
    autoplay: true,
    path: "assets/js/DNA_loop_animation.json",
    rendererSettings: {
      preserveAspectRatio: "xMidYMid meet",
    },
  });
  dnaAnim.setSpeed(0.6);
}

async function loadScriptOnce(src) {
  return new Promise((resolve, reject) => {
    const existing = document.querySelector(`script[data-dynamic-src="${src}"]`);
    if (existing) {
      existing.addEventListener("load", () => resolve(true), { once: true });
      existing.addEventListener("error", () => reject(new Error(`Failed to load ${src}`)), { once: true });
      return;
    }
    const s = document.createElement("script");
    s.src = src;
    s.async = true;
    s.dataset.dynamicSrc = src;
    s.onload = () => resolve(true);
    s.onerror = () => reject(new Error(`Failed to load ${src}`));
    document.head.appendChild(s);
  });
}

async function ensureLottieLoaded() {
  if (window.lottie) return true;
  const sources = [
    "https://cdnjs.cloudflare.com/ajax/libs/bodymovin/5.12.2/lottie.min.js",
    "https://cdn.jsdelivr.net/npm/lottie-web@5.12.2/build/player/lottie.min.js",
  ];
  for (const src of sources) {
    try {
      await loadScriptOnce(src);
      if (window.lottie) return true;
    } catch (_) {
      // Try next CDN source.
    }
  }
  return Boolean(window.lottie);
}

window.addEventListener("resize", () => {
  drawBarChart();
  drawScatter();
  drawGenerationCharts();
});

async function bootstrap() {
  try {
    await ensureLottieLoaded();
    initHeroDnaAnimation();
    await loadData();
    initStats();
    renderTaskCards();
    initNav();
    initTaskBrowser();
    initGenerationPills();
    initScatterControls();
    initScatterHover();
    renderAll();
    drawGenerationCharts();
  } catch (err) {
    renderLoadError(err);
  }
}

bootstrap();
