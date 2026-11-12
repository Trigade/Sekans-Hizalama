import init, { align, substitutionScore, alignGlobalScore, version } from "./pkg/seq_align_core.js";
import { MatrixView, DIR } from "./animate.js";
import { runBenchmark, randomDna } from "./bench.js";

// ---------------------------------------------------------------------------
// Sabitler
// ---------------------------------------------------------------------------

/** Bu boyutun üzerindeki matrisler DOM'a çizilmez (yalnızca sonuç gösterilir). */
const MAX_VIS_CELLS = 120 * 120;

/** Hizalama çıktısında satır başına karakter. */
const ALN_BLOCK = 60;

const BENCH_SIZES = [500, 1000, 2000, 4000];

/** Gotoh katmanlarında "erişilemez" anlamına gelen değerler (Rust tarafında i32::MIN / 4). */
const NEG_LIMIT = -1e8;

const ALGORITHMS = {
  "needleman-wunsch": {
    name: "Needleman–Wunsch",
    subtitle: "Needleman–Wunsch global hizalama · DP matrisi adım adım",
  },
  "smith-waterman": {
    name: "Smith–Waterman",
    subtitle: "Smith–Waterman lokal hizalama · en benzer alt bölgeyi bulur",
  },
  "gotoh-global": {
    name: "Gotoh · global",
    subtitle: "Gotoh global hizalama · afin boşluk cezası, üç matris (M, X, Y)",
  },
  "gotoh-local": {
    name: "Gotoh · lokal",
    subtitle: "Gotoh lokal hizalama · afin boşluk cezası, üç matris (M, X, Y)",
  },
  hirschberg: {
    name: "Hirschberg",
    subtitle: "Hirschberg global hizalama · böl ve yönet, lineer bellek",
  },
};

/** Deterministik mutasyon: yaklaşık `rate` oranında ikame, ekleme ve silme. */
function mutate(seq, rate, seed) {
  let x = seed >>> 0;
  const rnd = () => ((x = (Math.imul(x, 1664525) + 1013904223) >>> 0) / 4294967296);
  const bases = "ACGT";
  let out = "";
  for (const c of seq) {
    const r = rnd();
    if (r < rate / 3) continue;
    if (r < (2 * rate) / 3) out += bases[(rnd() * 4) | 0];
    else if (r < rate) out += c + bases[(rnd() * 4) | 0];
    else out += c;
  }
  return out;
}

const LONG_A = randomDna(3000, 7);
const LONG_B = mutate(LONG_A, 0.12, 11);

const SIMPLE = (match, mismatch) => ({ substitution: "simple", match, mismatch });
const LINEAR = (gapExtend) => ({ gapModel: "linear", gapOpen: 0, gapExtend });
const AFFINE = (gapOpen, gapExtend) => ({ gapModel: "affine", gapOpen, gapExtend });

const SAMPLES = [
  { label: "Wikipedia örneği", a: "GATTACA", b: "GCATGCU",
    opts: { mode: "global", ...SIMPLE(1, -1), ...LINEAR(-1) } },
  { label: "Lokal (Wikipedia SW)", a: "TGTTACGG", b: "GGTTGACTA",
    opts: { mode: "local", ...SIMPLE(3, -3), ...LINEAR(-2) } },
  { label: "Gömülü motif", a: "TTTTTGAATTCAAAAAA", b: "CCGAATTCGG",
    opts: { mode: "local", ...SIMPLE(2, -1), ...LINEAR(-2) } },
  { label: "Afin boşluk (Gotoh)", a: "ACGTTTTACG", b: "ACGACG",
    opts: { mode: "global", ...SIMPLE(2, -1), ...AFFINE(-4, -1) } },
  { label: "Protein · BLOSUM62", a: "HEAGAWGHEE", b: "PAWHEAE",
    opts: { mode: "global", substitution: "blosum62", match: 1, mismatch: -1, ...LINEAR(-4) } },
  { label: "DNA mutasyonu", a: "ACGTTGCAATGCCGTAGCTA", b: "ACGTGCAATGCCCGTAGTA",
    opts: { mode: "global", ...SIMPLE(1, -1), ...LINEAR(-1) } },
  {
    label: "Hemoglobin α / β · BLOSUM62 + afin",
    a: "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK",
    b: "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQ",
    opts: { mode: "global", substitution: "blosum62", match: 1, mismatch: -1, ...AFFINE(-11, -1) },
  },
  { label: "Uzun DNA · 3.000 harf (Hirschberg)", a: LONG_A, b: LONG_B,
    opts: { mode: "global", ...SIMPLE(1, -1), ...LINEAR(-1) } },
];

const prefersReducedMotion = matchMedia("(prefers-reduced-motion: reduce)").matches;

// ---------------------------------------------------------------------------
// DOM
// ---------------------------------------------------------------------------

const $ = (id) => document.getElementById(id);

const ui = {
  seqA: $("seqA"), seqB: $("seqB"), lenA: $("lenA"), lenB: $("lenB"),
  samples: $("samples"),
  match: $("match"), mismatch: $("mismatch"), gap: $("gap"), gapOpen: $("gapOpen"), speed: $("speed"),
  matchOut: $("matchOut"), mismatchOut: $("mismatchOut"), gapOut: $("gapOut"),
  gapOpenOut: $("gapOpenOut"), speedOut: $("speedOut"),
  gapLabel: $("gapLabel"), gapOpenRow: $("gapOpenRow"), simpleRows: $("simpleRows"),
  scoringHint: $("scoringHint"),
  arrows: $("arrows"),
  run: $("run"), skip: $("skip"), error: $("error"),
  stScore: $("stScore"), stIdentity: $("stIdentity"), stMatches: $("stMatches"),
  stMismatches: $("stMismatches"), stGaps: $("stGaps"), stTime: $("stTime"), stAlgo: $("stAlgo"),
  alignment: $("alignment"),
  matrix: $("matrix"), cellInfo: $("cellInfo"), matrixNote: $("matrixNote"),
  benchRun: $("benchRun"), benchOut: $("benchOut"),
  version: $("version"),
  subtitle: $("subtitle"),
  modes: document.querySelectorAll('input[name="mode"]'),
  substs: document.querySelectorAll('input[name="subst"]'),
  gapModels: document.querySelectorAll('input[name="gapModel"]'),
};

const view = new MatrixView(ui.matrix);

/** Uygulama durumu: tek bir aktif çalıştırma ve onun iptal kontrolcüsü. */
const state = {
  result: null,
  a: "",
  b: "",
  opts: null,
  mode: "global",
  controller: null,
};

function readMode() {
  return document.querySelector('input[name="mode"]:checked')?.value === "local" ? "local" : "global";
}

function setMode(mode) {
  for (const r of ui.modes) r.checked = r.value === mode;
}

const checkedValue = (name, fallback) =>
  document.querySelector(`input[name="${name}"]:checked`)?.value ?? fallback;

function setRadio(list, value) {
  for (const r of list) r.checked = r.value === value;
}

/** Bir örneğin ayarlarını kontrollere uygular. */
function applyOptions(o) {
  setMode(o.mode);
  setRadio(ui.substs, o.substitution);
  setRadio(ui.gapModels, o.gapModel);
  ui.match.value = String(o.match);
  ui.mismatch.value = String(o.mismatch);
  if (o.gapModel === "affine") ui.gapOpen.value = String(o.gapOpen);
  ui.gap.value = String(o.gapExtend);
  updateScoringUI();
}

/** Seçilen modele göre kaydırıcıları açar/kapatır ve açıklamayı günceller. */
function updateScoringUI() {
  const blosum = checkedValue("subst", "simple") === "blosum62";
  const affine = checkedValue("gapModel", "linear") === "affine";

  ui.simpleRows.classList.toggle("disabled", blosum);
  ui.match.disabled = blosum;
  ui.mismatch.disabled = blosum;
  ui.gapOpenRow.hidden = !affine;
  ui.gapLabel.textContent = affine ? "Uzatma" : "Boşluk";

  const hints = [];
  if (blosum) hints.push("BLOSUM62: amino asit çiftlerinin evrimsel değişim sıklığından türetilmiş puan tablosu. Protein dizileri için.");
  if (affine) hints.push("Afin boşluk: k uzunluklu boşluk = açma + k × uzatma. Tek uzun boşluk, birçok kısa boşluktan ucuzdur.");
  ui.scoringHint.textContent = hints.join(" ");
  ui.scoringHint.hidden = hints.length === 0;
  updateOutputs();
}

// ---------------------------------------------------------------------------
// Yardımcılar
// ---------------------------------------------------------------------------

/** FASTA başlıklarını (">" ile başlayan satırlar) ve boşlukları atar, büyük harfe çevirir. */
function cleanSequence(raw) {
  return raw
    .split(/\r?\n/)
    .filter((line) => !line.trimStart().startsWith(">"))
    .join("")
    .replace(/\s+/g, "")
    .toUpperCase();
}

const signed = (n) => (n > 0 ? `+${n}` : n < 0 ? `−${Math.abs(n)}` : "0");

/** Kontrollerden WASM'a gönderilecek seçenek nesnesini üretir. */
function readOptions() {
  const affine = checkedValue("gapModel", "linear") === "affine";
  return {
    mode: readMode(),
    substitution: checkedValue("subst", "simple"),
    match: Number(ui.match.value),
    mismatch: Number(ui.mismatch.value),
    gapOpen: affine ? Number(ui.gapOpen.value) : 0,
    gapExtend: Number(ui.gap.value),
  };
}

/** Hız kaydırıcısı (1–10) → saniyedeki hücre sayısı (üstel ölçek). */
function cellsPerSecond(cells) {
  const level = Number(ui.speed.value);
  const base = 12 * Math.pow(1.85, level - 1); // 12 … ~3.1k hücre/sn
  // Büyük matrislerde animasyon makul sürede bitsin (en fazla ~12 sn).
  return Math.max(base, cells / 12);
}

function showError(message) {
  ui.error.textContent = message;
  ui.error.hidden = !message;
}

function setBusy(busy) {
  ui.skip.disabled = !busy;
}

const isAbort = (e) => e?.name === "AbortError";

// ---------------------------------------------------------------------------
// Render: istatistikler ve hizalama çıktısı
// ---------------------------------------------------------------------------

function renderStats(r, ms) {
  ui.stScore.textContent = signed(r.score);
  ui.stIdentity.textContent = `%${r.stats.identity.toFixed(1)}`;
  ui.stMatches.textContent = r.stats.matches;
  ui.stMismatches.textContent = r.stats.mismatches;
  ui.stGaps.textContent = r.stats.gaps;
  ui.stTime.textContent = ms < 1 ? `${(ms * 1000).toFixed(0)} µs` : `${ms.toFixed(2)} ms`;
  ui.stAlgo.textContent = ALGORITHMS[r.algorithm]?.name ?? r.algorithm;
}

function clearStats() {
  for (const el of [ui.stScore, ui.stIdentity, ui.stMatches, ui.stMismatches, ui.stGaps, ui.stTime, ui.stAlgo]) {
    el.textContent = "–";
  }
}

function renderAlignment(r, { animate }) {
  const A = r.aligned_a;
  const B = r.aligned_b;

  if (A.length === 0) {
    const p = document.createElement("p");
    p.className = "empty";
    p.textContent = state.mode === "local"
      ? "Pozitif skorlu bir lokal hizalama bulunamadı: diziler arasında yeterince benzer bir bölge yok."
      : "Her iki dizi de boş.";
    ui.alignment.classList.remove("reveal");
    ui.alignment.replaceChildren(p);
    return;
  }

  const frag = document.createDocumentFragment();
  let k = 0;
  let posA = r.start[0];
  let posB = r.start[1];

  for (let start = 0; start < A.length; start += ALN_BLOCK) {
    const end = Math.min(A.length, start + ALN_BLOCK);
    const block = document.createElement("div");
    block.className = "aln-block";

    const rowA = document.createElement("div");
    const rowM = document.createElement("div");
    const rowB = document.createElement("div");
    rowA.className = "aln-row";
    rowM.className = "aln-row aln-mid";
    rowB.className = "aln-row";

    const labelA = posA + 1;
    const labelB = posB + 1;

    for (let p = start; p < end; p++) {
      const x = A[p];
      const y = B[p];
      const cls = x === "-" || y === "-" ? "g" : x === y ? "m" : "x";
      const mark = cls === "m" ? "|" : cls === "x" ? "·" : " ";
      if (x !== "-") posA++;
      if (y !== "-") posB++;

      for (const [row, ch] of [[rowA, x], [rowM, mark], [rowB, y]]) {
        const s = document.createElement("span");
        s.className = cls;
        s.textContent = ch;
        if (animate) s.style.setProperty("--k", String(k));
        row.append(s);
      }
      k++;
    }

    block.append(label(labelA), rowA, label(""), rowM, label(labelB), rowB);
    frag.append(block);
  }

  ui.alignment.classList.toggle("reveal", animate);
  ui.alignment.replaceChildren(frag);

  function label(text) {
    const el = document.createElement("span");
    el.className = "aln-label";
    el.textContent = String(text);
    return el;
  }
}

// ---------------------------------------------------------------------------
// Hücre açıklaması (hover)
// ---------------------------------------------------------------------------

const fmtScore = (v) => (v <= NEG_LIMIT ? "−∞" : String(v));

/** "taban ± fark" ifadesini okunur hâle getirir; taban −∞ ise fark yazılmaz. */
function plus(base, delta) {
  if (base <= NEG_LIMIT) return "−∞";
  return `${base} ${delta < 0 ? "−" : "+"} ${Math.abs(delta)}`;
}

function substitutionLabel(ca, cb, sub) {
  if (state.opts.substitution === "blosum62") return `BLOSUM62(${ca}, ${cb}) = ${sub}`;
  return `${ca} vs ${cb} → ${ca === cb ? "eşleşme" : "uyuşmazlık"} (${signed(sub)})`;
}

function explainCell(i, j) {
  return state.result.layers ? explainAffine(i, j) : explainLinear(i, j);
}

function explainLinear(i, j) {
  const { result: r, a, b, opts, mode } = state;
  const { cols, matrix, dirs } = r;
  const v = matrix[i * cols + j];
  const local = mode === "local";
  const gap = opts.gapExtend;

  if (local && (i === 0 || j === 0)) {
    return `<b>H(${i},${j}) = 0</b> — lokal hizalamada kenarlar her zaman 0: hizalama herhangi bir yerden başlayabilir`;
  }
  if (i === 0 && j === 0) return `<b>H(0,0) = 0</b> — başlangıç hücresi`;
  if (i === 0) return `<b>H(0,${j}) = ${j} × (${gap}) = ${v}</b> — B'nin ilk ${j} harfi boşlukla eşleşir`;
  if (j === 0) return `<b>H(${i},0) = ${i} × (${gap}) = ${v}</b> — A'nın ilk ${i} harfi boşlukla eşleşir`;

  const ca = a[i - 1];
  const cb = b[j - 1];
  const sub = substitutionScore(ca, cb, opts);
  const chosen = dirs[i * cols + j];

  const term = (dir, sym, base, delta, note) => {
    const txt = `${sym} ${plus(base, delta)} = ${base + delta} <small>(${note})</small>`;
    return dir === chosen ? `<span class="win">${txt}</span>` : txt;
  };

  const zero = local
    ? (chosen === DIR.NONE ? `<span class="win">0 <small>(yeniden başla)</small></span>, ` : `0, `)
    : "";
  const peak = local && r.end[0] === i && r.end[1] === j ? ` &nbsp;← <b>matrisin maksimumu, traceback buradan başlar</b>` : "";

  return (
    `<b>H(${i},${j})</b> &nbsp;${substitutionLabel(ca, cb, sub)}${peak}<br>` +
    `max( ${zero}${term(DIR.DIAG, "↖", matrix[(i - 1) * cols + (j - 1)], sub, "ikame")}, ` +
    `${term(DIR.UP, "↑", matrix[(i - 1) * cols + j], gap, "A'da harf, B'de boşluk")}, ` +
    `${term(DIR.LEFT, "←", matrix[i * cols + (j - 1)], gap, "A'da boşluk, B'de harf")} ) = <b>${v}</b>`
  );
}

function explainAffine(i, j) {
  const { result: r, a, b, opts, mode } = state;
  const { cols, matrix } = r;
  const { m, x, y } = r.layers;
  const local = mode === "local";
  const k = i * cols + j;
  const v = matrix[k];
  const o = opts.gapOpen;
  const e = opts.gapExtend;
  const oe = o + e;

  const best = Math.max(m[k], x[k], y[k]);
  const win = (val, txt) => (val === best && (!local || best > 0) ? `<span class="win">${txt}</span>` : txt);
  const cellLine =
    `Hücre = max( ${local ? "0, " : ""}${win(m[k], `M ${fmtScore(m[k])}`)}, ` +
    `${win(x[k], `X ${fmtScore(x[k])}`)}, ${win(y[k], `Y ${fmtScore(y[k])}`)} ) = <b>${v}</b>`;

  if (i === 0 || j === 0) {
    const note = local
      ? "lokal hizalamada kenarlar 0"
      : i === 0 && j === 0
        ? "başlangıç hücresi"
        : `ilk ${i || j} harf tek bir boşlukla eşleşir: ${o} + ${i || j} × (${e}) = ${o + (i || j) * e}`;
    return `<b>H(${i},${j})</b> — ${note}<br>${cellLine}`;
  }

  const ca = a[i - 1];
  const cb = b[j - 1];
  const sub = substitutionScore(ca, cb, opts);
  const d = (i - 1) * cols + (j - 1);
  const u = (i - 1) * cols + j;
  const l = i * cols + (j - 1);
  const peak = local && r.end[0] === i && r.end[1] === j ? ` &nbsp;← <b>en yüksek M, traceback buradan başlar</b>` : "";

  return (
    `<b>H(${i},${j})</b> &nbsp;${substitutionLabel(ca, cb, sub)} · açma ${o}, uzatma ${e}${peak}<br>` +
    `M = ${signed(sub)} + max( ${local ? "0, " : ""}M↖ ${fmtScore(m[d])}, X↖ ${fmtScore(x[d])}, Y↖ ${fmtScore(y[d])} ) = ${fmtScore(m[k])} <small>(harf–harf)</small><br>` +
    `X = max( M↑ ${plus(m[u], oe)}, X↑ ${plus(x[u], e)}, Y↑ ${plus(y[u], oe)} ) = ${fmtScore(x[k])} <small>(B'de boşluk: açılır ya da uzar)</small><br>` +
    `Y = max( M← ${plus(m[l], oe)}, X← ${plus(x[l], oe)}, Y← ${plus(y[l], e)} ) = ${fmtScore(y[k])} <small>(A'da boşluk)</small><br>` +
    cellLine
  );
}

function onCellHover(ev) {
  const cell = ev.target.closest(".cell.on");
  if (!cell || !state.result) return;
  const i = Number(cell.dataset.i);
  const j = Number(cell.dataset.j);
  view.highlight(i, j);
  ui.cellInfo.innerHTML = explainCell(i, j);
}

// ---------------------------------------------------------------------------
// Ana akış
// ---------------------------------------------------------------------------

/**
 * Hizalamayı çalıştırır.
 * @param {{ animate: boolean }} opts
 */
async function run({ animate }) {
  state.controller?.abort();
  const controller = new AbortController();
  state.controller = controller;
  const { signal } = controller;

  showError("");
  view.unhighlight();

  const a = cleanSequence(ui.seqA.value);
  const b = cleanSequence(ui.seqB.value);
  const opts = readOptions();
  const mode = opts.mode;

  let result;
  let ms;
  try {
    const t0 = performance.now();
    result = align(a, b, opts);
    ms = performance.now() - t0;
  } catch (e) {
    showError(e?.message ?? String(e));
    state.result = null;
    view.clear();
    clearStats();
    ui.alignment.replaceChildren();
    setBusy(false);
    return;
  }

  Object.assign(state, { result, a, b, opts, mode });
  renderStats(result, ms);
  ui.subtitle.textContent = ALGORITHMS[result.algorithm]?.subtitle ?? "";

  const cells = result.rows * result.cols;
  const drawable = cells <= MAX_VIS_CELLS && result.matrix.length === cells;
  ui.matrixNote.hidden = drawable;

  if (!drawable) {
    view.clear();
    const size =
      `${result.rows.toLocaleString("tr-TR")} × ${result.cols.toLocaleString("tr-TR")} ` +
      `(${cells.toLocaleString("tr-TR")} hücre)`;
    ui.matrixNote.textContent = result.algorithm === "hirschberg"
      ? `Matris ${size} görselleştirme için çok büyük. Hirschberg algoritması tam matrisi hiç oluşturmadan, ` +
        `yalnızca iki satır tutarak böl ve yönet ile hizaladı. Sonuç yukarıda.`
      : `Matris ${size} görselleştirme için çok büyük. Hizalama sonucu yukarıda.`;
    renderAlignment(result, { animate: false });
    setBusy(false);
    return;
  }

  view.render(result, a, b);
  ui.cellInfo.textContent = "Bir hücrenin üzerine gel ya da dokun: değerinin nasıl hesaplandığını gör.";

  if (!animate) {
    view.revealAll();
    renderAlignment(result, { animate: false });
    setBusy(false);
    return;
  }

  // Kullanıcı animasyonu açıkça istedi (buton, örnek, mod): sistemdeki
  // "animasyonları azalt" ayarı açık olsa bile bu çalıştırmada oynat.
  document.documentElement.classList.add("motion-ok");
  ui.alignment.replaceChildren();
  setBusy(true);
  try {
    await view.fill({ cellsPerSecond: cellsPerSecond(cells), signal });
    const stepMs = Math.max(15, Math.min(90, 2400 / Math.max(1, result.path.length)));
    await view.trace({ stepMs, signal });
    renderAlignment(result, { animate: true });
  } catch (e) {
    if (!isAbort(e)) throw e;
  } finally {
    if (state.controller === controller) setBusy(false);
  }
}

/** Animasyonu iptal edip sonucu anında gösterir. */
function skip() {
  if (!state.result) return;
  state.controller?.abort();
  view.revealAll();
  renderAlignment(state.result, { animate: false });
  setBusy(false);
}

// ---------------------------------------------------------------------------
// Benchmark
// ---------------------------------------------------------------------------

async function bench() {
  ui.benchRun.disabled = true;
  const progress = document.createElement("p");
  progress.className = "bench-progress";
  ui.benchOut.replaceChildren(progress);

  try {
    const results = await runBenchmark({
      wasmScore: alignGlobalScore,
      sizes: BENCH_SIZES,
      reps: 5,
      onProgress: (msg) => { progress.textContent = msg; },
    });
    renderBench(results);
  } catch (e) {
    progress.textContent = `Benchmark başarısız: ${e?.message ?? e}`;
  } finally {
    ui.benchRun.disabled = false;
  }
}

function renderBench(results) {
  const fmtMs = (x) => (x < 10 ? x.toFixed(2) : x.toFixed(1));
  const geo = Math.exp(results.reduce((acc, r) => acc + Math.log(r.speedup), 0) / results.length);
  const maxMs = Math.max(...results.flatMap((r) => [r.jsMs, r.wasmMs]));
  const largest = results.at(-1);

  const summary = document.createElement("div");
  summary.className = "bench-summary";
  summary.innerHTML =
    `<span>Ortalama hızlanma (geometrik)</span><strong>${geo.toFixed(2)}×</strong>` +
    `<span>En büyük test: ${largest.n.toLocaleString("tr-TR")} × ${largest.n.toLocaleString("tr-TR")} · ` +
    `WASM ${largest.wasmMcups.toFixed(0)} MCUPS, JS ${largest.jsMcups.toFixed(0)} MCUPS</span>`;

  const legend = document.createElement("div");
  legend.className = "bench-legend";
  legend.innerHTML = `<span class="l-js">JavaScript</span><span class="l-wasm">Rust/WASM</span>`;

  const rows = results
    .map((r) => {
      const wj = (r.jsMs / maxMs) * 100;
      const ww = (r.wasmMs / maxMs) * 100;
      return `<tr>
        <td class="num">${r.n.toLocaleString("tr-TR")} × ${r.n.toLocaleString("tr-TR")}</td>
        <td class="num">${(r.cells / 1e6).toFixed(2)} M</td>
        <td class="num">${fmtMs(r.jsMs)} ms</td>
        <td class="num">${fmtMs(r.wasmMs)} ms</td>
        <td class="num"><b>${r.speedup.toFixed(2)}×</b></td>
        <td><div class="bars">
          <div class="bar js" style="width:${wj}%"></div>
          <div class="bar wasm" style="width:${ww}%"></div>
        </div></td>
      </tr>`;
    })
    .join("");

  const wrap = document.createElement("div");
  wrap.className = "bench-wrap";
  wrap.innerHTML = `<table class="bench-table">
    <thead><tr>
      <th>Boyut</th><th>Hücre</th><th>JS</th><th>WASM</th><th>Hızlanma</th><th>Süre</th>
    </tr></thead>
    <tbody>${rows}</tbody>
  </table>`;

  ui.benchOut.replaceChildren(summary, legend, wrap);
}

// ---------------------------------------------------------------------------
// Olaylar
// ---------------------------------------------------------------------------

function updateCounts() {
  ui.lenA.textContent = cleanSequence(ui.seqA.value).length;
  ui.lenB.textContent = cleanSequence(ui.seqB.value).length;
}

function updateOutputs() {
  ui.matchOut.textContent = signed(Number(ui.match.value));
  ui.mismatchOut.textContent = signed(Number(ui.mismatch.value));
  ui.gapOut.textContent = signed(Number(ui.gap.value));
  ui.gapOpenOut.textContent = signed(Number(ui.gapOpen.value));
  ui.speedOut.textContent = ui.speed.value;
}

function debounce(fn, ms) {
  let t;
  return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
}

function loadSample(sample, { animate = true } = {}) {
  ui.seqA.value = sample.a;
  ui.seqB.value = sample.b;
  if (sample.opts) applyOptions(sample.opts);
  updateCounts();
  run({ animate });
}

function bindEvents() {
  for (const sample of SAMPLES) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "chip";
    btn.textContent = sample.label;
    btn.addEventListener("click", () => loadSample(sample));
    ui.samples.append(btn);
  }

  ui.run.addEventListener("click", () => run({ animate: true }));
  ui.skip.addEventListener("click", skip);

  // Skor parametreleri değişince animasyonsuz anında yeniden hesapla.
  const rerun = debounce(() => run({ animate: false }), 120);
  for (const el of [ui.match, ui.mismatch, ui.gap, ui.gapOpen]) {
    el.addEventListener("input", () => { updateOutputs(); rerun(); });
  }
  for (const el of [...ui.substs, ...ui.gapModels]) {
    el.addEventListener("change", () => { updateScoringUI(); rerun(); });
  }
  ui.speed.addEventListener("input", updateOutputs);

  for (const r of ui.modes) {
    r.addEventListener("change", () => run({ animate: true }));
  }

  ui.arrows.addEventListener("change", () => {
    ui.matrix.classList.toggle("show-arrows", ui.arrows.checked);
  });

  for (const el of [ui.seqA, ui.seqB]) {
    el.addEventListener("input", updateCounts);
    el.addEventListener("keydown", (e) => {
      if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) { e.preventDefault(); run({ animate: true }); }
    });
  }

  ui.matrix.addEventListener("pointerover", onCellHover);
  ui.matrix.addEventListener("pointerdown", onCellHover);
  ui.matrix.addEventListener("pointerleave", () => view.unhighlight());

  ui.benchRun.addEventListener("click", bench);
}

// ---------------------------------------------------------------------------
// Başlatma
// ---------------------------------------------------------------------------

async function main() {
  try {
    await init();
  } catch (e) {
    showError(
      "WASM modülü yüklenemedi. Sayfayı bir HTTP sunucusu üzerinden açtığından emin ol " +
      "(ör. web klasöründe: python -m http.server 8000). Ayrıntı: " + (e?.message ?? e),
    );
    ui.run.disabled = true;
    ui.benchRun.disabled = true;
    return;
  }

  ui.version.textContent = `v${version()}`;
  bindEvents();
  updateOutputs();
  // İlk açılış otomatik: sistem ayarı animasyonu kısıtlıyorsa saygı göster.
  loadSample(SAMPLES[0], { animate: !prefersReducedMotion });
}

main();