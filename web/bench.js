/**
 * JavaScript vs Rust/WASM karşılaştırması.
 *
 * Adil karşılaştırma için:
 *  - JS tarafı Rust'taki `needleman_wunsch_score` ile birebir aynı algoritmadır
 *    (iki satırlık lineer bellek, Int32Array, aynı yön önceliği).
 *  - Diziler sabit tohumlu PRNG ile üretilir → her çalıştırmada aynı girdiler.
 *  - Her boyut için ısınma turu (JIT/WASM derlemesi) + N ölçüm, medyan alınır.
 *  - Her ölçümde iki tarafın skoru karşılaştırılır; uyuşmazlık hata sayılır.
 */

/** Saf JS referans uygulaması: yalnızca skor, O(m) bellek. */
export function nwScoreJS(a, b, match, mismatch, gap) {
  const n = a.length;
  const m = b.length;
  let prev = new Int32Array(m + 1);
  let curr = new Int32Array(m + 1);

  for (let j = 0; j <= m; j++) prev[j] = j * gap;

  for (let i = 1; i <= n; i++) {
    const ca = a[i - 1];
    curr[0] = i * gap;
    for (let j = 1; j <= m; j++) {
      const diag = prev[j - 1] + (ca === b[j - 1] ? match : mismatch);
      const up = prev[j] + gap;
      const left = curr[j - 1] + gap;
      let best = diag > up ? diag : up;
      if (left > best) best = left;
      curr[j] = best;
    }
    const t = prev; prev = curr; curr = t;
  }
  return prev[m];
}

/** Mulberry32: küçük, hızlı, deterministik PRNG. */
function mulberry32(seed) {
  return () => {
    seed |= 0;
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function randomDna(length, seed) {
  const rnd = mulberry32(seed);
  const alphabet = "ACGT";
  let s = "";
  for (let i = 0; i < length; i++) s += alphabet[(rnd() * 4) | 0];
  return s;
}

const median = (xs) => {
  const s = [...xs].sort((x, y) => x - y);
  const mid = s.length >> 1;
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
};

/** UI'ın nefes alabilmesi için olay döngüsüne bir tur verir. */
const yieldToUI = () => new Promise((r) => setTimeout(r, 0));

/**
 * @param {object}   opts
 * @param {Function} opts.wasmScore   (a: string, b: string, m, mm, g) => number
 * @param {number[]} opts.sizes       dizi uzunlukları
 * @param {number}   opts.reps        ölçüm tekrarı
 * @param {Function} opts.onProgress  (mesaj: string) => void
 */
export async function runBenchmark({ wasmScore, sizes, reps = 5, onProgress = () => {} }) {
  const S = { match: 1, mismatch: -1, gap: -1 };
  const enc = new TextEncoder();
  const results = [];

  for (const [idx, n] of sizes.entries()) {
    const aStr = randomDna(n, 1000 + idx);
    const bStr = randomDna(n, 2000 + idx);
    const aBytes = enc.encode(aStr);
    const bBytes = enc.encode(bStr);

    const timeJS = () => {
      const t0 = performance.now();
      const score = nwScoreJS(aBytes, bBytes, S.match, S.mismatch, S.gap);
      return [performance.now() - t0, score];
    };
    const timeWasm = () => {
      const t0 = performance.now();
      const score = wasmScore(aStr, bStr, S.match, S.mismatch, S.gap);
      return [performance.now() - t0, score];
    };

    onProgress(`${n.toLocaleString("tr-TR")} × ${n.toLocaleString("tr-TR")}: ısınma…`);
    await yieldToUI();
    timeJS();
    timeWasm();

    const js = [];
    const wasm = [];
    for (let r = 0; r < reps; r++) {
      onProgress(`${n.toLocaleString("tr-TR")} × ${n.toLocaleString("tr-TR")}: ölçüm ${r + 1}/${reps}`);
      await yieldToUI();

      // Sırayı her turda değiştir: önbellek/termal etkiler bir tarafı kayırmasın.
      const order = r % 2 === 0 ? [timeJS, timeWasm] : [timeWasm, timeJS];
      const out = order.map((f) => f());
      const [tj, sj] = r % 2 === 0 ? out[0] : out[1];
      const [tw, sw] = r % 2 === 0 ? out[1] : out[0];

      if (sj !== sw) throw new Error(`Skor uyuşmazlığı (n=${n}): JS=${sj}, WASM=${sw}`);
      js.push(tj);
      wasm.push(tw);
    }

    const jsMs = median(js);
    const wasmMs = median(wasm);
    results.push({
      n,
      cells: n * n,
      jsMs,
      wasmMs,
      speedup: jsMs / wasmMs,
      // Milyon hücre güncellemesi / saniye (MCUPS): biyoenformatikte standart metrik.
      jsMcups: (n * n) / (jsMs * 1000),
      wasmMcups: (n * n) / (wasmMs * 1000),
    });
  }

  return results;
}