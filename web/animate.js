/**
 * DP matrisinin DOM'a çizilmesi ve animasyonu.
 *
 * Tasarım notları:
 *  - Hücreler bir kez oluşturulur, animasyon yalnızca CSS sınıfı ekler
 *    (layout'u tetiklemeyen opacity/transform geçişleri → akıcı 60 FPS).
 *  - Doldurma sırası "anti-diyagonal dalga cephesi"dir: (i + j) sabit olan
 *    hücreler birbirinden bağımsızdır; bu, algoritmanın gerçek veri
 *    bağımlılığını görsel olarak gösterir.
 *  - Tüm uzun süreli işlemler AbortSignal ile iptal edilebilir.
 */

export const DIR = Object.freeze({ NONE: 0, DIAG: 1, UP: 2, LEFT: 3 });
const ARROW = ["", "↖", "↑", "←"];

/** Bir sonraki animasyon karesini bekler; iptal edilirse AbortError fırlatır. */
function nextFrame(signal) {
  return new Promise((resolve, reject) => {
    if (signal?.aborted) return reject(abortError());
    const id = requestAnimationFrame(resolve);
    signal?.addEventListener("abort", () => { cancelAnimationFrame(id); reject(abortError()); }, { once: true });
  });
}

function abortError() {
  return new DOMException("Animasyon iptal edildi", "AbortError");
}

export class MatrixView {
  /** @param {HTMLElement} root */
  constructor(root) {
    this.root = root;
    this.cells = [];
    this.rowHeads = [];
    this.colHeads = [];
    this.result = null;
    this.highlighted = [];
  }

  /**
   * Matrisi sıfırdan oluşturur (tüm hücreler gizli başlar).
   * @param {object} result  WASM'dan dönen AlignmentResult
   * @param {string} a       normalize edilmiş dizi A (satırlar)
   * @param {string} b       normalize edilmiş dizi B (sütunlar)
   */
  render(result, a, b) {
    this.result = result;
    this.a = a;
    this.b = b;
    const { rows, cols, matrix, dirs } = result;

    this.root.style.setProperty("--cols", String(cols + 1));
    this.root.style.setProperty("--cell", `${this.#cellSize(cols + 1)}px`);

    const frag = document.createDocumentFragment();
    this.cells = new Array(rows * cols);
    this.rowHeads = new Array(rows);
    this.colHeads = new Array(cols);

    // Üst başlık satırı: köşe, boşluk sütunu, sonra B'nin harfleri
    frag.append(this.#head("", "hd corner"));
    for (let j = 0; j < cols; j++) {
      const h = this.#head(j === 0 ? "–" : b[j - 1], "hd");
      this.colHeads[j] = h;
      frag.append(h);
    }

    for (let i = 0; i < rows; i++) {
      const h = this.#head(i === 0 ? "–" : a[i - 1], "hd");
      this.rowHeads[i] = h;
      frag.append(h);

      for (let j = 0; j < cols; j++) {
        const k = i * cols + j;
        const el = document.createElement("div");
        el.className = "cell";
        el.dataset.i = String(i);
        el.dataset.j = String(j);
        el.textContent = String(matrix[k]);

        const d = dirs[k];
        if (d !== DIR.NONE) {
          const arr = document.createElement("span");
          arr.className = "arr";
          arr.setAttribute("aria-hidden", "true");
          arr.textContent = ARROW[d];
          el.prepend(arr);
        }
        this.cells[k] = el;
        frag.append(el);
      }
    }

    this.root.replaceChildren(frag);
  }

  clear() {
    this.result = null;
    this.cells = [];
    this.root.replaceChildren();
  }

  /** Hücreleri dalga cephesi sırasıyla açar. */
  async fill({ cellsPerSecond, signal }) {
    const { rows, cols } = this.result;
    const order = waveOrder(rows, cols);
    this.root.classList.remove("instant");

    let shown = 0;
    let budget = 0;
    let last = performance.now();

    while (shown < order.length) {
      const now = await nextFrame(signal);
      budget += ((now - last) / 1000) * cellsPerSecond;
      last = now;

      // En az 1 hücre aç ki çok yavaş hızlarda da ilerleme görünsün.
      const n = Math.max(1, Math.floor(budget));
      budget -= n;
      const end = Math.min(order.length, shown + n);
      for (; shown < end; shown++) this.cells[order[shown]].classList.add("on");
    }
  }

  /** Traceback yolunu (sağ alttan başlayarak) sırayla vurgular. */
  async trace({ stepMs, signal }) {
    const { path, cols } = this.result;
    let last = performance.now();
    let acc = stepMs; // ilk adım hemen

    for (let p = path.length - 1; p >= 0; ) {
      const now = await nextFrame(signal);
      acc += now - last;
      last = now;
      while (acc >= stepMs && p >= 0) {
        const [i, j] = path[p--];
        this.cells[i * cols + j].classList.add("path");
        acc -= stepMs;
      }
    }
  }

  /** Animasyonsuz anında tamamlama (atla / yeniden hesaplama). */
  revealAll({ withPath = true } = {}) {
    if (!this.result) return;
    this.root.classList.add("instant");
    for (const el of this.cells) el.classList.add("on");
    if (withPath) {
      const { path, cols } = this.result;
      for (const [i, j] of path) this.cells[i * cols + j].classList.add("path");
    }
  }

  /** Hover: hücrenin hesaplanmasına katılan komşuları ve başlıkları vurgular. */
  highlight(i, j) {
    this.unhighlight();
    if (!this.result) return;
    const { cols, dirs } = this.result;
    const add = (el, cls) => { el.classList.add(cls); this.highlighted.push([el, cls]); };

    add(this.cells[i * cols + j], "cur");
    add(this.rowHeads[i], "lit");
    add(this.colHeads[j], "lit");

    const d = dirs[i * cols + j];
    const nbs = [
      [DIR.DIAG, i - 1, j - 1],
      [DIR.UP, i - 1, j],
      [DIR.LEFT, i, j - 1],
    ];
    for (const [dir, ni, nj] of nbs) {
      if (ni < 0 || nj < 0) continue;
      add(this.cells[ni * cols + nj], dir === d ? "src" : "nb");
    }
  }

  unhighlight() {
    for (const [el, cls] of this.highlighted) el.classList.remove(cls);
    this.highlighted = [];
  }

  #head(text, cls) {
    const el = document.createElement("div");
    el.className = cls;
    el.textContent = text;
    return el;
  }

  /** Hücre boyutunu panel genişliğine göre ayarlar (22–40 px arası). */
  #cellSize(n) {
    const avail = (this.root.parentElement?.clientWidth ?? 800) - 24;
    const fit = Math.floor(avail / n) - 2;
    return Math.max(22, Math.min(40, fit));
  }
}

/** (i + j) artan sırada hücre indeksleri: anti-diyagonal dalga cephesi. */
function waveOrder(rows, cols) {
  const order = new Uint32Array(rows * cols);
  let k = 0;
  for (let d = 0; d <= rows + cols - 2; d++) {
    const iStart = Math.max(0, d - cols + 1);
    const iEnd = Math.min(d, rows - 1);
    for (let i = iStart; i <= iEnd; i++) order[k++] = i * cols + (d - i);
  }
  return order;
}