const STORAGE_KEY = "seqalign.introSeen.v1";

const root = document.getElementById("intro");
const steps = [...root.querySelectorAll(".intro-step")];
const dotsBox = document.getElementById("introDots");
const prevBtn = root.querySelector('[data-intro="prev"]');
const nextBtn = root.querySelector('[data-intro="next"]');
const openBtn = document.getElementById("introOpen");

let current = 0;
let lastFocus = null;

const MINI = {
  a: "GCAT",
  b: "GAT",
  values: [
    [0, -1, -2, -3],
    [-1, 1, 0, -1],
    [-2, 0, 0, -1],
    [-3, -1, 1, 0],
    [-4, -2, 0, 2],
  ],
  arrows: [
    ["", "←", "←", "←"],
    ["↑", "↖", "←", "←"],
    ["↑", "↑", "↖", "↖"],
    ["↑", "↑", "↖", "←"],
    ["↑", "↑", "↑", "↖"],
  ],
  path: [[4, 3], [3, 2], [2, 1], [1, 1], [0, 0]],
};

function buildMiniGrid() {
  const grid = document.getElementById("miniGrid");
  const rows = MINI.values.length;
  const cols = MINI.values[0].length;
  const pathOrder = new Map(MINI.path.map(([i, j], k) => [`${i},${j}`, k]));
  const frag = document.createDocumentFragment();

  const head = (text) => {
    const el = document.createElement("div");
    el.className = "mg-h";
    el.textContent = text;
    return el;
  };

  frag.append(head(""));
  for (let j = 0; j < cols; j++) frag.append(head(j === 0 ? "–" : MINI.b[j - 1]));

  for (let i = 0; i < rows; i++) {
    frag.append(head(i === 0 ? "–" : MINI.a[i - 1]));
    for (let j = 0; j < cols; j++) {
      const el = document.createElement("div");
      el.className = "mg-c";
      el.style.setProperty("--d", String(i + j));
      const p = pathOrder.get(`${i},${j}`);
      if (p !== undefined) {
        el.classList.add("on-path");
        el.style.setProperty("--p", String(p));
      }
      const arr = MINI.arrows[i][j];
      el.innerHTML = `${arr ? `<span class="mg-a">${arr}</span>` : ""}${MINI.values[i][j]}`;
      frag.append(el);
    }
  }
  grid.style.setProperty("--cols", String(cols + 1));
  grid.replaceChildren(frag);
}

function playMini() {
  const grid = document.getElementById("miniGrid");
  grid.classList.remove("play");
  void grid.offsetWidth;
  grid.classList.add("play");
}

function buildDots() {
  steps.forEach((_, k) => {
    const b = document.createElement("button");
    b.type = "button";
    b.className = "dot";
    b.setAttribute("aria-label", `${k + 1}. adım`);
    b.addEventListener("click", () => go(k));
    dotsBox.append(b);
  });
}

function go(k) {
  current = Math.max(0, Math.min(steps.length - 1, k));
  steps.forEach((s, idx) => { s.hidden = idx !== current; });
  [...dotsBox.children].forEach((d, idx) => d.setAttribute("aria-current", idx === current ? "step" : "false"));

  prevBtn.disabled = current === 0;
  const last = current === steps.length - 1;
  nextBtn.textContent = last ? "Uygulamaya geç" : "İleri";

  root.querySelector(".intro-body").scrollTo({ top: 0 });
  steps[current].classList.remove("enter");
  void steps[current].offsetWidth;
  steps[current].classList.add("enter");

  if (steps[current].dataset.step === "2") playMini();
}

function open(fromUser = false) {
  lastFocus = document.activeElement;
  if (fromUser) document.documentElement.classList.add("motion-ok");
  root.hidden = false;
  document.body.classList.add("intro-open");
  go(0);
  nextBtn.focus();
}

function close() {
  root.hidden = true;
  document.body.classList.remove("intro-open");
  try { localStorage.setItem(STORAGE_KEY, "1"); } catch { }
  (lastFocus instanceof HTMLElement ? lastFocus : document.getElementById("run"))?.focus();
}

function seenBefore() {
  try { return localStorage.getItem(STORAGE_KEY) === "1"; } catch { return false; }
}

root.addEventListener("click", (e) => {
  const action = e.target.closest("[data-intro]")?.dataset.intro;
  if (action === "close") close();
  if (action === "prev") go(current - 1);
  if (action === "next") (current === steps.length - 1 ? close() : go(current + 1));
});

document.addEventListener("keydown", (e) => {
  if (root.hidden) return;
  if (e.key === "Escape") close();
  else if (e.key === "ArrowRight") go(current + 1);
  else if (e.key === "ArrowLeft") go(current - 1);
});

openBtn.addEventListener("click", () => open(true));
document.getElementById("miniReplay").addEventListener("click", () => {
  document.documentElement.classList.add("motion-ok");
  playMini();
});

buildMiniGrid();
buildDots();
if (!seenBefore()) open(false);