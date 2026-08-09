// Widget builders for the Boolean Gallery's banded sidebar ("layout 2a").
//
// These are the pieces the generic controls.ts does not cover: segmented
// controls, Venn-glyph operation buttons, toggle pills, chips, a disclosure
// row that prints its own summary while closed, and the A/B operand card.
// They are deliberately dumb — every one takes its current value and a change
// callback and returns a handle whose `set()` re-renders it — so the gallery's
// state stays in boolean-gallery.ts and the panel assembly in
// boolean-panel.ts.
//
// Styling lives in styles/boolean-panel.css (all class names start `bgw-`).

export interface Segment<T extends string> {
  value: T;
  text: string;
  /** Full phrase behind an abbreviated label (title tooltip). */
  title?: string;
}

export interface Toggleable<T> {
  el: HTMLElement;
  /** Reflect a value chosen elsewhere; null selects no segment. */
  set(value: T): void;
}

function button(cls: string, text: string, title: string | undefined, onClick: () => void): HTMLButtonElement {
  const b = document.createElement('button');
  b.type = 'button';
  b.className = cls;
  b.textContent = text;
  if (title) b.title = title;
  b.addEventListener('click', onClick);
  return b;
}

/// N-up segmented control. `compact` is the inline variant used for Engine
/// (smaller, neutral selected fill); the default fills the panel width.
export function createSegmented<T extends string>(
  segments: Segment<T>[],
  value: T | null,
  onChange: (v: T) => void,
  variant: 'full' | 'compact' = 'full',
): Toggleable<T | null> {
  const el = document.createElement('div');
  el.className = variant === 'compact' ? 'bgw-seg bgw-seg-compact' : 'bgw-seg';
  el.setAttribute('role', 'group');
  const buttons: HTMLButtonElement[] = segments.map(seg => {
    const b = button('bgw-seg-btn', seg.text, seg.title, () => {
      set(seg.value);
      onChange(seg.value);
    });
    el.appendChild(b);
    return b;
  });
  function set(v: T | null) {
    segments.forEach((seg, i) => {
      const on = seg.value === v;
      buttons[i].classList.toggle('is-on', on);
      buttons[i].setAttribute('aria-pressed', String(on));
    });
  }
  set(value);
  return { el, set };
}

/// Venn glyphs, viewBox 0 0 30 18: two r=7.5 circles at cx 11 / 19, cy 9.
/// `vs` is the solid part (accent when the cell is selected), `vo` the thin
/// outline. The same family ships as SVG icons in MatterCAD's toolbar
/// (StaticData/Icons/{combine,subtract,intersect,subtract_and_replace}.svg),
/// which is why Difference draws the kept crescent directly instead of
/// knocking the right circle out with a background-colored disc: a painted
/// knockout is only right while the cell background is exactly that color
/// (it was already wrong on the selected cell's tint, and there is no one
/// right color at all in an app with a dark theme).
const VENN: Record<string, string> = {
  union:
    '<circle class="vs" cx="11" cy="9" r="7.5"/><circle class="vs" cx="19" cy="9" r="7.5"/>',
  intersection:
    '<circle class="vo" cx="11" cy="9" r="7.5"/><circle class="vo" cx="19" cy="9" r="7.5"/>' +
    '<path class="vs" d="M15 2.6557A7.5 7.5 0 0 1 15 15.3443A7.5 7.5 0 0 1 15 2.6557Z"/>',
  difference:
    '<path class="vs" d="M15 2.6557A7.5 7.5 0 1 0 15 15.3443A7.5 7.5 0 0 1 15 2.6557Z"/>' +
    '<circle class="vo" cx="11" cy="9" r="7.5"/><circle class="vo" cx="19" cy="9" r="7.5"/>',
};

export interface OpChoice {
  /** Index handed back to the caller (matches the wasm op enum). */
  value: number;
  /** Key into the glyph table. */
  glyph: keyof typeof VENN | string;
  text: string;
}

/// 3-column grid of Venn-glyph operation buttons (replaces the old dropdown).
export function createOpGrid(
  choices: OpChoice[],
  value: number,
  onChange: (v: number) => void,
): Toggleable<number> {
  const el = document.createElement('div');
  el.className = 'bgw-ops';
  const cells: HTMLButtonElement[] = choices.map(choice => {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'bgw-op';
    b.title = choice.text;
    const svg = `<svg class="bgw-op-glyph" viewBox="0 0 30 18" width="32" height="19" ` +
      `aria-hidden="true">${VENN[choice.glyph] ?? ''}</svg>`;
    b.innerHTML = `${svg}<span class="bgw-op-label"></span>`;
    (b.querySelector('.bgw-op-label') as HTMLElement).textContent = choice.text;
    b.addEventListener('click', () => { set(choice.value); onChange(choice.value); });
    el.appendChild(b);
    return b;
  });
  function set(v: number) {
    choices.forEach((choice, i) => {
      const on = choice.value === v;
      cells[i].classList.toggle('is-on', on);
      cells[i].setAttribute('aria-pressed', String(on));
    });
  }
  set(value);
  return { el, set };
}

/// One rounded on/off pill (View band).
export function createPill(
  label: string, title: string, on: boolean, onChange: (v: boolean) => void,
): Toggleable<boolean> {
  let state = on;
  const el = button('bgw-pill', label, title, () => { set(!state); onChange(state); });
  function set(v: boolean) {
    state = v;
    el.classList.toggle('is-on', v);
    el.setAttribute('aria-pressed', String(v));
  }
  set(on);
  return { el, set };
}

/// Small bordered text button (the Result band's debug actions).
export function createChip(label: string, title: string, onClick: () => void): HTMLButtonElement {
  return button('bgw-chip', label, title, onClick);
}

/// Square icon button (Load Random Pair's ⟳).
export function createIconButton(glyph: string, title: string, onClick: () => void): HTMLButtonElement {
  const b = button('bgw-iconbtn', glyph, title, onClick);
  b.setAttribute('aria-label', title);
  return b;
}

export interface CheckRow {
  el: HTMLElement;
  input: HTMLInputElement;
  /** Grey out (and block) the row when the choice cannot apply. */
  setEnabled(enabled: boolean): void;
  set(checked: boolean): void;
}

export function createCheckRow(
  label: string, title: string, checked: boolean, onChange: (v: boolean) => void,
): CheckRow {
  const el = document.createElement('label');
  el.className = 'bgw-check';
  el.title = title;
  const input = document.createElement('input');
  input.type = 'checkbox';
  input.checked = checked;
  input.addEventListener('change', () => onChange(input.checked));
  const span = document.createElement('span');
  span.textContent = label;
  el.appendChild(input);
  el.appendChild(span);
  return {
    el, input,
    setEnabled(enabled: boolean) {
      input.disabled = !enabled;
      el.classList.toggle('is-disabled', !enabled);
    },
    set(v: boolean) { input.checked = v; },
  };
}

/// Bare select with the panel's caret styling; no label of its own (the row
/// it sits in explains it).
export function createSelect(
  options: { value: string; text: string }[],
  value: string,
  onChange: (v: string) => void,
  title?: string,
): { el: HTMLElement; select: HTMLSelectElement } {
  const el = document.createElement('div');
  el.className = 'bgw-select';
  const select = document.createElement('select');
  for (const opt of options) {
    const o = document.createElement('option');
    o.value = opt.value;
    o.textContent = opt.text;
    select.appendChild(o);
  }
  select.value = value;
  if (title) select.title = title;
  select.addEventListener('change', () => onChange(select.value));
  el.appendChild(select);
  return { el, select };
}

export interface Disclosure {
  el: HTMLElement;
  body: HTMLElement;
  /** Live summary printed on the closed row. */
  setSummary(text: string, muted: boolean): void;
  setOpen(open: boolean): void;
  isOpen(): boolean;
}

/// Collapsed-by-default row that keeps showing its values while closed, so
/// folding the offsets away costs no information.
export function createDisclosure(
  title: string, open: boolean, onToggle: (open: boolean) => void,
): Disclosure {
  const el = document.createElement('div');
  el.className = 'bgw-disc';
  const head = document.createElement('button');
  head.type = 'button';
  head.className = 'bgw-disc-head';
  const left = document.createElement('span');
  left.className = 'bgw-disc-left';
  const caret = document.createElement('span');
  caret.className = 'bgw-caret';
  caret.textContent = '▸';
  const titleEl = document.createElement('span');
  titleEl.className = 'bgw-disc-title';
  titleEl.textContent = title;
  left.appendChild(caret);
  left.appendChild(titleEl);
  const summary = document.createElement('span');
  summary.className = 'bgw-disc-sum bgw-mono';
  head.appendChild(left);
  head.appendChild(summary);
  const body = document.createElement('div');
  body.className = 'bgw-disc-body';
  el.appendChild(head);
  el.appendChild(body);

  let state = open;
  function setOpen(v: boolean) {
    state = v;
    el.classList.toggle('is-open', v);
    head.setAttribute('aria-expanded', String(v));
    body.style.display = v ? '' : 'none';
  }
  head.addEventListener('click', () => { setOpen(!state); onToggle(state); });
  setOpen(open);
  return {
    el, body,
    setSummary(text: string, muted: boolean) {
      summary.textContent = text;
      summary.classList.toggle('is-muted', muted);
    },
    setOpen,
    isOpen: () => state,
  };
}

export interface OperandLine {
  /** Swatch colour key: which viewport mesh this row describes. */
  slot: 'a' | 'b';
  title: string;
  /** `#849728 · 702 tris · manifold` */
  meta: string;
  /** Colour the meta line as a caution (non-manifold input). */
  caution: boolean;
}

export interface OperandCard {
  el: HTMLElement;
  /** Replace the mesh rows and the rotation line. */
  setLines(lines: OperandLine[], rotation: string | null): void;
  /** Transient status ("Fetching #123…", a load failure) under the rows. */
  setMessage(text: string | null): void;
  setVisible(visible: boolean): void;
}

/// Inset card describing the two operands — plain information, not an alert.
export function createOperandCard(): OperandCard {
  const el = document.createElement('div');
  el.className = 'bgw-operands';
  const rows = document.createElement('div');
  rows.className = 'bgw-op-rows';
  const rot = document.createElement('div');
  rot.className = 'bgw-op-rot bgw-mono';
  const msg = document.createElement('div');
  msg.className = 'bgw-op-msg';
  msg.style.display = 'none';
  el.appendChild(rows);
  el.appendChild(rot);
  el.appendChild(msg);
  return {
    el,
    setLines(lines: OperandLine[], rotation: string | null) {
      rows.replaceChildren();
      for (const line of lines) {
        const row = document.createElement('div');
        row.className = 'bgw-op-row';
        const swatch = document.createElement('span');
        swatch.className = `bgw-swatch bgw-swatch-${line.slot}`;
        const text = document.createElement('div');
        text.className = 'bgw-op-text';
        const t = document.createElement('div');
        t.className = 'bgw-op-title';
        t.textContent = line.title;
        t.title = line.title;
        const m = document.createElement('div');
        m.className = line.caution ? 'bgw-op-meta bgw-mono is-caution' : 'bgw-op-meta bgw-mono';
        m.textContent = line.meta;
        text.appendChild(t);
        text.appendChild(m);
        row.appendChild(swatch);
        row.appendChild(text);
        rows.appendChild(row);
      }
      rot.textContent = rotation ?? '';
      rot.style.display = rotation ? '' : 'none';
    },
    setMessage(text: string | null) {
      msg.textContent = text ?? '';
      msg.style.display = text ? '' : 'none';
    },
    setVisible(visible: boolean) { el.style.display = visible ? '' : 'none'; },
  };
}
