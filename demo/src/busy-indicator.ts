// Lightweight "computing…" overlay shown while a long-running operation is in
// flight (today: boolean evaluations dispatched by boolean-runner.ts).
//
// The wasm boolean cannot report progress yet, so the card shows a CSS spinner,
// a title, a ticking elapsed-seconds counter and an optional cancel button.
// The layout deliberately reserves the two rows a real progress API will need —
// a phase label and a progress track (indeterminate until a fraction arrives) —
// so `setPhase()` can fill them in later without the card resizing or moving.
//
// Dependency-free: styles live in styles/main.css (.busy-indicator*).

export class BusyIndicator {
  private root: HTMLElement | null = null;
  private titleEl!: HTMLElement;
  private phaseEl!: HTMLElement;
  private elapsedEl!: HTMLElement;
  private barEl!: HTMLElement;
  private trackEl!: HTMLElement;
  private cancelBtn!: HTMLButtonElement;

  private tickTimer = 0;
  private showTimer = 0;
  private startMs = 0;
  private visible = false;
  private onCancel: (() => void) | null = null;

  /** Delay before the card appears, so sub-second runs never flash it. */
  constructor(private readonly showDelayMs = 250) {}

  private ensureDom(): HTMLElement {
    if (this.root) return this.root;

    const root = document.createElement('div');
    root.className = 'busy-indicator';
    root.setAttribute('role', 'status');
    root.setAttribute('aria-live', 'polite');

    const head = document.createElement('div');
    head.className = 'busy-head';

    const spinner = document.createElement('div');
    spinner.className = 'busy-spinner';
    head.appendChild(spinner);

    const text = document.createElement('div');
    text.className = 'busy-text';
    this.titleEl = document.createElement('div');
    this.titleEl.className = 'busy-title';
    this.phaseEl = document.createElement('div');
    this.phaseEl.className = 'busy-phase';
    text.appendChild(this.titleEl);
    text.appendChild(this.phaseEl);
    head.appendChild(text);

    this.elapsedEl = document.createElement('div');
    this.elapsedEl.className = 'busy-elapsed';
    head.appendChild(this.elapsedEl);

    this.trackEl = document.createElement('div');
    this.trackEl.className = 'busy-track indeterminate';
    this.barEl = document.createElement('div');
    this.barEl.className = 'busy-bar';
    this.trackEl.appendChild(this.barEl);

    this.cancelBtn = document.createElement('button');
    this.cancelBtn.className = 'busy-cancel';
    this.cancelBtn.type = 'button';
    this.cancelBtn.textContent = 'Cancel';
    this.cancelBtn.addEventListener('click', () => this.onCancel?.());

    root.appendChild(head);
    root.appendChild(this.trackEl);
    root.appendChild(this.cancelBtn);
    document.body.appendChild(root);
    this.root = root;
    return root;
  }

  /** Cancel affordance: pass null (or omit) to hide the button entirely. */
  setCancel(handler: (() => void) | null) {
    this.onCancel = handler;
    if (this.root) this.cancelBtn.style.display = handler ? '' : 'none';
  }

  /**
   * Start (or retarget) the indicator. Calling show() again while visible
   * updates the title and keeps the elapsed clock running, which is what a
   * coalesced back-to-back run should look like.
   */
  show(title: string, phase = 'Working…') {
    const root = this.ensureDom();
    this.titleEl.textContent = title;
    this.phaseEl.textContent = phase;
    this.cancelBtn.style.display = this.onCancel ? '' : 'none';
    if (this.visible) return;

    this.visible = true;
    this.startMs = performance.now();
    this.setProgressFraction(null);
    this.renderElapsed();
    window.clearTimeout(this.showTimer);
    this.showTimer = window.setTimeout(() => root.classList.add('visible'), this.showDelayMs);
    window.clearInterval(this.tickTimer);
    this.tickTimer = window.setInterval(() => this.renderElapsed(), 100);
  }

  /**
   * Forward-looking hook for the per-phase progress API: updates the phase
   * label and the track. `fraction` in [0,1] switches the bar from its
   * indeterminate sweep to a determinate fill; null keeps it indeterminate.
   */
  setPhase(label: string, fraction: number | null = null) {
    if (!this.root) return;
    this.phaseEl.textContent = label;
    this.setProgressFraction(fraction);
  }

  private setProgressFraction(fraction: number | null) {
    if (!this.root) return;
    if (fraction === null || !Number.isFinite(fraction)) {
      this.trackEl.classList.add('indeterminate');
      this.barEl.style.width = '';
    } else {
      const clamped = Math.max(0, Math.min(1, fraction));
      this.trackEl.classList.remove('indeterminate');
      this.barEl.style.width = `${(clamped * 100).toFixed(1)}%`;
    }
  }

  private renderElapsed() {
    if (!this.root) return;
    const secs = (performance.now() - this.startMs) / 1000;
    this.elapsedEl.textContent = `${secs.toFixed(1)} s`;
  }

  hide() {
    window.clearTimeout(this.showTimer);
    window.clearInterval(this.tickTimer);
    this.showTimer = 0;
    this.tickTimer = 0;
    this.visible = false;
    this.root?.classList.remove('visible');
  }

  dispose() {
    this.hide();
    this.root?.remove();
    this.root = null;
  }
}
