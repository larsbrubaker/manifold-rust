// Reusable UI control builders

export function createSlider(
  label: string, min: number, max: number, value: number, step: number,
  onChange: (val: number) => void,
): HTMLElement {
  const group = document.createElement('div');
  group.className = 'control-group';
  const lbl = document.createElement('label');
  const labelText = document.createTextNode(label);
  const valSpan = document.createElement('span');
  valSpan.className = 'slider-value';
  valSpan.textContent = String(value);
  lbl.appendChild(labelText);
  lbl.appendChild(valSpan);

  const input = document.createElement('input');
  input.type = 'range';
  input.min = String(min);
  input.max = String(max);
  input.step = String(step);
  input.value = String(value);
  input.addEventListener('input', () => {
    const v = parseFloat(input.value);
    valSpan.textContent = String(v);
    onChange(v);
  });

  group.appendChild(lbl);
  group.appendChild(input);
  return group;
}

export function createDropdown(
  label: string,
  options: { value: string; text: string }[],
  selectedValue: string,
  onChange: (val: string) => void,
): HTMLElement {
  const group = document.createElement('div');
  group.className = 'control-group';
  const lbl = document.createElement('label');
  lbl.textContent = label;

  const select = document.createElement('select');
  for (const opt of options) {
    const o = document.createElement('option');
    o.value = opt.value;
    o.textContent = opt.text;
    select.appendChild(o);
  }
  select.value = selectedValue;
  select.addEventListener('change', () => onChange(select.value));

  group.appendChild(lbl);
  group.appendChild(select);
  return group;
}

export function createCheckbox(
  label: string, checked: boolean,
  onChange: (val: boolean) => void,
): HTMLElement {
  const wrapper = document.createElement('label');
  wrapper.className = 'control-checkbox';
  const input = document.createElement('input');
  input.type = 'checkbox';
  input.checked = checked;
  input.addEventListener('change', () => onChange(input.checked));
  const span = document.createElement('span');
  span.textContent = label;
  wrapper.appendChild(input);
  wrapper.appendChild(span);
  return wrapper;
}

export function createButton(label: string, onClick: () => void): HTMLButtonElement {
  const btn = document.createElement('button');
  btn.className = 'control-btn';
  btn.textContent = label;
  btn.addEventListener('click', onClick);
  return btn;
}

export function createReadout(): HTMLElement {
  const box = document.createElement('div');
  box.className = 'info-readout';
  return box;
}

export function updateReadout(el: HTMLElement, entries: { label: string; value: string }[]) {
  el.innerHTML = entries.map(e =>
    `<span class="label">${e.label}:</span> <span class="value">${e.value}</span>`
  ).join('<br>');
}
