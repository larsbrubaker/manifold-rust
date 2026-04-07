// Menger Sponge demo: recursive boolean subtraction fractal

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { mengerSpongeMesh, type MeshData } from '../wasm.ts';
import { loadSetting, saveSetting } from '../settings.ts';

const DEMO = 'menger-sponge';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Menger Sponge</h2>
        <p>Recursive boolean subtraction creates this classic fractal. Each depth level removes cross-shaped holes from every sub-cube.</p>
      </div>
      <div class="demo-layout">
        <div class="demo-canvas-area" id="viewer-container"></div>
        <div class="demo-controls" id="controls"></div>
      </div>
    </div>
  `;

  const viewerEl = document.getElementById('viewer-container')!;
  const controlsEl = document.getElementById('controls')!;
  const viewer = new ThreeViewer(viewerEl);

  let depth = loadSetting(DEMO, 'depth', 0);
  const readout = createReadout();

  function showReadout(data: MeshData) {
    updateReadout(readout, [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ]);
  }

  function update() {
    const data = mengerSpongeMesh(depth);
    viewer.setMesh(data);
    viewer.setColor(0xcc8844);
    showReadout(data);
  }

  controlsEl.appendChild(createSlider('Depth ', 0, 2, depth, 1, v => { saveSetting(DEMO, 'depth', v); depth = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', loadSetting(DEMO, 'wireframe', false), v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  controlsEl.appendChild(readout);

  update();
  return () => viewer.dispose();
}
