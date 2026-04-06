// Extrude & Revolve demo

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { extrudeMesh, revolveMesh, type MeshData } from '../wasm.ts';

type Mode = 'extrude' | 'revolve';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Extrude &amp; Revolve</h2>
        <p>Turn 2D cross-sections into 3D solids.</p>
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

  let mode: Mode = 'extrude';
  let radius = 0.5;
  let segments = 32;
  let extrudeHeight = 2;
  let revolveDegrees = 360;

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
    let data: MeshData;
    if (mode === 'extrude') {
      data = extrudeMesh(radius, segments, extrudeHeight);
    } else {
      data = revolveMesh(radius, segments, revolveDegrees);
    }
    viewer.setMesh(data);
    showReadout(data);
  }

  function buildControls() {
    controlsEl.innerHTML = '';

    controlsEl.appendChild(createDropdown('Mode', [
      { value: 'extrude', text: 'Extrude' },
      { value: 'revolve', text: 'Revolve' },
    ], mode, (v) => { mode = v as Mode; buildControls(); update(); }));

    controlsEl.appendChild(createSlider('Radius ', 0.1, 2, radius, 0.1, v => { radius = v; update(); }));
    controlsEl.appendChild(createSlider('Segments ', 4, 64, segments, 4, v => { segments = v; update(); }));

    if (mode === 'extrude') {
      controlsEl.appendChild(createSlider('Height ', 0.1, 5, extrudeHeight, 0.1, v => { extrudeHeight = v; update(); }));
    } else {
      controlsEl.appendChild(createSlider('Degrees ', 30, 360, revolveDegrees, 30, v => { revolveDegrees = v; update(); }));
    }

    controlsEl.appendChild(createCheckbox('Wireframe', false, (v) => viewer.setWireframe(v)));
    controlsEl.appendChild(readout);
  }

  buildControls();
  update();

  return () => viewer.dispose();
}
