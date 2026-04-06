// Extrude with twist and taper demo

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { extrudeTwistMesh, type MeshData } from '../wasm.ts';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Extrude with Twist</h2>
        <p>Extrude a circular cross-section with twist rotation and top-scaling (taper). Adjust divisions for smoother twists.</p>
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

  let radius = 0.5;
  let segments = 32;
  let height = 2.0;
  let twistDegrees = 180;
  let nDivisions = 20;
  let scaleTop = 0.5;

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
    const data = extrudeTwistMesh(radius, segments, height, twistDegrees, nDivisions, scaleTop);
    viewer.setMesh(data);
    viewer.setColor(0x88aa44);
    showReadout(data);
  }

  controlsEl.appendChild(createSlider('Radius ', 0.1, 1.5, radius, 0.1, v => { radius = v; update(); }));
  controlsEl.appendChild(createSlider('Segments ', 4, 64, segments, 4, v => { segments = v; update(); }));
  controlsEl.appendChild(createSlider('Height ', 0.5, 5.0, height, 0.5, v => { height = v; update(); }));
  controlsEl.appendChild(createSlider('Twist (deg) ', 0, 720, twistDegrees, 10, v => { twistDegrees = v; update(); }));
  controlsEl.appendChild(createSlider('Divisions ', 1, 50, nDivisions, 1, v => { nDivisions = v; update(); }));
  controlsEl.appendChild(createSlider('Scale Top ', 0.1, 2.0, scaleTop, 0.1, v => { scaleTop = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', false, v => viewer.setWireframe(v)));
  controlsEl.appendChild(readout);

  update();
  return () => viewer.dispose();
}
