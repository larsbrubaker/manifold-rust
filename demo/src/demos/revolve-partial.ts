// Partial revolve demo: revolve a 2D profile less than 360 degrees

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { revolvePartialMesh, type MeshData } from '../wasm.ts';

const PROFILES = [
  { value: '0', text: 'Circle' },
  { value: '1', text: 'Rectangle' },
];

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Partial Revolve</h2>
        <p>Revolve a 2D profile (offset from the Y axis) by a partial angle to create torus-like shapes, arcs, and rings.</p>
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

  let profile = 0;
  let segments = 32;
  let degrees = 270;

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
    const data = revolvePartialMesh(profile, segments, degrees);
    viewer.setMesh(data);
    viewer.setColor(0xaa6688);
    showReadout(data);
  }

  controlsEl.appendChild(createDropdown('Profile', PROFILES, String(profile), v => { profile = parseInt(v); update(); }));
  controlsEl.appendChild(createSlider('Segments ', 4, 64, segments, 4, v => { segments = v; update(); }));
  controlsEl.appendChild(createSlider('Degrees ', 10, 360, degrees, 10, v => { degrees = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', false, v => viewer.setWireframe(v)));
  controlsEl.appendChild(readout);

  update();
  return () => viewer.dispose();
}
