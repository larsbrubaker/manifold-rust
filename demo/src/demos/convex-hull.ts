// Convex Hull demo

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createButton, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { hullMesh, type MeshData } from '../wasm.ts';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Convex Hull</h2>
        <p>Generate convex hulls from random 3D point clouds using QuickHull.</p>
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

  let numPoints = 50;
  let points = new Float32Array(0);

  const readout = createReadout();

  function showReadout(data: MeshData) {
    updateReadout(readout, [
      { label: 'Input Points', value: String(points.length / 3) },
      { label: 'Hull Vertices', value: String(data.num_vert) },
      { label: 'Hull Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ]);
  }

  function generatePoints() {
    points = new Float32Array(numPoints * 3);
    for (let i = 0; i < points.length; i++) {
      points[i] = (Math.random() - 0.5) * 2;
    }
  }

  function update() {
    if (points.length < 12) return; // need at least 4 points
    const data = hullMesh(points);
    viewer.setMesh(data);
    viewer.setColor(0x8866cc);
    showReadout(data);
  }

  controlsEl.appendChild(createSlider('Point Count ', 10, 500, numPoints, 10, v => {
    numPoints = v;
    generatePoints();
    update();
  }));
  controlsEl.appendChild(createButton('Regenerate', () => { generatePoints(); update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', false, (v) => viewer.setWireframe(v)));
  controlsEl.appendChild(readout);

  generatePoints();
  update();

  return () => viewer.dispose();
}
