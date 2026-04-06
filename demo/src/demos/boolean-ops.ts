// Boolean Operations demo

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { unionMesh, intersectMesh, differenceMesh, type MeshData } from '../wasm.ts';

type BoolOp = 'union' | 'intersection' | 'difference';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Boolean Operations</h2>
        <p>Union, intersection, and difference of two unit cubes with adjustable X offset.</p>
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

  let op: BoolOp = 'union';
  let offset = 0.5; // overlapping by default to showcase boolean operations

  const readout = createReadout();
  const errorBox = document.createElement('div');
  errorBox.className = 'demo-note';
  errorBox.style.display = 'none';

  function showReadout(data: MeshData) {
    errorBox.style.display = 'none';
    updateReadout(readout, [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ]);
  }

  function showError() {
    errorBox.style.display = 'block';
    errorBox.innerHTML = '<strong>Boolean operation failed.</strong><br>Try adjusting the offset.';
    updateReadout(readout, []);
  }

  function update() {
    try {
      let data: MeshData;
      switch (op) {
        case 'union':
          data = unionMesh(offset);
          break;
        case 'intersection':
          data = intersectMesh(offset);
          break;
        case 'difference':
          data = differenceMesh(offset);
          break;
      }
      viewer.setMesh(data);
      viewer.setColor(op === 'union' ? 0x4488cc : op === 'intersection' ? 0x44aa44 : 0xcc4444);
      showReadout(data);
    } catch {
      showError();
    }
  }

  controlsEl.appendChild(createDropdown('Operation', [
    { value: 'union', text: 'Union' },
    { value: 'intersection', text: 'Intersection' },
    { value: 'difference', text: 'Difference' },
  ], op, (v) => { op = v as BoolOp; update(); }));

  controlsEl.appendChild(createSlider('Offset X ', -0.5, 3, offset, 0.1, v => { offset = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', false, (v) => viewer.setWireframe(v)));
  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  update();

  return () => viewer.dispose();
}
