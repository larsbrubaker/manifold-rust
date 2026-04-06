// Boolean Gallery: sphere-sphere, sphere-cube, cylinder-sphere combinations

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { booleanGalleryMesh, type MeshData } from '../wasm.ts';

const SHAPES = [
  { value: '0', text: 'Cube' },
  { value: '1', text: 'Sphere' },
  { value: '2', text: 'Cylinder' },
];

const OPS = [
  { value: '0', text: 'Union' },
  { value: '1', text: 'Intersection' },
  { value: '2', text: 'Difference' },
];

const OP_COLORS = [0x4488cc, 0x44aa44, 0xcc4444];

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Boolean Gallery</h2>
        <p>Combine cubes, spheres, and cylinders with union, intersection, and difference. Adjust 3D offset to explore how shapes interact.</p>
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

  let shapeA = 1; // sphere
  let shapeB = 0; // cube
  let op = 2;     // difference
  let offsetX = 0.3;
  let offsetY = 0.0;
  let offsetZ = 0.0;

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

  function update() {
    try {
      const data = booleanGalleryMesh(shapeA, shapeB, op, offsetX, offsetY, offsetZ);
      viewer.setMesh(data);
      viewer.setColor(OP_COLORS[op] || 0x4488cc);
      showReadout(data);
    } catch {
      errorBox.style.display = 'block';
      errorBox.innerHTML = '<strong>Boolean operation failed.</strong><br>Try adjusting the offset.';
      updateReadout(readout, []);
    }
  }

  controlsEl.appendChild(createDropdown('Shape A', SHAPES, String(shapeA), v => { shapeA = parseInt(v); update(); }));
  controlsEl.appendChild(createDropdown('Shape B', SHAPES, String(shapeB), v => { shapeB = parseInt(v); update(); }));
  controlsEl.appendChild(createDropdown('Operation', OPS, String(op), v => { op = parseInt(v); update(); }));
  controlsEl.appendChild(createSlider('Offset X ', -1.5, 1.5, offsetX, 0.1, v => { offsetX = v; update(); }));
  controlsEl.appendChild(createSlider('Offset Y ', -1.5, 1.5, offsetY, 0.1, v => { offsetY = v; update(); }));
  controlsEl.appendChild(createSlider('Offset Z ', -1.5, 1.5, offsetZ, 0.1, v => { offsetZ = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', false, v => viewer.setWireframe(v)));
  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  update();
  return () => viewer.dispose();
}
