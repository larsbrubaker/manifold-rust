// Properties demo: display mesh metrics for various shapes

import { ThreeViewer } from '../three-viewer.ts';
import { createDropdown, createSlider, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { cubeMesh, sphereMesh, cylinderMesh, tetrahedronMesh, extrudeMesh, type MeshData } from '../wasm.ts';

type ShapeChoice = 'cube' | 'sphere' | 'cylinder' | 'tetrahedron' | 'extruded-circle';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Mesh Properties</h2>
        <p>Inspect volume, surface area, vertex count, and triangle count for various meshes.</p>
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

  let shape: ShapeChoice = 'sphere';
  let segments = 32;

  const readout = createReadout();

  function showReadout(data: MeshData) {
    updateReadout(readout, [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(6) },
      { label: 'Surface Area', value: data.surface_area.toFixed(6) },
    ]);
  }

  function update() {
    let data: MeshData;
    switch (shape) {
      case 'cube':
        data = cubeMesh(1, 1, 1, true);
        break;
      case 'sphere':
        data = sphereMesh(1, segments);
        break;
      case 'cylinder':
        data = cylinderMesh(2, 1, 1, segments);
        break;
      case 'tetrahedron':
        data = tetrahedronMesh();
        break;
      case 'extruded-circle':
        data = extrudeMesh(1, segments, 2);
        break;
    }
    viewer.setMesh(data);
    viewer.setColor(0x44aa88);
    showReadout(data);
  }

  function buildControls() {
    controlsEl.innerHTML = '';

    controlsEl.appendChild(createDropdown('Shape', [
      { value: 'cube', text: 'Unit Cube' },
      { value: 'sphere', text: 'Unit Sphere' },
      { value: 'cylinder', text: 'Cylinder (h=2, r=1)' },
      { value: 'tetrahedron', text: 'Tetrahedron' },
      { value: 'extruded-circle', text: 'Extruded Circle' },
    ], shape, (v) => { shape = v as ShapeChoice; buildControls(); update(); }));

    if (shape !== 'cube' && shape !== 'tetrahedron') {
      controlsEl.appendChild(createSlider('Segments ', 4, 128, segments, 4, v => { segments = v; update(); }));
    }

    controlsEl.appendChild(createCheckbox('Wireframe', false, (v) => viewer.setWireframe(v)));
    controlsEl.appendChild(readout);
  }

  buildControls();
  update();

  return () => viewer.dispose();
}
