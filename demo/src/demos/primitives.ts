// Primitives demo: cube, sphere, cylinder, tetrahedron with parameter sliders

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { cubeMesh, sphereMesh, cylinderMesh, tetrahedronMesh, type MeshData } from '../wasm.ts';

type Shape = 'cube' | 'sphere' | 'cylinder' | 'tetrahedron';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Primitives</h2>
        <p>Interactive 3D mesh primitives with adjustable parameters.</p>
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

  let shape: Shape = 'cube';
  let cubeSize = { x: 1, y: 1, z: 1 };
  let sphereRadius = 1;
  let sphereSegments = 16;
  let cylHeight = 2;
  let cylRadiusLow = 1;
  let cylRadiusHigh = 1;
  let cylSegments = 32;

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
    switch (shape) {
      case 'cube':
        data = cubeMesh(cubeSize.x, cubeSize.y, cubeSize.z, true);
        break;
      case 'sphere':
        data = sphereMesh(sphereRadius, sphereSegments);
        break;
      case 'cylinder':
        data = cylinderMesh(cylHeight, cylRadiusLow, cylRadiusHigh, cylSegments);
        break;
      case 'tetrahedron':
        data = tetrahedronMesh();
        break;
    }
    viewer.setMesh(data);
    showReadout(data);
  }

  function buildControls() {
    controlsEl.innerHTML = '';

    controlsEl.appendChild(createDropdown('Shape', [
      { value: 'cube', text: 'Cube' },
      { value: 'sphere', text: 'Sphere' },
      { value: 'cylinder', text: 'Cylinder' },
      { value: 'tetrahedron', text: 'Tetrahedron' },
    ], shape, (v) => { shape = v as Shape; buildControls(); update(); }));

    switch (shape) {
      case 'cube':
        controlsEl.appendChild(createSlider('Size X ', 0.1, 3, cubeSize.x, 0.1, v => { cubeSize.x = v; update(); }));
        controlsEl.appendChild(createSlider('Size Y ', 0.1, 3, cubeSize.y, 0.1, v => { cubeSize.y = v; update(); }));
        controlsEl.appendChild(createSlider('Size Z ', 0.1, 3, cubeSize.z, 0.1, v => { cubeSize.z = v; update(); }));
        break;
      case 'sphere':
        controlsEl.appendChild(createSlider('Radius ', 0.1, 3, sphereRadius, 0.1, v => { sphereRadius = v; update(); }));
        controlsEl.appendChild(createSlider('Segments ', 4, 64, sphereSegments, 4, v => { sphereSegments = v; update(); }));
        break;
      case 'cylinder':
        controlsEl.appendChild(createSlider('Height ', 0.1, 4, cylHeight, 0.1, v => { cylHeight = v; update(); }));
        controlsEl.appendChild(createSlider('Radius Low ', 0.1, 3, cylRadiusLow, 0.1, v => { cylRadiusLow = v; update(); }));
        controlsEl.appendChild(createSlider('Radius High ', 0.1, 3, cylRadiusHigh, 0.1, v => { cylRadiusHigh = v; update(); }));
        controlsEl.appendChild(createSlider('Segments ', 4, 64, cylSegments, 4, v => { cylSegments = v; update(); }));
        break;
      case 'tetrahedron':
        // No parameters
        break;
    }

    controlsEl.appendChild(createCheckbox('Wireframe', false, (v) => viewer.setWireframe(v)));
    controlsEl.appendChild(readout);
  }

  buildControls();
  update();

  return () => viewer.dispose();
}
