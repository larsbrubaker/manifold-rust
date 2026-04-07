// Boolean Gallery: sphere-sphere, sphere-cube, cylinder-sphere combinations

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { booleanGalleryMesh, booleanGalleryMeshRotated, type MeshData } from '../wasm.ts';
import { loadSetting, saveSetting } from '../settings.ts';

const DEMO = 'boolean-gallery';

const SHAPES = [
  { value: '0', text: 'Cube' },
  { value: '1', text: 'Sphere' },
  { value: '2', text: 'Cylinder' },
  { value: '3', text: 'Spiky Dodecahedron' },
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
        <p>Combine cubes, spheres, and cylinders with union, intersection, and difference. Toggle Animate to see Shape B rotate continuously.</p>
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

  let shapeA = loadSetting(DEMO, 'shapeA', 1);
  let shapeB = loadSetting(DEMO, 'shapeB', 0);
  let op = loadSetting(DEMO, 'op', 2);
  let offsetX = loadSetting(DEMO, 'offsetX', 0.3);
  let offsetY = loadSetting(DEMO, 'offsetY', 0.0);
  let offsetZ = loadSetting(DEMO, 'offsetZ', 0.0);
  let wireframe = loadSetting(DEMO, 'wireframe', false);
  let animate = loadSetting(DEMO, 'animate', true);
  let animating = false;
  let animId = 0;
  let rotX = 0, rotY = 0, rotZ = 0;
  const ROT_SPEED_X = 0.7;  // degrees per frame
  const ROT_SPEED_Y = 1.5;
  const ROT_SPEED_Z = 0.3;

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

  function update(silent = false) {
    try {
      let data: MeshData;
      if (rotX !== 0 || rotY !== 0 || rotZ !== 0) {
        data = booleanGalleryMeshRotated(shapeA, shapeB, op, offsetX, offsetY, offsetZ, rotX, rotY, rotZ);
      } else {
        data = booleanGalleryMesh(shapeA, shapeB, op, offsetX, offsetY, offsetZ);
      }
      viewer.setMesh(data);
      viewer.setColor(OP_COLORS[op] || 0x4488cc);
      errorBox.style.display = 'none';
      showReadout(data);
    } catch {
      if (!silent) {
        errorBox.style.display = 'block';
        errorBox.innerHTML = '<strong>Boolean operation failed.</strong><br>Try adjusting the offset.';
        updateReadout(readout, []);
      }
    }
  }

  function animateStep() {
    if (!animating) return;
    rotX = (rotX + ROT_SPEED_X) % 360;
    rotY = (rotY + ROT_SPEED_Y) % 360;
    rotZ = (rotZ + ROT_SPEED_Z) % 360;
    update(true);
    animId = requestAnimationFrame(animateStep);
  }

  function toggleAnimate(on: boolean) {
    saveSetting(DEMO, 'animate', on);
    if (on) {
      animating = true;
      animateStep();
    } else {
      animating = false;
      cancelAnimationFrame(animId);
      // Stay at current rotation — don't reset
    }
  }

  controlsEl.appendChild(createDropdown('Shape A', SHAPES, String(shapeA), v => { shapeA = parseInt(v); saveSetting(DEMO, 'shapeA', shapeA); update(); }));
  controlsEl.appendChild(createDropdown('Shape B', SHAPES, String(shapeB), v => { shapeB = parseInt(v); saveSetting(DEMO, 'shapeB', shapeB); update(); }));
  controlsEl.appendChild(createDropdown('Operation', OPS, String(op), v => { op = parseInt(v); saveSetting(DEMO, 'op', op); update(); }));
  controlsEl.appendChild(createSlider('Offset X ', -1.5, 1.5, offsetX, 0.1, v => { offsetX = v; saveSetting(DEMO, 'offsetX', v); update(); }));
  controlsEl.appendChild(createSlider('Offset Y ', -1.5, 1.5, offsetY, 0.1, v => { offsetY = v; saveSetting(DEMO, 'offsetY', v); update(); }));
  controlsEl.appendChild(createSlider('Offset Z ', -1.5, 1.5, offsetZ, 0.1, v => { offsetZ = v; saveSetting(DEMO, 'offsetZ', v); update(); }));
  controlsEl.appendChild(createCheckbox('Animate', animate, toggleAnimate));
  controlsEl.appendChild(createCheckbox('Wireframe', wireframe, v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  update();
  if (wireframe) viewer.setWireframe(true);
  toggleAnimate(animate);

  return () => {
    animating = false;
    cancelAnimationFrame(animId);
    viewer.dispose();
  };
}
