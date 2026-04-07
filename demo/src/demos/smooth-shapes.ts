// Smooth/Refined shapes demo: subdivision refinement of base shapes

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { refinedShapeMesh, type MeshData } from '../wasm.ts';
import { loadSetting, saveSetting } from '../settings.ts';

const DEMO = 'smooth-shapes';

const SHAPES = [
  { value: '0', text: 'Tetrahedron' },
  { value: '1', text: 'Cube' },
  { value: '2', text: 'Cylinder' },
];

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Subdivision Refinement</h2>
        <p>Increase the refinement level to subdivide base shapes, adding more triangles for smoother geometry.</p>
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

  let shape = loadSetting(DEMO, 'shape', 0);
  let refineLevel = loadSetting(DEMO, 'refineLevel', 1);
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
    const data = refinedShapeMesh(shape, refineLevel);
    viewer.setMesh(data);
    viewer.setColor(0x6688aa);
    showReadout(data);
  }

  controlsEl.appendChild(createDropdown('Base Shape', SHAPES, String(shape), v => { saveSetting(DEMO, 'shape', v); shape = parseInt(v); update(); }));
  controlsEl.appendChild(createSlider('Refine Level ', 1, 5, refineLevel, 1, v => { saveSetting(DEMO, 'refineLevel', v); refineLevel = v; update(); }));
  controlsEl.appendChild(createCheckbox('Wireframe', loadSetting(DEMO, 'wireframe', false), v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  controlsEl.appendChild(readout);

  update();
  return () => viewer.dispose();
}
