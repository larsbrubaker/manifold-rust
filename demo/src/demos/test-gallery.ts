// Test Gallery: WASM visualizations of key test cases

import { ThreeViewer } from '../three-viewer.ts';
import { createReadout, updateReadout } from '../controls.ts';
import { loadSetting, saveSetting } from '../settings.ts';
import {
  testMirrorUnionMesh, testSplitByPlaneMesh, testVugMesh, testWarpMesh,
  testSpiralMesh, testSphereDiffMesh, testCubesUnionMesh, testBatchSubtractMesh,
  type MeshData,
} from '../wasm.ts';

interface TestCase {
  name: string;
  description: string;
  color: number;
  generate: () => MeshData;
}

const TEST_CASES: TestCase[] = [
  {
    name: 'Mirror Union',
    description: 'Three cubes: one centered, one offset, one mirrored across the (1,1,0) plane, all unioned. Tests mirror() + union.',
    color: 0x4488cc,
    generate: () => testMirrorUnionMesh(),
  },
  {
    name: 'Split by Plane (Top)',
    description: 'A rotated cube split by a plane at z=1. Showing the top half. Tests split_by_plane().',
    color: 0x44aa88,
    generate: () => testSplitByPlaneMesh(0),
  },
  {
    name: 'Split by Plane (Bottom)',
    description: 'A rotated cube split by a plane at z=1. Showing the bottom half. Tests split_by_plane().',
    color: 0x88aa44,
    generate: () => testSplitByPlaneMesh(1),
  },
  {
    name: 'Vug (Cavity)',
    description: 'Outer 4x4x4 cube with inner 1x1x1 cube subtracted, then split to reveal the cavity inside. Tests difference() + split_by_plane().',
    color: 0xcc8844,
    generate: () => testVugMesh(),
  },
  {
    name: 'Warp (Parabolic)',
    description: 'Extruded square with warp function x += z\u00B2 applied. The shape bends parabolically. Tests warp().',
    color: 0xaa4488,
    generate: () => testWarpMesh(),
  },
  {
    name: 'Spiral',
    description: 'Recursive spiral of 11 unit cubes placed along an Archimedean spiral path, all unioned. Tests deep recursive boolean.',
    color: 0x6644cc,
    generate: () => testSpiralMesh(),
  },
  {
    name: 'Sphere Difference',
    description: 'Unit sphere with an offset sphere subtracted, leaving a spherical cap cutout. Tests sphere boolean precision.',
    color: 0xcc4444,
    generate: () => testSphereDiffMesh(),
  },
  {
    name: 'Overlapping Cubes',
    description: 'Three unit cubes at different offsets, all unioned. Tests multi-body union with partial overlap.',
    color: 0x44ccaa,
    generate: () => testCubesUnionMesh(),
  },
  {
    name: 'Batch Subtract',
    description: 'Flat slab with three cylinders punched through. Tests batch boolean subtraction.',
    color: 0x888888,
    generate: () => testBatchSubtractMesh(),
  },
];

const DEMO = 'test-gallery';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Test Gallery</h2>
        <p>WASM visualizations of key test operations: mirror, split, warp, spiral, boolean precision, and more.</p>
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

  let currentIndex = loadSetting(DEMO, 'currentIndex', 0);
  const readout = createReadout();
  const descEl = document.createElement('div');
  descEl.style.cssText = 'padding:8px 0;color:var(--text-muted);font-size:0.9em;line-height:1.4;';

  function showReadout(data: MeshData) {
    updateReadout(readout, [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ]);
  }

  function update() {
    const tc = TEST_CASES[currentIndex];
    descEl.textContent = tc.description;
    viewer.resetFrame();
    try {
      const data = tc.generate();
      viewer.setMesh(data);
      viewer.setColor(tc.color);
      showReadout(data);
    } catch (e) {
      updateReadout(readout, [{ label: 'Error', value: String(e) }]);
    }
  }

  // Test case selector
  const selectLabel = document.createElement('div');
  selectLabel.textContent = 'Test Case';
  selectLabel.style.cssText = 'font-weight:600;margin-bottom:4px;';
  const select = document.createElement('select');
  select.style.cssText = 'width:100%;padding:6px 8px;border-radius:6px;border:1px solid var(--border);background:var(--surface);color:var(--text);font-size:0.95em;';
  TEST_CASES.forEach((tc, i) => {
    const opt = document.createElement('option');
    opt.value = String(i);
    opt.textContent = tc.name;
    select.appendChild(opt);
  });
  select.value = String(currentIndex);
  select.addEventListener('change', () => {
    currentIndex = parseInt(select.value);
    saveSetting(DEMO, 'currentIndex', currentIndex);
    update();
  });

  // Prev/Next buttons
  const navRow = document.createElement('div');
  navRow.style.cssText = 'display:flex;gap:8px;margin-top:8px;';
  const prevBtn = document.createElement('button');
  prevBtn.textContent = '\u25C0 Prev';
  prevBtn.style.cssText = 'flex:1;padding:8px;border-radius:6px;border:1px solid var(--border);background:var(--surface);color:var(--text);cursor:pointer;font-size:0.9em;';
  prevBtn.addEventListener('click', () => {
    currentIndex = (currentIndex - 1 + TEST_CASES.length) % TEST_CASES.length;
    saveSetting(DEMO, 'currentIndex', currentIndex);
    select.value = String(currentIndex);
    update();
  });
  const nextBtn = document.createElement('button');
  nextBtn.textContent = 'Next \u25B6';
  nextBtn.style.cssText = 'flex:1;padding:8px;border-radius:6px;border:1px solid var(--border);background:var(--surface);color:var(--text);cursor:pointer;font-size:0.9em;';
  nextBtn.addEventListener('click', () => {
    currentIndex = (currentIndex + 1) % TEST_CASES.length;
    saveSetting(DEMO, 'currentIndex', currentIndex);
    select.value = String(currentIndex);
    update();
  });
  navRow.appendChild(prevBtn);
  navRow.appendChild(nextBtn);

  // Wireframe toggle
  const wireLabel = document.createElement('label');
  wireLabel.style.cssText = 'display:flex;align-items:center;gap:8px;margin-top:12px;cursor:pointer;';
  const wireCb = document.createElement('input');
  wireCb.type = 'checkbox';
  wireCb.checked = loadSetting(DEMO, 'wireframe', false);
  wireCb.addEventListener('change', () => { saveSetting(DEMO, 'wireframe', wireCb.checked); viewer.setWireframe(wireCb.checked); });
  wireLabel.appendChild(wireCb);
  wireLabel.appendChild(document.createTextNode('Wireframe'));

  controlsEl.appendChild(selectLabel);
  controlsEl.appendChild(select);
  controlsEl.appendChild(navRow);
  controlsEl.appendChild(descEl);
  controlsEl.appendChild(wireLabel);
  controlsEl.appendChild(readout);

  if (wireCb.checked) viewer.setWireframe(true);
  update();
  return () => viewer.dispose();
}
