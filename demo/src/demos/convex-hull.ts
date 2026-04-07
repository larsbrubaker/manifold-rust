// Convex Hull demo

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createButton, createCheckbox, createReadout, updateReadout } from '../controls.ts';
import { hullMesh, type MeshData } from '../wasm.ts';
import { loadSetting, saveSetting } from '../settings.ts';

const DEMO = 'convex-hull';

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Convex Hull</h2>
        <p>Generate convex hulls from random 3D point clouds using QuickHull. Toggle Animate to see the hull update in real-time as points bounce around.</p>
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

  let numPoints = loadSetting(DEMO, 'numPoints', 50);
  let points = new Float32Array(0);
  let velocities = new Float32Array(0);
  let initialPoints = new Float32Array(0);
  let animating = loadSetting(DEMO, 'animating', false);
  let animId = 0;
  const BOUNDS = 1.0;
  const SPEED = 0.008;

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
    velocities = new Float32Array(numPoints * 3);
    for (let i = 0; i < points.length; i++) {
      points[i] = (Math.random() - 0.5) * 2 * BOUNDS;
      velocities[i] = (Math.random() - 0.5) * 2 * SPEED;
    }
    initialPoints = new Float32Array(points);
  }

  function update() {
    if (points.length < 12) return; // need at least 4 points
    const data = hullMesh(points);
    viewer.setMesh(data);
    viewer.setColor(0x8866cc);
    showReadout(data);
  }

  function animateStep() {
    if (!animating) return;
    for (let i = 0; i < points.length; i++) {
      points[i] += velocities[i];
      if (points[i] > BOUNDS) {
        points[i] = BOUNDS;
        velocities[i] = -Math.abs(velocities[i]);
      } else if (points[i] < -BOUNDS) {
        points[i] = -BOUNDS;
        velocities[i] = Math.abs(velocities[i]);
      }
    }
    update();
    animId = requestAnimationFrame(animateStep);
  }

  function toggleAnimate(on: boolean) {
    if (on) {
      animating = true;
      if (velocities.length !== points.length) {
        velocities = new Float32Array(points.length);
        for (let i = 0; i < velocities.length; i++) {
          velocities[i] = (Math.random() - 0.5) * 2 * SPEED;
        }
      }
      animateStep();
    } else {
      animating = false;
      cancelAnimationFrame(animId);
      // Stay at current point positions — don't reset
    }
  }

  controlsEl.appendChild(createSlider('Point Count ', 10, 500, numPoints, 10, v => {
    saveSetting(DEMO, 'numPoints', v);
    numPoints = v;
    generatePoints();
    update();
  }));
  controlsEl.appendChild(createButton('Regenerate', () => { generatePoints(); update(); }));
  controlsEl.appendChild(createCheckbox('Animate', animating, (v) => { saveSetting(DEMO, 'animating', v); toggleAnimate(v); }));
  controlsEl.appendChild(createCheckbox('Wireframe', loadSetting(DEMO, 'wireframe', false), (v) => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  controlsEl.appendChild(readout);

  generatePoints();
  update();
  toggleAnimate(animating);

  return () => {
    animating = false;
    cancelAnimationFrame(animId);
    viewer.dispose();
  };
}
