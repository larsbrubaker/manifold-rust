// Reusable Three.js 3D viewer component for mesh display

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import type { MeshData } from './wasm.ts';

export class ThreeViewer {
  private renderer: THREE.WebGLRenderer;
  private scene: THREE.Scene;
  private camera: THREE.PerspectiveCamera;
  private controls: OrbitControls;
  private solidMesh: THREE.Mesh | null = null;
  private wireMesh: THREE.LineSegments | null = null;
  private showWireframe = false;
  private animId = 0;
  private disposed = false;

  constructor(container: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setClearColor(0xf0f0f0);
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();

    this.camera = new THREE.PerspectiveCamera(50, 1, 0.01, 100);
    this.camera.position.set(2, 2, 3);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.1;

    // Lighting
    const ambient = new THREE.AmbientLight(0xffffff, 0.5);
    this.scene.add(ambient);
    const dir = new THREE.DirectionalLight(0xffffff, 0.8);
    dir.position.set(3, 5, 4);
    this.scene.add(dir);
    const dir2 = new THREE.DirectionalLight(0xffffff, 0.3);
    dir2.position.set(-3, -2, -4);
    this.scene.add(dir2);

    // Grid helper
    const grid = new THREE.GridHelper(4, 20, 0xcccccc, 0xe0e0e0);
    this.scene.add(grid);

    // Resize handling
    const ro = new ResizeObserver(() => this.resize());
    ro.observe(container);

    this.resize();
    this.animate();
  }

  private resize() {
    const parent = this.renderer.domElement.parentElement;
    if (!parent) return;
    const w = parent.clientWidth;
    const h = parent.clientHeight;
    this.renderer.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  private animate() {
    if (this.disposed) return;
    this.animId = requestAnimationFrame(() => this.animate());
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  setMesh(data: MeshData) {
    this.clearMesh();

    const geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
    geom.setAttribute('normal', new THREE.BufferAttribute(data.normals, 3));
    geom.setIndex(new THREE.BufferAttribute(data.indices, 1));

    const mat = new THREE.MeshStandardMaterial({
      color: 0x4488cc,
      roughness: 0.5,
      metalness: 0.1,
      side: THREE.DoubleSide,
    });
    this.solidMesh = new THREE.Mesh(geom, mat);
    this.scene.add(this.solidMesh);

    if (this.showWireframe) {
      this.addWireframe(geom);
    }

    this.frameCamera(geom);
  }

  private addWireframe(geom: THREE.BufferGeometry) {
    const wireGeom = new THREE.WireframeGeometry(geom);
    const wireMat = new THREE.LineBasicMaterial({ color: 0x000000, opacity: 0.15, transparent: true });
    this.wireMesh = new THREE.LineSegments(wireGeom, wireMat);
    this.scene.add(this.wireMesh);
  }

  setWireframe(show: boolean) {
    this.showWireframe = show;
    if (this.wireMesh) {
      this.scene.remove(this.wireMesh);
      this.wireMesh.geometry.dispose();
      this.wireMesh = null;
    }
    if (show && this.solidMesh) {
      this.addWireframe(this.solidMesh.geometry);
    }
  }

  setColor(color: number) {
    if (this.solidMesh) {
      (this.solidMesh.material as THREE.MeshStandardMaterial).color.set(color);
    }
  }

  private clearMesh() {
    if (this.solidMesh) {
      this.scene.remove(this.solidMesh);
      this.solidMesh.geometry.dispose();
      (this.solidMesh.material as THREE.Material).dispose();
      this.solidMesh = null;
    }
    if (this.wireMesh) {
      this.scene.remove(this.wireMesh);
      this.wireMesh.geometry.dispose();
      this.wireMesh = null;
    }
  }

  private frameCamera(geom: THREE.BufferGeometry) {
    geom.computeBoundingSphere();
    const sphere = geom.boundingSphere;
    if (!sphere) return;
    const center = sphere.center;
    const radius = sphere.radius || 1;
    this.controls.target.copy(center);
    const dist = radius / Math.tan((this.camera.fov / 2) * Math.PI / 180) * 1.5;
    this.camera.position.set(center.x + dist * 0.5, center.y + dist * 0.5, center.z + dist);
    this.controls.update();
  }

  dispose() {
    this.disposed = true;
    cancelAnimationFrame(this.animId);
    this.clearMesh();
    this.controls.dispose();
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
