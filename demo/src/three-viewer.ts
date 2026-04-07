// Reusable Three.js 3D viewer component for mesh display

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import type { MeshData } from './wasm.ts';

export class ThreeViewer {
  private renderer: THREE.WebGLRenderer;
  private scene: THREE.Scene;
  private camera: THREE.PerspectiveCamera;
  private controls: OrbitControls;
  private backMesh: THREE.Mesh | null = null;
  private frontMesh: THREE.Mesh | null = null;
  private wireMesh: THREE.LineSegments | null = null;
  private showWireframe = false;
  private animId = 0;
  private disposed = false;
  private hasFramed = false;

  constructor(container: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
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

    const hasColors = data.has_colors && data.colors;
    let isTransparent = false;

    if (hasColors) {
      // Colors are RGBA (4 floats per vertex), split into RGB + alpha
      const colorsRGBA = data.colors!;
      const vertCount = colorsRGBA.length / 4;
      const rgb = new Float32Array(vertCount * 3);
      for (let i = 0; i < vertCount; i++) {
        rgb[i * 3] = colorsRGBA[i * 4];
        rgb[i * 3 + 1] = colorsRGBA[i * 4 + 1];
        rgb[i * 3 + 2] = colorsRGBA[i * 4 + 2];
        if (colorsRGBA[i * 4 + 3] < 0.999) isTransparent = true;
      }
      geom.setAttribute('color', new THREE.BufferAttribute(rgb, 3));
    }

    if (isTransparent) {
      // Per-vertex alpha via separate attribute + shader patch
      const colorsRGBA = data.colors!;
      const vertCount = colorsRGBA.length / 4;
      const alphas = new Float32Array(vertCount);
      for (let i = 0; i < vertCount; i++) {
        alphas[i] = colorsRGBA[i * 4 + 3];
      }
      geom.setAttribute('aAlpha', new THREE.BufferAttribute(alphas, 1));

      // Shader patch: pass per-vertex alpha through to fragment, multiply gl_FragColor.a
      const patchAlpha = (mat: THREE.MeshStandardMaterial) => {
        mat.onBeforeCompile = (shader) => {
          // Vertex: declare attribute + varying, assign in main
          shader.vertexShader = 'attribute float aAlpha;\nvarying float vAlpha;\n' + shader.vertexShader;
          shader.vertexShader = shader.vertexShader.replace(
            '#include <begin_vertex>',
            '#include <begin_vertex>\nvAlpha = aAlpha;'
          );
          // Fragment: declare varying, multiply final alpha
          shader.fragmentShader = 'varying float vAlpha;\n' + shader.fragmentShader;
          shader.fragmentShader = shader.fragmentShader.replace(
            '#include <dithering_fragment>',
            '#include <dithering_fragment>\ngl_FragColor.a *= vAlpha;'
          );
        };
      };

      // Two-pass rendering for correct transparency:
      // Pass 1 (renderOrder 0): back faces, no depth write — so back faces show behind transparent front faces
      const backMat = new THREE.MeshStandardMaterial({
        color: 0xffffff,
        vertexColors: true,
        roughness: 0.5,
        metalness: 0.1,
        side: THREE.BackSide,
        flatShading: true,
        transparent: true,
        depthWrite: false,
      });
      patchAlpha(backMat);
      this.backMesh = new THREE.Mesh(geom, backMat);
      this.backMesh.renderOrder = 0;
      this.scene.add(this.backMesh);

      // Pass 2 (renderOrder 1): front faces, depth write on — opaque fragments block, transparent blend
      const frontMat = new THREE.MeshStandardMaterial({
        color: 0xffffff,
        vertexColors: true,
        roughness: 0.5,
        metalness: 0.1,
        side: THREE.FrontSide,
        flatShading: true,
        transparent: true,
        depthWrite: true,
      });
      patchAlpha(frontMat);
      this.frontMesh = new THREE.Mesh(geom, frontMat);
      this.frontMesh.renderOrder = 1;
      this.scene.add(this.frontMesh);
    } else {
      // Opaque single-pass
      const mat = new THREE.MeshStandardMaterial({
        color: hasColors ? 0xffffff : 0x4488cc,
        vertexColors: !!hasColors,
        roughness: 0.5,
        metalness: 0.1,
        side: THREE.DoubleSide,
        flatShading: true,
      });
      this.frontMesh = new THREE.Mesh(geom, mat);
      this.scene.add(this.frontMesh);
    }

    if (this.showWireframe) {
      this.addWireframe(geom);
    }

    if (!this.hasFramed) {
      this.frameCamera(geom);
      this.hasFramed = true;
    }
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
    const mesh = this.frontMesh || this.backMesh;
    if (show && mesh) {
      this.addWireframe(mesh.geometry);
    }
  }

  setColor(color: number) {
    if (this.frontMesh) {
      const mat = this.frontMesh.material as THREE.MeshStandardMaterial;
      if (!mat.vertexColors) {
        mat.color.set(color);
      }
    }
  }

  private clearMesh() {
    const sharedGeom = this.backMesh && this.frontMesh && this.backMesh.geometry === this.frontMesh.geometry;
    if (this.backMesh) {
      this.scene.remove(this.backMesh);
      if (!sharedGeom) this.backMesh.geometry.dispose();
      (this.backMesh.material as THREE.Material).dispose();
      this.backMesh = null;
    }
    if (this.frontMesh) {
      this.scene.remove(this.frontMesh);
      this.frontMesh.geometry.dispose();
      (this.frontMesh.material as THREE.Material).dispose();
      this.frontMesh = null;
    }
    if (this.wireMesh) {
      this.scene.remove(this.wireMesh);
      this.wireMesh.geometry.dispose();
      this.wireMesh = null;
    }
  }

  resetFrame() {
    this.hasFramed = false;
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
