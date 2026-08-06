// Depth peeling for order-independent transparency, ported (trimmed) from
// NodeDesigner's rendering/depth-peeling.js — itself based on gkjohnson's
// depth-peeling-demo. The gallery's X-Ray mode uses it so boolean results
// render with *correct* see-through layering (every internal wall visible in
// the right order), which two-pass alpha blending cannot guarantee for
// self-overlapping or non-manifold outputs.
//
// Differences from the source: instance-based (no module globals) so each
// viewer owns its resources, no SSAA coupling (callers pass drawing-buffer
// pixels), and only the MeshStandardMaterial mixin path — the demo has no
// custom shader pipeline.

import * as THREE from 'three';

/** Fullscreen quad (from three's postprocessing Pass). */
class FullScreenQuad {
  private mesh: THREE.Mesh;
  private cam = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
  private scene = new THREE.Scene();

  constructor(material: THREE.Material) {
    this.mesh = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material);
    this.mesh.frustumCulled = false;
    this.scene.add(this.mesh);
  }

  get material(): THREE.MeshBasicMaterial {
    return this.mesh.material as THREE.MeshBasicMaterial;
  }

  render(renderer: THREE.WebGLRenderer) {
    renderer.render(this.scene, this.cam);
  }

  dispose() {
    this.mesh.geometry.dispose();
    (this.mesh.material as THREE.Material).dispose();
  }
}

/**
 * MeshStandardMaterial extended with the peeling discards: fragments behind
 * the opaque pass or at/in front of the previous peel layer are dropped.
 * Exact depth equality (no epsilon) — float32 depth + invariant
 * rasterization make equality sufficient, and any bias eats back faces near
 * silhouettes (see the NodeDesigner source for the war story).
 */
export class DepthPeelStandardMaterial extends THREE.MeshStandardMaterial {
  isDepthPeelMaterial = true;
  private _enablePeel = false;
  private _uniforms = {
    nearDepth: { value: null as THREE.Texture | null },
    opaqueDepth: { value: null as THREE.Texture | null },
    resolution: { value: new THREE.Vector2() },
  };

  get enableDepthPeeling() {
    return this._enablePeel;
  }

  set enableDepthPeeling(v: boolean) {
    if (this._enablePeel !== v) {
      this._enablePeel = v;
      this.needsUpdate = true;
    }
  }

  get nearDepth() {
    return this._uniforms.nearDepth.value;
  }

  set nearDepth(v: THREE.Texture | null) {
    if ((this._uniforms.nearDepth.value === null) !== (v === null)) {
      this.needsUpdate = true;
    }
    this._uniforms.nearDepth.value = v;
  }

  set opaqueDepth(v: THREE.Texture | null) {
    this._uniforms.opaqueDepth.value = v;
  }

  get resolution(): THREE.Vector2 {
    return this._uniforms.resolution.value;
  }

  customProgramCacheKey(): string {
    return `peel${Number(this._enablePeel)}|${Number(this._uniforms.nearDepth.value !== null)}`;
  }

  onBeforeCompile(shader: THREE.WebGLProgramParametersWithUniforms) {
    Object.assign(shader.uniforms, this._uniforms);
    shader.fragmentShader = /* glsl */ `
      #define DEPTH_PEELING ${Number(this._enablePeel)}
      #define FIRST_PASS ${Number(this._uniforms.nearDepth.value === null)}
      #if DEPTH_PEELING
      uniform sampler2D nearDepth;
      uniform sampler2D opaqueDepth;
      uniform vec2 resolution;
      #endif
      ${shader.fragmentShader}
    `.replace(
      'void main() {',
      /* glsl */ `
      void main() {
        #if DEPTH_PEELING
        vec2 screenUV = gl_FragCoord.xy / resolution;
        // Behind the opaque pass (1.0 when nothing opaque rendered).
        if ( texture2D( opaqueDepth, screenUV ).r < gl_FragCoord.z ) discard;
        #if ! FIRST_PASS
        // At or in front of the previously peeled layer.
        if ( texture2D( nearDepth, screenUV ).r >= gl_FragCoord.z ) discard;
        #endif
        #endif
      `,
    );
  }
}

const DEPTH_TYPE = THREE.FloatType;

/** Per-viewer depth-peeling renderer. */
export class DepthPeeler {
  private depthA: THREE.DepthTexture | null = null;
  private depthB: THREE.DepthTexture | null = null;
  private opaqueDepth: THREE.DepthTexture | null = null;
  private opaqueTarget: THREE.WebGLRenderTarget | null = null;
  private layers: THREE.WebGLRenderTarget[] = [];
  private quad: FullScreenQuad | null = null;
  private w = 0;
  private h = 0;
  private supported = true;
  private clearColor = new THREE.Color();

  setSize(w: number, h: number) {
    if (w < 10 || h < 10 || (w === this.w && h === this.h)) return;
    this.w = w;
    this.h = h;
    this.disposeTargets();
  }

  private ensure() {
    if (!this.supported || this.opaqueTarget || this.w < 10) return;
    try {
      this.depthA = new THREE.DepthTexture(this.w, this.h, DEPTH_TYPE);
      this.depthB = new THREE.DepthTexture(this.w, this.h, DEPTH_TYPE);
      this.opaqueDepth = new THREE.DepthTexture(this.w, this.h, DEPTH_TYPE);
      this.opaqueTarget = new THREE.WebGLRenderTarget(this.w, this.h, {
        colorSpace: THREE.NoColorSpace,
        depthBuffer: true,
      });
      this.quad = new FullScreenQuad(new THREE.MeshBasicMaterial());
    } catch (e) {
      console.warn('Depth peeling unavailable:', e);
      this.supported = false;
    }
  }

  /**
   * Render `scene` with `peelGroup` depth-peeled over everything else.
   * Both groups must be children of `scene`; lights stay wherever they are.
   * Falls back to a plain render when unsupported.
   */
  render(
    renderer: THREE.WebGLRenderer,
    scene: THREE.Scene,
    camera: THREE.Camera,
    peelGroup: THREE.Object3D,
    numLayers = 5,
  ) {
    this.ensure();
    if (!this.supported || !this.opaqueTarget || !this.quad) {
      renderer.render(scene, camera);
      return;
    }

    while (this.layers.length < numLayers) {
      this.layers.push(new THREE.WebGLRenderTarget(this.w, this.h, {
        colorSpace: THREE.NoColorSpace,
        depthBuffer: true,
      }));
    }
    while (this.layers.length > numLayers) {
      this.layers.pop()!.dispose();
    }
    for (const l of this.layers) {
      if (l.width !== this.w || l.height !== this.h) l.setSize(this.w, this.h);
    }

    const peelVisible = peelGroup.visible;

    // Pass 1: everything except the peel group, capturing opaque depth.
    peelGroup.visible = false;
    this.opaqueTarget.depthTexture = this.opaqueDepth;
    renderer.setRenderTarget(this.opaqueTarget);
    renderer.render(scene, camera);
    renderer.setRenderTarget(null);
    const qm = this.quad.material;
    qm.map = this.opaqueTarget.texture;
    qm.blending = THREE.NoBlending;
    qm.transparent = false;
    qm.depthTest = false;
    qm.depthWrite = false;
    qm.needsUpdate = true;
    this.quad.render(renderer);
    this.opaqueTarget.depthTexture = null as any;

    const clearAlpha = renderer.getClearAlpha();
    renderer.getClearColor(this.clearColor);
    // Peel layers must hold ONLY peeled fragments over transparent black —
    // a scene background would repaint over the opaque pass on composite.
    const savedBackground = scene.background;
    scene.background = null;

    // Passes 2..N: peel front-to-back with everything else hidden (lights
    // stay active — they are not "rendered").
    peelGroup.visible = true;
    const others: THREE.Object3D[] = [];
    for (const child of scene.children) {
      if (child !== peelGroup && !(child as any).isLight) {
        if (child.visible) {
          others.push(child);
          child.visible = false;
        }
      }
    }
    const res = new THREE.Vector2(this.w, this.h);
    for (let i = 0; i < numLayers; i++) {
      const writeDepth = i % 2 === 0 ? this.depthA! : this.depthB!;
      const nearDepth = i === 0 ? null : (i % 2 === 0 ? this.depthB! : this.depthA!);
      this.configureMaterials(peelGroup, true, nearDepth, res);
      const target = this.layers[i];
      target.depthTexture = writeDepth;
      renderer.setRenderTarget(target);
      renderer.setClearColor(0, 0);
      renderer.render(scene, camera);
      renderer.setRenderTarget(null);
      target.depthTexture = null as any;
    }
    for (const child of others) child.visible = true;

    renderer.setClearColor(this.clearColor, clearAlpha);
    scene.background = savedBackground;

    // Composite back-to-front over the opaque pass.
    renderer.autoClear = false;
    for (let i = numLayers - 1; i >= 0; i--) {
      qm.map = this.layers[i].texture;
      qm.blending = THREE.NormalBlending;
      qm.transparent = true;
      qm.depthTest = false;
      qm.depthWrite = false;
      qm.needsUpdate = true;
      this.quad.render(renderer);
    }
    renderer.autoClear = true;

    this.configureMaterials(peelGroup, false, null, res);
    peelGroup.visible = peelVisible;
  }

  private configureMaterials(
    obj: THREE.Object3D,
    peeling: boolean,
    nearDepth: THREE.Texture | null,
    res: THREE.Vector2,
  ) {
    obj.traverse(o => {
      const mat = (o as THREE.Mesh).material as DepthPeelStandardMaterial | undefined;
      if (!mat || !(mat as any).isDepthPeelMaterial) return;
      mat.enableDepthPeeling = peeling;
      mat.nearDepth = nearDepth;
      mat.opaqueDepth = peeling ? this.opaqueDepth : null;
      mat.resolution.copy(res);
      if (peeling) {
        // Straight replace into the layer target; blending happens in the
        // composite pass.
        mat.blending = THREE.CustomBlending;
        mat.blendSrc = THREE.OneFactor;
        mat.blendDst = THREE.ZeroFactor;
        mat.depthWrite = true;
        (mat as any).forceSinglePass = true;
      } else {
        mat.blending = THREE.NormalBlending;
        mat.depthWrite = false;
        (mat as any).forceSinglePass = false;
      }
    });
  }

  private disposeTargets() {
    this.depthA?.dispose();
    this.depthB?.dispose();
    this.opaqueDepth?.dispose();
    this.opaqueTarget?.dispose();
    for (const l of this.layers) l.dispose();
    this.layers = [];
    this.quad?.dispose();
    this.depthA = this.depthB = this.opaqueDepth = null;
    this.opaqueTarget = null;
    this.quad = null;
  }

  dispose() {
    this.disposeTargets();
  }
}
