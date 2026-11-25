# aortamesh design notes

Decisions, status and open questions, collected from planning discussions (Nov 2026).
Update this file as decisions change.

## 1. Goal

Produce a labeled tetrahedral **fluid** mesh of the aorta, with an optional boundary layer
(split into tets), from a Slicer segmentation and a set of cut planes. One script, run on a laptop or
Great Lakes in seconds to minutes, with no license or remote desktop.

The workflow it replaces:

1. Slicer segmentation → STL export.
2. SimModeler (remote, U-M-managed Windows desktop): cut planes, caps, mesh. Crashes often.
3. Export Abaqus 2D *and* 3D: the 3D `.inp` export is unreliable, and the 2D export gives
   explicit boundary triangles.
4. MATLAB script: labels each cap uniquely and converts tets for CHeart.

SimModeler's one real advantage is that every cap gets a unique label. In aortamesh the
**labels come from the plane names by construction**, so no face classification or relabeling
step is needed.

Wall/FSI meshes (TBAD flap, non-manifold geometry) are out of scope for the first version. A
normal aorta *fluid* domain is currently desired.

## 2. Pipeline

| Stage          | What                                                             | Status   |
| -------------- | ---------------------------------------------------------------- | -------- |
| Slicer         | segment, export closed surface (VTP), place planes (`.mrk.json`) | manual   |
| Load           | `io.load_case` → surface + ordered cuts, LPS, float64            | **done** |
| Cross-sections | ray-cast each plane, re-center, clip polygon                     | **done** |
| Cut            | split surface with each clip polygon, remove the piece beyond it | **next** |
| Mesh           | gmsh: reparametrize, boundary layer, cap, fill, label            | todo     |
| Neutral mesh   | `VolumeMesh` dataclass                                           | todo     |
| Writers        | CHeart X/T/B, XDMF, VTU                                          | todo     |
| Checks         | quality, label coverage, cap areas                               | todo     |

## 3. Inputs from Slicer

- **Surface:** smooth the segmentation deliberately (Segmentations → Representations → Closed
  surface → smoothing factor; this changes geometry), export the segment to a model, and save
  it as **VTP, not STL**. STL has 32-bit coordinates and no labels.
- **Coordinates:** Slicer saves models and markups in LPS (the scene itself is RAS). The
  `.mrk.json` has `coordinateSystem` and `coordinateUnits` fields; the VTP has a `SPACE` field
  array. The readers convert both to LPS.
- **Planes:** read `center` and `normal` directly; Slicer computes them from `baseToNode` and
  `objectToBase` when it saves, so they already include handle rotations. The matrices can be
  ignored.
- **Placing planes:** in the 3D view, points snap to the visible surface. Place the center in a
  2D slice view, or drag it with the interaction handles, to get it inside the lumen. It
  doesn't need to be centered: the code re-centers.

### Surface undulations (deferred)

The Slicer surface shows smooth ripples left over from voxel staircasing. Slicer runs
discrete flying edges on the binary labelmap, then `vtkWindowedSincPolyDataFilter` (the
smoothing factor sets its passband; Taubin's λ|μ is the first-order version). The binary mask
fixes the wall only to within a voxel. Where the wall is nearly parallel to the grid, the
staircase error has a wavelength of several voxels, which a low-pass filter can't remove
without also flattening real anatomy. With anisotropic voxels the ripples usually follow the
axial slices, one per slice thickness.

It matters for WSS (curvature-sensitive, shows as stripes), much less for pressure drop. gmsh
remeshing won't remove it: `classifySurfaces` reparametrizes the triangulation it is given.

Options, best first:
1. Isosurface of the CT intensity (partial-volume gives sub-voxel wall position), restricted
   to a dilated mask, at the lumen/wall midpoint threshold. Recovers information the binary
   mask discarded.
2. Anti-alias the binary mask: ITK `AntiAliasBinaryImageFilter` (Whitaker 1998; in
   SimpleITK), the smoothest level set consistent with the inside/outside voxels, then
   extract the zero isosurface. Doesn't shrink small vessels. Cruder version: Gaussian-blur
   the mask (σ ≈ 1 voxel), take the 0.5 isosurface.
3. In Slicer: oversample the segmentation (Segment Editor → Specify geometry) and smooth the
   labelmap before surface generation. Reduces the ripples, adds no information.
4. Smooth the surface harder: rounds ostia and narrows small vessels. Least attractive.

1 and 2 would add an optional stage upstream of `surface.vtp` (reading `.seg.nrrd`, plus
the CT for 1) and likely a SimpleITK dependency. To confirm the cause first, compare the
ripple wavelength and direction with the voxel spacing. Not pursuing for now (Nov 2026).

## 4. Cross-sections and clip polygons (`clip.py`, done)

For each plane:

1. Cast 72 rays from the center, **within the plane**, and record the first wall hit along
   each (and the second surface hit beyond it, if any).
2. Move the center to the area centroid of the hit polygon and recast. Repeat until the step
   is < 1e-3 × equivalent radius (2 iterations on the example case).
3. Clip polygon: star-shaped polygon through the hits, pushed out to `r·(1 + margin) + pad`
   (defaults 0.25 and 1 mm), but never more than halfway to the next surface along that ray.
   Fan-triangulated from the center. The polygon spans this lumen and never reaches a
   neighbouring vessel, unlike an infinite plane.

Errors and warnings:
- Any ray missing the wall → error (center outside the lumen).
- Any ray > 2× the median radius → warning (ray escaped into a branch ostium, or the plane is
  badly oblique; move the plane).

Example case (`with_arch_branches`):

| Plane   | r_eq (mm) | Center moved (mm) | Nearest other vessel (mm) |
| ------- | --------- | ----------------- | ------------------------- |
| inlet   | 14.4      | 2.4               | 35                        |
| branch1 | 4.6       | 0.4               | 3.2                       |
| branch2 | 2.6       | 0.3               | 3.8                       |
| branch3 | 4.8       | 0.3               | none                      |
| outlet  | 10.4      | 0.6               | none                      |

Verified: hits lie on the wall (< 1e-4 mm); no polygon vertex is inside the lumen; every
surface point inside a polygon belongs to that polygon's own vessel (< 0.08 mm). Before the
halfway limit, the branch1 polygon came within ~1 mm of the neighbouring branch.

## 5. Next: cutting

Plan:
- Cut the surface along each clip polygon (only where the polygon is, not the whole plane).
- Of the two pieces, **discard the smaller by area**: the stub beyond the cut (branch tip,
  root/LV below the inlet). The main aorta is always the larger piece.
- Log the area ratio; a ratio near 1 means the plane is in a strange place (e.g. cutting the
  arch in half). Stop and report rather than guess.
- Flip each stored normal to point into the discarded piece, giving consistent outward cap
  normals (also needed for inflow profiles).
- **Don't cap in pyvista** if a boundary layer will be used (see §6). Output the open wall
  surface plus, for each opening, its plane and label name.
- Cut planes should be close to perpendicular to the vessel, or the boundary layer at the rim
  comes out skewed. The centerline (`centerline.vtp`) could be used to check or set this.

## 6. Meshing in gmsh

- Merge the open wall surface; `classifySurfaces` + `createGeometry` so gmsh can
  reparametrize and remesh it. Feed it a reasonably smooth, well-shaped triangulation.
- **Boundary layer:** `gmsh.model.geo.extrudeBoundaryLayer`, extruding **from the wall only**,
  then cap. If caps exist when extruding, rim nodes average the wall and cap normals and the
  layer collapses at every inlet/outlet. This is why gmsh's `aneurysm.py` example (copy in
  `mesh.py`) extrudes first and caps after.
- `recombine=False` → prisms split into tets (CHeart and dolfinx need linear tets only).
- Geometric growth: one element per layer, **cumulative** heights, sign on every entry
  (negative = inward for that example's normals). Heights are absolute distances, unlike
  `geo.extrude` where they're fractions.
  ```python
  heights = list(h1 * np.cumsum(r ** np.arange(n)))
  extrudeBoundaryLayer(walls, numElements=[1]*n, heights=[-h for h in heights],
                       recombine=False, second=...)
  ```
- Keep total layer thickness < ~1/3 of the smallest local radius (converging normals in
  concave corners otherwise fold over). Match the core element size to the last layer size,
  `h1·r^(n-1)`. `viewIndex` can scale the thickness per node if needed.
- gmsh does not detect folded layers. Check `getElementQualities` (min SICN > 0).
- **Caps and labels:** find the open curve loops on the inner layer surface, match each loop
  to the nearest plane → that plane's label, then `addPlaneSurface`. The rims are planar by
  construction.
- **Outlet labels include the boundary layer's annular side faces**, not only the core's flat
  cap. Easy to miss.
- Size guide: one tet ≈ 0.12·h³ of volume. For a ~160 mL aorta, h = 2 mm gives roughly
  150–250k tets. Curvature refinement (`h_min`/`h_max`) helps the small arch branches.

Open: remesh the surface in gmsh (size control, curvature refinement), or keep the smoothed
pyvista surface and only fill the volume? The second is simpler but fixes surface element
size. Leaning toward gmsh remeshing.

## 7. Neutral mesh

```python
@dataclass(frozen=True)
class VolumeMesh:
    points:     np.ndarray  # (N, 3) float64, LPS, mm
    tets:       np.ndarray  # (M, 4) int64, 0-based
    bnd_tris:   np.ndarray  # (K, 3) int64
    bnd_labels: np.ndarray  # (K,)   int32
    bnd_tet:    np.ndarray  # (K,)   int64, parent tet (needed for .B)
    label_names: dict[int, str]
    # validate(), to_pyvista()
```

Every writer reads only this. Parent tets: build all tet faces as sorted triples and match
(`np.unique` or a dict).

## 8. Writers

**CHeart** (`.X` nodes, `.T` linear tets, `.B` labeled boundary triangles):
- `.B` rows are `element n1 n2 n3 label`, 1-based. Each triangle must be a face of the tet
  in its first column.
- `.T` holds linear tets only.
- **TODO:** confirm the exact header lines from a real SimModeler-derived file set.
- **TODO:** copy label numbering from the existing MATLAB script, so existing CHeart input
  files keep working.
- Don't depend on Will's `cheartpy`: it needs physical groups named `Surface1…`, and its
  `get_gmsh_space` calls `get_nodes(dim=3)`, which silently drops surface-classified nodes
  and corrupts gmsh-native meshes. A small writer from `VolumeMesh` is simpler.

**XDMF** (for FEniCSx) via meshio/h5py, so dolfinx isn't a dependency. **VTU** for ParaView.

## 9. Units

CHeart has no units: mesh length unit + time unit + the density and viscosity entered define
them, and outputs follow (mm mesh + seconds → velocity in mm/s).

| Density entered | System         | Velocity | Pressure |
| --------------- | -------------- | -------- | -------- |
| 1060            | SI (mesh in m) | m/s      | Pa       |
| 1.06e-3         | mm–g–s         | mm/s     | Pa       |
| 1.06e-6         | mm–kg–s        | mm/s     | kPa      |
| 1.06            | cm–g–s         | cm/s     | 0.1 Pa   |

mm–g–s is convenient: viscosity stays ~0.004 and pressure comes out in Pa. A mm mesh with
density 1060 silently solves a different problem. The lab's actual convention is unconfirmed:
check the density in an existing CHeart input file. Store the mesh unit and give writers a
`scale` (CHeart in mm; FEniCSx solver in m, `1e-3`).

## 10. Checks (run every time)

- Zero untagged boundary faces.
- Every `.B` triangle is a face of its listed tet.
- Min SICN > 0.
- Each cap's area ≈ the cross-section area from §4.
- Validation vs SimModeler: one geometry through both pipelines; compare time per
  iteration, element quality, and pressure drop.

## 11. Open questions

1. Label numbering expected by existing CHeart inputs (from the MATLAB script).
2. Exact `.X/.T/.B` header format.
3. Lab's unit convention (density in an existing input file).
4. gmsh surface remeshing vs. keeping the pyvista surface (§6).
5. Boundary-layer parameters (first-layer height from Reynolds number, growth, count).

Decided: `requires-python` stays `>=3.14` (Nov 2026).
