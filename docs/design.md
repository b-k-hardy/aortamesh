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

| Stage          | What                                                                | Status   |
| -------------- | ------------------------------------------------------------------- | -------- |
| Slicer         | segment, export closed surface (VTP), place planes (`.mrk.json`)    | manual   |
| Load           | `io.load_case` → surface + ordered cuts, LPS, float64               | **done** |
| Cross-sections | ray-cast each plane, re-center, clip polygon                        | **done** |
| Cut            | split surface with each clip polygon, remove the piece beyond it    | **done** |
| Mesh           | gmsh: reparametrize, boundary layer, cap, fill, label               | **done** |
| Neutral mesh   | `VolumeMesh` dataclass                                              | **done** |
| Config         | TOML → `Config` / `MeshOptions` (`config.py`, `sample_config.toml`) | **done** |
| Writers        | CHeart X/T/B (**done**, `cheart.py`); VTU via `to_pyvista`; XDMF    | partial  |
| Checks         | quality, label coverage, cap areas                                  | partial  |

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

## 5. Cutting (`cut.py`, done)

`cut_surface(surface, sections) -> (wall, openings)`. Cuts run in order, each on the
already-cut surface. For each section, only the plane inside its clip polygon is used:

1. Snap vertices within `snap` × median edge length (default 0.2) of the plane, and inside
   the polygon, onto the plane. Without this the cut leaves slivers next to vertices.
2. The triangles inside the polygon that cross or touch the plane form a closed strip around
   the vessel. Without the strip, the rest of the surface splits into one piece per side.
   Anything other than exactly two pieces → error.
3. **Discard the smaller piece by area** (branch tip, root below the inlet). Ratio >
   `max_area_ratio` (0.5) → error: the plane probably cuts across the main vessel.
4. Clip the strip triangles exactly at the plane, keeping the kept side. Points on crossed
   edges are shared, so the wall stays watertight apart from the openings, and each rim
   lies in its plane.

Each `Opening` has the name, section center, **outward normal (into the discarded piece)**,
ordered rim points, cap area, discarded area and area ratio. **No caps** (boundary layer first,
§6). The wall has no `label` cell data yet (numbering is open, §11). `boundary_loops(tris)`
orders the open edges into loops.

Example case (`with_arch_branches`), about 1.2 s for 5 cuts:

| Plane   | Cap area (mm²) | Cap / section | Discarded (mm²) | Area ratio |
| ------- | -------------- | ------------- | --------------- | ---------- |
| inlet   | 649.1          | 1.0017        | 2655            | 0.096      |
| branch1 | 67.3           | 1.0018        | 121             | 0.004      |
| branch2 | 21.1           | 1.0006        | 48              | 0.002      |
| branch3 | 71.3           | 1.0000        | 187             | 0.006      |
| outlet  | 339.3          | 1.0012        | 1437            | 0.055      |

Verified:
- Topology: exactly one boundary loop per plane, no non-manifold edges.
- Rims: on their planes to < 6e-14 mm, with radii matching the ray-cast radii.
- Cap areas: within 0.2% of the section area. Before snapping they were all 1.0013, the
  expected ratio for a 72-ray polygon inscribed in a much finer rim.
- Normals: no kept wall vertex near a rim lies on the outward side, and 1 mm out along the
  normal is inside the original lumen.
- Unchanged surface: only the 184 snapped vertices moved (max 0.09 mm, 0.14 of an edge);
  every other point is exactly on the original surface. Kept plus discarded area matches the
  original to 0.18 mm² (the snapping).
- Mesh quality: worst triangle quality 0.21 (original 0.44). With `snap` = 0.05, 17
  triangles fell below 0.1.
- Rejection: a plane across the mid-aorta (from the centerline) is rejected (ratio 0.79).

Still open: cut planes should be close to perpendicular to the vessel, or the boundary layer
at the rim comes out skewed. The centerline could be used to check or set this.
`centerline.vtp` holds 4 polylines, each a path from the inlet to one end (outlet,
branch3, branch2, branch1). They overlap along the ascending aorta and arch, so a polyline
is not one vessel.

## 6. Meshing in gmsh (`mesh.py`, done)

`mesh_volume(wall, openings, MeshOptions, centerline) -> VolumeMesh`. The wall goes in through the API as a
discrete entity (no temp files), gmsh remeshes it, and the boundary layer is extruded inward
(height sign **negative** for our outward-oriented walls). Each hole in the inner surface gets
a plane cap, and the core volume fills the inside. Boundary labels are assigned after meshing:
wall entities → 1, and every other boundary surface (caps plus the layer's side faces) goes to
the nearest opening's cap disk → 2, 3, … in cut order. **This numbering is provisional** (§11).

Options (`MeshOptions`, `BoundaryLayer`; all also in the TOML): `size_max` 2.0,
`size_min` 0.5, `radius_elements` 16, `curvature_elements` 20, `optimize_netgen` on,
`min_quality` 0.05 (warn), and a boundary layer of 3 layers, first height 0.15 mm, growth
1.5 (0.71 mm total). `sample_config.toml` currently uses first height 0.2 (0.95 mm total,
which triggers the thickness warning). Set `boundary_layer = None` / `enabled = false` to
turn the layer off.

**Sizing.** `radius_elements` N sets the size to 2πr/N, with r the `Radius` of the nearest
`centerline.vtp` point, through gmsh's size callback. On **curves** (patch seams, rims) this
size replaces gmsh's. On **surfaces and in the volume** the callback takes the smaller of the
two, so curvature sizing still refines tight spots. Both are needed:
- Curvature alone refines every patch seam into a band of tiny elements (next bullet list).
- Radius alone is too coarse where the branches join the aorta: those wall points are
  nearest the aorta centerline (2 mm) but the wall bends sharply. That gave 5 tets with
  SICN ≈ 0, and without a boundary layer gmsh failed ("segment and facet intersect").
- Combined, on the example case: wall 19k triangles (26k with curvature alone). Within
  0.5 mm of a seam: 0.7% of triangles (was 16%). Min SICN 0.096–0.126 with a layer, 0.26
  without one. In branches (r < 8 mm) curvature is finer and wins: edges ≈ 0.7 × 2πr/16.

Example case. With curvature sizing only and the 0.71 mm layer: 331,573 tets (236k in the
layer, 96k core), 63k nodes, 23 s. With combined sizing and the 0.95 mm layer from
`sample_config.toml`: 256,624 tets (174k in the layer), 48k nodes, 26 s, min SICN 0.126, and
the layer at exactly 0.2/0.5/0.95 mm. The checks below were run on the first and re-run on the
second.
Verified:
- `validate()` passes: no untagged faces, each boundary triangle is a face of its tet, the
  boundary is closed, and every tet has positive volume. All boundary faces point outward.
- Wall: remeshed wall nodes lie exactly on the cut wall (max distance 0.0000 mm).
- Volume: 0.16% below the volume enclosed by the cut wall (chordal error from remeshing).
- Boundary layer: every one of the 13,206 wall nodes has a copy at 0.150, 0.375 and
  0.7125 mm from the wall, exactly the cumulative heights.
- Cap areas projected onto each plane are within 1% of the opening area.
- Quality: min SICN 0.107 with Netgen (0.018 without). Without the boundary layer, one
  sliver (SICN 0.004) survives, and the warning reports it.

Findings:
- **Sizing trap:** with gmsh's defaults `Mesh.MeshSizeFromPoints = 1` and
  `MeshSizeExtendFromBoundary = 1`, sizes come from the input triangulation's nodes, so the
  wall keeps its 0.65 mm edges whatever `size_max` says (1.87M tets). Both are set to 0.
- **Patch seams:** `classifySurfaces` splits the wall into 15 patches (29 seam curves) for
  reparametrization. The seams follow the input's zigzag edges, so curvature sizing on curves
  made a dense ring at every seam: 69% of the smallest 5% of triangles were within 0.5 mm of
  a seam. Fixed by radius sizing on curves (above).
- **Caps are not flat with a boundary layer.** The layer's inner rim follows the wall
  normal, and where the wall meets the plane obliquely the rim leaves the plane: up to
  0.59 mm at the inlet (the wall there is up to 60° from perpendicular to the plane),
  0.19–0.29 mm elsewhere (0.75 mm at the inlet with the 0.95 mm layer). The offset is
  ≈ thickness × sin(tilt). gmsh's `aneurysm.py` has the same effect but only ≤ 0.15 (tilt
  ≤ 18°, and its STL openings aren't flat to begin with), and it never checks. Without the
  layer, caps are exactly flat. More precisely (0.95 mm layer): the cap's centre is flat
  (residual 0.0000 mm) but lies on gmsh's best-fit plane through the non-flat inner rim. That
  plane is tilted 0.5° (inlet, outlet) to 3.9° (branch1) from the cut plane and shifted up
  to 0.32 mm along it. The layer's ring-shaped side faces, labelled as cap, tilt with the
  wall (up to 50° at the inlet).
  Fixes (open question):
  1. Our own extrusion normals through `extrudeBoundaryLayer(..., viewIndex)`. **Prototyped
     and working: see "Potential direction" below.**
  2. Build the layer ourselves in numpy and let gmsh fill only the core.
  3. Flow extensions (cylinders along the plane normal), so the wall meets the cap at 90°.
     gmsh's `naca_boundary_layer_3d.py` gets a flat symmetry plane this way: the wing meets
     it at 90°.
  With 1 or 2, a steeply tilted rim gives skewed layer cells (perpendicular thickness about
  cos(tilt) at the rim), so more perpendicular planes still help.
- **Most tets are in the boundary layer** (3 layers × 3 tets per prism × wall triangles).
  The last layer (0.34 mm) is much thinner than the core edges (median 2.4 mm), a 7× jump.
  The total thickness is capped by the smallest branch (r ≈ 2.6 mm → < 0.86 mm), so a smooth
  transition needs a finer core near the wall (a size field) or a thinner wall mesh.

### Potential direction: custom extrusion normals (prototyped Nov 2026, not built)

`extrudeBoundaryLayer(..., viewIndex=v)` with a **vector** view uses the view as the extrusion
normals, and the **length of each vector multiplies the layer heights** at that point. The
number of layers stays the same everywhere. (A scalar view only scales the thickness.) That
gives two features with no change to the rest of the pipeline.

**How the prototype did it.** All of this happens before meshing, and there's still one
`generate(3)`.
- Outward area-weighted vertex normals on the *input* cut wall.
- Within `blend` mm of a rim, remove each normal's component along the cap normal `c`:
  `n − (1 − w)(n·c)c`, with `w` a smoothstep from 0 at the rim to 1 at `blend`, then
  renormalize.
- Pass them as a list-based view: `gmsh.view.addListData(v, "VT", ...)`, a vector at each
  corner of every input triangle.
- gmsh looks a list view up **by position**, and the remeshed wall nodes lie on the input
  surface. So the view doesn't have to wait for the wall to be remeshed, which was the
  worry before prototyping.
- Height sign is the same as without a view: negative heights grow the layer against the
  vectors, i.e. inward for outward normals.

**Flat caps** (unit-length normals, `blend` = 3 mm, 0.83 mm layer):

| | Gmsh normals (current) | Bent normals |
| --- | --- | --- |
| Cap nodes off the cut plane, max | 0.70 mm inlet, 0.22–0.37 mm others | ≤ 2e-13 mm, every cap incl. side ring |
| Layer heights away from rims | 0.175 / 0.437 / 0.831 | unchanged |
| `validate`, volume vs enclosed | ok, −0.18% | ok, −0.27% |
| Worst tet (min SICN) | 0.112 | 0.089 at the inlet rim (2 tets < 0.1); other rims ≥ 0.14 |

`blend` = 6 mm made almost no difference (worst 0.091). The trade-off is the expected one:
the layer leans where the wall is tilted, and its thickness perpendicular to the wall at the
rim is cos(tilt). That's 0.51× at the inlet (tilt 59.6°), 0.90–0.93× at the branches and
0.96× at the outlet.

**Vector length.** The layer tops land at local length × configured height:
- uniform 0.5: distance/length = 0.1749 / 0.4374 / 0.8310 in the aorta, the configured
  heights exactly
- local radius / 10 mm (0.23–1.59): the same ratios in every vessel size (77–98% of layer
  nodes within 2% of a layer top; the spread comes from the nearest-centerline lookup near
  junctions)
- uniform 2: gmsh fails ("segment and facet intersect"). 1.66 mm is too thick for branch2.

**What it enables:**
1. **Full thickness at tilted rims:** scale the bent rim normals by 1/cos(tilt) to keep the
   thickness perpendicular to the wall (about 2× at the inlet's steepest point). The layer
   then reaches further across the cap plane.
2. **Thickness that follows the vessel** with a fixed number of layers, e.g. ∝ local
   radius. The aorta gets a thick layer and the branches a thin one, and the one-third limit
   holds everywhere, instead of the smallest branch setting one thickness for all (§6
   Findings, open question 5).

**Before building it:**
- Smooth the length field over the wall. The nearest centerline radius jumps at branch
  junctions.
- Thin layers mean lower quality. The uniform 0.5 run had 43k tets with SICN < 0.1 (worst
  0.048), and the radius/10 run had 1.9k (worst 0.055). That's expected for thin layers,
  but `min_quality` (0.05) would start warning.
- With a fixed layer count, the first-layer height scales too (about 0.04 mm in branch2 at
  radius/10). Whether thickness should go with r, r/√Re or similar is a modelling choice.
- Integration sketch: build the normals (and lengths) in `mesh_volume`, add the view before
  `_extrude_boundary_layer`, and pass its index. New `BoundaryLayer` options might be
  `cap_blend` (mm; 0 = current behaviour) and a thickness mode (`uniform` / `radius`).
- The prototype script was a scratchpad file, not kept in the repo. Everything it did is
  described above.

Background and gmsh notes:

- Merge the open wall surface; `classifySurfaces` + `createGeometry` so gmsh can
  reparametrize and remesh it. Feed it a reasonably smooth, well-shaped triangulation.
- **Boundary layer:** `gmsh.model.geo.extrudeBoundaryLayer`, extruding **from the wall only**,
  then cap. If caps exist when extruding, rim nodes average the wall and cap normals and the
  layer collapses at every inlet/outlet. This is why gmsh's `aneurysm.py` example (in git at
  `3bb40e8`) extrudes first and caps after.
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
  `h1·r^(n-1)`. `viewIndex` can scale the thickness per node (scalar view) or set
  the normals (vector view).
- gmsh does not detect folded layers. Check `getElementQualities` (min SICN > 0).
- **Caps and labels:** find the open curve loops on the inner layer surface, then
  `addPlaneSurface`. Labels are matched to planes after meshing, by geometry. The outer
  rims are planar by construction; the inner rims are not (see Findings).
- **Outlet labels include the boundary layer's annular side faces**, not only the core's flat
  cap. Easy to miss.
- Size guide: one tet ≈ 0.12·h³ of volume. For a ~160 mL aorta, h = 2 mm gives roughly
  150–250k tets. Curvature refinement (`h_min`/`h_max`) helps the small arch branches.

Built: gmsh remeshing of the surface. Keeping the pyvista surface as-is (no remesh) isn't
implemented.

## 7. Neutral mesh (`volume.py`, done)

As below, plus `from_tets` (finds the boundary and parent tets, makes every tet positively
oriented, winds boundary triangles outward), `volume`, `label_areas`, `summary`,
`to_pyvista` (tets) and `boundary_pyvista` (labeled boundary). The CLI writes `mesh.vtu` and
`boundary.vtp` with these.

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

**CHeart** (`cheart.py`, done): `write_cheart(vm, dir, CheartOptions(name, scale))` and
`read_cheart`. Whitespace-separated text, one header line, every index 1-based:

| File | Header | Rows (one per)                                                |
| ---- | ------ | ------------------------------------------------------------- |
| `.X` | `N 3`  | node: `x y z` (`%.17g`, round-trips float64)                  |
| `.T` | `M N`  | tet: `n1 n2 n3 n4` (rows of `.X`)                             |
| `.B` | `K`    | boundary triangle: `tet n1 n2 n3 label` (`tet` = row of `.T`) |

- Nodes and tets are written in `VolumeMesh` order, so `.B`'s indices agree with `.T` and
  `.X` by construction. `write_cheart` runs `validate()` first, which now also rejects nodes
  no tet uses.
- `scale` multiplies coordinates only (1.0 = mm, 1e-3 = m; §9). `.T`/`.B` don't change.
- `normals = true` also writes `normals.D`, structured like `.X` (header `N 3`, row i = node
  i). Each row is the outward unit normal from `VolumeMesh.node_normals()`: the area-weighted
  average of the adjacent boundary triangles' normals, and zero at interior nodes. It's
  computed in numpy, not pyvista's `compute_normals`, which can split or renumber points.
  Not scaled. Also in `boundary.vtp` as `point_data["normal"]`.
  - Verified (0.95 mm layer): 48,424 rows match `.X`, and the non-zero rows are exactly the
    10,451 nodes in `.B`, all unit length. A step of 0.05 or 0.01 mm along +n leaves the
    domain and along −n enters it, at every boundary node. Wall nodes agree with pyvista
    to a median of 0.24° (p99 2.1°; the gap is area weighting).
  - Rim nodes (136, where the wall meets a cap) average both surfaces. The thin side-ring
    triangles count for little, so the normal leans toward the wall: median 81° from the
    cap normal (VTK's unweighted average would give about 45°). There's no single correct
    normal at a corner. If a solver needs cap-only or wall-only normals there, they'd have
    to be computed per label.
- Tets are positively oriented and `.B` triangles wound outward. Whether CHeart needs a
  particular node order within a tet is unconfirmed.
- Verified on the example case, parsing the files independently of the package:
  - headers match the row counts; all indices fall in 1..N and 1..M, and every node is used
  - all 20,898 `.B` triangles are faces of their listed tet, and the set of `.B` triangles is
    exactly the set of boundary faces of `.T`
  - every tet read from the files has positive volume, and all `.B` triangles are wound
    outward
  - coordinates round-trip exactly; label counts match
  - `read_cheart` gives back identical arrays, and `validate()` passes on them
- **TODO:** confirm against a real SimModeler-derived file set: header lines, any trailing
  metadata, tet node order.
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
4. Non-flat caps with a boundary layer (§6 Findings). Custom extrusion normals fix it
   (prototyped, §6 "Potential direction"); not yet built. Alternatives: more perpendicular
   planes, or flow extensions.
5. Boundary-layer parameters (first-layer height from Reynolds number, growth, count), and
   the layer-to-core size jump (§6 Findings). Variable thickness with a fixed layer count is
   possible through the view's vector lengths (§6 "Potential direction").
6. Cut options (`margin`, `pad`, `snap`, `max_area_ratio`) aren't in the TOML yet.

Decided: `requires-python` stays `>=3.14` (Nov 2026).
