# aortamesh

Python package that turns a 3D Slicer aortic segmentation plus cut planes into a labeled
tetrahedral fluid mesh for CFD (CHeart `.X/.T/.B`, plus XDMF/VTU). It replaces the lab's
Slicer → STL → SimModeler → Abaqus 2D/3D → MATLAB relabeling workflow.

**Read `docs/design.md` before starting work.** It holds the design decisions, what's done,
what's next, and the open questions.

## How we work

- Brandon designs and implements the structure; Claude fills in details when asked. Propose
  before making large structural changes.
- Verify geometry, don't just run it: check results against the surface (distances, areas,
  label coverage) and say what was checked. A run that doesn't crash isn't a pass.
- Brandon edits files himself. Re-read a file before editing it, and never overwrite his
  changes. Match the surrounding formatting; don't run `ruff format` over files he has
  edited (his style keeps trailing commas). `ruff check` is fine.

## Setup and commands

```bash
uv sync                                                  # Python 3.14 (.python-version)
uv run aortamesh -c sample_config.toml                   # ~25 s on the example case
uv run aortamesh -cd examples/with_arch_branches       # defaults; or: uv run python -m aortamesh
uv run ruff check src
uv run pytest                                            # no tests yet
```

The CLI (`cli.py`) prints one line per plane and per opening, then a mesh summary. It writes
`<case>/preview/*.vtp` (surface, clip polygons, outlines, cut wall) and `<output>/mesh.vtu` +
`boundary.vtp` + CHeart `mesh.X/.T/.B` (default output `<case>/mesh/`). `-cd`/`-o` override the config. Lint config is `ruff.toml`; `.vscode/settings.json`
makes the VS Code extension use it over user settings.

## Layout

- `src/aortamesh/io.py`: `read_surface`, `read_cut`, `read_cuts`, `read_centerline` → `Centerline(points, radius)`, `load_case` → `Case(surface, cuts, centerline)` (centerline is None if there's no `centerline.vtp`)
- `src/aortamesh/clip.py`: `RayCaster`, `fit_cross_section` → `CrossSection`, `clip_polygon` (array), `clipping_surface` / `section_outline` (PolyData for ParaView)
- `src/aortamesh/geometry.py`: shared helpers `plane_basis`, `polygon_centroid_2d`, `triangle_areas`
- `src/aortamesh/cut.py`: `cut_surface(surface, sections)` → open wall `pv.PolyData` + `list[Opening]`; `boundary_loops`
- `src/aortamesh/cli.py` (+ `__main__.py`): `main(argv=None)`, entry point `aortamesh`
- `src/aortamesh/mesh.py`: `mesh_volume(wall, openings, MeshOptions, centerline)` → `VolumeMesh`; `BoundaryLayer`, `enclosed_volume`
- `src/aortamesh/volume.py`: `VolumeMesh` (`from_tets`, `validate`, `summary`, `to_pyvista`, `boundary_pyvista`)
- `src/aortamesh/cheart.py`: `write_cheart(vm, dir, CheartOptions(name, scale))` → `<name>.X/.T/.B`; `read_cheart` (round trip, header checks)
- `src/aortamesh/config.py`: `load_config(toml)` → `Config(directory, output, mesh, cheart)`; unknown keys are errors
- `sample_config.toml`: every option with its default
- `examples/with_arch_branches/`: `surface.vtp`, `centerline.vtp`, `inlet/branch1..3/outlet.mrk.json`
- `examples/no_arch_branches/`: `surface.vtp`, `centerline.vtp` (no planes yet)

## Conventions

- Surfaces are `pv.PolyData` with `cell_data["label"]`; no custom surface class.
- numpy arrays, float64, 0-based indices everywhere. Convert to 1-based only inside the CHeart writer.
- Coordinates in LPS, lengths in mm, throughout the package. Convert only at the readers and writers.
- Plain functions plus small frozen dataclasses. No abstract base classes.
- No direct `vtk` imports; use pyvista. (pyvista still depends on vtk.)
- gmsh state stays inside one context manager; convert to numpy before returning, so gmsh tags never leak out.

## Input conventions (a case folder)

- `surface.vtp`: closed lumen surface exported from Slicer. The `SPACE` field gives the
  coordinate system (LPS). The readers reject surfaces with open edges.
- One `.mrk.json` plane per cut, named `inlet`, `branch1`…`branchN` (in flow order along the
  arch), `outlet`. The file name becomes the cap label name.
- Plane centers must be inside the lumen (`load_case` checks this) but needn't be centered.
- Plane normal sign is arbitrary; the code decides which side to discard.
- `centerline.vtp`: from Slicer's Extract Centerline, with a `Radius` point array. Optional,
  but radius sizing (`radius_elements > 0`, the default) needs it.

## Known traps

- `pv.PolyData.ray_trace(..., first_point=True)` can return the *farthest* hit (its de-dup
  sorts by coordinate). Use `RayCaster`, which sorts by distance along the ray.
- A ray through a shared edge reports the same hit once per cell; `RayCaster` merges these.
  Without that, a duplicate looks like a neighbouring vessel at 0 mm.
- gmsh sizes the remeshed wall from the input triangulation unless `Mesh.MeshSizeFromPoints`
  and `Mesh.MeshSizeExtendFromBoundary` are 0 (`_set_mesh_options` sets both). Otherwise
  `size_max` is silently ignored.
- Don't use curvature sizing on its own: gmsh applies it to the patch seams too, which
  zigzag along the input edges and become bands of tiny elements. And don't use radius
  sizing on its own: it's too coarse at the branch ostia and produces flat tets. The size
  callback in `_set_radius_sizing` combines them (radius on curves, the smaller of the two
  elsewhere).
- `.mrk.json` files don't include parent transforms from the Slicer scene. Planes and model
  must not sit under a transform when saved.
