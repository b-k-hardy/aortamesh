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
uv run python examples/preview_cuts.py examples/with_arch_branches
uv run ruff check src examples
uv run pytest                                            # no tests yet
```

`preview_cuts.py` prints one line per plane and writes `<case>/preview/*.vtp` for ParaView.

## Layout

- `src/aortamesh/io.py`: `read_surface`, `read_cut`, `read_cuts`, `load_case` → `Case(surface, cuts)`
- `src/aortamesh/clip.py`: `RayCaster`, `fit_cross_section` → `CrossSection`, `clipping_surface`, `section_outline`
- `examples/with_arch_branches/`: `surface.vtp`, `centerline.vtp`, `inlet/branch1..3/outlet.mrk.json`
- `examples/no_arch_branches/`: `surface.vtp`, `centerline.vtp` (no planes yet)
- `mesh.py` (repo root): scratch copy of gmsh's `aneurysm.py` boundary-layer example; not part of the package

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

## Known traps

- `pv.PolyData.ray_trace(..., first_point=True)` can return the *farthest* hit (its de-dup
  sorts by coordinate). Use `RayCaster`, which sorts by distance along the ray.
- A ray through a shared edge reports the same hit once per cell; `RayCaster` merges these.
  Without that, a duplicate looks like a neighbouring vessel at 0 mm.
- `.mrk.json` files don't include parent transforms from the Slicer scene. Planes and model
  must not sit under a transform when saved.
