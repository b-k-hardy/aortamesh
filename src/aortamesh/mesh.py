"""Tetrahedral fluid mesh of the cut wall, with an optional geometric boundary layer (gmsh).

Follows gmsh's ``aneurysm.py`` example:

1. Load the open wall as a discrete surface, then ``classifySurfaces`` + ``createGeometry`` so
   gmsh can reparametrize and remesh it at the requested size.
2. Extrude the boundary layer inward **from the wall only**. If the caps existed already, rim
   nodes would average wall and cap normals and the layer would collapse at every opening.
3. Fill each hole in the inner surface with a plane cap, then the core volume inside.
4. Mesh. Prisms are split into tets (``recombine=False``): CHeart needs linear tets only.

All gmsh state stays inside ``_gmsh_session``. Everything that leaves this module is numpy,
and boundary labels are assigned by geometry (each cap entity goes to the nearest opening),
so no gmsh tags leak out.

Labels are provisional until the MATLAB script's numbering is known: 1 = wall, then 2, 3, ...
for the openings in cut order (inlet, branch1..N, outlet).
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field

import gmsh
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from aortamesh.cut import Opening
from aortamesh.io import Centerline
from aortamesh.volume import VolumeMesh

WALL_LABEL = 1


# =============================================================================
# Options
# =============================================================================
@dataclass(frozen=True)
class BoundaryLayer:
    """Prism layers grown inward from the wall, split into tets.

    Layer ``k`` (0-based) is ``first_height * growth**k`` thick. Defaults are for a rough mesh:
    0.15 + 0.225 + 0.34 = 0.71 mm in total.
    """

    layers: int = 3
    first_height: float = 0.15  # mm
    growth: float = 1.5

    def __post_init__(self) -> None:
        if self.layers < 1 or self.first_height <= 0 or self.growth <= 0:
            raise ValueError(f"invalid boundary layer: {self}")

    @property
    def heights(self) -> np.ndarray:
        """Cumulative distance from the wall to the top of each layer (what gmsh expects)."""
        return self.first_height * np.cumsum(self.growth ** np.arange(self.layers))

    @property
    def thickness(self) -> float:
        return float(self.heights[-1])


@dataclass(frozen=True)
class MeshOptions:
    """Defaults give a rough mesh (~300k tets for an adult thoracic aorta)."""

    size_max: float = 2.0  # mm
    size_min: float = 0.5  # mm
    # Elements around a vessel's circumference: size 2*pi*r / N, with r the radius at the
    # nearest centerline point. Needs a centerline. 0 turns this off.
    radius_elements: int = 16
    # Elements per 2*pi of wall curvature; 0 turns this off. With radius sizing on, it's
    # applied to surfaces only (see _set_radius_sizing). Without it, gmsh also applies it to
    # the seams between the wall's patches, and every seam becomes a band of tiny elements.
    curvature_elements: int = 20
    optimize_netgen: bool = True
    min_quality: float = (
        0.05  # warn below this min SICN (gmsh's quality: 1 = regular, <= 0 = inverted)
    )
    boundary_layer: BoundaryLayer | None = field(default_factory=BoundaryLayer)
    verbose: bool = False  # print gmsh's log

    def __post_init__(self) -> None:
        if not 0 < self.size_min <= self.size_max:
            raise ValueError(f"need 0 < size_min <= size_max, got {self.size_min}, {self.size_max}")


# =============================================================================
# Meshing
# =============================================================================
def mesh_volume(
    wall: pv.PolyData,
    openings: list[Opening],
    options: MeshOptions | None = None,
    centerline: Centerline | None = None,
) -> VolumeMesh:
    """Tet-mesh the region enclosed by ``wall`` and plane caps over its ``openings``.

    ``centerline`` is needed when ``options.radius_elements`` is non-zero.
    """
    options = options or MeshOptions()
    if options.radius_elements and centerline is None:
        raise ValueError(
            "radius_elements needs a centerline (centerline.vtp in the case folder); "
            "or set radius_elements = 0",
        )
    bl = options.boundary_layer
    if bl is not None:
        r_min = min(math.sqrt(o.area / math.pi) for o in openings)
        if bl.thickness > r_min / 3:
            warnings.warn(
                f"boundary layer is {bl.thickness:.2f} mm thick, more than a third of the "
                f"smallest opening radius ({r_min:.2f} mm); layers may fold in tight spots",
                stacklevel=2,
            )

    pts = np.asarray(wall.points, dtype=np.float64)
    tris = np.asarray(wall.regular_faces, dtype=np.int64)

    with _gmsh_session(options.verbose):
        walls = _add_wall(pts, tris)
        inner = _extrude_boundary_layer(walls, bl) if bl is not None else walls
        _add_caps_and_core(inner, len(openings))
        _set_mesh_options(options)
        if options.radius_elements:
            _set_radius_sizing(centerline, options)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.renumberNodes()

        points, tets, tet_tags = _get_tets()
        sicn = gmsh.model.mesh.getElementQualities(tet_tags, "minSICN")
        surfaces = _get_boundary_surfaces()

    if sicn.min() <= 0:
        raise RuntimeError(
            f"{(sicn <= 0).sum()} inverted or flat tets (min SICN {sicn.min():.3f}); the "
            "boundary layer has probably folded: make it thinner or coarsen the wall",
        )
    if sicn.min() < options.min_quality:
        warnings.warn(
            f"{(sicn < options.min_quality).sum()} tets with min SICN < {options.min_quality} "
            f"(worst {sicn.min():.3f})",
            stacklevel=2,
        )

    wall_tags = {tag for _, tag in walls}
    tol = max(2 * (bl.thickness if bl else 0.0), options.size_max)
    labeled, labels = [], []
    for tag, tri in surfaces.items():
        lab = WALL_LABEL if tag in wall_tags else 2 + _nearest_opening(points, tri, openings, tol)
        labeled.append(tri)
        labels.append(np.full(len(tri), lab))
    names = {WALL_LABEL: "wall"} | {2 + i: o.name for i, o in enumerate(openings)}

    vm = VolumeMesh.from_tets(points, tets, np.vstack(labeled), np.concatenate(labels), names)
    vm.validate()

    # Wrong extrusion direction (outward) would add about wall area x thickness of volume.
    expected = enclosed_volume(wall, openings)
    allowed = 0.01 * expected + (0.25 * wall.area * bl.thickness if bl else 0.0)
    if abs(vm.volume - expected) > allowed:
        raise RuntimeError(
            f"mesh volume {vm.volume:.0f} mm^3 vs {expected:.0f} mm^3 enclosed by the wall; "
            "did the boundary layer grow outward?",
        )
    return vm


def enclosed_volume(wall: pv.PolyData, openings: list[Opening]) -> float:
    """Volume inside ``wall`` closed by flat fans from each opening's center to its rim."""
    p = np.asarray(wall.points)
    t = wall.regular_faces
    vol = np.einsum("ij,ij->i", p[t[:, 0]], np.cross(p[t[:, 1]], p[t[:, 2]])).sum()
    for o in openings:
        # The rim runs along the wall's boundary edges, so the cap uses them reversed.
        rim, nxt = o.rim, np.roll(o.rim, -1, axis=0)
        vol += np.einsum("j,ij->i", o.center, np.cross(nxt, rim)).sum()
    return float(vol / 6.0)


# =============================================================================
# gmsh steps (only called inside _gmsh_session)
# =============================================================================
@contextmanager
def _gmsh_session(verbose: bool) -> Iterator[None]:
    if gmsh.isInitialized():
        raise RuntimeError("gmsh is already initialized elsewhere in this process")
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add("aortamesh")
        yield
    finally:
        gmsh.finalize()


def _add_wall(pts: np.ndarray, tris: np.ndarray) -> list[tuple[int, int]]:
    """Load the wall triangulation and turn it into remeshable surface patches."""
    tag = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, tag, np.arange(1, len(pts) + 1), pts.ravel())
    gmsh.model.mesh.addElementsByType(tag, 2, [], (tris + 1).ravel())
    # angle = pi: no splitting at feature edges (the wall is smooth); patches are only split
    # as needed for reparametrization.
    gmsh.model.mesh.classifySurfaces(math.pi, True, True, math.pi)
    gmsh.model.mesh.createGeometry()
    return gmsh.model.getEntities(2)


def _extrude_boundary_layer(
    walls: list[tuple[int, int]],
    bl: BoundaryLayer,
) -> list[tuple[int, int]]:
    """Grow the layers inward; returns the inner ("top") surfaces."""
    gmsh.option.setNumber("Geometry.ExtrudeReturnLateralEntities", 0)
    # One element per layer, cumulative heights. Negative = against the surface normals, which
    # point out of the lumen (Slicer's closed surfaces are outward oriented).
    out = gmsh.model.geo.extrudeBoundaryLayer(
        walls,
        [1] * bl.layers,
        (-bl.heights).tolist(),
        False,
    )
    gmsh.model.geo.synchronize()
    return [e for e in out if e[0] == 2]


def _add_caps_and_core(inner: list[tuple[int, int]], n_openings: int) -> None:
    curves = [c[1] for c in gmsh.model.getBoundary(inner)]
    loops = gmsh.model.geo.addCurveLoops(curves)
    if len(loops) != n_openings:
        raise RuntimeError(f"found {len(loops)} holes to cap, expected {n_openings}")
    caps = [gmsh.model.geo.addPlaneSurface([loop]) for loop in loops]
    shell = gmsh.model.geo.addSurfaceLoop([s for _, s in inner] + caps)
    gmsh.model.geo.addVolume([shell])
    gmsh.model.geo.synchronize()


def _set_mesh_options(o: MeshOptions) -> None:
    gmsh.option.setNumber("Mesh.Algorithm", 1)  # MeshAdapt: robust on reparametrized patches
    gmsh.option.setNumber("Mesh.MeshSizeMin", o.size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", o.size_max)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", o.curvature_elements)
    # Otherwise sizes come from the input triangulation's nodes and spread from the rims, and
    # the remeshed wall keeps the input's edge length whatever size_max says.
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", int(o.optimize_netgen))


def _set_radius_sizing(centerline: Centerline, o: MeshOptions) -> None:
    """Size by the nearest centerline radius, 2*pi*r / radius_elements.

    On curves (patch seams, opening rims) this size *replaces* gmsh's: the curves follow the
    input triangulation's zigzag edges, so curvature sizing would refine every seam into a
    band of tiny elements. On surfaces and in the volume it only lowers gmsh's size, so
    curvature sizing still refines tight spots such as the branch ostia, where the wall is
    nearest the aorta's centerline but bends sharply. gmsh then clamps to
    [size_min, size_max], and holds the callback until ``finalize``.
    """
    tree = cKDTree(centerline.points)
    target = 2 * np.pi * centerline.radius / o.radius_elements

    def size(dim: int, _tag: int, x: float, y: float, z: float, lc: float) -> float:
        _, i = tree.query((x, y, z))
        h = float(target[i])
        return h if dim == 1 else min(lc, h)

    gmsh.model.mesh.setSizeCallback(size)


def _get_tets() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """All nodes (0-based, in tag order) and the linear tets."""
    types = gmsh.model.mesh.getElementTypes(3)
    if list(types) != [4]:
        raise RuntimeError(f"expected only linear tets, got gmsh element types {list(types)}")
    # Every node, whatever entity it is classified on (dim=3 alone would drop surface nodes).
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(node_tags)
    node_tags = node_tags[order].astype(np.int64)
    if not np.array_equal(node_tags, np.arange(1, len(node_tags) + 1)):
        raise RuntimeError("gmsh node tags are not contiguous after renumbering")
    points = coords.reshape(-1, 3)[order]

    tet_tags, tet_nodes = gmsh.model.mesh.getElementsByType(4)
    tets = tet_nodes.reshape(-1, 4).astype(np.int64) - 1
    return points, tets, tet_tags


def _get_boundary_surfaces() -> dict[int, np.ndarray]:
    """Triangles (0-based) of each surface entity on the domain boundary (one volume beside it)."""
    out = {}
    for _, tag in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, tag)
        if len(up) == 1:
            _, nodes = gmsh.model.mesh.getElementsByType(2, tag)
            out[tag] = nodes.reshape(-1, 3).astype(np.int64) - 1
    return out


# =============================================================================
# Labels
# =============================================================================
def _nearest_opening(
    points: np.ndarray,
    tris: np.ndarray,
    openings: list[Opening],
    tol: float,
) -> int:
    """Index of the opening whose cap disk is nearest to this surface's nodes."""
    x = points[np.unique(tris)]
    best, best_d = -1, np.inf
    for i, o in enumerate(openings):
        s = (x - o.center) @ o.normal
        radial = np.linalg.norm((x - o.center) - np.outer(s, o.normal), axis=1)
        r_rim = np.linalg.norm(o.rim - o.center, axis=1).max()
        d = float(np.sqrt(s**2 + np.maximum(radial - r_rim, 0.0) ** 2).max())
        if d < best_d:
            best, best_d = i, d
    if best_d > tol:
        raise RuntimeError(
            f"a cap surface is {best_d:.2f} mm from the nearest opening ({openings[best].name}); "
            f"expected < {tol:.2f} mm",
        )
    return best
