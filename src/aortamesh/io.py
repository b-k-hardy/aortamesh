"""Readers for Slicer exports: the lumen surface (.vtp/.vtk/.stl) and cut planes (.mrk.json).

Everything is returned in one coordinate system (LPS by default, which is what Slicer writes
to disk) and in float64, so downstream geometry never has to think about either.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pyvista as pv

# RAS <-> LPS is a sign flip on x and y; the transform is its own inverse.
_RAS_LPS = np.array([-1.0, -1.0, 1.0])
_SYSTEMS = ("LPS", "RAS")


def _convert(xyz: np.ndarray, src: str, dst: str) -> np.ndarray:
    src, dst = src.upper(), dst.upper()
    if src not in _SYSTEMS or dst not in _SYSTEMS:
        raise ValueError(f"unknown coordinate system: {src!r} -> {dst!r}")
    return xyz if src == dst else xyz * _RAS_LPS


# =============================================================================
# Surface
# =============================================================================
def read_surface(path: str | Path, coordinate_system: str = "LPS") -> pv.PolyData:
    """Read a closed lumen surface and return it as a clean, all-triangle float64 PolyData.

    The source coordinate system is taken from the ``SPACE`` field array that Slicer writes
    into .vtp files. Files without it (e.g. STL) are assumed to be LPS, Slicer's default.
    """
    path = Path(path)
    surf = pv.read(path)
    if not isinstance(surf, pv.PolyData):
        surf = surf.extract_surface()

    src = "LPS"
    if "SPACE" in surf.field_data:
        src = str(np.asarray(surf.field_data["SPACE"]).ravel()[0])

    surf = surf.triangulate().clean()
    surf.points = _convert(
        np.asarray(surf.points, dtype=np.float64), src, coordinate_system
    )
    # Slicer's normals are float32 and would be stale after any flip; recompute on demand.
    surf.point_data.pop("Normals", None)
    surf.field_data["SPACE"] = [coordinate_system.upper()]

    if surf.n_open_edges:
        raise ValueError(
            f"{path.name}: surface has {surf.n_open_edges} open edges; expected a closed "
            "segmentation surface (cuts are made later, by aortamesh)."
        )
    return surf


# =============================================================================
# Cut planes
# =============================================================================
@dataclass(frozen=True)
class Cut:
    """One Slicer plane markup: where to cut the vessel and what to call the cap."""

    name: str  # from the file stem: "inlet", "branch1", ..., "outlet"
    center: (
        np.ndarray
    )  # (3,) a point on the plane, inside the lumen (not necessarily centered)
    normal: np.ndarray  # (3,) unit normal; its sign is arbitrary and is fixed later
    source: Path | None = field(default=None, compare=False)


def read_cut(path: str | Path, coordinate_system: str = "LPS") -> Cut:
    """Read a single-plane ``.mrk.json`` file written by Slicer."""
    path = Path(path)
    with path.open() as f:
        doc = json.load(f)

    markups = doc.get("markups", [])
    if len(markups) != 1 or markups[0].get("type") != "Plane":
        raise ValueError(f"{path.name}: expected exactly one Plane markup")
    m = markups[0]

    units = m.get("coordinateUnits", "mm")
    if units != "mm":
        raise ValueError(f"{path.name}: coordinateUnits is {units!r}, expected 'mm'")

    src = m.get("coordinateSystem", "LPS")
    center = _convert(np.asarray(m["center"], dtype=np.float64), src, coordinate_system)
    normal = _convert(np.asarray(m["normal"], dtype=np.float64), src, coordinate_system)
    normal /= np.linalg.norm(normal)

    name = path.name.removesuffix(".mrk.json")
    return Cut(name=name, center=center, normal=normal, source=path)


def _cut_order(name: str) -> tuple[int, int, str]:
    """inlet first, then branch1, branch2, ... in flow order, then outlet, then anything else."""
    if name == "inlet":
        return (0, 0, name)
    if m := re.fullmatch(r"branch(\d+)", name):
        return (1, int(m.group(1)), name)
    if name == "outlet":
        return (2, 0, name)
    return (3, 0, name)


def read_cuts(directory: str | Path, coordinate_system: str = "LPS") -> list[Cut]:
    """Read every ``*.mrk.json`` in a folder, ordered inlet, branch1..N, outlet."""
    paths = sorted(Path(directory).glob("*.mrk.json"))
    if not paths:
        raise FileNotFoundError(f"no .mrk.json files in {directory}")
    cuts = [read_cut(p, coordinate_system) for p in paths]
    return sorted(cuts, key=lambda c: _cut_order(c.name))


# =============================================================================
# A whole case
# =============================================================================
@dataclass(frozen=True)
class Case:
    surface: pv.PolyData
    cuts: list[Cut]
    coordinate_system: str = "LPS"


def load_case(
    directory: str | Path,
    surface_name: str = "surface.vtp",
    coordinate_system: str = "LPS",
) -> Case:
    """Load ``<directory>/surface.vtp`` and all ``*.mrk.json`` planes beside it."""
    directory = Path(directory)
    surface = read_surface(directory / surface_name, coordinate_system)
    cuts = read_cuts(directory, coordinate_system)

    # Every plane center must sit inside the lumen, or rays cast from it are meaningless.
    probe = pv.PolyData(np.array([c.center for c in cuts]))
    inside = probe.select_interior_points(surface)["selected_points"].astype(bool)
    if not inside.all():
        bad = [c.name for c, ok in zip(cuts, inside) if not ok]
        raise ValueError(f"plane centers outside the lumen: {bad}")

    return Case(surface=surface, cuts=cuts, coordinate_system=coordinate_system)
