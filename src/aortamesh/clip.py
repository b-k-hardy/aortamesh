"""Ray-cast cross-sections at each cut plane and build local clipping surfaces from them.

For each plane we cast rays outward from the plane center, *within* the plane (perpendicular
to the normal), and keep the first wall hit on each. Those hits trace the lumen boundary in
that plane, wherever the user happened to click. Then:

1. Re-center: the area centroid of the hit polygon becomes the new center, and the rays are
   recast from there, repeated until the center stops moving. The center Slicer gives us is
   only "somewhere inside the lumen".
2. Clipping surface: a star-shaped polygon through the hits, pushed outward by a margin. It
   spans the lumen and slightly beyond it, but (unlike an infinite plane) not the rest of the
   aorta, so a cut at one vessel can't slice through a neighbouring one.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pyvista as pv

from aortamesh.geometry import plane_basis, polygon_centroid_2d
from aortamesh.io import Cut


# =============================================================================
# Ray casting
# =============================================================================
class RayCaster:
    """Ray queries against a fixed surface.

    Uses ``PolyData.intersect_with_line``, which caches a cell locator on the mesh, so
    only the first query pays to build it. Deliberately not
    ``ray_trace(first_point=True)``: its de-duplication sorts hits by coordinate rather
    than by distance, so the "first" point it returns can be the farthest one. Here hits
    are sorted by distance along the ray.
    """

    def __init__(self, surface: pv.PolyData, merge_tol: float = 1e-6) -> None:
        self.surface = surface
        # Bounding-box diagonal: no ray needs to be longer.
        self.max_length = float(surface.length)
        # A ray through a shared edge or vertex reports the same hit once per cell;
        # hits closer than this (relative to the surface size) count as one.
        self._merge = merge_tol * self.max_length

    def distances(self, origin: np.ndarray, direction: np.ndarray) -> np.ndarray:
        """Sorted, de-duplicated distances to every surface crossing along one ray."""
        end = origin + self.max_length * direction
        pts, _ = self.surface.intersect_with_line(origin, end, deduplicate_points=False)
        if pts.size == 0:
            return np.empty(0)
        t = np.sort((np.asarray(pts, dtype=np.float64) - origin) @ direction)
        t = t[t > self._merge]
        if t.size > 1:
            t = t[np.r_[True, np.diff(t) > self._merge]]
        return t

    def hits(self, origin: np.ndarray, directions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Distances to the first and second surface hits along each unit direction.

        From inside the lumen, the first hit is the wall of this vessel and the second
        (if any) is where the ray enters a neighbouring one. NaN where there is no hit.
        """
        first = np.full(len(directions), np.nan)
        second = np.full(len(directions), np.nan)
        for i, d in enumerate(directions):
            t = self.distances(origin, d)
            if t.size:
                first[i] = t[0]
            if t.size > 1:
                second[i] = t[1]
        return first, second


# =============================================================================
# Cross-sections
# =============================================================================
@dataclass(frozen=True)
class CrossSection:
    """The lumen boundary in one cut plane, as seen from a re-centered origin."""

    cut: Cut
    center: np.ndarray  # (3,) area centroid of the hit polygon (refined)
    normal: np.ndarray  # (3,) unit normal (same line as cut.normal; sign not yet fixed)
    directions: np.ndarray  # (n, 3) in-plane unit ray directions
    radii: np.ndarray  # (n,) distance from center to the wall along each ray
    clearance: np.ndarray  # (n,) distance to the next surface beyond the wall (inf if none)
    area: float  # area of the hit polygon, mm^2
    shift: float  # how far the center moved from the user's click, mm
    iterations: int

    @property
    def hits(self) -> np.ndarray:
        """(n, 3) wall points, in ray order around the normal."""
        return self.center + self.radii[:, None] * self.directions

    @property
    def equivalent_radius(self) -> float:
        """Radius of the circle with the same area."""
        return float(np.sqrt(self.area / np.pi))

    def summary(self) -> str:
        r = self.radii
        gap = self.clearance.min()
        gap_s = "none in plane" if np.isinf(gap) else f"{gap:.2f} mm"
        return (
            f"{self.cut.name:>8}: r_eq {self.equivalent_radius:5.2f} mm  "
            f"r [{r.min():5.2f}, {r.max():5.2f}]  area {self.area:7.1f} mm^2  "
            f"center moved {self.shift:5.2f} mm in {self.iterations} it  "
            f"nearest other vessel: {gap_s}"
        )


def fit_cross_section(
    caster: RayCaster,
    cut: Cut,
    n_rays: int = 72,
    max_iter: int = 20,
    tol: float = 1e-3,
    outlier_factor: float = 2.0,
) -> CrossSection:
    """Cast ``n_rays`` in-plane rays from the cut center and iteratively re-center them.

    ``tol`` is relative to the equivalent radius. Raises if any ray misses the wall (the
    center is outside the lumen or the plane leaves the vessel). Warns if a ray is more than
    ``outlier_factor`` times the median radius, which usually means it escaped through a branch
    ostium or the plane is badly oblique: move the plane rather than trusting the section.
    """
    n = cut.normal / np.linalg.norm(cut.normal)
    u, v = plane_basis(n)
    theta = np.linspace(0.0, 2 * np.pi, n_rays, endpoint=False)
    dirs = np.cos(theta)[:, None] * u + np.sin(theta)[:, None] * v

    center = cut.center.copy()
    for it in range(1, max_iter + 1):
        radii, _ = caster.hits(center, dirs)
        if np.isnan(radii).any():
            raise RuntimeError(
                f"{cut.name}: {np.isnan(radii).sum()}/{n_rays} rays missed the wall; "
                "is the plane center inside the lumen?",
            )
        xy = radii[:, None] * np.c_[np.cos(theta), np.sin(theta)]  # hits in (u, v) coords
        c2, area = polygon_centroid_2d(xy)
        step = float(np.linalg.norm(c2))
        center = center + c2[0] * u + c2[1] * v
        if step < tol * np.sqrt(area / np.pi):
            break
    else:
        warnings.warn(
            f"{cut.name}: re-centering did not converge in {max_iter} iterations",
        )

    # Final rays from the converged center, so radii and center agree exactly.
    radii, beyond = caster.hits(center, dirs)
    clearance = np.where(np.isnan(beyond), np.inf, beyond - radii)
    xy = radii[:, None] * np.c_[np.cos(theta), np.sin(theta)]
    _, area = polygon_centroid_2d(xy)

    med = np.median(radii)
    if (radii > outlier_factor * med).any():
        k = int((radii > outlier_factor * med).sum())
        warnings.warn(
            f"{cut.name}: {k} ray(s) exceed {outlier_factor}x the median radius "
            f"({radii.max():.1f} vs {med:.1f} mm); the plane may cut into a branch ostium",
        )

    return CrossSection(
        cut=cut,
        center=center,
        normal=n,
        directions=dirs,
        radii=radii,
        clearance=clearance,
        area=float(area),
        shift=float(np.linalg.norm(center - cut.center)),
        iterations=it,
    )


# =============================================================================
# Clipping surfaces
# =============================================================================
def clip_polygon(
    section: CrossSection,
    margin: float = 0.25,
    pad: float = 1.0,
) -> np.ndarray:
    """(n, 3) vertices of the star-shaped clipping polygon, one per ray, in ray order.

    Each vertex sits at ``r * (1 + margin) + pad`` along its ray (a relative margin plus an
    absolute one in mm), but never more than halfway across the gap to the next surface the
    ray meets. So it always reaches past this vessel's wall and never into a neighbour's.
    """
    r = section.radii * (1.0 + margin) + pad
    r = np.minimum(r, section.radii + 0.5 * section.clearance)
    return section.center + r[:, None] * section.directions


def clipping_surface(
    section: CrossSection,
    margin: float = 0.25,
    pad: float = 1.0,
) -> pv.PolyData:
    """The clipping polygon as a surface, for ParaView.

    Triangulated as a fan from the center, which is valid because the polygon is star-shaped
    about it by construction.
    """
    rim = clip_polygon(section, margin, pad)
    pts = np.vstack([section.center, rim])
    n = len(rim)
    i = np.arange(n)
    faces = np.c_[np.full(n, 3), np.zeros(n, int), 1 + i, 1 + (i + 1) % n].ravel()
    disc = pv.PolyData(pts, faces)
    disc.field_data["name"] = [section.cut.name]
    return disc


def section_outline(section: CrossSection) -> pv.PolyData:
    """Closed polyline through the wall hits, for checking a section in ParaView."""
    hits = section.hits
    n = len(hits)
    line = pv.PolyData(hits, lines=np.r_[n + 1, np.arange(n), 0])
    line.field_data["name"] = [section.cut.name]
    return line
