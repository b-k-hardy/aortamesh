"""Cut the closed lumen surface open at each cross-section, keeping the main vessel.

For each section, only the part of the plane inside its clipping polygon (see clip.py) is
used, so the cut follows this vessel's wall and can't touch a neighbouring one:

1. Snap vertices lying almost on the plane onto it, so the cut leaves no sliver triangles.
2. The triangles that cross or touch the plane inside the polygon form a closed strip around
   the vessel. Without the strip, the rest of the surface falls into two pieces: one on
   each side of the cut.
3. Discard the smaller piece by area (a branch tip, or the root below the inlet). A ratio
   near 1 means the plane cuts across the main vessel, so stop rather than guess.
4. Clip the strip triangles exactly at the plane and keep the part on the kept side. The rim
   of the opening then lies in the plane, which the caps rely on.

No caps are added: the boundary layer must be extruded from the open wall before capping.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from aortamesh.clip import CrossSection, clip_polygon
from aortamesh.geometry import plane_basis, polygon_centroid_2d, triangle_areas


@dataclass(frozen=True)
class Opening:
    """One open end of the cut wall surface, and the cap that will close it."""

    name: str  # cap label name, from the plane file
    center: np.ndarray  # (3,) section centroid, on the plane
    normal: np.ndarray  # (3,) unit, outward: points into the discarded piece
    rim: np.ndarray  # (k, 3) boundary loop, in the plane, ordered (first point not repeated)
    area: float  # area enclosed by the rim, mm^2
    discarded_area: float  # surface area removed by this cut, mm^2
    area_ratio: float  # discarded piece / kept piece, by surface area

    def summary(self) -> str:
        return (
            f"{self.name:>8}: rim {len(self.rim):4d} pts  cap area {self.area:7.1f} mm^2  "
            f"discarded {self.discarded_area:8.1f} mm^2  ratio {self.area_ratio:.4f}"
        )


def cut_surface(
    surface: pv.PolyData,
    sections: list[CrossSection],
    margin: float = 0.25,
    pad: float = 1.0,
    snap: float = 0.2,
    max_area_ratio: float = 0.5,
) -> tuple[pv.PolyData, list[Opening]]:
    """Cut ``surface`` open at every section, in order. Returns the open wall and its openings.

    ``margin`` and ``pad`` size the clipping polygons (see ``clip_polygon``). ``snap`` is
    relative to the median edge length: vertices within that distance of a plane are moved
    onto it. Raises if a cut doesn't split the surface in two, or if the smaller piece is
    more than ``max_area_ratio`` of the larger.
    """
    pts = np.asarray(surface.points, dtype=np.float64).copy()
    tris = np.asarray(surface.regular_faces, dtype=np.int64).copy()
    h = float(np.median(_edge_lengths(pts, tris)))

    openings = []
    for section in sections:
        pts, tris, opening = _cut_one(
            pts,
            tris,
            section,
            margin,
            pad,
            snap * h,
            max_area_ratio,
        )
        openings.append(opening)

    wall = pv.PolyData.from_regular_faces(pts, tris)
    if "SPACE" in surface.field_data:
        wall.field_data["SPACE"] = surface.field_data["SPACE"]
    return wall, openings


# =============================================================================
# One cut
# =============================================================================
def _cut_one(
    pts: np.ndarray,
    tris: np.ndarray,
    section: CrossSection,
    margin: float,
    pad: float,
    snap: float,
    max_area_ratio: float,
) -> tuple[np.ndarray, np.ndarray, Opening]:
    name = section.cut.name
    c, n = section.center, section.normal
    u, v = plane_basis(n)
    polygon = clip_polygon(section, margin, pad)
    inside = _in_star_polygon(pts, c, u, v, polygon)

    # 1. Snap near-plane vertices onto the plane (only inside the polygon).
    pts = pts.copy()
    s = (pts - c) @ n
    near = inside & (np.abs(s) < snap)
    pts[near] -= s[near, None] * n
    s[near] = 0.0

    # 2. The strip: triangles inside the polygon that cross or touch the plane.
    pos, neg, on = s > 0, s < 0, s == 0
    crossing = (pos[tris].any(1) & neg[tris].any(1)) | on[tris].any(1)
    strip = crossing & inside[tris].any(1)
    if not strip.any():
        raise RuntimeError(f"{name}: the plane doesn't cross the surface inside its polygon")

    rest = tris[~strip]
    comp = _vertex_components(len(pts), rest)
    rest_comp = comp[rest[:, 0]]

    # The rest triangles touching the strip on each side should be one piece per side.
    used = np.zeros(len(pts), dtype=bool)
    used[rest.ravel()] = True
    sv = np.unique(tris[strip])
    plus = np.unique(comp[sv[pos[sv] & used[sv]]])
    minus = np.unique(comp[sv[neg[sv] & used[sv]]])
    if len(plus) != 1 or len(minus) != 1 or plus[0] == minus[0]:
        raise RuntimeError(
            f"{name}: the cut doesn't split the surface into two pieces "
            f"({len(plus)} piece(s) on one side, {len(minus)} on the other). Does the plane "
            "cross the vessel cleanly, inside its clipping polygon?",
        )

    # 3. Discard the smaller piece.
    rest_area = triangle_areas(pts, rest)
    a_plus = float(rest_area[rest_comp == plus[0]].sum())
    a_minus = float(rest_area[rest_comp == minus[0]].sum())
    keep = 1.0 if a_plus >= a_minus else -1.0
    ratio = min(a_plus, a_minus) / max(a_plus, a_minus)
    if ratio > max_area_ratio:
        raise ValueError(
            f"{name}: the cut leaves pieces of {a_plus:.0f} and {a_minus:.0f} mm^2 "
            f"(ratio {ratio:.2f} > {max_area_ratio}); the plane may cut across the main "
            "vessel rather than off a branch or end",
        )
    discard = minus[0] if keep > 0 else plus[0]

    # 4. Clip the strip at the plane, keeping the kept side.
    clipped, new_pts = _clip_strip(pts, tris[strip], keep * s, c, n)
    new_tris = np.vstack([rest[rest_comp != discard], clipped])
    all_pts = np.vstack([pts, new_pts]) if len(new_pts) else pts
    all_pts, new_tris = _drop_unused(all_pts, new_tris)

    area_before = float(triangle_areas(pts, tris).sum())
    area_after = float(triangle_areas(all_pts, new_tris).sum())

    rim = _find_rim(all_pts, new_tris, c, n, u, v, polygon, snap, name)
    _, rim_area = polygon_centroid_2d(np.c_[(rim - c) @ u, (rim - c) @ v])

    opening = Opening(
        name=name,
        center=c,
        normal=-keep * n,
        rim=rim,
        area=rim_area,
        discarded_area=area_before - area_after,
        area_ratio=ratio,
    )
    return all_pts, new_tris, opening


def _clip_strip(
    pts: np.ndarray,
    strip: np.ndarray,
    s: np.ndarray,
    c: np.ndarray,
    n: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Keep the ``s >= 0`` part of each strip triangle, fan-triangulated.

    New points go on crossed edges and are shared between the two triangles on each edge,
    so the result stays watertight. They are indexed after the existing ``pts``.
    """
    edge_point: dict[tuple[int, int], int] = {}
    new_pts: list[np.ndarray] = []

    def crossing(a: int, b: int) -> int:
        key = (a, b) if a < b else (b, a)
        if key not in edge_point:
            t = s[a] / (s[a] - s[b])
            p = pts[a] + t * (pts[b] - pts[a])
            new_pts.append(p - ((p - c) @ n) * n)  # exactly on the plane
            edge_point[key] = len(pts) + len(new_pts) - 1
        return edge_point[key]

    out = []
    for tri in strip:
        if (s[tri] == 0).all():
            continue  # lies in the plane: that's cap, not wall
        poly = []
        for i in range(3):
            a, b = int(tri[i]), int(tri[(i + 1) % 3])
            if s[a] >= 0:
                poly.append(a)
            if s[a] * s[b] < 0:
                poly.append(crossing(a, b))
        out.extend((poly[0], poly[j], poly[j + 1]) for j in range(1, len(poly) - 1))

    clipped = np.array(out, dtype=np.int64).reshape(-1, 3)
    return clipped, np.array(new_pts).reshape(-1, 3)


def _find_rim(
    pts: np.ndarray,
    tris: np.ndarray,
    c: np.ndarray,
    n: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    polygon: np.ndarray,
    tol: float,
    name: str,
) -> np.ndarray:
    """The boundary loop this cut created: in the plane and inside the clipping polygon."""
    matches = []
    for loop in boundary_loops(tris):
        p = pts[loop]
        if np.abs((p - c) @ n).max() < tol and _in_star_polygon(p, c, u, v, polygon).all():
            matches.append(p)
    if len(matches) != 1:
        raise RuntimeError(f"{name}: expected one new boundary loop, found {len(matches)}")
    return matches[0]


# =============================================================================
# Mesh helpers
# =============================================================================
def boundary_loops(tris: np.ndarray) -> list[np.ndarray]:
    """Ordered vertex loops along the open edges of a consistently oriented triangle mesh."""
    edges = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    key = np.sort(edges, axis=1)
    _, inv, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    open_edges = edges[counts[inv] == 1]

    nxt = dict(zip(open_edges[:, 0].tolist(), open_edges[:, 1].tolist(), strict=True))
    if len(nxt) != len(open_edges):
        raise RuntimeError("non-manifold boundary: a vertex starts two open edges")

    loops = []
    while nxt:
        start, b = nxt.popitem()
        loop = [start]
        while b != start:
            loop.append(b)
            b = nxt.pop(b)
        loops.append(np.array(loop))
    return loops


def _vertex_components(n_pts: int, tris: np.ndarray) -> np.ndarray:
    """Connected-component label of each vertex, linking vertices that share a triangle."""
    a = np.concatenate([tris[:, 0], tris[:, 1]])
    b = np.concatenate([tris[:, 1], tris[:, 2]])
    graph = coo_matrix((np.ones(len(a)), (a, b)), shape=(n_pts, n_pts))
    _, labels = connected_components(graph, directed=False)
    return labels


def _drop_unused(pts: np.ndarray, tris: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    used, inv = np.unique(tris, return_inverse=True)
    return pts[used], inv.reshape(tris.shape)


def _edge_lengths(pts: np.ndarray, tris: np.ndarray) -> np.ndarray:
    p = pts[tris]
    return np.linalg.norm(p[:, [1, 2, 0]] - p, axis=2).ravel()


def _in_star_polygon(
    p: np.ndarray,
    c: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    polygon: np.ndarray,
) -> np.ndarray:
    """Whether each point's projection onto the plane lies inside a polygon star-shaped about c.

    ``polygon`` vertices must be in increasing angle about the normal, as clip.py builds them.
    """
    q = np.c_[(polygon - c) @ u, (polygon - c) @ v]
    xy = np.c_[(p - c) @ u, (p - c) @ v]
    ang_q = np.mod(np.arctan2(q[:, 1], q[:, 0]), 2 * np.pi)
    ang_p = np.mod(np.arctan2(xy[:, 1], xy[:, 0]), 2 * np.pi)

    order = np.argsort(ang_q)
    q, ang_q = q[order], ang_q[order]
    i = np.searchsorted(ang_q, ang_p) - 1  # wedge [i, i+1); -1 wraps to the last wedge
    a, b = q[i % len(q)], q[(i + 1) % len(q)]
    cross = (b[:, 0] - a[:, 0]) * (xy[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (xy[:, 0] - a[:, 0])
    return cross >= 0
