"""Small geometry helpers shared across modules."""

from __future__ import annotations

import numpy as np


def plane_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two unit vectors (u, v) spanning the plane, with (u, v, normal) right-handed."""
    n = normal / np.linalg.norm(normal)
    helper = np.eye(3)[np.argmin(np.abs(n))]  # the axis least aligned with n
    u = np.cross(n, helper)
    u /= np.linalg.norm(u)
    return u, np.cross(n, u)


def polygon_centroid_2d(xy: np.ndarray) -> tuple[np.ndarray, float]:
    """Area centroid and (unsigned) area of a simple polygon (vertices in order, not repeated)."""
    x, y = xy[:, 0], xy[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cross = x * yn - xn * y
    area = 0.5 * cross.sum()
    cx = ((x + xn) * cross).sum() / (6 * area)
    cy = ((y + yn) * cross).sum() / (6 * area)
    return np.array([cx, cy]), abs(float(area))


def triangle_areas(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """(k,) areas of the triangles ``tris`` (k, 3) indexing ``points``."""
    a, b, c = (points[tris[:, i]] for i in range(3))
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
