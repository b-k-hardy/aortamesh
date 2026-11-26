"""The neutral volume mesh: linear tets plus labeled boundary triangles. Every writer reads this."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyvista as pv

# Faces of a positively oriented tet (v0..v3), each wound so its normal points out of the tet.
_TET_FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])


@dataclass(frozen=True)
class VolumeMesh:
    points: np.ndarray  # (N, 3) float64, LPS, mm
    tets: np.ndarray  # (M, 4) int64, 0-based, positively oriented
    bnd_tris: np.ndarray  # (K, 3) int64, wound with the normal pointing out of the domain
    bnd_labels: np.ndarray  # (K,) int32
    bnd_tet: np.ndarray  # (K,) int64, parent tet of each boundary triangle (needed for .B)
    label_names: dict[int, str]

    @classmethod
    def from_tets(
        cls,
        points: np.ndarray,
        tets: np.ndarray,
        labeled_tris: np.ndarray,
        tri_labels: np.ndarray,
        label_names: dict[int, str],
    ) -> VolumeMesh:
        """Find the boundary of ``tets`` and label it from ``labeled_tris`` (any winding).

        Boundary faces not among ``labeled_tris`` get label 0, which ``validate`` rejects.
        """
        points = np.asarray(points, dtype=np.float64)
        tets = np.asarray(tets, dtype=np.int64).copy()
        neg = _signed_volumes(points, tets) < 0
        tets[neg] = tets[neg][:, [0, 1, 3, 2]]

        faces = tets[:, _TET_FACES].reshape(-1, 3)
        key = np.sort(faces, axis=1)
        _, inv, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        on_bnd = counts[inv] == 1

        lookup = {
            tuple(t): int(lab)
            for t, lab in zip(np.sort(labeled_tris, axis=1).tolist(), tri_labels, strict=True)
        }
        labels = np.array([lookup.get(tuple(k), 0) for k in key[on_bnd].tolist()], dtype=np.int32)

        return cls(
            points=points,
            tets=tets,
            bnd_tris=faces[on_bnd],
            bnd_labels=labels,
            bnd_tet=np.flatnonzero(on_bnd) // 4,
            label_names=dict(label_names),
        )

    # -------------------------------------------------------------------------
    def validate(self) -> None:
        """Raise ValueError listing every structural problem found."""
        problems = []
        vol = _signed_volumes(self.points, self.tets)
        if (vol <= 0).any():
            problems.append(f"{(vol <= 0).sum()} tets with non-positive volume")

        used = np.zeros(len(self.points), dtype=bool)
        used[self.tets.ravel()] = True
        if not used.all():
            problems.append(f"{(~used).sum()} nodes not used by any tet")

        untagged = int((self.bnd_labels == 0).sum())
        if untagged:
            problems.append(f"{untagged} untagged boundary triangles")
        unknown = set(np.unique(self.bnd_labels).tolist()) - set(self.label_names) - {0}
        if unknown:
            problems.append(f"boundary labels without names: {sorted(unknown)}")

        parent = self.tets[self.bnd_tet]
        in_parent = (self.bnd_tris[:, :, None] == parent[:, None, :]).any(2).all(1)
        if not in_parent.all():
            problems.append(f"{(~in_parent).sum()} boundary triangles not a face of their tet")

        t = self.bnd_tris
        edges = np.sort(np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        if (counts != 2).any():
            problems.append(f"boundary not closed: {(counts != 2).sum()} edges not shared by 2")

        if problems:
            raise ValueError("invalid VolumeMesh: " + "; ".join(problems))

    # -------------------------------------------------------------------------
    @property
    def volume(self) -> float:
        return float(_signed_volumes(self.points, self.tets).sum())

    def label_areas(self) -> dict[str, float]:
        a = _triangle_areas(self.points, self.bnd_tris)
        return {
            name: float(a[self.bnd_labels == lab].sum()) for lab, name in self.label_names.items()
        }

    def summary(self) -> str:
        lines = [
            (
                f"{len(self.points)} nodes, {len(self.tets)} tets, {len(self.bnd_tris)} "
                f"boundary triangles, volume {self.volume:.0f} mm^3"
            ),
        ]
        areas = self.label_areas()
        for lab, name in self.label_names.items():
            n = int((self.bnd_labels == lab).sum())
            lines.append(f"  label {lab:2d} {name:>8}: {n:6d} tris  {areas[name]:9.1f} mm^2")
        return "\n".join(lines)

    def to_pyvista(self) -> pv.UnstructuredGrid:
        cells = np.c_[np.full(len(self.tets), 4), self.tets].ravel()
        types = np.full(len(self.tets), pv.CellType.TETRA, dtype=np.uint8)
        return pv.UnstructuredGrid(cells, types, self.points)

    def boundary_pyvista(self) -> pv.PolyData:
        """Boundary triangles with ``cell_data["label"]``, for checking labels in ParaView."""
        surf = pv.PolyData.from_regular_faces(self.points, self.bnd_tris)
        surf.cell_data["label"] = self.bnd_labels
        surf.field_data["label_names"] = [
            f"{lab}:{name}" for lab, name in self.label_names.items()
        ]
        return surf.clean(point_merging=False)  # drops interior points


def _signed_volumes(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    a, b, c, d = (points[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def _triangle_areas(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    a, b, c = (points[tris[:, i]] for i in range(3))
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
