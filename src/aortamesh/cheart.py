"""CHeart mesh files: ``.X`` nodes, ``.T`` linear tets, ``.B`` labeled boundary triangles.

All three are whitespace-separated text with a one-line header, and every index is 1-based:

    .X   N 3               then N rows:  x y z
    .T   M N               then M rows:  n1 n2 n3 n4         (rows of .X)
    .B   K                 then K rows:  tet n1 n2 n3 label  (tet = row of .T, n = rows of .X)

Each ``.B`` triangle is a face of the tet in its first column. Nodes are written in
``VolumeMesh`` order and tets in ``VolumeMesh`` order, so the indices in all three files agree
by construction. This is the only place indices become 1-based.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from aortamesh.volume import VolumeMesh


@dataclass(frozen=True)
class CheartOptions:
    name: str = "mesh"  # files are <output>/<name>.X, .T, .B
    # Multiplies coordinates on output: 1.0 keeps mm, 1e-3 gives m. The density and viscosity
    # entered in CHeart must match (design.md §9).
    scale: float = 1.0

    def __post_init__(self) -> None:
        if not self.name or "/" in self.name:
            raise ValueError(f"invalid CHeart file name: {self.name!r}")
        if self.scale <= 0:
            raise ValueError(f"scale must be positive, got {self.scale}")


def write_cheart(
    vm: VolumeMesh,
    directory: str | Path,
    options: CheartOptions | None = None,
) -> tuple[Path, Path, Path]:
    """Write ``<directory>/<name>.X/.T/.B``. Validates the mesh first. Returns the three paths."""
    options = options or CheartOptions()
    vm.validate()
    stem = Path(directory) / options.name
    x_path, t_path, b_path = (Path(f"{stem}{s}") for s in (".X", ".T", ".B"))

    n, m, k = len(vm.points), len(vm.tets), len(vm.bnd_tris)
    # %.17g round-trips float64 exactly.
    np.savetxt(x_path, vm.points * options.scale, fmt="%.17g", header=f"{n} 3", comments="")
    np.savetxt(t_path, vm.tets + 1, fmt="%d", header=f"{m} {n}", comments="")
    rows = np.c_[vm.bnd_tet + 1, vm.bnd_tris + 1, vm.bnd_labels]
    np.savetxt(b_path, rows, fmt="%d", header=f"{k}", comments="")
    return x_path, t_path, b_path


def read_cheart(
    directory: str | Path,
    name: str = "mesh",
    scale: float = 1.0,
) -> VolumeMesh:
    """Read ``<directory>/<name>.X/.T/.B`` back into a VolumeMesh (0-based; coords / scale).

    Checks each header against its file. Label names aren't stored in CHeart files, so they
    come back as ``"label<k>"``.
    """
    stem = Path(directory) / name
    x_head, x = _read(Path(f"{stem}.X"), float)
    t_head, t = _read(Path(f"{stem}.T"), np.int64)
    b_head, b = _read(Path(f"{stem}.B"), np.int64)
    x, t, b = x.reshape(-1, 3), t.reshape(-1, 4), b.reshape(-1, 5)

    problems = []
    if x_head != [len(x), 3]:
        problems.append(f".X header {x_head}, expected [{len(x)}, 3]")
    if t_head != [len(t), len(x)]:
        problems.append(f".T header {t_head}, expected [{len(t)}, {len(x)}]")
    if b_head != [len(b)]:
        problems.append(f".B header {b_head}, expected [{len(b)}]")
    if t.min() < 1 or t.max() > len(x):
        problems.append(f".T node indices span {t.min()}..{t.max()}, expected 1..{len(x)}")
    if b[:, 0].min() < 1 or b[:, 0].max() > len(t):
        problems.append(
            f".B tet indices span {b[:, 0].min()}..{b[:, 0].max()}, expected 1..{len(t)}"
        )
    if b[:, 1:4].min() < 1 or b[:, 1:4].max() > len(x):
        problems.append(
            f".B node indices span {b[:, 1:4].min()}..{b[:, 1:4].max()}, expected 1..{len(x)}"
        )
    if problems:
        raise ValueError(f"{stem}: " + "; ".join(problems))

    labels = b[:, 4].astype(np.int32)
    return VolumeMesh(
        points=x / scale,
        tets=t - 1,
        bnd_tris=b[:, 1:4] - 1,
        bnd_labels=labels,
        bnd_tet=b[:, 0] - 1,
        label_names={int(lab): f"label{lab}" for lab in np.unique(labels)},
    )


def _read(path: Path, dtype: type) -> tuple[list[int], np.ndarray]:
    with path.open() as f:
        header = [int(v) for v in f.readline().split()]
        data = np.loadtxt(f, dtype=dtype, ndmin=2)
    return header, data
