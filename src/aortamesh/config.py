"""Run configuration from a TOML file. Missing keys take the defaults in the option classes."""

from __future__ import annotations

import dataclasses
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from aortamesh.cheart import CheartOptions
from aortamesh.mesh import BoundaryLayer, MeshOptions


@dataclass(frozen=True)
class Config:
    directory: Path | None = None  # case folder: surface.vtp and *.mrk.json
    output: Path | None = None  # where mesh files go; default <directory>/mesh
    mesh: MeshOptions = field(default_factory=MeshOptions)
    cheart: CheartOptions = field(default_factory=CheartOptions)

    @property
    def output_dir(self) -> Path:
        if self.output is not None:
            return self.output
        if self.directory is None:
            raise ValueError("no case directory given")
        return self.directory / "mesh"


def load_config(path: str | Path) -> Config:
    """Read a config file. Relative paths in it are relative to the file, not the shell."""
    path = Path(path)
    with path.open("rb") as f:
        doc = tomllib.load(f)
    base = path.parent

    mesh_tbl = dict(doc.pop("mesh", {}))
    bl_tbl = dict(doc.pop("boundary_layer", {}))
    cheart_tbl = dict(doc.pop("cheart", {}))
    top = _take(doc, Config, "top level", exclude={"mesh", "cheart"})

    bl = None
    if bl_tbl.pop("enabled", True):
        bl = BoundaryLayer(**_take(bl_tbl, BoundaryLayer, "[boundary_layer]"))
    mesh = MeshOptions(
        **_take(mesh_tbl, MeshOptions, "[mesh]", exclude={"boundary_layer"}),
        boundary_layer=bl,
    )

    def resolve(p: str | None) -> Path | None:
        return base / p if p else None  # "" means unset

    return Config(
        directory=resolve(top.get("directory")),
        output=resolve(top.get("output")),
        mesh=mesh,
        cheart=CheartOptions(**_take(cheart_tbl, CheartOptions, "[cheart]")),
    )


def _take(
    table: dict[str, Any],
    cls: type,
    where: str,
    exclude: frozenset[str] | set[str] = frozenset(),
) -> dict[str, Any]:
    """The entries of ``table`` that are fields of ``cls``; anything else is an error."""
    allowed = {f.name for f in dataclasses.fields(cls)} - set(exclude)
    unknown = set(table) - allowed
    if unknown:
        raise ValueError(
            f"unknown key(s) {sorted(unknown)} in {where}; expected some of {sorted(allowed)}",
        )
    return dict(table)
