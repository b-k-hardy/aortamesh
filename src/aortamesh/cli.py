"""Cut and mesh one case, writing intermediate and final results for ParaView.

    uv run aortamesh -c sample_config.toml
    uv run aortamesh -cd examples/with_arch_branches       # all defaults

Writes <case>/preview/: the surface, one clip_<name>.vtp and outline_<name>.vtp per plane, and
the cut wall.vtp. Writes the output folder (default <case>/mesh/): mesh.vtu (tets) and
boundary.vtp (boundary triangles with cell_data["label"]), plus the CHeart mesh <name>.X/.T/.B
(see [cheart] in the config).
"""

import argparse
import dataclasses
import time
from pathlib import Path

from aortamesh.cheart import write_cheart
from aortamesh.clip import (
    RayCaster,
    clipping_surface,
    fit_cross_section,
    section_outline,
)
from aortamesh.config import Config, load_config
from aortamesh.cut import cut_surface
from aortamesh.io import load_case
from aortamesh.mesh import mesh_volume


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Create a non-dissected thoracic aorta fluid mesh")
    p.add_argument("-c", "--config", type=Path, help="TOML config file (see sample_config.toml)")
    p.add_argument(
        "-cd", "--case_dir", type=Path, help="Directory with 3D Slicer exports (overrides config)"
    )
    p.add_argument("-o", "--output", type=Path, help="Output directory (overrides config)")

    args = p.parse_args(argv)

    config = load_config(args.config) if args.config else Config()
    if args.case_dir:
        config = dataclasses.replace(config, directory=args.case_dir)
    if args.output:
        config = dataclasses.replace(config, output=args.output)
    if config.directory is None:
        p.error("give a case directory with -cd or in the config file")

    case_dir = config.directory
    print(case_dir)
    case = load_case(case_dir)
    caster = RayCaster(case.surface)

    out = case_dir / "preview"
    out.mkdir(exist_ok=True)
    case.surface.save(out / "surface.vtp")

    sections = []
    for cut in case.cuts:
        section = fit_cross_section(caster, cut)
        print(section.summary())
        clipping_surface(section).save(out / f"clip_{cut.name}.vtp")
        section_outline(section).save(out / f"outline_{cut.name}.vtp")
        sections.append(section)

    wall, openings = cut_surface(case.surface, sections)
    for o in openings:
        print(o.summary())
    wall.save(out / "wall.vtp")

    print(f"wrote {out}")

    t0 = time.perf_counter()
    vm = mesh_volume(wall, openings, config.mesh, case.centerline)
    print(f"meshed in {time.perf_counter() - t0:.1f} s: {vm.summary()}")

    mesh_out = config.output_dir
    mesh_out.mkdir(parents=True, exist_ok=True)
    vm.to_pyvista().save(mesh_out / "mesh.vtu")
    vm.boundary_pyvista().save(mesh_out / "boundary.vtp")
    write_cheart(vm, mesh_out, config.cheart)
    print(f"wrote {mesh_out} (CHeart: {config.cheart.name}.X/.T/.B, scale {config.cheart.scale})")


if __name__ == "__main__":
    main()
