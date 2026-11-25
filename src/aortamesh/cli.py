"""Fit every cut plane in a case and write the results for inspection in ParaView.

    uv run python examples/preview_cuts.py examples/with_arch_branches

Writes <case>/preview/: one clip_<name>.vtp and outline_<name>.vtp per plane, plus the
surface, so you can load the whole folder and check each section against the wall.
"""

import argparse
import sys
from pathlib import Path

from aortamesh.clip import (
    RayCaster,
    clipping_surface,
    fit_cross_section,
    section_outline,
)
from aortamesh.io import load_case


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Create a non-dissected thoracic aorta fluid mesh")
    p.add_argument(
        "-cd", "--case_dir", type=str, required=True, help="Directory with 3D Slicer exports"
    )

    args = p.parse_args(argv)

    print(args.case_dir)

    case_dir = Path(args.case_dir)
    case = load_case(case_dir)
    caster = RayCaster(case.surface)

    out = case_dir / "preview"
    out.mkdir(exist_ok=True)
    case.surface.save(out / "surface.vtp")

    for cut in case.cuts:
        section = fit_cross_section(caster, cut)
        print(section.summary())
        clipping_surface(section).save(out / f"clip_{cut.name}.vtp")
        section_outline(section).save(out / f"outline_{cut.name}.vtp")

    print(f"wrote {out}")


if __name__ == "__main__":
    main()
