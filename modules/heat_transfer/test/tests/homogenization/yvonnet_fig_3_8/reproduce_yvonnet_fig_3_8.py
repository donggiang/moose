#!/usr/bin/env python3
"""Reproduce the conductivity sweep in Yvonnet Fig. 3.8 with MOOSE."""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


K_MATRIX = 1.0
K_FIBER = 5.0
FRACTIONS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# Table 3.1 values from Yvonnet, Computational Homogenization of Heterogeneous Materials.
YVONNET_FINE_PER = {
    0.0: 1.0000,
    0.05: 1.0672,
    0.1: 1.1419,
    0.2: 1.3055,
    0.3: 1.4989,
    0.4: 1.7279,
    0.5: 2.0104,
    0.6: 2.3704,
    0.7: 2.8643,
}
YVONNET_FINE_UTG = {
    0.0: 1.0000,
    0.05: 1.0674,
    0.1: 1.1428,
    0.2: 1.3099,
    0.3: 1.5103,
    0.4: 1.7508,
    0.5: 2.0500,
    0.6: 2.4312,
    0.7: 2.9486,
}


@dataclass(frozen=True)
class MooseResult:
    fraction: float
    mode: str
    k11: float
    actual_fraction: float


def parse_fractions(raw: list[str]) -> list[float]:
    if not raw:
        return list(FRACTIONS)

    values: list[float] = []
    for item in raw:
        for token in item.replace(",", " ").split():
            values.append(float(token))
    return values


def polygon_radius_for_area(area: float, n_segments: int) -> float:
    polygon_factor = 0.5 * n_segments * math.sin(2.0 * math.pi / n_segments)
    return math.sqrt(area / polygon_factor)


def point_list(points: list[tuple[float, float]]) -> str:
    return "\n              ".join(f"{x:.14g} {y:.14g} 0" for x, y in points)


def circle_points(fraction: float, n_segments: int) -> list[tuple[float, float]]:
    radius = polygon_radius_for_area(fraction, n_segments)
    return [
        (
            0.5 + radius * math.cos(2.0 * math.pi * i / n_segments),
            0.5 + radius * math.sin(2.0 * math.pi * i / n_segments),
        )
        for i in range(n_segments)
    ]


def mesh_block(fraction: float, args: argparse.Namespace) -> str:
    square = point_list([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    circle = point_list(circle_points(fraction, args.circle_segments))
    return f"""
[Mesh]
  [square_boundary]
    type = PolyLineMeshGenerator
    points = '{square}'
    loop = true
    nums_edges_between_points = {args.boundary_edges}
  []
  [fiber_boundary]
    type = PolyLineMeshGenerator
    points = '{circle}'
    loop = true
    nums_edges_between_points = 1
  []
  [fiber]
    type = XYDelaunayGenerator
    boundary = fiber_boundary
    refine_boundary = false
    desired_area = {args.desired_area:g}
    tri_element_type = TRI6
    interior_points = '0.5 0.5 0'
    output_subdomain_id = 2
    output_subdomain_name = fiber
    output_boundary = fiber_boundary
  []
  [cell]
    type = XYDelaunayGenerator
    boundary = square_boundary
    holes = fiber
    stitch_holes = true
    refine_holes = false
    verify_holes = false
    refine_boundary = false
    desired_area = {args.desired_area:g}
    tri_element_type = TRI6
    output_subdomain_id = 1
    output_subdomain_name = matrix
    output_boundary = outer
  []
  [left]
    type = ParsedGenerateSideset
    input = cell
    new_sideset_name = left
    included_boundaries = outer
    combinatorial_geometry = 'x < 1e-10'
  []
  [right]
    type = ParsedGenerateSideset
    input = left
    new_sideset_name = right
    included_boundaries = outer
    combinatorial_geometry = 'x > 1 - 1e-10'
  []
  [bottom]
    type = ParsedGenerateSideset
    input = right
    new_sideset_name = bottom
    included_boundaries = outer
    combinatorial_geometry = 'y < 1e-10'
  []
  [top]
    type = ParsedGenerateSideset
    input = bottom
    new_sideset_name = top
    included_boundaries = outer
    combinatorial_geometry = 'y > 1 - 1e-10'
  []
  [pin]
    type = ExtraNodesetGenerator
    input = top
    new_boundary = pin
    coord = '0.5 0.5 0'
    use_closest_node = true
  []
[]
"""


def material_block() -> str:
    return f"""
[Materials]
  [matrix]
    type = HeatConductionMaterial
    block = matrix
    specific_heat = 1
    thermal_conductivity = {K_MATRIX:g}
  []
  [fiber]
    type = HeatConductionMaterial
    block = fiber
    specific_heat = 1
    thermal_conductivity = {K_FIBER:g}
  []
[]
"""


def common_postprocessors() -> str:
    return """
  [total_area]
    type = VolumePostprocessor
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [fiber_area]
    type = VolumePostprocessor
    block = fiber
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [actual_fraction]
    type = ParsedPostprocessor
    pp_names = 'fiber_area total_area'
    pp_symbols = 'Af A'
    expression = 'Af/A'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
"""


def executioner_block() -> str:
    return """
[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-10
  l_max_its = 50
[]
"""


def output_block(file_base: str) -> str:
    return f"""
[Outputs]
  exodus = false
  file_base = '{file_base}'
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
"""


def periodic_input(fraction: float, file_base: str, args: argparse.Namespace) -> str:
    return f"""{mesh_block(fraction, args)}
[Variables]
  [chi_x]
    order = FIRST
    family = LAGRANGE
  []
  [chi_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [heat_x]
    type = HeatConduction
    variable = chi_x
  []
  [heat_y]
    type = HeatConduction
    variable = chi_y
  []
  [rhs_x]
    type = HomogenizedHeatConduction
    variable = chi_x
    component = 0
  []
  [rhs_y]
    type = HomogenizedHeatConduction
    variable = chi_y
    component = 1
  []
[]

[BCs]
  [Periodic]
    [left_right]
      primary = left
      secondary = right
      translation = '1 0 0'
    []
    [bottom_top]
      primary = bottom
      secondary = top
      translation = '0 1 0'
    []
  []
  [pin_x]
    type = DirichletBC
    variable = chi_x
    boundary = pin
    value = 0
  []
  [pin_y]
    type = DirichletBC
    variable = chi_y
    boundary = pin
    value = 0
  []
[]
{material_block()}
{executioner_block()}
[Postprocessors]
{common_postprocessors()}
  [k11]
    type = HomogenizedThermalConductivity
    chi = 'chi_x chi_y'
    row = 0
    col = 0
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
{output_block(file_base)}
"""


def utg_input(fraction: float, file_base: str, args: argparse.Namespace) -> str:
    return f"""{mesh_block(fraction, args)}
[Variables]
  [temp]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [flux_x]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [linear_x]
    type = ParsedFunction
    expression = 'x'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = temp
  []
[]

[AuxKernels]
  [flux_x]
    type = DiffusionFluxAux
    variable = flux_x
    diffusion_variable = temp
    diffusivity = thermal_conductivity
    component = x
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]

[BCs]
  [linear_temperature]
    type = FunctionDirichletBC
    variable = temp
    boundary = 'left right bottom top'
    function = linear_x
  []
[]
{material_block()}
{executioner_block()}
[Postprocessors]
{common_postprocessors()}
  [average_flux_x]
    type = ElementAverageValue
    variable = flux_x
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [k11]
    type = ParsedPostprocessor
    pp_names = 'average_flux_x'
    pp_symbols = 'qx'
    expression = '-qx'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
{output_block(file_base)}
"""


def run_moose(exe: Path, input_file: Path) -> None:
    cmd = [str(exe), "-i", input_file.name, "--color", "off"]
    result = subprocess.run(
        cmd,
        cwd=input_file.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        tail = result.stdout[-5000:]
        raise RuntimeError(
            f"MOOSE failed for {input_file.name} with exit code {result.returncode}\n{tail}"
        )


def read_csv_value(csv_file: Path) -> tuple[float, float]:
    with csv_file.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No data rows found in {csv_file}")

    row = rows[-1]
    return float(row["k11"]), float(row["actual_fraction"])


def analytic_values(fraction: float) -> tuple[float, float, float]:
    voigt = (1.0 - fraction) * K_MATRIX + fraction * K_FIBER
    reuss = 1.0 / ((1.0 - fraction) / K_MATRIX + fraction / K_FIBER)
    maxwell = K_MATRIX + (
        2.0 * K_MATRIX * fraction * (K_FIBER - K_MATRIX)
        / (2.0 * K_MATRIX + (1.0 - fraction) * (K_FIBER - K_MATRIX))
    )
    return voigt, reuss, maxwell


def write_results(results: list[MooseResult], out_csv: Path) -> None:
    rows = []
    result_map = {(r.fraction, r.mode): r for r in results}
    for fraction in sorted({r.fraction for r in results}):
        voigt, reuss, maxwell = analytic_values(fraction)
        rows.append(
            {
                "target_fraction": fraction,
                "actual_fraction_per": result_map.get((fraction, "per"), None).actual_fraction
                if (fraction, "per") in result_map
                else "",
                "k11_moose_per": result_map.get((fraction, "per"), None).k11
                if (fraction, "per") in result_map
                else "",
                "actual_fraction_utg": result_map.get((fraction, "utg"), None).actual_fraction
                if (fraction, "utg") in result_map
                else "",
                "k11_moose_utg": result_map.get((fraction, "utg"), None).k11
                if (fraction, "utg") in result_map
                else "",
                "k11_yvonnet_fine_per": YVONNET_FINE_PER.get(fraction, ""),
                "k11_yvonnet_fine_utg": YVONNET_FINE_UTG.get(fraction, ""),
                "k11_voigt": voigt,
                "k11_reuss": reuss,
                "k11_maxwell": maxwell,
            }
        )

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(results_csv: Path, png_file: Path, mplconfigdir: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))
    mplconfigdir.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel

    with results_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    fractions = [float(row["target_fraction"]) for row in rows]
    per = [float(row["k11_moose_per"]) if row["k11_moose_per"] else math.nan for row in rows]
    utg = [float(row["k11_moose_utg"]) if row["k11_moose_utg"] else math.nan for row in rows]
    ref_per = [
        float(row["k11_yvonnet_fine_per"]) if row["k11_yvonnet_fine_per"] else math.nan
        for row in rows
    ]
    ref_utg = [
        float(row["k11_yvonnet_fine_utg"]) if row["k11_yvonnet_fine_utg"] else math.nan
        for row in rows
    ]

    curve_f = [i / 200.0 for i in range(201)]
    voigt = [analytic_values(f)[0] for f in curve_f]
    reuss = [analytic_values(f)[1] for f in curve_f]
    maxwell = [analytic_values(f)[2] for f in curve_f]

    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    ax.plot(curve_f, voigt, color="0.25", linewidth=1.4, label="Voigt")
    ax.plot(curve_f, reuss, color="0.25", linewidth=1.4, linestyle="--", label="Reuss")
    ax.plot(curve_f, maxwell, color="#2c6b9f", linewidth=1.6, label="Maxwell")
    ax.plot(fractions, utg, marker="o", color="#b33f2e", linewidth=1.2, label="MOOSE, UTG")
    ax.plot(fractions, per, marker="s", color="#20734d", linewidth=1.2, label="MOOSE, PER")
    ax.plot(
        fractions,
        ref_utg,
        marker="o",
        linestyle="-",
        linewidth=2.4,
        markerfacecolor="none",
        markeredgecolor="#b33f2e",
        markeredgewidth=1.6,
        alpha=0.75,
        label="Yvonnet fine, UTG",
    )
    ax.plot(
        fractions,
        ref_per,
        marker="s",
        linestyle="-",
        linewidth=2.4,
        markerfacecolor="none",
        markeredgecolor="#20734d",
        markeredgewidth=1.6,
        alpha=0.75,
        label="Yvonnet fine, PER",
    )
    ax.set_xlabel("fiber volume fraction")
    ax.set_ylabel("$k_{11}^{hom}$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1.0, 5.05)
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.legend(frameon=False, fontsize=8, ncols=2)
    fig.savefig(png_file, dpi=200)
    plt.close(fig)


def make_result_for_zero(modes: list[str]) -> list[MooseResult]:
    return [MooseResult(0.0, mode, K_MATRIX, 0.0) for mode in modes]


def main() -> int:
    repo_heat_transfer = Path(__file__).resolve().parents[4]
    script_dir = Path(__file__).resolve().parent
    default_exe = repo_heat_transfer / "heat_transfer-opt"
    default_workdir = Path(tempfile.gettempdir()) / "moose_yvonnet_fig_3_8"
    default_png_file = script_dir / "yvonnet_fig_3_8.png"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=default_exe, help="Path to heat_transfer-opt.")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=default_workdir,
        help="Directory for generated inputs and outputs.",
    )
    parser.add_argument(
        "--fractions",
        nargs="*",
        default=[],
        help="Fractions to run. Accepts space- or comma-separated values.",
    )
    parser.add_argument("--modes", nargs="+", choices=("per", "utg"), default=["per", "utg"])
    parser.add_argument("--boundary-edges", type=int, default=64)
    parser.add_argument("--circle-segments", type=int, default=128)
    parser.add_argument("--desired-area", type=float, default=1.0e-3)
    parser.add_argument(
        "--png-file",
        type=Path,
        default=default_png_file,
        help="Path for the reproduced figure PNG.",
    )
    parser.add_argument("--no-run", action="store_true", help="Only write generated input files.")
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG plot generation.")
    args = parser.parse_args()

    if not args.exe.exists():
        raise FileNotFoundError(
            f"{args.exe} does not exist. Build it with make -C modules/heat_transfer."
        )
    if args.circle_segments < 16:
        raise ValueError("--circle-segments should be at least 16")
    if args.boundary_edges < 4:
        raise ValueError("--boundary-edges should be at least 4")

    fractions = parse_fractions(args.fractions)
    if any(f < 0.0 or f >= math.pi * 0.5**2 for f in fractions):
        raise ValueError("Fractions must be in [0, pi/4) for a centered circular fiber.")

    args.workdir.mkdir(parents=True, exist_ok=True)
    results: list[MooseResult] = []

    if 0.0 in fractions:
        results.extend(make_result_for_zero(args.modes))

    for fraction in fractions:
        if fraction == 0.0:
            continue
        tag = f"f{fraction:.4f}".replace(".", "p")
        for mode in args.modes:
            file_base = f"yvonnet_{mode}_{tag}"
            input_file = args.workdir / f"{file_base}.i"
            csv_file = args.workdir / f"{file_base}.csv"
            text = (
                periodic_input(fraction, file_base, args)
                if mode == "per"
                else utg_input(fraction, file_base, args)
            )
            input_file.write_text(text)
            if args.no_run:
                continue

            print(f"Running {mode.upper()} f={fraction:g}")
            run_moose(args.exe, input_file)
            k11, actual_fraction = read_csv_value(csv_file)
            results.append(MooseResult(fraction, mode, k11, actual_fraction))

    if args.no_run:
        print(f"Wrote input files to {args.workdir}")
        return 0

    out_csv = args.workdir / "yvonnet_fig_3_8_results.csv"
    write_results(results, out_csv)
    print(f"Wrote {out_csv}")

    if not args.no_plot:
        try:
            out_png = args.png_file
            out_png.parent.mkdir(parents=True, exist_ok=True)
            plot_results(out_csv, out_png, args.workdir / "mplconfig")
            print(f"Wrote {out_png}")
        except ImportError as err:
            print(f"Skipping plot because matplotlib could not be imported: {err}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
