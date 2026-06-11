#!/usr/bin/env python3
"""Run a small FE2-style heat conduction demonstration based on Yvonnet's RVE."""

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
DEFAULT_FRACTION = 0.30


@dataclass(frozen=True)
class MacroResult:
    case: str
    elemental_error: float
    nodal_error: float
    center_temperature: float
    average_k_eff: float | None = None
    min_k_eff: float | None = None
    max_k_eff: float | None = None
    average_macro_grad_T_x: float | None = None
    average_macro_grad_T_y: float | None = None
    average_received_macro_grad_T_x: float | None = None
    average_received_macro_grad_T_y: float | None = None
    average_qbar_x: float | None = None
    average_qbar_y: float | None = None


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


def rve_mesh_block(fraction: float, args: argparse.Namespace) -> str:
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
    desired_area = {args.rve_desired_area:g}
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
    desired_area = {args.rve_desired_area:g}
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


def rve_materials() -> str:
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


def rve_executioner() -> str:
    # QuadraturePointMultiApp derives from TransientMultiApp, so the RVE input is
    # written as a one-step quasi-static transient solve.
    return """
[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-10
  l_max_its = 50
[]
"""


def rve_postprocessors() -> str:
    return """
[Postprocessors]
  [macro_T]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END FINAL'
  []
  [macro_grad_T_x]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END FINAL'
  []
  [macro_grad_T_y]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END FINAL'
  []
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
  [k11]
    type = HomogenizedThermalConductivity
    chi = 'chi_x chi_y'
    row = 0
    col = 0
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [k22]
    type = HomogenizedThermalConductivity
    chi = 'chi_x chi_y'
    row = 1
    col = 1
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [qbar_x]
    type = ParsedPostprocessor
    pp_names = 'k11 macro_grad_T_x'
    pp_symbols = 'k gx'
    expression = '-k*gx'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [qbar_y]
    type = ParsedPostprocessor
    pp_names = 'k22 macro_grad_T_y'
    pp_symbols = 'k gy'
    expression = '-k*gy'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]
"""


def rve_output(file_base: str, csv_output: bool) -> str:
    if csv_output:
        file_base_line = f"  file_base = '{file_base}'"
        csv_block = """
  [csv]
    type = CSV
    execute_on = FINAL
  []
"""
    else:
        # MultiApp instances must not share an explicit file_base.
        file_base_line = ""
        csv_block = ""

    return f"""
[Outputs]
  exodus = false
{file_base_line}
{csv_block}[]
"""


def rve_input(
    fraction: float, file_base: str, args: argparse.Namespace, csv_output: bool = False
) -> str:
    return f"""{rve_mesh_block(fraction, args)}
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
{rve_materials()}
{rve_executioner()}
{rve_postprocessors()}
{rve_output(file_base, csv_output)}
"""


def macro_common(k_ref: float, file_base: str, args: argparse.Namespace) -> str:
    source_scale = 2.0 * k_ref * math.pi * math.pi
    return f"""
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = {args.macro_n}
  ny = {args.macro_n}
  xmax = 1
  ymax = 1
  elem_type = QUAD4
[]

[Variables]
  [T]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [exact_ic]
    type = FunctionIC
    variable = T
    function = exact
  []
[]

[Functions]
  [exact]
    type = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []
  [source_shape]
    type = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = source_shape
    value = {source_scale:.16g}
  []
[]

[BCs]
  [exact_boundary]
    type = FunctionDirichletBC
    variable = T
    boundary = 'left right bottom top'
    function = exact
  []
[]

[Postprocessors]
  [nodal_error]
    type = NodalL2Error
    variable = T
    function = exact
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [elemental_error]
    type = ElementL2Error
    variable = T
    function = exact
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [center_temperature]
    type = PointValue
    variable = T
    point = '0.5 0.5 0'
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]

[VectorPostprocessors]
  [centerline]
    type = LineValueSampler
    variable = T
    start_point = '0 0.5 0'
    end_point = '1 0.5 0'
    num_points = 41
    sort_by = x
    execute_on = FINAL
  []
[]

[Outputs]
  exodus = false
  file_base = '{file_base}'
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
"""


def macro_reference_input(k_ref: float, file_base: str, args: argparse.Namespace) -> str:
    return f"""{macro_common(k_ref, file_base, args)}
[Materials]
  [effective_conductivity]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat density'
    prop_values = '{k_ref:.16g} 1 1'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-11
  l_max_its = 50
[]
"""


def macro_fe2_input(k_ref: float, file_base: str, rve_file_name: str, args: argparse.Namespace) -> str:
    return f"""{macro_common(k_ref, file_base, args)}
[AuxVariables]
  [k_eff]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = {k_ref:.16g}
  []
  [grad_T_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [grad_T_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [received_macro_T]
    family = MONOMIAL
    order = CONSTANT
  []
  [received_macro_grad_T_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [received_macro_grad_T_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [qbar_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [qbar_y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [grad_T_x]
    type = VariableGradientComponent
    variable = grad_T_x
    gradient_variable = T
    component = x
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END FINAL'
  []
  [grad_T_y]
    type = VariableGradientComponent
    variable = grad_T_y
    gradient_variable = T
    component = y
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END FINAL'
  []
[]

[Materials]
  [effective_conductivity_from_rve]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = k_eff
    expression = 'k_eff'
  []
  [thermal_defaults]
    type = GenericConstantMaterial
    prop_names = 'specific_heat density'
    prop_values = '1 1'
  []
[]

[Postprocessors]
  [average_k_eff]
    type = ElementAverageValue
    variable = k_eff
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [min_k_eff]
    type = ElementExtremeValue
    variable = k_eff
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [max_k_eff]
    type = ElementExtremeValue
    variable = k_eff
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_macro_grad_T_x]
    type = ElementAverageValue
    variable = grad_T_x
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_macro_grad_T_y]
    type = ElementAverageValue
    variable = grad_T_y
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_received_macro_grad_T_x]
    type = ElementAverageValue
    variable = received_macro_grad_T_x
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_received_macro_grad_T_y]
    type = ElementAverageValue
    variable = received_macro_grad_T_y
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_qbar_x]
    type = ElementAverageValue
    variable = qbar_x
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
  [average_qbar_y]
    type = ElementAverageValue
    variable = qbar_y
    execute_on = 'INITIAL TIMESTEP_END FINAL'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu       101'
  line_search = none
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-11
  l_max_its = 50
[]

[MultiApps]
  [rves]
    type = QuadraturePointMultiApp
    input_files = '{rve_file_name}'
    execute_on = TIMESTEP_BEGIN
  []
[]

[Transfers]
  [macro_T_to_rve]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    to_multi_app = rves
    source_variable = T
    postprocessor = macro_T
    execute_on = TIMESTEP_BEGIN
  []
  [macro_grad_T_x_to_rve]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    to_multi_app = rves
    source_variable = grad_T_x
    postprocessor = macro_grad_T_x
    execute_on = TIMESTEP_BEGIN
  []
  [macro_grad_T_y_to_rve]
    type = MultiAppVariableValueSamplePostprocessorTransfer
    to_multi_app = rves
    source_variable = grad_T_y
    postprocessor = macro_grad_T_y
    execute_on = TIMESTEP_BEGIN
  []
  [rve_k11_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = k11
    variable = k_eff
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
  [rve_macro_T_echo_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = macro_T
    variable = received_macro_T
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
  [rve_macro_grad_T_x_echo_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = macro_grad_T_x
    variable = received_macro_grad_T_x
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
  [rve_macro_grad_T_y_echo_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = macro_grad_T_y
    variable = received_macro_grad_T_y
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
  [rve_qbar_x_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = qbar_x
    variable = qbar_x
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
  [rve_qbar_y_to_macro]
    type = MultiAppPostprocessorInterpolationTransfer
    from_multi_app = rves
    postprocessor = qbar_y
    variable = qbar_y
    num_points = 4
    execute_on = TIMESTEP_BEGIN
  []
[]
"""


def run_moose(exe: Path, input_file: Path) -> str:
    cmd = [str(exe), "-i", input_file.name, "--color", "off"]
    result = subprocess.run(
        cmd,
        cwd=input_file.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        tail = result.stdout[-6000:]
        raise RuntimeError(
            f"MOOSE failed for {input_file.name} with exit code {result.returncode}\n{tail}"
        )
    return result.stdout


def read_last_row(csv_file: Path) -> dict[str, str]:
    with csv_file.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No data rows found in {csv_file}")
    return rows[-1]


def read_rve_result(csv_file: Path) -> tuple[float, float]:
    row = read_last_row(csv_file)
    return float(row["k11"]), float(row["actual_fraction"])


def read_macro_result(case: str, csv_file: Path) -> MacroResult:
    row = read_last_row(csv_file)
    return MacroResult(
        case=case,
        nodal_error=float(row["nodal_error"]),
        elemental_error=float(row["elemental_error"]),
        center_temperature=float(row["center_temperature"]),
        average_k_eff=float(row["average_k_eff"]) if row.get("average_k_eff") else None,
        min_k_eff=float(row["min_k_eff"]) if row.get("min_k_eff") else None,
        max_k_eff=float(row["max_k_eff"]) if row.get("max_k_eff") else None,
        average_macro_grad_T_x=float(row["average_macro_grad_T_x"])
        if row.get("average_macro_grad_T_x")
        else None,
        average_macro_grad_T_y=float(row["average_macro_grad_T_y"])
        if row.get("average_macro_grad_T_y")
        else None,
        average_received_macro_grad_T_x=float(row["average_received_macro_grad_T_x"])
        if row.get("average_received_macro_grad_T_x")
        else None,
        average_received_macro_grad_T_y=float(row["average_received_macro_grad_T_y"])
        if row.get("average_received_macro_grad_T_y")
        else None,
        average_qbar_x=float(row["average_qbar_x"]) if row.get("average_qbar_x") else None,
        average_qbar_y=float(row["average_qbar_y"]) if row.get("average_qbar_y") else None,
    )


def newest_centerline_csv(workdir: Path, file_base: str) -> Path:
    matches = sorted(workdir.glob(f"{file_base}_centerline_*.csv"))
    if not matches:
        raise RuntimeError(f"No centerline CSV found for {file_base} in {workdir}")
    return matches[-1]


def read_centerline(csv_file: Path) -> list[dict[str, float]]:
    with csv_file.open(newline="") as f:
        rows = list(csv.DictReader(f))
    values: list[dict[str, float]] = []
    for row in rows:
        values.append({key: float(value) for key, value in row.items() if value != ""})
    return values


def centerline_l2_difference(
    left: list[dict[str, float]], right: list[dict[str, float]], variable: str = "T"
) -> float:
    if len(left) != len(right):
        raise RuntimeError("Cannot compare centerlines with different lengths")
    return math.sqrt(sum((a[variable] - b[variable]) ** 2 for a, b in zip(left, right)) / len(left))


def centerline_l2_error(line: list[dict[str, float]]) -> float:
    return math.sqrt(
        sum((row["T"] - math.sin(math.pi * row["x"]) * math.sin(math.pi * row["y"])) ** 2 for row in line)
        / len(line)
    )


def write_summary(
    summary_csv: Path,
    fraction: float,
    actual_fraction: float,
    k_ref: float,
    reference: MacroResult,
    fe2: MacroResult,
    reference_centerline_error: float,
    fe2_centerline_error: float,
    centerline_difference: float,
) -> None:
    rows = [
        {
            "case": "rve",
            "target_fraction": fraction,
            "actual_fraction": actual_fraction,
            "k_eff": k_ref,
            "elemental_error": "",
            "nodal_error": "",
            "center_temperature": "",
            "average_k_eff": "",
            "min_k_eff": "",
            "max_k_eff": "",
            "centerline_l2_vs_exact": "",
            "centerline_l2_vs_reference": "",
            "average_macro_grad_T_x": "",
            "average_macro_grad_T_y": "",
            "average_received_macro_grad_T_x": "",
            "average_received_macro_grad_T_y": "",
            "average_qbar_x": "",
            "average_qbar_y": "",
        },
        {
            "case": "macro_reference",
            "target_fraction": fraction,
            "actual_fraction": actual_fraction,
            "k_eff": k_ref,
            "elemental_error": reference.elemental_error,
            "nodal_error": reference.nodal_error,
            "center_temperature": reference.center_temperature,
            "average_k_eff": k_ref,
            "min_k_eff": k_ref,
            "max_k_eff": k_ref,
            "centerline_l2_vs_exact": reference_centerline_error,
            "centerline_l2_vs_reference": 0.0,
            "average_macro_grad_T_x": "",
            "average_macro_grad_T_y": "",
            "average_received_macro_grad_T_x": "",
            "average_received_macro_grad_T_y": "",
            "average_qbar_x": "",
            "average_qbar_y": "",
        },
        {
            "case": "fe2_qp_rves",
            "target_fraction": fraction,
            "actual_fraction": actual_fraction,
            "k_eff": k_ref,
            "elemental_error": fe2.elemental_error,
            "nodal_error": fe2.nodal_error,
            "center_temperature": fe2.center_temperature,
            "average_k_eff": fe2.average_k_eff,
            "min_k_eff": fe2.min_k_eff,
            "max_k_eff": fe2.max_k_eff,
            "centerline_l2_vs_exact": fe2_centerline_error,
            "centerline_l2_vs_reference": centerline_difference,
            "average_macro_grad_T_x": fe2.average_macro_grad_T_x,
            "average_macro_grad_T_y": fe2.average_macro_grad_T_y,
            "average_received_macro_grad_T_x": fe2.average_received_macro_grad_T_x,
            "average_received_macro_grad_T_y": fe2.average_received_macro_grad_T_y,
            "average_qbar_x": fe2.average_qbar_x,
            "average_qbar_y": fe2.average_qbar_y,
        },
    ]

    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(
    png_file: Path,
    reference_line: list[dict[str, float]],
    fe2_line: list[dict[str, float]],
    k_ref: float,
    reference: MacroResult,
    fe2: MacroResult,
    mplconfigdir: Path,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))
    mplconfigdir.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel

    xs = [row["x"] for row in reference_line]
    exact = [math.sin(math.pi * row["x"]) * math.sin(math.pi * row["y"]) for row in reference_line]
    ref = [row["T"] for row in reference_line]
    fe2_values = [row["T"] for row in fe2_line]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
    ax0.plot(xs, exact, color="0.15", linewidth=2.1, label="analytical")
    ax0.plot(xs, ref, color="#2c6b9f", linewidth=1.8, marker="o", markersize=8, label="macro reference")
    ax0.plot(xs, fe2_values, color="#b33f2e", linewidth=1.8, marker="s", markersize=3, label="FE2 QP RVEs")
    ax0.set_xlabel("x at y = 0.5")
    ax0.set_ylabel("temperature")
    ax0.grid(True, color="0.88", linewidth=0.8)
    ax0.legend(frameon=False, fontsize=8)

    labels = ["macro ref", "FE2"]
    elemental_errors = [reference.elemental_error, fe2.elemental_error]
    ax1.bar(labels, elemental_errors, color=["#2c6b9f", "#b33f2e"], width=0.55)
    ax1.set_ylabel("elemental L2 error vs exact")
    ax1.grid(True, axis="y", color="0.88", linewidth=0.8)
    ax1.set_title(f"RVE k11 = {k_ref:.6g}", fontsize=10)

    fig.savefig(png_file, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    heat_transfer_dir = Path(__file__).resolve().parents[4]

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=heat_transfer_dir / "heat_transfer-opt")
    parser.add_argument("--fraction", type=float, default=DEFAULT_FRACTION)
    parser.add_argument("--macro-n", type=int, default=4, help="Number of macro elements per side.")
    parser.add_argument("--boundary-edges", type=int, default=24)
    parser.add_argument("--circle-segments", type=int, default=48)
    parser.add_argument("--rve-desired-area", type=float, default=2.0e-2)
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "moose_yvonnet_fe2_heat",
        help="Directory for generated MOOSE inputs and raw outputs.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=script_dir / "yvonnet_fe2_heat_summary.csv",
        help="Summary CSV written in the MOOSE source tree by default.",
    )
    parser.add_argument(
        "--png-file",
        type=Path,
        default=script_dir / "yvonnet_fe2_heat_comparison.png",
        help="Comparison PNG written in the MOOSE source tree by default.",
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.exe.exists():
        raise FileNotFoundError(f"{args.exe} does not exist. Build modules/heat_transfer first.")
    if args.macro_n < 1:
        raise ValueError("--macro-n must be positive")
    if args.circle_segments < 16:
        raise ValueError("--circle-segments should be at least 16")
    if args.boundary_edges < 4:
        raise ValueError("--boundary-edges should be at least 4")
    if args.fraction <= 0.0 or args.fraction >= math.pi * 0.5**2:
        raise ValueError("The centered circular fiber fraction must be in (0, pi/4).")

    args.workdir.mkdir(parents=True, exist_ok=True)

    rve_reference_base = "rve_reference"
    rve_reference_file = args.workdir / f"{rve_reference_base}.i"
    rve_reference_file.write_text(rve_input(args.fraction, rve_reference_base, args, csv_output=True))

    rve_subapp_base = "rve_qp_subapp"
    rve_subapp_file = args.workdir / f"{rve_subapp_base}.i"
    rve_subapp_file.write_text(rve_input(args.fraction, rve_subapp_base, args, csv_output=False))

    print(f"Running standalone RVE in {args.workdir}")
    run_moose(args.exe, rve_reference_file)
    k_ref, actual_fraction = read_rve_result(args.workdir / f"{rve_reference_base}.csv")
    print(f"RVE k11 = {k_ref:.10g} at actual area fraction {actual_fraction:.10g}")

    reference_base = "macro_reference"
    reference_file = args.workdir / f"{reference_base}.i"
    reference_file.write_text(macro_reference_input(k_ref, reference_base, args))

    fe2_base = "macro_fe2_qp"
    fe2_file = args.workdir / f"{fe2_base}.i"
    fe2_file.write_text(macro_fe2_input(k_ref, fe2_base, rve_subapp_file.name, args))

    print("Running one-scale macro reference")
    run_moose(args.exe, reference_file)
    print("Running FE2-style macro solve with one RVE per macro quadrature point")
    run_moose(args.exe, fe2_file)

    reference = read_macro_result("macro_reference", args.workdir / f"{reference_base}.csv")
    fe2 = read_macro_result("fe2_qp_rves", args.workdir / f"{fe2_base}.csv")

    reference_line = read_centerline(newest_centerline_csv(args.workdir, reference_base))
    fe2_line = read_centerline(newest_centerline_csv(args.workdir, fe2_base))
    reference_centerline_error = centerline_l2_error(reference_line)
    fe2_centerline_error = centerline_l2_error(fe2_line)
    centerline_difference = centerline_l2_difference(reference_line, fe2_line)

    write_summary(
        args.summary_csv,
        args.fraction,
        actual_fraction,
        k_ref,
        reference,
        fe2,
        reference_centerline_error,
        fe2_centerline_error,
        centerline_difference,
    )
    print(f"Wrote {args.summary_csv}")

    if not args.no_plot:
        try:
            args.png_file.parent.mkdir(parents=True, exist_ok=True)
            plot_results(
                args.png_file,
                reference_line,
                fe2_line,
                k_ref,
                reference,
                fe2,
                args.workdir / "mplconfig",
            )
            print(f"Wrote {args.png_file}")
        except ImportError as err:
            print(f"Skipping plot because matplotlib could not be imported: {err}", file=sys.stderr)

    print(
        "FE2 centerline RMS difference from one-scale macro = "
        f"{centerline_difference:.6e}; average transferred k_eff = {fe2.average_k_eff:.10g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
