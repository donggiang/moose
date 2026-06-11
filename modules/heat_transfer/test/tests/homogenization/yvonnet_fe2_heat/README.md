# Yvonnet FE2 Heat Conduction Demo

This directory contains a runnable FE2-style heat conduction demonstration built
from the Yvonnet circular-fiber RVE used in `../yvonnet_fig_3_8`.

The script performs three solves:

1. A standalone periodic RVE solve at fiber area fraction `f = 0.3` to compute
   `k11` with `HomogenizedHeatConduction` and `HomogenizedThermalConductivity`.
2. A one-scale 2D macro heat equation using that RVE-computed conductivity as a
   constant `thermal_conductivity`.
3. A multiscale macro heat equation using `QuadraturePointMultiApp`, with one
   Yvonnet-style RVE subapp at each macro quadrature point. The RVE `k11`
   postprocessor is transferred back to a macro elemental AuxVariable `k_eff`,
   and `ParsedMaterial` exposes `k_eff` as the macro `thermal_conductivity`.

The macro manufactured solution is

```text
T(x,y) = sin(pi*x) * sin(pi*y)
```

on the unit square with Dirichlet boundary conditions from the exact solution.
The source term is `2*k_ref*pi^2*sin(pi*x)*sin(pi*y)`, where `k_ref` is the
standalone RVE conductivity.

Run with the MOOSE conda environment:

```bash
conda run -n moose python modules/heat_transfer/test/tests/homogenization/yvonnet_fe2_heat/run_yvonnet_fe2_heat.py \
  --workdir /private/tmp/moose_yvonnet_fe2_heat
```

The script writes generated MOOSE inputs and raw outputs to `--workdir`; if this
argument is omitted, it uses Python's system temporary directory. It writes the
summary CSV and comparison PNG into this directory:

- `yvonnet_fe2_heat_summary.csv`
- `yvonnet_fe2_heat_comparison.png`

The `extracted_inputs/` directory contains the generated `.i` files from the
8x8 macro-mesh FE2 run:

- `rve_reference.i`: standalone RVE solve used to compute `k_ref`
- `macro_reference.i`: one-scale macro reference solve
- `rve_qp_subapp.i`: RVE input instantiated at each macro quadrature point
- `macro_fe2_qp.i`: macro FE2 solve with `QuadraturePointMultiApp`

This is a loose FE2/MultiApp prototype for linear heat conduction. Because the
RVE is linear and homogeneous at the macro scale, the effective conductivity is
constant; the point of this case is to verify the macro/RVE execution and data
transfer path before extending it to gradient-dependent or nonlinear RVEs.
