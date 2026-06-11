# Yvonnet Figure 3.8 Reproduction

This workflow reproduces the heat-conduction homogenization sweep from
Yvonnet, *Computational Homogenization of Heterogeneous Materials*, Fig. 3.8.
The model is a 2D unit cell with a centered circular fiber. The matrix
conductivity is `k1 = 1 W/m/K`, the fiber conductivity is `k2 = 5 W/m/K`, and
the plotted response is the effective conductivity component `k11` versus fiber
area fraction. The generated polygonal fiber radius is chosen so the polygon
area equals the requested target fraction.

The periodic calculation uses the heat-transfer homogenization objects:

- `HomogenizedHeatConduction`
- `HomogenizedThermalConductivity`

The uniform-temperature-gradient calculation applies `T = x` on the outer
boundary and evaluates `k11` from the volume average of `k dT/dx`.

Run from the repository root with the MOOSE conda environment:

```bash
conda run -n moose make -C modules/heat_transfer -j 8
conda run -n moose python modules/heat_transfer/test/tests/homogenization/yvonnet_fig_3_8/reproduce_yvonnet_fig_3_8.py
```

By default the script writes generated inputs, per-run CSV files, and
`yvonnet_fig_3_8_results.csv` under the system temporary directory:

```text
<tempdir>/moose_yvonnet_fig_3_8
```

Use `--workdir /path/to/output` to choose a fixed output directory.
The reproduced PNG is written beside this README as `yvonnet_fig_3_8.png`.
Use `--png-file /path/to/figure.png` to choose a different figure path.

For a faster smoke test:

```bash
conda run -n moose python modules/heat_transfer/test/tests/homogenization/yvonnet_fig_3_8/reproduce_yvonnet_fig_3_8.py \
  --fractions 0 0.3 \
  --boundary-edges 16 \
  --circle-segments 48 \
  --desired-area 0.01
```

For a finer mesh, reduce `--desired-area` and increase
`--circle-segments`/`--boundary-edges`. The analytic curves in the output are:

- Voigt: `k = (1 - f) k1 + f k2`
- Reuss: `k = 1 / ((1 - f) / k1 + f / k2)`
- Maxwell, with `D = 2`:
  `k = k1 + D k1 f (k2 - k1) / (D k1 + (1 - f)(k2 - k1))`
