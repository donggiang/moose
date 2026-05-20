#!/usr/bin/env python3
"""CG wedge-mesh DCB CZM result vs experimental data (Davila DCB) and
LEFM/beam-theory landmarks.

Inputs:
  tier2_1_dcb_cg_davila_wedge_out.csv  -- CG point-pin (y=+/-1.98) on
                                          msh_dcb_h0p5_2blocks_rf3.e
  dg_notes/dcb_exp.csv                 -- experimental P vs opening
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Davila DCB parameters
E       = 150e3
nu      = 0.25
h       = 1.98
a0      = 55.0
width   = 20.0
GIc     = 0.268

E_eff   = E / (1.0 - nu * nu)
P_c     = (E_eff * h ** 3 * GIc / (12.0 * a0 ** 2)) ** 0.5
P_c_N   = P_c * width                                # ~60.6 N
slope_per_width = E_eff * h ** 3 / (8.0 * a0 ** 3)
slope           = slope_per_width * width            # ~18.66 N/mm of opening


def read_csv(path, header=True):
    with path.open() as f:
        rows = list(csv.reader(f))
    if header:
        keys = [k.strip() for k in rows[0]]
        return {k: np.array([float(r[i]) for r in rows[1:]])
                for i, k in enumerate(keys)}
    arr = np.array([[float(c) for c in r] for r in rows])
    return arr


def main():
    here = Path(__file__).parent

    sim = read_csv(here / "tier2_1_dcb_cg_davila_wedge_out.csv")
    sim_open = sim['opening']
    sim_P    = sim['P_top_N']
    sim_peak = sim_P.max()
    sim_peak_open = sim_open[sim_P.argmax()]

    exp = read_csv(here / "dg_notes" / "dcb_exp.csv", header=False)
    exp_open = exp[:, 0]
    exp_P    = exp[:, 1]
    exp_peak = exp_P.max()
    exp_peak_open = exp_open[exp_P.argmax()]

    fig, ax = plt.subplots(figsize=(8.5, 6))

    open_max = max(sim_open.max(), exp_open.max())
    line = np.linspace(0, open_max, 50)
    ax.plot(line, slope * line, "-", color="black", lw=1.4,
            label=fr"Beam theory linear elastic, $P/\delta$ = {slope:.2f} N/mm")
    ax.axhline(P_c_N, color="tab:green", lw=1.6, ls="--",
               label=fr"LEFM peak ($a_0$ = {a0:.0f} mm): {P_c_N:.1f} N")

    ax.plot(exp_open, exp_P, "o", color="tab:red", ms=6, mfc='white', mew=1.2,
            label=fr"Experiment (Davila DCB) -- peak {exp_peak:.1f} N "
                  fr"at {exp_peak_open:.2f} mm")
    ax.plot(sim_open, sim_P, "-^", color="tab:purple", ms=5, lw=1.6,
            label=fr"CG wedge mesh, outer pin -- peak {sim_peak:.1f} N "
                  fr"at {sim_peak_open:.2f} mm")

    ax.set_xlabel("end opening [mm]")
    ax.set_ylabel("reaction P [N]  (NodalSum $\\cdot$ width)")
    ax.set_title("Tier 2.1 DCB CZM (Davila Fig 8 material) -- "
                 "CG wedge mesh vs experiment")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, open_max * 1.02)
    ax.set_ylim(0, max(sim_P.max(), exp_P.max(), P_c_N) * 1.20)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    out = here / "tier2_1_wedge_vs_exp.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    print(f"Experiment peak  : {exp_peak:.2f} N at opening {exp_peak_open:.2f} mm")
    print(f"CG wedge peak    : {sim_peak:.2f} N at opening {sim_peak_open:.2f} mm "
          f"({100.0 * sim_peak / exp_peak:.1f} % of experiment)")
    print(f"LEFM peak        : {P_c_N:.1f} N "
          f"({100.0 * sim_peak / P_c_N:.1f} % of LEFM)")
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
