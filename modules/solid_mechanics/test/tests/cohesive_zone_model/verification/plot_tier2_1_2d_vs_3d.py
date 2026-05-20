#!/usr/bin/env python3
"""Tier 2.1 DCB CZM: 2D wedge vs 2D MOOSE 2-block vs 3D MOOSE 2-block,
overlaid on Davila experiment and LEFM/beam-theory landmarks.

In 3D the pin becomes a LINE load across the full width (z direction),
so the NodalSum reaction is the TOTAL force; no per-width multiplication.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
slope           = slope_per_width * width            # ~18.66 N/mm


def read_csv(path, header=True):
    with path.open() as f:
        rows = list(csv.reader(f))
    if header:
        keys = [k.strip() for k in rows[0]]
        return {k: np.array([float(r[i]) for r in rows[1:]])
                for i, k in enumerate(keys)}
    return np.array([[float(c) for c in r] for r in rows])


def main():
    here = Path(__file__).parent

    runs = [
        ("2D Cubit wedge (h=0.5)",  "tier2_1_dcb_cg_davila_wedge_h0p5_out.csv",
         "P_top_N", {"color": "tab:purple", "marker": "^"}),
        ("2D MOOSE 2-block (h=0.5)", "tier2_1_dcb_cg_davila_2blk_out.csv",
         "P_top_N", {"color": "tab:blue",   "marker": "D"}),
        ("3D MOOSE 2-block (nx=150, ny=4, nz=2)", "tier2_1_dcb_cg_davila_3d_2blk_out.csv",
         "P_top_abs", {"color": "tab:green", "marker": "o"}),
    ]
    sims = []
    for label, fname, col, style in runs:
        d = read_csv(here / fname)
        opening = d['opening']
        P       = d[col]
        order   = np.argsort(opening)
        opening = opening[order]
        P       = P[order]
        peak    = P.max()
        peak_open = opening[P.argmax()]
        sims.append({'label': label, 'style': style, 'opening': opening, 'P': P,
                     'peak': peak, 'peak_open': peak_open,
                     'jump_final': d['max_normal_jump'][np.argmax(d['opening'])]})
        print(f"{label:<40} peak {peak:>5.2f} N at opening {peak_open:>4.2f} mm "
              f"(final opening {opening[-1]:.2f} mm, final P {P[-1]:.2f} N)")

    exp = read_csv(here / "dg_notes" / "dcb_exp.csv", header=False)
    exp_open = exp[:, 0]
    exp_P    = exp[:, 1]
    exp_peak = exp_P.max()
    exp_peak_open = exp_open[exp_P.argmax()]
    print(f"Experiment                               peak {exp_peak:.2f} N at "
          f"opening {exp_peak_open:.2f} mm")

    fig, ax = plt.subplots(figsize=(10, 6))

    open_max = max(max(s['opening'].max() for s in sims), exp_open.max())
    line = np.linspace(0, open_max, 50)
    ax.plot(line, slope * line, "-", color="black", lw=1.4,
            label=fr"Beam theory linear elastic, $P/\delta$ = {slope:.2f} N/mm")
    ax.axhline(P_c_N, color="tab:gray", lw=1.6, ls="--",
               label=fr"LEFM peak ($a_0$ = {a0:.0f} mm): {P_c_N:.1f} N")
    ax.plot(exp_open, exp_P, "o", color="tab:red", ms=6, mfc='white', mew=1.2,
            label=fr"Experiment (Davila DCB) -- peak {exp_peak:.1f} N "
                  fr"at {exp_peak_open:.2f} mm")
    for s in sims:
        ax.plot(s['opening'], s['P'], "-", **s['style'], ms=4, lw=1.4,
                label=fr"{s['label']} -- peak {s['peak']:.1f} N "
                      fr"at {s['peak_open']:.2f} mm")

    ax.set_xlabel("end opening [mm]")
    ax.set_ylabel("reaction P [N]")
    ax.set_title("Tier 2.1 DCB CZM CG -- 2D vs 3D 2-block vs experiment")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, open_max * 1.02)
    y_top = max(max(s['P'].max() for s in sims), exp_P.max(), P_c_N) * 1.20
    ax.set_ylim(0, y_top)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    out = here / "tier2_1_2d_vs_3d.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
