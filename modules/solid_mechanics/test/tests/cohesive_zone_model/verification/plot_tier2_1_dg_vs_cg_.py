#!/usr/bin/env python3
"""DG cohesive vs CG cohesive on Davila DCB."""
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


def read(path, header=True):
    with path.open() as f:
        rows = list(csv.reader(f))
    if header:
        keys = [k.strip() for k in rows[0]]
        return {k: np.array([float(r[i]) for r in rows[1:]])
                for i, k in enumerate(keys)}
    return np.array([[float(c) for c in r] for r in rows])


def main():
    here = Path(__file__).parent

    cg2d = read(here / "tier2_1_dcb_cg_davila_2blk_out.csv")
    dg     = read(here / "tier2_1_dcb_dg_davila_out.csv")
    #dg_e5  = read(here / "tier2_1_dcb_dg_davila_eta5_out.csv")
    dg_e5v = read(here / "tier2_1_dcb_dg_davila_out.csv")
    #dg_al  = read(here / "tier2_1_dcb_dg_davila_AL_out.csv")
    #cg_al  = read(here / "tier2_1_dcb_cg_davila_AL_out.csv")
    #cg3d_al = read(here / "tier2_1_dcb_cg_davila_3d_AL_out.csv")
    exp  = read(here /"dcb_exp.csv", header=False)

    cg2d_open = cg2d['opening']
    cg2d_P    = cg2d['P_top_N']
    cg2d_peak = cg2d_P.max()
    cg2d_peak_open = cg2d_open[cg2d_P.argmax()]

    dg_open = dg['opening']
    dg_P    = dg['P_top_N']
    dg_peak = dg_P.max()
    dg_peak_open = dg_open[dg_P.argmax()]

#     dg5_open = dg_e5['opening']
#     dg5_P    = dg_e5['P_top_N']
#     dg5_peak = dg5_P.max()
#     dg5_peak_open = dg5_open[dg5_P.argmax()]

    dgv_open = dg_e5v['opening']
    dgv_P    = dg_e5v['P_top_N']
    dgv_peak = dgv_P.max()
    dgv_peak_open = dgv_open[dgv_P.argmax()]

#     al_open = dg_al['opening']
#     al_P    = dg_al['P_AL_N']    # AL-correct reaction = lambda * f_ref * area
#     al_peak = al_P.max()
#     al_peak_open = al_open[al_P.argmax()]

    # CG AL: predictor sign went negative on first step (closing branch);
    # plot |opening| vs |P_AL_N| to compare on the same axes with the others.
#     cgal_open = np.abs(cg_al['opening'])
#     cgal_P    = np.abs(cg_al['P_AL_N'])
#     cgal_peak = cgal_P.max()
#     cgal_peak_open = cgal_open[cgal_P.argmax()]

#     # CG 3D AL: same direction issue; |abs| compare.
#     cg3d_open = np.abs(cg3d_al['opening'])
#     cg3d_P    = np.abs(cg3d_al['P_AL_N'])
#     cg3d_peak = cg3d_P.max()
#     cg3d_peak_open = cg3d_open[cg3d_P.argmax()]

    exp_open = exp[:, 0]
    exp_P    = exp[:, 1]
    exp_peak = exp_P.max()
    exp_peak_open = exp_open[exp_P.argmax()]

#     print(f"CG 2D 2-block:  peak {cg2d_peak:.2f} N at opening {cg2d_peak_open:.2f} mm, "
#           f"final P={cg2d_P[-1]:.2f} at opening={cg2d_open[-1]:.2f}")
#     print(f"DG cohesive eta=5:  peak {dg5_peak:.2f} N at opening {dg5_peak_open:.2f} mm, "
#           f"final P={dg5_P[-1]:.2f} at opening={dg5_open[-1]:.2f}")

   # print(f"Experiment:     peak {exp_peak:.2f} N at opening {exp_peak_open:.2f} mm")

    fig, ax = plt.subplots(figsize=(9.5, 6))

    open_max = max(cg2d_open.max(), exp_open.max(), dg_open.max())
    line = np.linspace(0, open_max, 50)
    ax.plot(line, slope * line, "-", color="black", lw=1.4,
            label=fr"Beam theory linear")
    ax.axhline(P_c_N, color="tab:gray", lw=1.6, ls="--",
               label=fr"LEFM peak")
    ax.plot(exp_open, exp_P, "o", color="tab:red", ms=6, mfc='white', mew=1.2,
            label=f"Experiment")
    ax.plot(cg2d_open, cg2d_P, "D-", color="tab:blue", ms=4, lw=1.4,
            label=f"CG")

    ax.plot(dgv_open, dgv_P, "v-", color="tab:purple", ms=4, lw=1.4,
            label=fr"DG")

    ax.set_xlabel("end opening [mm]")
    ax.set_ylabel("reaction P [N]")
    #ax.set_title("Tier 2.1 DCB CZM CG (BreakMeshByBlock + CZM block) vs DG cohesive")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, open_max * 1.02)
    ax.set_ylim(0, max(cg2d_P.max(), exp_P.max(), dg_P.max(), P_c_N) * 1.20)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    out = here / "tier2_1_dg_vs_cg.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
