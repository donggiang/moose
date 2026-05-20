#!/usr/bin/env python3
"""Tier 2.1 DCB CZM, displacement-controlled point pin:
CG (tier2_1_dcb_cg_davila_out.csv) vs DG-CZM (dg_phase5_dcb_davila_dgczm_out.csv)
vs LEFM/Davila Figure 8 expected peak.

P_top_N is the total upper-arm reaction (per-width x 20 mm), opening is the
relative y-displacement between the upper and lower pins.

LEFM peak (plane strain, a0 = 55 mm):
  P_c = sqrt(E_eff h^3 G_Ic / (12 a^2)) per unit width  -> 3.03 N/mm
  P_c x width  = 60.6 N   (Davila Figure 8 measured ~63 N)
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Davila DCB parameters (Fig 8 / Table 1)
E       = 150e3
nu      = 0.25
h       = 1.98
a0      = 55.0
width   = 20.0
GIc     = 0.268

E_eff   = E / (1.0 - nu * nu)
I_arm   = h ** 3 / 12.0
P_c     = (E_eff * h ** 3 * GIc / (12.0 * a0 ** 2)) ** 0.5
P_c_N   = P_c * width

# Linear structural slope (no crack growth): two cantilever arms of length a0
# loaded equal-and-opposite at the tip, opening = 2 * tip_disp
slope_per_width = E_eff * h ** 3 / (8.0 * a0 ** 3)
slope           = slope_per_width * width    # [N / mm of end-opening]


def read(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}


def main():
    here = Path(__file__).parent

    cg       = read(here / "tier2_1_dcb_cg_davila_out.csv")
    cg_wedge = read(here / "tier2_1_dcb_cg_davila_wedge_out.csv")
    dg       = read(here / "dg_notes" / "dg_phase5_dcb_davila_dgczm_out.csv")

    cg_open   = cg['opening']
    cg_P      = cg['P_top_N']
    cgw_open  = cg_wedge['opening']
    cgw_P     = cg_wedge['P_top_N']
    cgw_jump  = cg_wedge['max_normal_jump']
    cgw_trac  = cg_wedge['max_normal_traction']
    dg_open   = dg['opening']
    dg_P      = dg['P_top_N']

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.5))

    # ---- Left: full range load-opening -----------------------------------
    open_max = max(cg_open.max(), dg_open.max())
    line = np.linspace(0, open_max, 50)
    axL.plot(line, slope * line, "-", color="black", lw=1.6,
             label=fr"Linear elastic (no crack), slope = {slope:.2f} N/mm")
    axL.axhline(P_c_N, color="tab:green", lw=1.6, ls="--",
                label=fr"LEFM peak ($a_0$=55): {P_c_N:.1f} N  (Davila Fig 8 ~ 63 N)")
    axL.plot(cg_open, cg_P, "s-", color="tab:red", ms=4, lw=1.4,
             label=f"CG point-pin generated mesh (y=$\\pm$0.99 pin)")
    cgw_peak = cgw_P.max()
    cgw_peak_open = cgw_open[cgw_P.argmax()]
    axL.plot(cgw_open, cgw_P, "^-", color="tab:purple", ms=4, lw=1.6,
             label=f"CG wedge mesh, outer pin (y=$\\pm$1.98) -- "
                   f"peak {cgw_peak:.1f} N at {cgw_peak_open:.2f} mm")
    axL.plot(dg_open, dg_P, "o-", color="tab:blue", ms=4, lw=1.4,
             label=f"DG-CZM face-pin")
    axL.set_xlabel("end opening [mm]  (= 2 * tip disp)")
    axL.set_ylabel("upper-arm reaction P [N]  (per-width x 20)")
    axL.set_title("Tier 2.1 DCB CZM (Davila Fig 8 material) -- displacement control")
    axL.grid(alpha=0.3)
    axL.legend(loc="upper left", fontsize=9)
    axL.set_xlim(0, open_max * 1.02)
    axL.set_ylim(0, max(cg_P.max(), dg_P.max()) * 1.08)

    # ---- Right: zoom near LEFM peak --------------------------------------
    zoom = 6.0   # mm
    axR.plot(line[line <= zoom], slope * line[line <= zoom], "-",
             color="black", lw=1.6, label=fr"Linear elastic, {slope:.2f} N/mm")
    axR.axhline(P_c_N, color="tab:green", lw=1.6, ls="--",
                label=fr"LEFM peak {P_c_N:.1f} N")
    axR.plot(cg_open, cg_P, "s-", color="tab:red", ms=5, lw=1.4,
             label="CG point-pin generated mesh")
    axR.plot(cgw_open, cgw_P, "^-", color="tab:purple", ms=5, lw=1.6,
             label=f"CG wedge mesh, outer pin (peak {cgw_peak:.1f} N)")
    axR.plot(dg_open, dg_P, "o-", color="tab:blue", ms=5, lw=1.4,
             label="DG-CZM face-pin")
    axR.set_xlabel("end opening [mm]")
    axR.set_ylabel("reaction P [N]")
    axR.set_title(fr"Zoom: opening $\leq$ {zoom:.0f} mm (LEFM-peak region)")
    axR.grid(alpha=0.3)
    axR.legend(loc="upper left", fontsize=9)
    axR.set_xlim(0, zoom)
    # y-limit: focus near LEFM peak band
    axR.set_ylim(0, max(P_c_N * 2.5,
                        cg_P[cg_open <= zoom].max() if (cg_open <= zoom).any() else P_c_N))

    fig.tight_layout()
    out = here / "tier2_1_cg_vs_dg_davila.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    # Brief summary
    print(f"Linear structural slope = {slope:.3f} N/mm of opening")
    print(f"LEFM peak (a0={a0}, GIc={GIc}) = {P_c_N:.1f} N  (Davila Fig 8 ~63 N)")
    print(f"CG (gen): final opening = {cg_open[-1]:.2f} mm, final P = {cg_P[-1]:.2f} N, "
          f"max_jump = {cg['max_normal_jump'][-1]:.2e} mm")
    print(f"CG (wedge): peak P = {cgw_peak:.2f} N at opening = {cgw_peak_open:.2f} mm, "
          f"final P = {cgw_P[-1]:.2f} N at opening = {cgw_open[-1]:.2f} mm, "
          f"final max_jump = {cgw_jump[-1]:.2e} mm")
    print(f"DG : final opening = {dg_open[-1]:.2f} mm, final P = {dg_P[-1]:.2f} N, "
          f"max_jump = {dg['max_normal_jump'][-1]:.2e} mm")
    print(f"Plot: {out}")


if __name__ == "__main__":
    main()
