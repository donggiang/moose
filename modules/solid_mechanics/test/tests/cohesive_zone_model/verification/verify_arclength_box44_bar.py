#!/usr/bin/env python3
"""
Box 4.4 arc-length verification on a 2D linear-elastic bar in plane strain.
Compares the (lambda, tip_disp_x) trajectory against the LEFM closed-form
relation tip_disp = lambda * sigma_ref * L / E_eff with
E_eff = E / (1 - nu**2) for plane strain.

Usage:
    python3 verify_arclength_box44_bar.py [path/to/arclength_box44_bar_out.csv]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference parameters (must match arclength_box44_bar.i).
E         = 1.0e5      # MPa
NU        = 0.3
SIGMA_REF = 1.0e3      # MPa, reference traction at lambda=1
L         = 10.0       # mm, bar length
DELTA_S   = 0.01       # arc-length step (input)

# Plane-strain effective modulus for uniaxial bar with sigma_yy=0:
E_EFF = E / (1.0 - NU**2)

# tip displacement at lambda=1 (LEFM):
TIP_AT_UNIT_LAMBDA = SIGMA_REF * L / E_EFF

# Pass tolerance — relative on slope of tip_disp vs lambda.
REL_TOL = 5e-3   # 0.5 % is comfortable; should be 0.2 % from observed run.


def lefm_tip_disp(lam: float) -> float:
    return lam * TIP_AT_UNIT_LAMBDA


def make_plot(times, lambdas, tips, errs, out_path: Path) -> None:
    lam_dense = np.linspace(0.0, max(lambdas) * 1.05, 200)
    tip_dense = np.array([lefm_tip_disp(l) for l in lam_dense])

    fig, (ax_curve, ax_err) = plt.subplots(
        2, 1, figsize=(7.5, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    ax_curve.plot(lam_dense, tip_dense, "-", color="tab:blue", lw=1.5,
                  label=r"LEFM: $u_{\rm tip} = \lambda\,\sigma_{\rm ref} L / E_{\rm eff}$")
    ax_curve.plot(lambdas, tips, "o", color="tab:red", ms=5, mfc="none",
                  label="Box 4.4 arc-length")
    ax_curve.set_ylabel(r"tip displacement $u_x$ [mm]")
    ax_curve.set_title(
        r"Box 4.4 arc-length on linear-elastic bar — "
        r"$\Delta s = " + f"{DELTA_S:g}" + r"$, $E_{\rm eff} = E/(1-\nu^2)$"
    )
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper left")

    floor = 1e-16
    ax_err.semilogy(lambdas, np.maximum(errs, floor), "o-", color="tab:purple",
                    ms=3, lw=0.8, label="absolute error")
    ax_err.axhline(REL_TOL * TIP_AT_UNIT_LAMBDA,
                   color="tab:red", lw=0.8, ls="--",
                   label=f"tolerance ({REL_TOL*100:g}% of $u(\\lambda=1)$)")
    ax_err.set_xlabel(r"load factor $\lambda$")
    ax_err.set_ylabel(r"$|u_{\rm sim} - u_{\rm LEFM}|$ [mm]")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "arclength_box44_bar_out.csv"
    )
    if not csv_path.is_file():
        print(f"FAIL: CSV not found at {csv_path}")
        return 2

    times, lambdas, tips = [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            lambdas.append(float(row["lambda"]))
            tips.append(float(row["tip_disp_x"]))

    times = np.asarray(times)
    lambdas = np.asarray(lambdas)
    tips = np.asarray(tips)

    tip_ref = np.array([lefm_tip_disp(l) for l in lambdas])
    errs = np.abs(tips - tip_ref)

    # Slope from least-squares (skip step 0 where both are zero).
    nonzero = lambdas > 0
    slope = float(np.polyfit(lambdas[nonzero], tips[nonzero], 1)[0])
    slope_ref = TIP_AT_UNIT_LAMBDA
    rel_slope_err = abs(slope - slope_ref) / slope_ref

    plot_path = csv_path.with_name("arclength_box44_bar_compare.png")
    make_plot(times, lambdas, tips, errs, plot_path)

    print("Box 4.4 arc-length — linear-elastic bar (plane strain)")
    print(f"  rows checked          : {len(times)}")
    print(f"  E_eff = E/(1-nu^2)    : {E_EFF:.6e}")
    print(f"  tip(lambda=1) ref     : {TIP_AT_UNIT_LAMBDA:.6e} mm")
    print(f"  measured slope        : {slope:.6e} mm/lambda")
    print(f"  expected slope        : {slope_ref:.6e} mm/lambda")
    print(f"  relative slope error  : {rel_slope_err*100:.3f} %")
    print(f"  max |u_sim - u_LEFM|  : {np.max(errs):.6e} mm")
    print(f"  comparison plot       : {plot_path}")

    if rel_slope_err < REL_TOL:
        print("PASS")
        return 0
    print(f"FAIL: slope error {rel_slope_err*100:.3f}% exceeds {REL_TOL*100:g}%")
    return 1


if __name__ == "__main__":
    sys.exit(main())
