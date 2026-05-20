#!/usr/bin/env python3
"""
Tier 1.2 verification: BiLinearMixedModeTraction under pure mode-II monotonic shear.

For pure mode II (delta_n = 0), the bilinear law reduces to the same shape
as Tier 1.1 with the shear-channel parameters:

  t_t(delta_s) = K * delta_s                          0 <= delta_s <= delta_s_0
               = S * (delta_s_f - delta_s)/(delta_s_f - delta_s_0)
                                                      delta_s_0 <= delta_s <= delta_s_f
               = 0                                    delta_s >= delta_s_f

with delta_s_0 = S/K, delta_s_f = 2 GIIc/S. We additionally require that
no spurious normal_traction develops since delta_n is pinned to zero by
the boundary conditions.

Usage:
    python3 verify_tier1_2.py [path/to/tier1_2_mode_II_monotonic_out.csv]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference law parameters (must match the input file).
K       = 5.0e3
S       = 30.0
GII_c   = 0.5

delta_0 = S / K              # 6.0e-3
delta_f = 2.0 * GII_c / S    # 3.333e-2

ABS_TOL_TRACTION = 1e-7
NORMAL_TOL       = 1e-9   # spurious mode-I traction must vanish


def analytical_traction(delta_s: float) -> float:
    if delta_s <= 0:
        return 0.0
    if delta_s <= delta_0:
        return K * delta_s
    if delta_s <= delta_f:
        return S * (delta_f - delta_s) / (delta_f - delta_0)
    return 0.0


def make_plot(deltas, tractions, ref, errs, out_path: Path) -> None:
    d_dense = np.linspace(0.0, max(delta_f, deltas.max()) * 1.05, 1000)
    t_dense = np.array([analytical_traction(d) for d in d_dense])

    fig, (ax_curve, ax_err) = plt.subplots(
        2, 1, figsize=(7.5, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    ax_curve.plot(d_dense, t_dense, "-", color="tab:blue", lw=1.5, label="analytical")
    ax_curve.plot(deltas, tractions, "o", color="tab:red", ms=4, mfc="none",
                  label="simulation")
    ax_curve.axvline(delta_0, color="0.6", lw=0.8, ls="--")
    ax_curve.axvline(delta_f, color="0.6", lw=0.8, ls="--")
    ax_curve.text(delta_0, S, r" $\delta_{s,0}$", va="bottom", color="0.4")
    ax_curve.text(delta_f, 0, r" $\delta_{s,f}$", va="bottom", color="0.4")
    ax_curve.set_ylabel(r"tangent traction $t_t$")
    ax_curve.set_title("Tier 1.2 — bilinear mode-II traction–separation")
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper right")

    floor = 1e-16
    ax_err.semilogy(deltas, np.maximum(errs, floor), "o-", color="tab:purple",
                    ms=3, lw=0.8)
    ax_err.axhline(ABS_TOL_TRACTION, color="tab:red", lw=0.8, ls="--",
                   label=f"tolerance {ABS_TOL_TRACTION:.0e}")
    err_top = max(ABS_TOL_TRACTION * 10,
                  errs.max() * 10 if errs.max() > 0 else ABS_TOL_TRACTION)
    ax_err.set_ylim(floor, err_top)
    ax_err.set_xlabel(r"tangential jump $\delta_s$")
    ax_err.set_ylabel(r"$|t_t^{\mathrm{sim}} - t_t^{\mathrm{ref}}|$")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "tier1_2_mode_II_monotonic_out.csv"
    )
    if not csv_path.is_file():
        print(f"FAIL: CSV not found at {csv_path}")
        return 2

    times, deltas, tractions, normals_t, normals_j = [], [], [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            deltas.append(float(row["tangent_jump"]))
            tractions.append(float(row["tangent_traction"]))
            normals_t.append(float(row["normal_traction"]))
            normals_j.append(float(row["normal_jump"]))

    times      = np.asarray(times)
    deltas     = np.asarray(deltas)
    tractions  = np.asarray(tractions)
    normals_t  = np.asarray(normals_t)
    normals_j  = np.asarray(normals_j)

    ref = np.array([analytical_traction(d) for d in deltas])
    errs = np.abs(tractions - ref)
    idx = int(np.argmax(errs))
    max_err = float(errs[idx])
    max_normal_traction = float(np.max(np.abs(normals_t)))
    max_normal_jump     = float(np.max(np.abs(normals_j)))

    plot_path = csv_path.with_name(
        csv_path.stem.replace("_out", "") + "_compare.png"
    )
    make_plot(deltas, tractions, ref, errs, plot_path)

    print("Tier 1.2 — mode-II monotonic shear")
    print(f"  rows checked          : {len(times)}")
    print(f"  delta_s_0 (analytical): {delta_0:.6e}")
    print(f"  delta_s_f (analytical): {delta_f:.6e}")
    print(f"  peak traction         : {S:.6e}")
    print(f"  max |tt - tt_ref|     : {max_err:.6e}   (tol {ABS_TOL_TRACTION:.0e})")
    t, dn, tt, tt_ref = times[idx], deltas[idx], tractions[idx], ref[idx]
    print(f"    at t = {t:.4f}, delta_s = {dn:.6e}, tt = {tt:.6e}, tt_ref = {tt_ref:.6e}")
    print(f"  max |normal_traction| : {max_normal_traction:.6e}   (tol {NORMAL_TOL:.0e})")
    print(f"  max |normal_jump|     : {max_normal_jump:.6e}")
    print(f"  comparison plot       : {plot_path}")

    if max_err < ABS_TOL_TRACTION and max_normal_traction < NORMAL_TOL:
        print("PASS")
        return 0
    if not max_err < ABS_TOL_TRACTION:
        print("FAIL: tangent traction error exceeds tolerance")
    if not max_normal_traction < NORMAL_TOL:
        print("FAIL: spurious normal traction in pure mode-II loading")
    return 1


if __name__ == "__main__":
    sys.exit(main())
