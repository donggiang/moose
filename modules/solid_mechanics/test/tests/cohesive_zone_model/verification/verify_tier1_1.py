#!/usr/bin/env python3
"""
Tier 1.1 verification: compare BiLinearMixedModeTraction against the closed-form
bilinear t_n(delta_n) curve.

Usage:
    python3 verify_tier1_1.py [path/to/tier1_1_mode_I_monotonic_out.csv]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference law parameters (must match the input file).
K       = 5.0e3   # penalty_stiffness — chosen so delta_0 = 0.01 (softening at u=0.01)
N       = 50.0    # normal_strength
GI_c    = 0.5     # mode-I fracture energy

# Analytical landmarks.
delta_0 = N / K                    # 1.0e-2 (softening onset)
delta_f = 2.0 * GI_c / N           # 2.0e-2 (full decohesion)

# Pass tolerance — absolute on traction.
ABS_TOL = 1e-7

# Tangent-traction tolerance — must remain numerical noise (no shear loading).
TANGENT_TOL = 1e-10


def analytical_traction(delta_n: float) -> float:
    """Closed-form bilinear t_n(delta_n) for monotonic loading."""
    if delta_n <= 0:
        return 0.0
    if delta_n <= delta_0:
        return K * delta_n
    if delta_n <= delta_f:
        return N * (delta_f - delta_n) / (delta_f - delta_0)
    return 0.0


def make_plot(deltas, tractions, tn_ref, errs, out_path: Path) -> None:
    # Dense analytical reference curve covering the full simulated range.
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
    ax_curve.text(delta_0, N, r" $\delta_0$", va="bottom", color="0.4")
    ax_curve.text(delta_f, 0,   r" $\delta_f$", va="bottom", color="0.4")
    ax_curve.set_ylabel(r"normal traction $t_n$")
    ax_curve.set_title("Tier 1.1 — bilinear mode-I traction–separation")
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper right")

    floor = 1e-16
    ax_err.semilogy(deltas, np.maximum(errs, floor), "o-", color="tab:purple",
                    ms=3, lw=0.8)
    ax_err.axhline(ABS_TOL, color="tab:red", lw=0.8, ls="--",
                   label=f"tolerance {ABS_TOL:.0e}")
    ax_err.set_ylim(floor, max(ABS_TOL * 10, errs.max() * 10 if errs.max() > 0 else 1e-10))
    ax_err.set_xlabel(r"normal jump $\delta_n$")
    ax_err.set_ylabel(r"$|t_n^{\mathrm{sim}} - t_n^{\mathrm{ref}}|$")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "tier1_1_mode_I_monotonic_out.csv"
    )
    if not csv_path.is_file():
        print(f"FAIL: CSV not found at {csv_path}")
        return 2

    times, deltas, tractions, tangents = [], [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            deltas.append(float(row["normal_jump"]))
            tractions.append(float(row["normal_traction"]))
            tangents.append(float(row["tangent_traction"]))

    times = np.asarray(times)
    deltas = np.asarray(deltas)
    tractions = np.asarray(tractions)
    tangents = np.asarray(tangents)

    tn_ref = np.array([analytical_traction(d) for d in deltas])
    errs = np.abs(tractions - tn_ref)
    idx = int(np.argmax(errs))
    max_err_traction = float(errs[idx])
    max_err_at = (times[idx], deltas[idx], tractions[idx], tn_ref[idx])
    max_tangent = float(np.max(np.abs(tangents)))
    rows_checked = len(times)

    plot_path = csv_path.with_name(csv_path.stem.replace("_out", "") + "_compare.png")
    make_plot(deltas, tractions, tn_ref, errs, plot_path)

    print("Tier 1.1 — mode-I monotonic pull")
    print(f"  rows checked         : {rows_checked}")
    print(f"  delta_0 (analytical) : {delta_0:.6e}")
    print(f"  delta_f (analytical) : {delta_f:.6e}")
    print(f"  peak traction        : {N:.6e}")
    print(f"  max |tn - tn_ref|    : {max_err_traction:.6e}   (tol {ABS_TOL:.0e})")
    if max_err_at:
        t, dn, tn, tn_ref = max_err_at
        print(f"    at t = {t:.4f}, delta_n = {dn:.6e}, tn = {tn:.6e}, tn_ref = {tn_ref:.6e}")
    print(f"  max |tangent_traction| : {max_tangent:.6e}   (tol {TANGENT_TOL:.0e})")
    print(f"  comparison plot      : {plot_path}")

    pass_traction = max_err_traction < ABS_TOL
    pass_tangent  = max_tangent      < TANGENT_TOL
    if pass_traction and pass_tangent:
        print("PASS")
        return 0
    if not pass_traction:
        print("FAIL: traction error exceeds tolerance")
    if not pass_tangent:
        print("FAIL: tangent traction is non-zero (mode-I should have no shear)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
