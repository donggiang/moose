#!/usr/bin/env python3
"""
Snap-back CZM verification.

Compares the (top_disp, F) trajectory traced by Box 4.4 arc-length against
the analytical bilinear envelope, and overlays the displacement-control
trajectory (which fails at the peak).

Usage:
    python3 verify_snap_back.py
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Material / cohesive parameters (must match snap_back_*.i).
N        = 100.0     # peak cohesive normal traction (MPa)
GIc      = 0.25      # mode-I fracture energy (N/mm)
K_PEN    = 1.0e5     # cohesive penalty stiffness (MPa/mm)

# Derived cohesive landmarks.
DELTA_0  = N / K_PEN          # 1e-3 mm — softening onset
DELTA_F  = 2.0 * GIc / N      # 5e-3 mm — full decohesion

# Reference (q^) traction — matches the AL `ref_load` in input.
SIGMA_REF = 100.0

# Tolerances.
# Pass criterion is "every AL point sits on the analytical envelope".
# Snap-back grid resolution (Δs) means we usually don't land exactly at
# the peak, so we evaluate envelope adherence at the points we did sample.
ABS_TOL_F = 1.0e-6   # MPa — tight, since AL converges to machine precision


def cohesive_traction(delta_n: float) -> float:
    """Bilinear cohesive law in the normal direction."""
    if delta_n <= 0.0:
        return K_PEN * delta_n      # compressive penalty (negative)
    if delta_n <= DELTA_0:
        return K_PEN * delta_n      # elastic penalty branch
    if delta_n <= DELTA_F:
        return N * (DELTA_F - delta_n) / (DELTA_F - DELTA_0)
    return 0.0


def analytical_curve(K_e_total: float):
    """Return arrays (delta_n, F, u_top, lam) tracing the analytical
    bilinear-cohesive + linear-elastic-bulk-in-series envelope."""
    deltas = np.concatenate([
        np.linspace(0.0, DELTA_0, 50),
        np.linspace(DELTA_0, DELTA_F, 200),
    ])
    F   = np.array([cohesive_traction(d) for d in deltas])
    u   = F / K_e_total + deltas
    lam = F / SIGMA_REF
    return deltas, F, u, lam


def read_csv(path: Path):
    if not path.is_file():
        return None
    rows = list(csv.DictReader(path.open()))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}


def main() -> int:
    here = Path(__file__).parent
    al   = read_csv(here / "snap_back_AL_out.csv")
    disp = read_csv(here / "snap_back_disp_out.csv")

    if al is None:
        print("FAIL: snap_back_AL_out.csv not found")
        return 2

    # Estimate K_e_total from the elastic phase of the AL data.
    #   u_top = F / K_e_total + delta_n  ⇒  K_e_total = F / (u_top - delta_n)
    # Use middle of elastic phase to avoid the small initial-condition offset.
    elastic_mask = (al["normal_jump"] > 0) & (al["normal_jump"] < 0.7 * DELTA_0)
    if elastic_mask.sum() < 2:
        # Use the smallest-jump non-zero point.
        elastic_mask = (al["normal_jump"] > 0)
    if elastic_mask.sum() == 0:
        print("FAIL: no elastic-phase points to estimate K_e_total")
        return 1
    F_elastic = al["normal_traction"][elastic_mask]
    elastic_strain = al["top_disp"][elastic_mask] - al["normal_jump"][elastic_mask]
    K_e_total = float(np.mean(F_elastic / elastic_strain))

    # Analytical curve using measured K_e_total.
    d_a, F_a, u_a, lam_a = analytical_curve(K_e_total)

    # Verification metrics.
    F_peak_sim = float(np.max(al["normal_traction"]))
    F_peak_err = abs(F_peak_sim - N) / N

    # Find AL points closest to the analytical envelope and compute mean error.
    al_idx_sort = np.argsort(al["normal_jump"])
    delta_sorted = al["normal_jump"][al_idx_sort]
    F_sorted = al["normal_traction"][al_idx_sort]
    F_ref = np.interp(delta_sorted, d_a, F_a)
    F_err_max = float(np.max(np.abs(F_sorted - F_ref)))

    # Plot.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Left panel: F vs u_top — the snap-back curve.
    ax1.plot(u_a, F_a, "-", color="tab:blue", lw=1.5,
             label="analytical (bilinear CZM + elastic bulk)")
    ax1.plot(al["top_disp"], al["normal_traction"], "o-", color="tab:red",
             ms=5, mfc="white", lw=1.0, label="Box 4.4 arc-length")
    if disp is not None:
        ax1.plot(disp["top_disp"], disp["normal_traction"], "s", color="tab:green",
                 ms=4, mfc="white", label="standard Newton (disp control)")
        # Mark where standard Newton stalled.
        last = -1
        ax1.annotate("Newton stalls at peak",
                     xy=(disp["top_disp"][last], disp["normal_traction"][last]),
                     xytext=(disp["top_disp"][last] + 0.005,
                             disp["normal_traction"][last] - 30),
                     fontsize=9, color="tab:green",
                     arrowprops=dict(arrowstyle="->", color="tab:green", lw=0.8))
    ax1.set_xlabel(r"top displacement $u_{\rm top}$ [mm]")
    ax1.set_ylabel(r"interface traction $t_n$ [MPa]")
    ax1.set_title("Snap-back load–displacement curve")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="upper left", fontsize=9)

    # Right panel: lambda vs delta_n — the AL parameterisation.
    ax2.plot(d_a, lam_a, "-", color="tab:blue", lw=1.5, label="analytical")
    ax2.plot(al["normal_jump"], al["lambda"], "o-", color="tab:red",
             ms=5, mfc="white", lw=1.0, label="Box 4.4 arc-length")
    ax2.axvline(DELTA_0, color="0.6", lw=0.8, ls="--")
    ax2.axvline(DELTA_F, color="0.6", lw=0.8, ls="--")
    ax2.text(DELTA_0, 0.05, r" $\delta_0$", color="0.4")
    ax2.text(DELTA_F, 0.05, r" $\delta_f$", color="0.4")
    ax2.set_xlabel(r"normal jump $\delta_n$ [mm]")
    ax2.set_ylabel(r"load factor $\lambda$")
    ax2.set_title(r"Arc-length parameterisation")
    ax2.grid(alpha=0.3)
    ax2.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    plot_path = here / "snap_back_compare.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    print("Snap-back CZM — Box 4.4 arc-length vs analytical")
    print(f"  rows AL                : {len(al['time'])}")
    print(f"  K_e_total (measured)   : {K_e_total:.4e} MPa/mm")
    print(f"  delta_0  (analytical)  : {DELTA_0:.4e} mm")
    print(f"  delta_f  (analytical)  : {DELTA_F:.4e} mm")
    print(f"  F_peak   (analytical)  : {N:.4e} MPa")
    print(f"  F_peak_sim             : {F_peak_sim:.4e} MPa")
    print(f"  rel(F_peak)            : {F_peak_err*100:.3f} %")
    print(f"  max |F_sim - F_ref|    : {F_err_max:.4e} MPa")
    print(f"  comparison plot        : {plot_path}")

    snap_back_seen = al["lambda"].size >= 7 and (
        al["lambda"][6] > al["lambda"][7] and
        al["top_disp"][6] > al["top_disp"][7]
    )
    print(f"  snap-back observed     : {snap_back_seen}")

    pass_envelope = F_err_max < ABS_TOL_F
    if pass_envelope and snap_back_seen:
        print("PASS")
        return 0
    if not pass_envelope:
        print(f"FAIL: max envelope error {F_err_max:.3e} MPa exceeds {ABS_TOL_F:.0e}")
    if not snap_back_seen:
        print("FAIL: snap-back (decreasing u and λ) not observed in trajectory")
    return 1


if __name__ == "__main__":
    sys.exit(main())
