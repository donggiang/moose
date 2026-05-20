#!/usr/bin/env python3
"""
Tier 1.4 verification: BiLinearMixedModeTraction at fixed mode-mix beta = 1.

Mixed-mode landmarks (Camanho-Davila with shear strength S = T):

  delta_m_0 = delta_n_0 * delta_s_0 * sqrt((1 + beta^2) /
                                           (delta_s_0^2 + (beta * delta_n_0)^2))

  delta_m_f = (2 / (K * delta_m_0))
              * [GI_c + (GII_c - GI_c) * (beta^2/(1+beta^2))^eta]   (B-K)

For monotonic loading the bilinear law in mixed-mode separation reduces to:

  t_m(delta_m) = K * delta_m                                   0 <= delta_m <= delta_m_0
               = (1 - d) * K * delta_m                         delta_m_0 <= delta_m <= delta_m_f
               = 0                                             delta_m >= delta_m_f
  d(delta_m_max) = delta_m_f * (delta_m_max - delta_m_0)
                 / (delta_m_max * (delta_m_f - delta_m_0))

We extract from the simulation:
  delta_m = sqrt(jump_x^2 + jump_y^2)
  t_m     = sqrt(traction_x^2 + traction_y^2)

and additionally verify that beta stays at 1 (i.e. jump_x ~= jump_y) and that
the per-component tractions are equal in magnitude (since K is the same in
normal and tangent and (1-d) multiplies both equally).

Usage:
    python3 verify_tier1_4.py [path/to/tier1_4_mixed_mode_beta1_out.csv]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference law parameters (must match input file).
K       = 5.0e3
N       = 50.0
S       = 30.0
GI_c    = 0.5
GII_c   = 0.5
eta     = 1.45
beta    = 1.0

delta_n_0 = N / K          # 1.0e-2
delta_s_0 = S / K          # 6.0e-3

# Mixed-mode landmarks.
delta_m_0 = (delta_n_0 * delta_s_0
             * np.sqrt((1.0 + beta**2)
                       / (delta_s_0**2 + (beta * delta_n_0)**2)))

GcMix = GI_c + (GII_c - GI_c) * (beta**2 / (1.0 + beta**2))**eta
delta_m_f = (2.0 / (K * delta_m_0)) * GcMix

t_m_peak = K * delta_m_0   # bilinear peak at delta_m_0

ABS_TOL_TRACTION = 1e-6     # mixed-mode geometry projection slightly looser than 1.1
BETA_TOL         = 1e-6
COMPONENT_TOL    = 1e-6


def damage(delta_max: float) -> float:
    if delta_max <= delta_m_0:
        return 0.0
    if delta_max >= delta_m_f:
        return 1.0
    return delta_m_f * (delta_max - delta_m_0) / (delta_max * (delta_m_f - delta_m_0))


def reference_t_m(delta_m: float) -> float:
    """Bilinear t_m(delta_m) for monotonic mixed-mode loading."""
    if delta_m <= 0:
        return 0.0
    if delta_m <= delta_m_0:
        return K * delta_m
    if delta_m <= delta_m_f:
        return (1.0 - damage(delta_m)) * K * delta_m
    return 0.0


def make_plot(delta_m, t_m, ref, errs, beta_obs, out_path: Path) -> None:
    d_dense = np.linspace(0.0, max(delta_m_f, delta_m.max()) * 1.05, 1000)
    t_dense = np.array([reference_t_m(d) for d in d_dense])

    fig, (ax_curve, ax_err) = plt.subplots(
        2, 1, figsize=(7.8, 7.2), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    ax_curve.plot(d_dense, t_dense, "-", color="tab:blue", lw=1.5,
                  label="analytical")
    ax_curve.plot(delta_m, t_m, "o", color="tab:red", ms=4, mfc="none",
                  label="simulation")
    ax_curve.axvline(delta_m_0, color="0.6", lw=0.8, ls="--")
    ax_curve.axvline(delta_m_f, color="0.6", lw=0.8, ls="--")
    ax_curve.text(delta_m_0, t_m_peak, r" $\delta_{m,0}$",
                  va="bottom", color="0.4")
    ax_curve.text(delta_m_f, 0, r" $\delta_{m,f}$",
                  va="bottom", color="0.4")
    ax_curve.set_ylabel(r"mixed-mode traction $t_m$")
    ax_curve.set_title(
        rf"Tier 1.4 — bilinear mixed-mode at $\beta = {beta:.0f}$"
    )
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper right")

    floor = 1e-16
    ax_err.semilogy(delta_m, np.maximum(errs, floor), "o-", color="tab:purple",
                    ms=3, lw=0.8)
    ax_err.axhline(ABS_TOL_TRACTION, color="tab:red", lw=0.8, ls="--",
                   label=f"tolerance {ABS_TOL_TRACTION:.0e}")
    err_top = max(ABS_TOL_TRACTION * 10,
                  errs.max() * 10 if errs.max() > 0 else ABS_TOL_TRACTION)
    ax_err.set_ylim(floor, err_top)
    ax_err.set_xlabel(r"mixed-mode jump $\delta_m$")
    ax_err.set_ylabel(r"$|t_m^{\mathrm{sim}} - t_m^{\mathrm{ref}}|$")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "tier1_4_mixed_mode_beta1_out.csv"
    )
    if not csv_path.is_file():
        print(f"FAIL: CSV not found at {csv_path}")
        return 2

    times, jx, jy, tx, ty = [], [], [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            jx.append(float(row["jump_x"]))
            jy.append(float(row["jump_y"]))
            tx.append(float(row["traction_x"]))
            ty.append(float(row["traction_y"]))

    times = np.asarray(times)
    jx = np.asarray(jx)
    jy = np.asarray(jy)
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    delta_m = np.sqrt(jx**2 + jy**2)
    t_m     = np.sqrt(tx**2 + ty**2)

    # Observed beta = jump_shear / jump_normal. With our BC choice,
    # tangential = jump_x, normal = jump_y.
    beta_obs = np.where(jy > 1e-14, jx / jy, np.nan)

    # Component magnitudes should be equal (K acts identically on normal and tangent,
    # (1-d) multiplies both equally).
    component_diff = np.abs(np.abs(tx) - np.abs(ty))

    ref = np.array([reference_t_m(d) for d in delta_m])
    errs = np.abs(t_m - ref)
    idx = int(np.argmax(errs))

    plot_path = csv_path.with_name(
        csv_path.stem.replace("_out", "") + "_compare.png"
    )
    make_plot(delta_m, t_m, ref, errs, beta_obs, plot_path)

    print("Tier 1.4 — mixed-mode at beta = 1")
    print(f"  rows checked          : {len(times)}")
    print(f"  delta_m_0 (analytical): {delta_m_0:.6e}")
    print(f"  delta_m_f (analytical): {delta_m_f:.6e}")
    print(f"  peak t_m (analytical) : {t_m_peak:.6e}")

    # Check beta in the meaningful range (skip the very first sample where everything is zero).
    mask = (jy > 1e-6)
    if mask.any():
        beta_min = float(np.nanmin(beta_obs[mask]))
        beta_max = float(np.nanmax(beta_obs[mask]))
        beta_err = max(abs(beta_min - beta), abs(beta_max - beta))
        print(f"  beta observed range   : [{beta_min:.6f}, {beta_max:.6f}]   "
              f"(target {beta:.0f}, max-err {beta_err:.2e}, tol {BETA_TOL:.0e})")
    else:
        beta_err = 0.0
        print("  beta observed range   : insufficient samples to evaluate")

    print(f"  max |tx| - |ty|       : {float(component_diff.max()):.6e}   "
          f"(tol {COMPONENT_TOL:.0e})")
    print(f"  max |tm - tm_ref|     : {float(errs[idx]):.6e}   "
          f"(tol {ABS_TOL_TRACTION:.0e})")
    print(f"    at t = {times[idx]:.4f}, delta_m = {delta_m[idx]:.6e}, "
          f"tm = {t_m[idx]:.6e}, ref = {ref[idx]:.6e}")
    print(f"  comparison plot       : {plot_path}")

    pass_traction  = errs.max() < ABS_TOL_TRACTION
    pass_beta      = beta_err < BETA_TOL
    pass_component = float(component_diff.max()) < COMPONENT_TOL
    if pass_traction and pass_beta and pass_component:
        print("PASS")
        return 0
    if not pass_traction:
        print("FAIL: |tm - tm_ref| exceeds tolerance")
    if not pass_beta:
        print("FAIL: beta drifted from 1 outside tolerance")
    if not pass_component:
        print("FAIL: |t_x| != |t_y| (mixed-mode symmetry violated)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
