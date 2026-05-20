#!/usr/bin/env python3
"""
Tier 1.3 verification: BiLinearMixedModeTraction under load / unload / reload.

Reference law for monotonically growing damage:

  d(delta) = delta_f * (delta - delta_0) / (delta * (delta_f - delta_0))    (clamped to [0,1])

For our load history we sweep delta_n^max over time. The closed-form traction
at sample k is:

  delta^max_k = max_{j<=k} delta_n_j     (irreversible internal variable)
  d_k         = d(delta^max_k)
  t_n_ref_k   = (1 - d_k) * K * delta_n_k          if delta_n_k >= 0

When delta_n_k = delta^max_k, this reduces to the bilinear envelope; otherwise
it is the secant unload/reload through the origin with frozen damage.

Usage:
    python3 verify_tier1_3.py [path/to/tier1_3_mode_I_load_unload_out.csv]
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
N       = 50.0
GI_c    = 0.5

delta_0 = N / K            # 1.0e-2
delta_f = 2.0 * GI_c / N   # 2.0e-2

# Pass tolerances.
ABS_TOL_TRACTION = 1e-7
TANGENT_TOL      = 1e-10


def damage(delta_max: float) -> float:
    if delta_max <= delta_0:
        return 0.0
    if delta_max >= delta_f:
        return 1.0
    return delta_f * (delta_max - delta_0) / (delta_max * (delta_f - delta_0))


def reference_traction(delta_n: float, delta_max: float) -> float:
    if delta_n <= 0.0:
        return 0.0
    d = damage(delta_max)
    return (1.0 - d) * K * delta_n


def make_plot(deltas, tractions, ref, errs, d_max_series, out_path: Path) -> None:
    # Bilinear envelope.
    d_dense = np.linspace(0.0, max(delta_f, deltas.max()) * 1.05, 1000)
    t_dense = np.array([reference_traction(d, d) for d in d_dense])

    fig, (ax_curve, ax_err) = plt.subplots(
        2, 1, figsize=(7.8, 7.2), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    # Envelope (solid blue), reference traction at every step (small dots), simulation (rings).
    ax_curve.plot(d_dense, t_dense, "-", color="tab:blue", lw=1.2,
                  label="bilinear envelope")
    ax_curve.plot(deltas, ref, ".", color="tab:green", ms=3,
                  label="analytical (with frozen $d$)")
    ax_curve.plot(deltas, tractions, "o", color="tab:red", ms=4, mfc="none",
                  label="simulation")

    ax_curve.axvline(delta_0, color="0.6", lw=0.8, ls="--")
    ax_curve.axvline(delta_f, color="0.6", lw=0.8, ls="--")
    ax_curve.text(delta_0, N, r" $\delta_0$", va="bottom", color="0.4")
    ax_curve.text(delta_f, 0, r" $\delta_f$", va="bottom", color="0.4")
    ax_curve.set_ylabel(r"normal traction $t_n$")
    ax_curve.set_title("Tier 1.3 — load / unload / reload")
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
    ax_err.set_xlabel(r"normal jump $\delta_n$")
    ax_err.set_ylabel(r"$|t_n^{\mathrm{sim}} - t_n^{\mathrm{ref}}|$")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    # Also a time-history plot so unload/reload is visually obvious.
    th_path = out_path.with_name(out_path.stem.replace("_compare", "") + "_history.png")
    fig2, (ax_d, ax_t) = plt.subplots(2, 1, figsize=(7.8, 6), sharex=True)
    times = np.linspace(0, 1, len(deltas))  # placeholder; caller may pass real times
    # We don't actually have time here without re-reading; skip in this helper.
    plt.close(fig2)


def make_history_plot(times, deltas, tractions, ref, d_max_series, out_path: Path) -> None:
    fig, (ax_d, ax_t) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    ax_d.plot(times, deltas, "-", color="tab:blue", lw=1.2, label=r"$\delta_n$ (sim)")
    ax_d.plot(times, d_max_series, "--", color="tab:orange", lw=1.0,
              label=r"$\delta_n^{\max}$ (running)")
    ax_d.set_ylabel("normal jump")
    ax_d.grid(alpha=0.3)
    ax_d.legend(loc="upper right", fontsize=9)
    ax_d.set_title("Tier 1.3 — time history")

    ax_t.plot(times, tractions, "-", color="tab:red", lw=1.2, label="simulation")
    ax_t.plot(times, ref, "--", color="tab:green", lw=1.0, label="analytical")
    ax_t.set_xlabel("time")
    ax_t.set_ylabel(r"normal traction $t_n$")
    ax_t.grid(alpha=0.3)
    ax_t.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "tier1_3_mode_I_load_unload_out.csv"
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

    # Build the running delta_n^max (irreversibility of damage).
    d_max_series = np.maximum.accumulate(np.maximum(deltas, 0.0))

    # Closed-form traction at every step using the running delta_max.
    ref = np.array([reference_traction(d, dmax)
                    for d, dmax in zip(deltas, d_max_series)])
    errs = np.abs(tractions - ref)
    idx = int(np.argmax(errs))
    max_err = float(errs[idx])
    max_tangent = float(np.max(np.abs(tangents)))

    # Sanity check: damage is non-decreasing (delta_max is monotone by construction;
    # confirm the sim didn't somehow produce a value outside the secant either).
    monotone = bool(np.all(np.diff(d_max_series) >= -1e-15))

    plot_path = csv_path.with_name(
        csv_path.stem.replace("_out", "") + "_compare.png"
    )
    history_path = csv_path.with_name(
        csv_path.stem.replace("_out", "") + "_history.png"
    )
    make_plot(deltas, tractions, ref, errs, d_max_series, plot_path)
    make_history_plot(times, deltas, tractions, ref, d_max_series, history_path)

    print("Tier 1.3 — mode-I load / unload / reload")
    print(f"  rows checked          : {len(times)}")
    print(f"  delta_0 (analytical)  : {delta_0:.6e}")
    print(f"  delta_f (analytical)  : {delta_f:.6e}")
    print(f"  peak delta_n reached  : {deltas.max():.6e}")
    print(f"  d at peak             : {damage(d_max_series.max()):.6f}")
    print(f"  max |tn - tn_ref|     : {max_err:.6e}   (tol {ABS_TOL_TRACTION:.0e})")
    t, dn, tn, tn_ref = times[idx], deltas[idx], tractions[idx], ref[idx]
    print(f"    at t = {t:.4f}, delta_n = {dn:.6e}, tn = {tn:.6e}, tn_ref = {tn_ref:.6e}")
    print(f"  max |tangent_traction| : {max_tangent:.6e}   (tol {TANGENT_TOL:.0e})")
    print(f"  delta_max monotone     : {monotone}")
    print(f"  comparison plot        : {plot_path}")
    print(f"  history plot           : {history_path}")

    if (max_err < ABS_TOL_TRACTION
            and max_tangent < TANGENT_TOL
            and monotone):
        print("PASS")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
