#!/usr/bin/env python3
"""
Tier 3.1 verification: J-integral / dissipated-energy conservation.

For a bilinear cohesive law driven monotonically past delta_n_f, the area
under the t_n(delta_n) curve equals the input fracture energy:

    integral_0^{delta_n_f} t_n(delta) ddelta  =  (1/2) * N * delta_n_f  =  G_Ic

Per Rice's J-integral identity in the small-strain limit, this is the
energy released per unit interface area at full decohesion.

Two estimators of the cumulative work per unit area are computed:

  1. Trapezoidal integral of (t_n) * d(delta_n) along the simulation
     trajectory.  Numerical, has O(dt) error at the bilinear kink.

  2. Closed-form W_diss(delta_max) at each step using the running
     delta_n^max history.  This must equal G_Ic exactly when
     delta_max >= delta_n_f.

Both are then multiplied by the interface measure (length in 2D) to give
total work, and compared to G_Ic * A_iface.

Usage:
    python3 verify_tier3_1.py [path/to/tier3_1_jintegral_out.csv]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference law parameters.
K       = 5.0e3
N       = 50.0
GI_c    = 0.5

delta_0 = N / K            # 1.0e-2
delta_f = 2.0 * GI_c / N   # 2.0e-2

# Acceptance tolerances.
ENERGY_TOL_EXACT = 1e-9    # closed-form must hit G_Ic exactly (machine precision)
ENERGY_TOL_TRAP  = 1e-3    # trap-rule allowed kink artifact ~O(dt)


def damage(delta_max: float) -> float:
    if delta_max <= delta_0:
        return 0.0
    if delta_max >= delta_f:
        return 1.0
    return delta_f * (delta_max - delta_0) / (delta_max * (delta_f - delta_0))


def analytical_dissipated(delta_max: float) -> float:
    """W_diss per unit interface area as a function of running delta_max."""
    if delta_max <= delta_0:
        return 0.0
    if delta_max >= delta_f:
        return GI_c
    peak = K * delta_0
    W_env = 0.5 * K * delta_0 * delta_0
    W_env += peak * (delta_f * (delta_max - delta_0)
                     - 0.5 * (delta_max**2 - delta_0**2)) / (delta_f - delta_0)
    d = damage(delta_max)
    U = 0.5 * (1.0 - d) * K * delta_max**2
    return W_env - U


def make_plot(times, deltas, tractions,
              W_trap_total, W_exact_total, A_iface,
              out_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=False)

    # Top: traction-separation curve with envelope (sanity).
    ax_curve = axes[0]
    d_dense = np.linspace(0.0, max(delta_f, deltas.max()) * 1.05, 500)
    t_dense = np.array([(1.0 - damage(d)) * K * d if d > 0 else 0.0 for d in d_dense])
    ax_curve.plot(d_dense, t_dense, "-", color="tab:blue", lw=1.2, label="bilinear envelope")
    ax_curve.plot(deltas, tractions, "o", color="tab:red", ms=4, mfc="none",
                  label="simulation")
    ax_curve.axvline(delta_0, color="0.6", lw=0.8, ls="--")
    ax_curve.axvline(delta_f, color="0.6", lw=0.8, ls="--")
    ax_curve.set_xlabel(r"normal jump $\delta_n$")
    ax_curve.set_ylabel(r"normal traction $t_n$")
    ax_curve.set_title("Tier 3.1 — bilinear traction (top), cumulative work (bottom)")
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper right")

    # Bottom: cumulative work vs time.
    ax_w = axes[1]
    target = GI_c * A_iface
    ax_w.plot(times, W_trap_total, "-", color="tab:purple", lw=1.4,
              label="W (trap rule)")
    ax_w.plot(times, W_exact_total, "--", color="tab:green", lw=1.4,
              label="W (closed form)")
    ax_w.axhline(target, color="tab:red", lw=0.8, ls="--",
                 label=f"$G_{{Ic}} \\cdot A = {target}$")
    ax_w.set_xlabel("time")
    ax_w.set_ylabel("cumulative work")
    ax_w.grid(alpha=0.3)
    ax_w.legend(loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).parent / "tier3_1_jintegral_out.csv"
    )
    if not csv_path.is_file():
        print(f"FAIL: CSV not found at {csv_path}")
        return 2

    times, deltas, tractions, iface_lengths = [], [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            deltas.append(float(row["normal_jump"]))
            tractions.append(float(row["normal_traction"]))
            iface_lengths.append(float(row["iface_length"]))

    times         = np.asarray(times)
    deltas        = np.asarray(deltas)
    tractions     = np.asarray(tractions)
    iface_lengths = np.asarray(iface_lengths)

    A_iface = float(iface_lengths[-1])     # constant in small strain

    # Trap-rule work per unit area.
    W_trap_per_area = np.zeros_like(times)
    for k in range(1, len(times)):
        W_trap_per_area[k] = (W_trap_per_area[k-1]
            + 0.5 * (tractions[k] + tractions[k-1])
                  * (deltas[k] - deltas[k-1]))

    # The traj is monotone here, so all elastic energy is recovered as the law
    # transitions to d = 1; trap rule integrates the area under the bilinear
    # curve directly.

    # Closed-form work per unit area at each step (uses running delta_max).
    d_max_series = np.maximum.accumulate(np.maximum(deltas, 0.0))
    W_exact_per_area = np.array([analytical_dissipated(dm) for dm in d_max_series])

    # Multiply by interface measure for total work.
    W_trap_total  = W_trap_per_area  * A_iface
    W_exact_total = W_exact_per_area * A_iface

    target = GI_c * A_iface

    final_trap  = float(W_trap_total[-1])
    final_exact = float(W_exact_total[-1])
    err_trap    = abs(final_trap  - target)
    err_exact   = abs(final_exact - target)

    plot_path = csv_path.with_name(
        csv_path.stem.replace("_out", "") + "_compare.png"
    )
    make_plot(times, deltas, tractions, W_trap_total, W_exact_total,
              A_iface, plot_path)

    print("Tier 3.1 — J-integral / dissipated-energy conservation")
    print(f"  rows checked          : {len(times)}")
    print(f"  interface length A    : {A_iface}")
    print(f"  target (G_Ic * A)     : {target}")
    print(f"  W (closed form, end)  : {final_exact:.10e}   |err| {err_exact:.3e}   "
          f"(tol {ENERGY_TOL_EXACT:.0e})")
    print(f"  W (trap rule, end)    : {final_trap:.10e}    |err| {err_trap:.3e}   "
          f"(tol {ENERGY_TOL_TRAP:.0e})")
    print(f"  comparison plot       : {plot_path}")

    pass_exact = err_exact < ENERGY_TOL_EXACT
    pass_trap  = err_trap  < ENERGY_TOL_TRAP
    if pass_exact and pass_trap:
        print("PASS")
        return 0
    print("FAIL")
    if not pass_exact:
        print(f"  - closed-form work != G_Ic*A within {ENERGY_TOL_EXACT:.0e}")
    if not pass_trap:
        print(f"  - trapezoidal work deviates more than {ENERGY_TOL_TRAP:.0e}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
