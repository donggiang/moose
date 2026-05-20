#!/usr/bin/env python3
"""
Tier 2.0 — penalty stiffness K conditioning sweep.

Runs `tier2_0_K_conditioning.i` with several values of penalty_stiffness on
AS4/PEEK material parameters (Turon 2006, Table 1). For each K it measures:

  1. whether the simulation converged at all
  2. max |t_n - bilinear(t_n; K, N, GIc)|   — pointwise law agreement
  3. spurious-oscillation amplitude near the kink, normalized by N
  4. peak normal traction reached (should be ~ N if cohesive activated)
  5. cohesive-vs-bulk compliance ratio (diagnostic only)

Then it picks the largest K that:
  - converged everywhere
  - reproduced the bilinear law to < 1 % of N (= 0.8 MPa)
  - had < 5 % oscillation near the kink

Usage:
    python3 run_tier2_0_sweep.py
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
INPUT = HERE / "tier2_0_K_conditioning.i"
SOLVER = HERE.parents[3] / "solid_mechanics-opt"

# AS4/PEEK material constants (Turon 2006 Table 1).
N    = 80.0          # normal_strength, MPa
GIc  = 0.969         # mode-I fracture energy, N/mm
delta_f = 2.0 * GIc / N    # = 0.024225 mm

K_MIN_VALID = N**2 / (2.0 * GIc)  # below this, delta_n_0 > delta_n_f → invalid law

# Sweep values (N/mm^3). All > K_MIN_VALID.
K_VALUES = [1e4, 3e4, 1e5, 3e5, 1e6]

# Acceptance thresholds.
ABS_TOL_TRACTION = 0.01 * N   # 1 % of peak strength
OSC_TOL          = 0.05 * N   # 5 % of peak strength
PEAK_FRAC_MIN    = 0.95       # measured peak must reach 95 % of N to count as activated


@dataclass
class KResult:
    K: float
    converged: bool
    delta_n_0: float
    rows: int
    max_traction_reached: float
    max_traction_err: float
    osc_amp_near_kink: float
    final_normal_jump: float
    csv_path: Path


def analytical_traction(delta_n: float, K: float) -> float:
    delta_0 = N / K
    if delta_n <= 0:
        return 0.0
    if delta_n <= delta_0:
        return K * delta_n
    if delta_n <= delta_f:
        return N * (delta_f - delta_n) / (delta_f - delta_0)
    return 0.0


def run_one_K(K: float) -> KResult:
    file_base = f"tier2_0_K_{K:.0e}".replace("+", "")
    # When Outputs/file_base is set explicitly, MOOSE writes <file_base>.csv
    # (not <file_base>_out.csv).
    csv_path = HERE / f"{file_base}.csv"
    # Clean any stale output.
    csv_path.unlink(missing_ok=True)

    cmd = [
        str(SOLVER),
        "-i", str(INPUT),
        f"Materials/czm/penalty_stiffness={K}",
        f"Outputs/file_base={file_base}",
    ]
    print(f"[K={K:.0e}] running...", flush=True)
    proc = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True, timeout=600)
    converged = (proc.returncode == 0) and csv_path.is_file()
    if not converged:
        print(f"[K={K:.0e}] FAILED — return code {proc.returncode}")
        print("--- stderr tail ---")
        print(proc.stderr[-1500:])
        return KResult(K, False, N/K, 0, 0.0, math.inf, math.inf, 0.0, csv_path)

    # Read CSV.
    times, deltas, tractions = [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["time"]))
            deltas.append(float(row["normal_jump"]))
            tractions.append(float(row["normal_traction"]))
    times     = np.asarray(times)
    deltas    = np.asarray(deltas)
    tractions = np.asarray(tractions)

    # Pointwise law agreement.
    refs  = np.array([analytical_traction(d, K) for d in deltas])
    errs  = np.abs(tractions - refs)
    max_err = float(errs.max())
    max_traction = float(tractions.max())
    final_jump = float(deltas[-1])

    # Oscillation amplitude near the elastic-softening kink.
    # Measure local peak-to-trough deviation from the analytical envelope
    # in a window of ±20 % of the bilinear range around delta_n_0.
    delta_0 = N / K
    window = 0.2 * (delta_f - delta_0)
    near_kink = np.abs(deltas - delta_0) < window
    if near_kink.any() and near_kink.sum() > 2:
        residuals = tractions[near_kink] - refs[near_kink]
        osc_amp = float(residuals.max() - residuals.min())
    else:
        osc_amp = 0.0

    return KResult(K, True, delta_0, len(times),
                   max_traction, max_err, osc_amp,
                   final_jump, csv_path)


def make_plot(results: List[KResult], out_path: Path) -> None:
    fig, (ax_curve, ax_err) = plt.subplots(
        2, 1, figsize=(8, 7.2), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    # Bilinear envelope is K-dependent (different δ_n^0 per K). Draw envelope
    # for each K. Use a colour gradient so K-evolution is visible.
    cmap = plt.colormaps["viridis"]
    for i, r in enumerate(results):
        if not r.converged:
            continue
        color = cmap(i / max(1, len(results) - 1))
        # Reload CSV for plotting.
        deltas, tractions = [], []
        with r.csv_path.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                deltas.append(float(row["normal_jump"]))
                tractions.append(float(row["normal_traction"]))
        deltas    = np.asarray(deltas)
        tractions = np.asarray(tractions)

        # Sim trace.
        ax_curve.plot(deltas, tractions, "-", color=color, lw=1.2,
                      label=f"K = {r.K:.0e}")
        # Per-K analytical envelope (sample at the same delta values for the
        # error subplot).
        refs = np.array([analytical_traction(d, r.K) for d in deltas])
        errs = np.abs(tractions - refs)
        floor = 1e-16
        ax_err.semilogy(deltas, np.maximum(errs, floor), "-",
                        color=color, lw=1.0)

    # Reference: any analytical envelope (use highest K so kink is leftmost).
    K_max = max(r.K for r in results if r.converged)
    d_dense = np.linspace(0.0, 1.05 * delta_f, 1000)
    t_dense = np.array([analytical_traction(d, K_max) for d in d_dense])
    ax_curve.plot(d_dense, t_dense, "k--", lw=0.8, alpha=0.5,
                  label=f"analytical (K={K_max:.0e})")
    ax_curve.axvline(delta_f, color="0.5", lw=0.8, ls=":")
    ax_curve.text(delta_f, 0, r" $\delta_f$", va="bottom", color="0.4")
    ax_curve.set_ylabel(r"normal traction $t_n$ (MPa)")
    ax_curve.set_title("Tier 2.0 — K-conditioning sweep on AS4/PEEK single element")
    ax_curve.grid(alpha=0.3)
    ax_curve.legend(loc="upper right", fontsize=8)
    ax_curve.set_ylim(-2, 1.25 * N)

    ax_err.axhline(ABS_TOL_TRACTION, color="tab:red", lw=0.8, ls="--",
                   label=f"tol {ABS_TOL_TRACTION:.2g}")
    ax_err.set_xlabel(r"normal jump $\delta_n$ (mm)")
    ax_err.set_ylabel(r"$|t_n^{\rm sim} - t_n^{\rm ref}|$")
    ax_err.grid(alpha=0.3, which="both")
    ax_err.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def pick_best_K(results: List[KResult]) -> Optional[float]:
    """Return the largest K that meets all acceptance criteria."""
    accepted = [r for r in results
                if r.converged
                and r.max_traction_reached >= PEAK_FRAC_MIN * N
                and r.max_traction_err < ABS_TOL_TRACTION
                and r.osc_amp_near_kink < OSC_TOL]
    if not accepted:
        return None
    return max(r.K for r in accepted)


def main() -> int:
    if not SOLVER.exists():
        print(f"FAIL: solid_mechanics-opt not found at {SOLVER}")
        return 2
    if not INPUT.is_file():
        print(f"FAIL: input not found at {INPUT}")
        return 2

    print(f"K_min_valid (K > N^2/(2 GIc)) = {K_MIN_VALID:.0f}\n")

    results: List[KResult] = []
    for K in K_VALUES:
        if K <= K_MIN_VALID:
            print(f"[K={K:.0e}] skipped — below validity threshold {K_MIN_VALID:.0f}")
            continue
        results.append(run_one_K(K))

    print()
    print(f"{'K':>10s} {'conv':>5s} {'delta_n_0':>11s} {'rows':>5s} "
          f"{'peak t_n':>10s} {'max err':>10s} {'osc/N':>8s}")
    for r in results:
        peak = r.max_traction_reached
        print(f"{r.K:>10.0e} {str(r.converged):>5s} {r.delta_n_0:>11.3e} "
              f"{r.rows:>5d} {peak:>10.3f} {r.max_traction_err:>10.3e} "
              f"{r.osc_amp_near_kink/N:>8.3f}")

    plot_path = HERE / "tier2_0_K_sweep.png"
    make_plot(results, plot_path)
    print(f"\nplot: {plot_path}")

    best = pick_best_K(results)
    print()
    if best is not None:
        print(f"RECOMMENDED K = {best:.0e}")
        print("(largest K that converged, reproduced bilinear within "
              f"{ABS_TOL_TRACTION:.2g}, and has oscillation < "
              f"{OSC_TOL/N*100:.0f}% of N near the kink)")
        return 0
    else:
        print("FAIL: no K met all acceptance criteria.")
        print("- check that at least one K activated softening (peak ~ N)")
        print("- consider relaxing tolerances or extending K range")
        return 1


if __name__ == "__main__":
    sys.exit(main())
