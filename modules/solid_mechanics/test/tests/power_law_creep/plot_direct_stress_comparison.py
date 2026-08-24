#!/usr/bin/env python3

import argparse
import csv
import statistics

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def read_history(filename, strain_column):
    with open(filename, newline="", encoding="utf-8") as csv_file:
        rows = list(csv.DictReader(csv_file))
    return {
        "time": [float(row["time"]) for row in rows],
        "strain": [float(row[strain_column]) for row in rows],
        "iterations": [float(row["local_iterations_average"]) for row in rows],
    }


def read_timings(filename):
    with open(filename, encoding="utf-8") as timing_file:
        return [float(value) for value in timing_file if value.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_csv")
    parser.add_argument("direct_csv")
    parser.add_argument("reference_timings")
    parser.add_argument("direct_timings")
    parser.add_argument("output_prefix")
    parser.add_argument("--strain-column", default="creep_strain_yy_average")
    args = parser.parse_args()

    reference = read_history(args.reference_csv, args.strain_column)
    direct = read_history(args.direct_csv, args.strain_column)
    if reference["time"] != direct["time"]:
        raise ValueError("The comparison histories use different time points")

    ratios = []
    relative_differences = []
    for reference_value, direct_value in zip(reference["strain"], direct["strain"]):
        if reference_value == 0.0:
            ratios.append(1.0 if direct_value == 0.0 else float("nan"))
            relative_differences.append(0.0 if direct_value == 0.0 else float("inf"))
        else:
            ratios.append(direct_value / reference_value)
            relative_differences.append((direct_value - reference_value) / reference_value)

    figure, axes = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True)
    axes[0].plot(reference["time"], reference["strain"], "o-", label="MOOSE radial return")
    axes[0].plot(direct["time"], direct["strain"], "s--", label="Direct stress")
    axes[0].set_ylabel(r"$\epsilon^{vp}_{yy}$")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(reference["time"], ratios, "o-")
    axes[1].axhline(1.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("Direct / reference")
    axes[1].yaxis.set_major_formatter(FormatStrFormatter("%.11f"))
    axes[1].grid(True, alpha=0.3)

    axes[2].semilogy(
        reference["time"],
        [max(abs(value), 1e-18) for value in relative_differences],
        "o-",
    )
    axes[2].set_xlabel("Time")
    axes[2].set_ylabel("Absolute relative difference")
    axes[2].grid(True, which="both", alpha=0.3)
    figure.tight_layout()
    figure.savefig(args.output_prefix + ".pdf")
    figure.savefig(args.output_prefix + ".png", dpi=200)

    with open(args.output_prefix + ".csv", "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                "time",
                "reference_creep_strain_yy",
                "direct_creep_strain_yy",
                "direct_over_reference",
                "relative_difference",
                "reference_local_iterations",
                "direct_local_iterations",
            ]
        )
        for values in zip(
            reference["time"],
            reference["strain"],
            direct["strain"],
            ratios,
            relative_differences,
            reference["iterations"],
            direct["iterations"],
        ):
            writer.writerow(values)

    reference_timings = read_timings(args.reference_timings)
    direct_timings = read_timings(args.direct_timings)
    reference_total_iterations = sum(reference["iterations"])
    direct_total_iterations = sum(direct["iterations"])
    maximum_relative_difference = max(abs(value) for value in relative_differences)
    with open(args.output_prefix + "_summary.txt", "w", encoding="utf-8") as summary:
        summary.write(f"Reference total converged-state local iterations: {reference_total_iterations:g}\n")
        summary.write(f"Direct total converged-state local iterations: {direct_total_iterations:g}\n")
        summary.write(
            f"Reference wall times (s): {', '.join(f'{value:.2f}' for value in reference_timings)}\n"
        )
        summary.write(
            f"Direct wall times (s): {', '.join(f'{value:.2f}' for value in direct_timings)}\n"
        )
        summary.write(f"Reference median wall time (s): {statistics.median(reference_timings):.3f}\n")
        summary.write(f"Direct median wall time (s): {statistics.median(direct_timings):.3f}\n")
        summary.write(
            "Direct/reference median wall-time ratio: "
            f"{statistics.median(direct_timings) / statistics.median(reference_timings):.3f}\n"
        )
        summary.write(f"Maximum absolute relative strain difference: {maximum_relative_difference:.6e}\n")


if __name__ == "__main__":
    main()
