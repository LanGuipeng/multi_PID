"""Render manuscript figures from recorded experiment results."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


def configure_style():
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = "Times New Roman" if "Times New Roman" in available else "STIXGeneral"
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [family],
            "mathtext.fontset": "stix",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(figure, output, name, dpi):
    figure.tight_layout()
    figure.savefig(output / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(output / f"{name}.png", dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_switching(record, output, dpi):
    figure, axes = plt.subplots(1, 2, figsize=(7.0, 2.55), sharey=True)
    states = [(row["state"], Fraction(row["probability"])) for row in record["table"]]
    costs = np.linspace(0.001, 5, 2500)
    for target, axis in enumerate(axes):
        prior = sum(p for state, p in states if state[0] == target)
        curves = []
        for source, color in [(1, "#2266a4"), (2, "#d07924")]:
            benefit = -np.maximum(1 - costs, 0)
            for observation in (0, 1):
                a = (
                    sum(
                        p
                        for state, p in states
                        if state[0] == target and state[source] == observation
                    )
                    / prior
                )
                b = sum(
                    p for state, p in states if state[0] != target and state[source] == observation
                ) / (1 - prior)
                benefit += np.maximum(float(a) - costs * float(b), 0)
            curves.append(benefit)
            axis.plot(costs, benefit, color=color, lw=1.5, label=rf"$h_{{{source},t}}$")
        axis.plot(
            costs, np.minimum(*curves), color="black", lw=1.7, ls="--", label=r"$r_{\alpha,t}$"
        )
        crossing = 4 / 3 if target == 0 else 3 / 4
        axis.axvline(crossing, color=".5", lw=0.7, ls=":")
        axis.text(crossing + 0.07, 0.47, r"$4/3$" if target == 0 else r"$3/4$")
        axis.set(
            xlim=(0, 5), ylim=(0, 0.64), xlabel=r"Cost $\lambda$", title=f"Target task $t={target}$"
        )
    axes[0].set_ylabel("Payoff improvement")
    axes[1].legend(frameon=False, loc="upper right")
    save_figure(figure, output, "envelope_switching", dpi)


def antichain_label(masks, n_sources=3):
    sets = [
        r"\{" + "".join(str(i + 1) for i in range(n_sources) if mask & (1 << i)) + r"\}"
        for mask in masks
    ]
    return "$" + r"\,".join(sets) + "$"


def plot_atoms(record, output, dpi):
    atoms = record["atoms"]
    figure, axis = plt.subplots(figsize=(7, 3.2))
    positions = np.arange(len(atoms))
    bottom = np.zeros(len(atoms))
    for target, color in enumerate(["#2266a4", "#d07924", "#508758"]):
        values = np.array([float(row["target_contributions_80d"][str(target)]) for row in atoms])
        axis.bar(positions, values, bottom=bottom, width=0.72, color=color, label=f"$t={target}$")
        bottom += values
    for index, row in enumerate(atoms):
        if row["structural_zero"]:
            axis.plot(index, 0, "x", color="black", ms=5)
    axis.set_xticks(positions)
    axis.set_xticklabels(
        [antichain_label(row["alpha"]) for row in atoms], rotation=75, ha="right", fontsize=8
    )
    axis.set_ylabel("Prior-weighted atom contribution (bit)")
    axis.legend(frameon=False, ncol=3, loc="upper right")
    axis.set_ylim(0, max(bottom) * 1.15)
    axis.text(0.99, 0.72, "x  exact zero", transform=axis.transAxes, ha="right", fontsize=9)
    save_figure(figure, output, "ternary_atoms", dpi)


def plot_precision(record, output, dpi):
    figure, axes = plt.subplots(1, 2, figsize=(7, 2.8))
    histories = [row for row in record["runs"] if row.get("history")]
    if histories:
        history = max(histories, key=lambda row: len(row["history"]))["history"]
        for key, label, color in [
            ("total_width", "Total width", "#222222"),
            ("tail", "Tail bound", "#2266a4"),
            ("endpoint", "Envelope allowance", "#d07924"),
        ]:
            axes[0].loglog(
                [row["cells"] for row in history],
                [float(Fraction(row["errors"][key])) for row in history],
                "-o",
                ms=3,
                label=label,
                color=color,
            )
        axes[0].legend(frameon=False, fontsize=8)
    else:
        axes[0].text(
            0.5, 0.5, "No completed endpoint levels", ha="center", transform=axes[0].transAxes
        )
    axes[0].set(xlabel="Cells in completed level", ylabel="Certified allowance (bit)")
    for method, label, color in [
        ("exact", "Exact segmentation", "#2266a4"),
        ("endpoint", "Endpoint refinement", "#d07924"),
    ]:
        rows = [
            row
            for row in record["runs"]
            if method in row and (method == "exact" or row[method]["tolerance_met"])
        ]
        axes[1].loglog(
            [float(row["epsilon"]) for row in rows],
            [row[method]["timings"]["total"] for row in rows],
            "-o",
            color=color,
            ms=4,
            label=label,
        )
    incomplete = [
        row for row in record["runs"] if "endpoint" in row and not row["endpoint"]["tolerance_met"]
    ]
    if incomplete:
        axes[1].plot(
            [float(row["epsilon"]) for row in incomplete],
            [row["endpoint"]["timings"]["total"] for row in incomplete],
            "x",
            color="#d07924",
            ms=7,
            label="Budget reached",
        )
    axes[1].set(xlabel=r"Requested width $\varepsilon$ (bit)", ylabel="Wall time (s)")
    axes[1].invert_xaxis()
    axes[1].legend(frameon=False, fontsize=8)
    save_figure(figure, output, "precision_cost", dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/reference"))
    parser.add_argument("--output", type=Path, default=Path("docs/figures"))
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    names = ["examples_results.json", "precision_results.json"]
    missing = [name for name in names if not (args.results / name).is_file()]
    if missing:
        parser.error("Missing result files: " + ", ".join(missing))
    examples = json.loads((args.results / names[0]).read_text(encoding="utf-8"))
    precision = json.loads((args.results / names[1]).read_text(encoding="utf-8"))
    args.output.mkdir(parents=True, exist_ok=True)
    configure_style()
    plot_switching(examples["switching"], args.output, args.dpi)
    plot_atoms(examples["ternary"], args.output, args.dpi)
    plot_precision(precision, args.output, args.dpi)
    print(f"Wrote three PDF/PNG figure pairs to {args.output}")


if __name__ == "__main__":
    main()
