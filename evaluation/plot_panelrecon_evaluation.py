#!/usr/bin/env python3
"""Plot PanelRecon evaluation results from TSV files."""

import argparse
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(
        description="Create plots from PanelRecon evaluation TSV outputs."
    )
    parser.add_argument(
        "--results_dir",
        type=Path,
        default=REPO_ROOT / "evaluation" / "read_depth_results",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
    )
    args = parser.parse_args()

    results_dir = args.results_dir
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = results_dir / "plots"

    run_metrics_path = results_dir / "run_metrics.tsv"
    sample_results_path = results_dir / "sample_results.tsv"

    missing_inputs = []
    if not run_metrics_path.exists():
        missing_inputs.append(str(run_metrics_path))
    if not sample_results_path.exists():
        missing_inputs.append(str(sample_results_path))

    if missing_inputs:
        print("ERROR: missing required input file(s):", file=sys.stderr)
        for missing_input in missing_inputs:
            print(f"  {missing_input}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_cache_dir = output_dir / "matplotlib_cache"
    matplotlib_cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache_dir))

    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        import pandas as pd
        import seaborn as sns
    except ImportError as error:
        print(f"ERROR: missing plotting dependency: {error}", file=sys.stderr)
        print("Install with: pip install pandas seaborn matplotlib", file=sys.stderr)
        return 1

    sns.set_theme(style="whitegrid")

    panel_plot_labels = {
        "3528311_Covered": "GenOncology-Dx_v1",
        "3528311_Covered.2bit": "GenOncology-Dx_v1",
        "3549481_Covered": "GenOncology-Dx_v1.6",
        "3549481_Covered.2bit": "GenOncology-Dx_v1.6",
    }

    run_metrics = pd.read_csv(run_metrics_path, sep="\t")
    sample_results = pd.read_csv(sample_results_path, sep="\t")

    for column in [
        "read_count",
        "accuracy",
        "overall_runtime_seconds",
        "index_runtime_seconds",
        "find_runtime_seconds",
        "peak_rss_mb",
    ]:
        if column in run_metrics.columns:
            run_metrics[column] = pd.to_numeric(run_metrics[column], errors="coerce")

    for column in [
        "read_count",
        "is_correct",
        "best_score",
        "best_panel_covered_kmers_pct",
    ]:
        if column in sample_results.columns:
            sample_results[column] = pd.to_numeric(sample_results[column], errors="coerce")

    run_metrics = run_metrics.sort_values("read_count")
    sample_results = sample_results.sort_values(["sample", "read_count"])

    accuracy_plot = output_dir / "accuracy_by_read_count.png"
    plt.figure(figsize=(7, 4.5))
    sns.lineplot(
        data=run_metrics,
        x="read_count",
        y="accuracy",
        marker="o",
        linewidth=2,
    )
    plt.title("Identification accuracy by N examined reads")
    plt.xlabel("N reads scanned")
    plt.ylabel("Accuracy")
    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig(accuracy_plot, dpi=200)
    plt.close()

    runtime_rows = []
    for _, row in run_metrics.iterrows():
        runtime_rows.append(
            {
                "read_count": row["read_count"],
                "phase": "overall",
                "seconds": row["overall_runtime_seconds"],
            }
        )
        runtime_rows.append(
            {
                "read_count": row["read_count"],
                "phase": "index loading",
                "seconds": row["index_runtime_seconds"],
            }
        )
        runtime_rows.append(
            {
                "read_count": row["read_count"],
                "phase": "find/scanning",
                "seconds": row["find_runtime_seconds"],
            }
        )

    runtime_plot = output_dir / "runtime_by_read_count.png"
    runtime_data = pd.DataFrame(runtime_rows)
    plt.figure(figsize=(8, 4.8))
    sns.lineplot(
        data=runtime_data,
        x="read_count",
        y="seconds",
        hue="phase",
        marker="o",
        linewidth=2,
    )
    plt.title("PanelRecon runtime by N examined reads")
    plt.xlabel("N reads scanned")
    plt.ylabel("Runtime (seconds)")
    plt.tight_layout()
    plt.savefig(runtime_plot, dpi=200)
    plt.close()

    memory_plot = output_dir / "memory_by_read_count.png"
    plt.figure(figsize=(7, 4.5))
    sns.barplot(
        data=run_metrics,
        x="read_count",
        y="peak_rss_mb",
        color="#4c72b0",
    )
    max_memory_mb = run_metrics["peak_rss_mb"].max()
    y_padding_mb = max(max_memory_mb * 0.05, 1.0)
    plt.ylim(0, max_memory_mb + y_padding_mb)
    plt.title("Peak memory by read count")
    plt.xlabel("N reads scanned")
    plt.ylabel("RAM Peak usage (Mb)")
    plt.tight_layout()
    plt.savefig(memory_plot, dpi=200)
    plt.close()

    heatmap_plot = output_dir / "sample_correctness_heatmap.png"
    sample_summary = (
        sample_results.groupby("sample", as_index=False)
        .agg(
            expected_panel=("expected_panel", "first"),
        )
    )
    sample_summary["is_negative_control"] = (
        sample_summary["expected_panel"].astype(str).str.lower() == "none"
    )
    sample_summary["sample_type_order"] = 1
    sample_summary.loc[sample_summary["is_negative_control"], "sample_type_order"] = 0
    sample_summary["expected_panel_for_sort"] = sample_summary["expected_panel"].astype(str)
    sample_summary.loc[
        sample_summary["is_negative_control"],
        "expected_panel_for_sort",
    ] = "none"
    sample_order = (
        sample_summary
        .sort_values(["sample_type_order", "expected_panel_for_sort", "sample"])
        ["sample"]
        .tolist()
    )
    correctness_matrix = sample_results.pivot(
        index="sample",
        columns="read_count",
        values="is_correct",
    )
    correctness_matrix = correctness_matrix.reindex(sample_order)
    correctness_matrix = correctness_matrix.reindex(
        sorted(correctness_matrix.columns),
        axis=1,
    )
    correctness_matrix = correctness_matrix.transpose()

    panel_by_sample = sample_summary.set_index("sample")["expected_panel"].astype(str)
    panel_names_for_samples = []
    for sample in correctness_matrix.columns:
        panel_name = panel_by_sample.get(sample, "unknown")
        if panel_name.lower() == "none":
            panel_name = "none"
        panel_names_for_samples.append(panel_name)

    expected_panel_names = []
    if "none" in panel_names_for_samples:
        expected_panel_names.append("none")
    for panel_name in sorted(set(panel_names_for_samples)):
        if panel_name != "none":
            expected_panel_names.append(panel_name)

    panel_colors = {}
    if "none" in expected_panel_names:
        panel_colors["none"] = "#bdbdbd"

    non_control_panels = []
    for panel_name in expected_panel_names:
        if panel_name != "none":
            non_control_panels.append(panel_name)

    panel_palette = sns.color_palette("tab10", n_colors=len(non_control_panels))
    for i in range(len(non_control_panels)):
        panel_colors[non_control_panels[i]] = panel_palette[i]

    panel_color_order = []
    panel_code_by_name = {}
    for i in range(len(expected_panel_names)):
        panel_name = expected_panel_names[i]
        panel_color_order.append(panel_colors[panel_name])
        panel_code_by_name[panel_name] = i

    panel_codes = []
    for panel_name in panel_names_for_samples:
        panel_codes.append(panel_code_by_name[panel_name])
    expected_panel_matrix = pd.DataFrame(
        [panel_codes],
        index=["Expected panel"],
        columns=correctness_matrix.columns,
    )

    panel_group_boundaries = []
    for i in range(1, len(panel_names_for_samples)):
        if panel_names_for_samples[i] != panel_names_for_samples[i - 1]:
            panel_group_boundaries.append(i)

    heatmap_width = max(8.0, 0.25 * len(correctness_matrix.columns))
    fig, (panel_axis, axis) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(heatmap_width, 6.4),
        gridspec_kw={"height_ratios": [0.35, 4.6], "hspace": 0.04},
        sharex=True,
    )
    sns.heatmap(
        expected_panel_matrix,
        ax=panel_axis,
        cmap=ListedColormap(panel_color_order),
        cbar=False,
        linewidths=0.3,
        linecolor="white",
        xticklabels=False,
    )
    panel_axis.set_xlabel("")
    panel_axis.set_ylabel("")
    panel_axis.set_yticklabels(["Expected panel"], rotation=0)
    panel_axis.tick_params(axis="x", bottom=False, labelbottom=False)

    axis = sns.heatmap(
        correctness_matrix,
        ax=axis,
        cmap=["#d95f02", "#1b9e77"],
        cbar=False,
        linewidths=0.3,
        linecolor="white",
    )
    legend_handles = [
        Patch(facecolor="#1b9e77", edgecolor="white", label="Correct"),
        Patch(facecolor="#d95f02", edgecolor="white", label="Incorrect"),
    ]
    panel_legend_handles = []
    for panel_name in expected_panel_names:
        panel_label = panel_name
        if panel_label == "none":
            panel_label = "Negative control"
        elif panel_label in panel_plot_labels:
            panel_label = panel_plot_labels[panel_label]
        elif panel_label.endswith(".2bit"):
            panel_label = panel_label[:-5]
        panel_legend_handles.append(
            Patch(facecolor=panel_colors[panel_name], edgecolor="white", label=panel_label)
        )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.24, 0.98),
        ncol=2,
        frameon=False,
        title="PanelRecon result",
    )
    fig.legend(
        handles=panel_legend_handles,
        title="Expected gene panel",
        loc="upper center",
        bbox_to_anchor=(0.70, 0.98),
        ncol=3,
        frameon=False,
    )
    axis.set_xlabel("Sample")
    axis.set_ylabel("N reads scanned")
    axis.tick_params(axis="x", rotation=90)
    axis.tick_params(axis="y", rotation=0)

    x_tick_labels = []
    for sample in correctness_matrix.columns:
        x_tick_labels.append(sample)
    axis.set_xticklabels(x_tick_labels, rotation=90)

    for panel_boundary in panel_group_boundaries:
        axis.axvline(panel_boundary, color="black", linewidth=0.8)
        panel_axis.axvline(panel_boundary, color="black", linewidth=0.8)

    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.42, top=0.76, hspace=0.04)
    fig.savefig(heatmap_plot, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"[plots] wrote {accuracy_plot}")
    print(f"[plots] wrote {runtime_plot}")
    print(f"[plots] wrote {memory_plot}")
    print(f"[plots] wrote {heatmap_plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
