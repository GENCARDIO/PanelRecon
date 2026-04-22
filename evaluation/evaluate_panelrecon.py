#!/usr/bin/env python3
"""Evaluate PanelRecon panel identification across different read limits."""

import argparse
import csv
import select
import shlex
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PANELRECON = REPO_ROOT / "PanelRecon"
DEFAULT_READ_COUNTS = "1000,5000,10000,30000,50000,100000"
POLL_INTERVAL_SECONDS = 0.05


def read_peak_rss_kb(pid):
    """Read peak resident memory from /proc when available."""
    status_path = Path("/proc") / str(pid) / "status"
    try:
        values = {}
        with status_path.open() as handle:
            for line in handle:
                if line.startswith(("VmHWM:", "VmRSS:")):
                    fields = line.split()
                    if len(fields) >= 2:
                        values[fields[0].rstrip(":")] = int(fields[1])
        return values.get("VmHWM", values.get("VmRSS", 0))
    except OSError:
        return 0


def main():
    default_fastq_list = REPO_ROOT / "evaluation" / "fastq_list.2.txt"
    if not default_fastq_list.exists():
        default_fastq_list = REPO_ROOT / "fastq_list.2.txt"

    default_known_panels = REPO_ROOT / "evaluation" / "known.panels.txt"
    if not default_known_panels.exists():
        default_known_panels = REPO_ROOT / "evaluation" / "known.panels.tsv"
    if not default_known_panels.exists():
        default_known_panels = REPO_ROOT / "evaluation" / "known.panes.tsv"

    parser = argparse.ArgumentParser(
        description="Evaluate PanelRecon find accuracy across multiple read limits."
    )
    parser.add_argument(
        "--index_dir",
        type=Path,
        default=REPO_ROOT / "test",
    )
    parser.add_argument(
        "--fastq_list",
        type=Path,
        default=default_fastq_list,
    )
    parser.add_argument(
        "--known_panels",
        type=Path,
        default=default_known_panels,
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=REPO_ROOT / "evaluation" / "read_depth_results",
    )
    parser.add_argument(
        "--read_counts",
        default=DEFAULT_READ_COUNTS,
    )
    args = parser.parse_args()

    missing_inputs = []
    for name, path in [
        ("PanelRecon", PANELRECON),
        ("index_dir", args.index_dir),
        ("fastq_list", args.fastq_list),
        ("known_panels", args.known_panels),
    ]:
        if not path.exists():
            missing_inputs.append(f"{name}: {path}")
    if missing_inputs:
        print("ERROR: missing required input(s):", file=sys.stderr)
        for missing_input in missing_inputs:
            print(f"  {missing_input}", file=sys.stderr)
        return 1

    read_counts = []
    for token in args.read_counts.replace(",", " ").split():
        try:
            read_count = int(token)
        except ValueError:
            print(f"ERROR: invalid read count: {token}", file=sys.stderr)
            return 1
        if read_count < 1:
            print(f"ERROR: read counts must be positive: {read_count}", file=sys.stderr)
            return 1
        read_counts.append(read_count)
    if not read_counts:
        print("ERROR: at least one read count is required", file=sys.stderr)
        return 1

    known_panels = {}
    with args.known_panels.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, start=1):
            is_empty_line = not fields or not fields[0].strip()
            is_comment_line = not is_empty_line and fields[0].startswith("#")
            if is_empty_line or is_comment_line:
                continue

            is_header = (
                line_number == 1
                and len(fields) >= 2
                and fields[0] == "sample"
                and fields[1] == "best_panel"
            )
            if is_header:
                continue

            has_missing_columns = len(fields) < 2
            has_empty_values = not has_missing_columns and (
                not fields[0].strip() or not fields[1].strip()
            )
            if has_missing_columns or has_empty_values:
                print(f"ERROR: invalid known panel line {line_number}", file=sys.stderr)
                return 1

            expected_panel = fields[1].strip()
            if expected_panel.lower() == "none":
                expected_panel = "none"

            known_panels[fields[0].strip()] = expected_panel

    if not known_panels:
        print(f"ERROR: no known panels found in {args.known_panels}", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = args.output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    sample_result_rows = []
    run_metric_rows = []
    sample_result_fields = [
        "read_count",
        "sample",
        "expected_panel",
        "predicted_panel",
        "is_correct",
        "scanned_reads",
        "best_score",
        "score_margin_vs_next",
        "best_panel_covered_kmers",
        "best_panel_covered_kmers_pct",
        "rank_file",
    ]
    run_metric_fields = [
        "read_count",
        "exit_code",
        "overall_runtime_seconds",
        "index_runtime_seconds",
        "find_runtime_seconds",
        "peak_rss_kb",
        "peak_rss_mb",
        "total_known_samples",
        "evaluated_samples",
        "missing_samples",
        "correct_samples",
        "accuracy",
        "rank_file",
        "stdout_log",
        "stderr_log",
        "command",
    ]

    for read_count in read_counts:
        rank_path = args.output_dir / f"panel_ranks.max_reads_{read_count}.tsv"
        stdout_log = log_dir / f"max_reads_{read_count}.stdout.log"
        stderr_log = log_dir / f"max_reads_{read_count}.stderr.log"
        command = [
            str(PANELRECON),
            "find",
            "--index_dir",
            str(args.index_dir),
            "--fastq_list",
            str(args.fastq_list),
            "--min_reads",
            str(read_count),
            "--max_reads",
            str(read_count),
            "--output",
            str(rank_path),
        ]

        print(f"[INFO] running max_reads={read_count}", flush=True)
        peak_rss_kb = 0
        index_start_time = None
        scan_start_time = None
        start_time = time.perf_counter()
        with stdout_log.open("w") as stdout_handle, stderr_log.open("w") as stderr_handle:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=stderr_handle,
                text=True,
                bufsize=1,
            )
            while process.poll() is None:
                peak_rss_kb = max(peak_rss_kb, read_peak_rss_kb(process.pid))

                readable_streams = []
                if process.stdout is not None:
                    readable_streams, _, _ = select.select(
                        [process.stdout],
                        [],
                        [],
                        POLL_INTERVAL_SECONDS,
                    )

                if readable_streams:
                    line = process.stdout.readline()
                    if line:
                        stdout_handle.write(line)
                        stdout_handle.flush()

                        now = time.perf_counter()
                        if "INFO: Reading panel index from" in line:
                            index_start_time = now

                        is_first_sample_line = line.startswith(" INFO: (")
                        if index_start_time is not None and scan_start_time is None:
                            if is_first_sample_line:
                                scan_start_time = now

            exit_code = process.wait()

            if process.stdout is not None:
                for line in process.stdout:
                    stdout_handle.write(line)
                    stdout_handle.flush()

                    now = time.perf_counter()
                    if "INFO: Reading panel index from" in line:
                        index_start_time = now

                    is_first_sample_line = line.startswith(" INFO: (")
                    if index_start_time is not None and scan_start_time is None:
                        if is_first_sample_line:
                            scan_start_time = now

        end_time = time.perf_counter()
        overall_runtime_seconds = end_time - start_time
        index_runtime_seconds = 0.0
        find_runtime_seconds = 0.0

        if index_start_time is not None and scan_start_time is not None:
            index_runtime_seconds = scan_start_time - index_start_time

        if scan_start_time is not None:
            find_runtime_seconds = end_time - scan_start_time

        rank_rows = {}
        if exit_code == 0:
            with rank_path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                has_required_columns = (
                    reader.fieldnames is not None
                    and "sample" in reader.fieldnames
                    and "best_panel" in reader.fieldnames
                )
                if not has_required_columns:
                    print(
                        f"ERROR: rank file is missing required columns: {rank_path}",
                        file=sys.stderr,
                    )
                    return 1
                for row in reader:
                    rank_rows[row["sample"]] = row

        correct_samples = 0
        missing_samples = 0
        for sample, expected_panel in sorted(known_panels.items()):
            rank_row = rank_rows.get(sample)
            if rank_row is None:
                missing_samples += 1
                predicted_panel = ""
                is_correct = False
                scanned_reads = ""
                best_score = ""
                score_margin_vs_next = ""
                best_panel_covered_kmers = ""
                best_panel_covered_kmers_pct = ""
            else:
                predicted_panel = rank_row.get("best_panel", "").strip()
                predicted_panel_for_check = predicted_panel
                if predicted_panel_for_check.lower() == "none":
                    predicted_panel_for_check = "none"
                is_correct = predicted_panel_for_check == expected_panel
                scanned_reads = rank_row.get("scanned_reads", "")
                best_score = rank_row.get("best_score", "")
                score_margin_vs_next = rank_row.get("score_margin_vs_next", "")
                best_panel_covered_kmers = rank_row.get("best_panel_covered_kmers", "")
                best_panel_covered_kmers_pct = rank_row.get("best_panel_covered_kmers_pct", "")

            if is_correct:
                correct_samples += 1

            sample_result_rows.append(
                {
                    "read_count": read_count,
                    "sample": sample,
                    "expected_panel": expected_panel,
                    "predicted_panel": predicted_panel,
                    "is_correct": int(is_correct),
                    "scanned_reads": scanned_reads,
                    "best_score": best_score,
                    "score_margin_vs_next": score_margin_vs_next,
                    "best_panel_covered_kmers": best_panel_covered_kmers,
                    "best_panel_covered_kmers_pct": best_panel_covered_kmers_pct,
                    "rank_file": str(rank_path),
                }
            )

        total_samples = len(known_panels)
        evaluated_samples = total_samples - missing_samples
        accuracy = 0.0
        if total_samples > 0:
            accuracy = correct_samples / total_samples

        run_metric_rows.append(
            {
                "read_count": read_count,
                "exit_code": exit_code,
                "overall_runtime_seconds": f"{overall_runtime_seconds:.3f}",
                "index_runtime_seconds": f"{index_runtime_seconds:.3f}",
                "find_runtime_seconds": f"{find_runtime_seconds:.3f}",
                "peak_rss_kb": peak_rss_kb,
                "peak_rss_mb": f"{peak_rss_kb / 1024.0:.2f}",
                "total_known_samples": total_samples,
                "evaluated_samples": evaluated_samples,
                "missing_samples": missing_samples,
                "correct_samples": correct_samples,
                "accuracy": f"{accuracy:.4f}",
                "rank_file": str(rank_path),
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "command": shlex.join(command),
            }
        )
        print(
            f"[INFO] max_reads={read_count} accuracy={accuracy:.4f} "
            f"overall={overall_runtime_seconds:.3f}s "
            f"index={index_runtime_seconds:.3f}s "
            f"find={find_runtime_seconds:.3f}s "
            f"peak_rss={peak_rss_kb / 1024.0:.2f}MB",
            flush=True,
        )

    sample_results_path = args.output_dir / "sample_results.tsv"
    with sample_results_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_result_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_result_rows)

    run_metrics_path = args.output_dir / "run_metrics.tsv"
    with run_metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=run_metric_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(run_metric_rows)

    print(f"[INFO] wrote {sample_results_path}")
    print(f"[INFO] wrote {run_metrics_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
