#!/usr/bin/env python3
"""Dataset-scale benchmark for the labeling workflow.

Runs ``analyze_labeling_example`` over a directory of structures (optionally an
evenly spaced, deterministic sample), serially and with a process pool, and
records per-structure wall time, peak RSS, label counts, and failures. Failures
are isolated per structure. Output (ignored by git) goes to ``benchmark_output/``.

    python benchmarks/dataset_benchmark.py --input DIR --cofactor HEM --limit 200 --workers 8
"""

import argparse
import json
import resource
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from modules.coordination_api import discover_structure_paths  # noqa: E402
from modules.ml_export import MixedCofactorError, analyze_labeling_example  # noqa: E402


def sample_paths(paths, limit):
    """Evenly spaced deterministic sample of a sorted path list."""
    paths = sorted(paths)
    if not limit or limit >= len(paths):
        return paths
    step = len(paths) / limit
    return [paths[int(i * step)] for i in range(limit)]


def run_one(args):
    path, cofactor = args
    start = time.perf_counter()
    row = {"structure": Path(path).name, "status": "ok", "error": ""}
    try:
        labels, contract = analyze_labeling_example(path, cofactor)
        row.update(
            n_labels=len(labels),
            n_sites=int(labels["site_id"].nunique()) if len(labels) else 0,
            n_primary=int(labels["is_primary"].sum()) if len(labels) else 0,
        )
        if not len(labels):
            row["status"] = "no_cofactor"
    except MixedCofactorError as exc:
        row.update(status="rejected", error=str(exc)[:200])
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}"[:200])
    row["wall_s"] = round(time.perf_counter() - start, 4)
    # ru_maxrss is bytes on macOS, KiB on Linux; it is the worker's high-water mark.
    scale = 1e6 if sys.platform == "darwin" else 1e3
    row["worker_peak_rss_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / scale, 1)
    return row


def run_pass(paths, cofactor, workers):
    jobs = [(str(p), cofactor) for p in paths]
    start = time.perf_counter()
    if workers == 1:
        rows = [run_one(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(run_one, jobs, chunksize=4))
    return pd.DataFrame(rows), time.perf_counter() - start


def summarize(frame, wall, workers):
    times = frame["wall_s"].tolist()
    q = statistics.quantiles(times, n=20) if len(times) >= 2 else times * 19
    return {
        "workers": workers,
        "structures": len(frame),
        "status_counts": frame["status"].value_counts().to_dict(),
        "wall_s": round(wall, 2),
        "structures_per_s": round(len(frame) / wall, 2) if wall else None,
        "per_structure_s": {"median": round(statistics.median(times), 3),
                            "p95": round(q[18], 3), "max": round(max(times), 3)},
        "max_worker_rss_mb": float(frame["worker_peak_rss_mb"].max()),
        "labels_total": int(frame["n_labels"].fillna(0).sum()),
        "error_types": frame.loc[frame["error"] != "", "error"].str.split(":").str[0]
                            .value_counts().to_dict(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--cofactor", required=True)
    parser.add_argument("--limit", type=int, default=200, help="0 = all structures")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--skip-serial", action="store_true")
    parser.add_argument("--serial-limit", type=int, default=50,
                        help="Serial pass uses only this many of the sampled structures")
    parser.add_argument("--output-dir", default=str(ROOT / "benchmark_output"))
    args = parser.parse_args()

    paths = sample_paths(discover_structure_paths(args.input), args.limit)
    report = {"input": str(args.input), "cofactor": args.cofactor, "sampled": len(paths), "passes": []}
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if not args.skip_serial:
        serial_paths = paths[:: max(1, len(paths) // args.serial_limit)][: args.serial_limit]
        frame, wall = run_pass(serial_paths, args.cofactor, 1)
        report["passes"].append(summarize(frame, wall, 1))
    frame, wall = run_pass(paths, args.cofactor, args.workers)
    report["passes"].append(summarize(frame, wall, args.workers))
    frame.to_csv(out / "dataset_benchmark_structures.csv", index=False)
    (out / "dataset_benchmark_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
