#!/usr/bin/env python3
"""Smoke benchmark over representative cofactor families.

Runs each case in a fresh call to ``analyze_structure`` and records wall time,
peak Python memory, output row counts, and per-case failures (a failing case
never stops the run). A final multiprocessing pass checks that ``batch_analyze``
matches serial results and isolates a deliberately missing input. Results go to
``benchmark_output/`` (ignored by git).
"""

import argparse
import json
import sys
import time
import tracemalloc
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from modules.coordination_api import analyze_structure, batch_analyze  # noqa: E402
from modules.ml_export import export_residue_labels  # noqa: E402

REF = ROOT / "reference_structures"
CASES = [
    {"name": "cu_plastocyanin", "path": REF / "0_Plastocyanin/1ag6.cif", "cofactor": "CU", "options": {}},
    {"name": "heme_myoglobin", "path": REF / "3_Myoglobin/1a6m.pdb", "cofactor": "HEM", "options": {}},
    {"name": "zn_carbonic_anhydrase", "path": REF / "3_Metals/1ca2.pdb", "cofactor": "ZN", "options": {}},
    {"name": "fes_ferredoxin_per_site", "path": REF / "3_Metals/2zvs.pdb", "cofactor": "SF4",
     "options": {"site_mode": "per-site"}},
    {"name": "nitrogenase_combinatorial", "path": REF / "2_Nitrogenase/3u7q_monomer.pdb",
     "cofactor": ["ICS", "CLF", "HCA"], "options": {"combinatorial": True}},
    {"name": "oec_cluster", "path": REF / "1_OEX/0_Aligned_Reduced/4ub6.pdb", "cofactor": "OEX",
     "options": {"site_mode": "per-site", "first_model_only": True}},
    {"name": "missing_input_isolated", "path": REF / "does_not_exist.pdb", "cofactor": "CU", "options": {}},
]


def run_case(case, shells):
    row = {"case": case["name"], "status": "ok", "error": ""}
    tracemalloc.start()
    start = time.perf_counter()
    try:
        tables = analyze_structure(case["path"], case["cofactor"], shells=shells, **case["options"])
        row["wall_s"] = round(time.perf_counter() - start, 3)
        for key, frame in tables.items():
            row[f"n_{key}"] = len(frame)
        row["n_sites"] = int(tables["atoms"]["site_id"].nunique()) if len(tables["atoms"]) else 0
        try:
            row["n_label_rows"] = len(export_residue_labels(tables))
        except ValueError as exc:  # documented policy rejection, not a failure
            row["n_label_rows"] = None
            row["label_export_note"] = str(exc)[:120]
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                   wall_s=round(time.perf_counter() - start, 3))
    row["peak_python_mb"] = round(tracemalloc.get_traced_memory()[1] / 1e6, 1)
    tracemalloc.stop()
    return row


def multiprocessing_check(workers):
    paths = [str(c["path"]) for c in CASES if c["name"] == "cu_plastocyanin"]
    paths += [str(REF / "3_Myoglobin/1a6m.pdb"), str(REF / "does_not_exist.pdb")]
    serial = batch_analyze(paths, ["CU", "HEM"], workers=1)
    parallel = batch_analyze(paths, ["CU", "HEM"], workers=workers)
    same = all(serial[k].reset_index(drop=True).equals(parallel[k].reset_index(drop=True))
               for k in ("residues", "atoms", "links", "contacts"))
    return {"serial_equals_parallel": bool(same), "workers": workers,
            "errors_isolated": len(parallel["errors"]), "structures_ok": int(parallel["atoms"]["structure_id"].nunique())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shells", type=int, default=3)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--output-dir", default=str(ROOT / "benchmark_output"))
    args = parser.parse_args()

    rows = [run_case(case, args.shells) for case in CASES]
    frame = pd.DataFrame(rows)
    mp = multiprocessing_check(args.workers)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "batch_smoke_benchmark.csv", index=False)
    (out / "batch_smoke_benchmark.json").write_text(
        json.dumps({"shells": args.shells, "cases": rows, "multiprocessing": mp}, indent=2, default=str))
    print(frame.to_string(index=False))
    print(json.dumps(mp))
    unexpected = frame[(frame["status"] == "failed") & (frame["case"] != "missing_input_isolated")]
    return 1 if len(unexpected) or not mp["serial_equals_parallel"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
