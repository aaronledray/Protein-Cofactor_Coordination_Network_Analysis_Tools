#!/usr/bin/env python3
"""
v0.0.1 (STMQNC): walks an AF3 result folder, grabs the best-ranked sample per
run (from ranking_scores.csv), skips the redundant run-level model, and calls
STSQNC v0.0.3 once per selected model.
"""
import argparse
import csv
import os
import subprocess
import sys
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


STSQC_SCRIPT = Path(__file__).with_name("2_Single_Template_Single_Query_Network_Comparison_STSQNC_v0.0.3.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-align a template against the best AF3 sample per run."
    )
    parser.add_argument("--template-csv", required=True, help="Coordination CSV produced by SSCNA.")
    parser.add_argument("--template-atoms-csv", help="Optional Coord_Breakdown_atoms CSV; inferred if omitted.")
    parser.add_argument("--af3-root", required=True, help="Root folder containing AF3 runs (e.g., .../11_5).")
    parser.add_argument("--job-prefix", default="stmqnc", help="Prefix for per-run job names and manifest.")
    parser.add_argument("--match-cutoff", type=float, default=2.0, help="Distance cutoff (Å) for STSQNC.")
    return parser.parse_args()


def find_run_dirs(root: Path) -> List[Path]:
    return sorted({p.parent for p in root.rglob("ranking_scores.csv") if p.is_file()})


def pick_best_sample(ranking_csv: Path) -> Optional[Tuple[str, str, float]]:
    with ranking_csv.open() as f:
        reader = csv.DictReader(f)
        best: Optional[Tuple[str, str, float]] = None
        for row in reader:
            seed = (row.get("seed") or "").strip()
            sample = (row.get("sample") or "").strip()
            score_val = row.get("ranking_score")
            try:
                score = float(score_val)
            except Exception:
                continue
            if not seed or not sample:
                continue
            if best is None or score > best[2]:
                best = (seed, sample, score)
    return best


def locate_model(run_dir: Path, seed: str, sample: str) -> Optional[Path]:
    sample_dir = f"seed-{seed}_sample-{sample}"
    direct = run_dir / sample_dir / "model.cif"
    if direct.is_file():
        return direct

    scores_csv = run_dir / "scores.csv"
    if scores_csv.is_file():
        with scores_csv.open() as f:
            for row in csv.DictReader(f):
                if (row.get("sample_dir") or "").strip() != sample_dir:
                    continue
                rel = (row.get("model_relpath") or "").strip()
                if rel:
                    candidate = (run_dir / rel).resolve()
                    if candidate.is_file():
                        return candidate
    return None


def load_confidences(run_dir: Path, seed: str, sample: str) -> Dict[str, Optional[float]]:
    """
    Load confidence scores for a specific sample. Falls back to run-level summary if missing.
    """
    fields = ["iptm", "ptm", "ranking_score", "fraction_disordered", "has_clash"]
    result = {k: None for k in fields}
    sample_conf = run_dir / f"seed-{seed}_sample-{sample}" / "summary_confidences.json"
    run_conf = run_dir / f"{run_dir.name}_summary_confidences.json"
    for path in (sample_conf, run_conf):
        if not path.is_file():
            continue
        try:
            data = path.read_text()
            obj = json.loads(data)
            for k in fields:
                if k in obj:
                    result[k] = obj[k]
        except Exception:
            continue
        break
    return result


def run_stsqnc(
    template_csv: Path,
    template_atoms_csv: Optional[Path],
    query_structure: Path,
    job_name: str,
    match_cutoff: float
) -> Tuple[bool, str]:
    cmd = [
        sys.executable,
        str(STSQC_SCRIPT),
        "--template-csv", str(template_csv),
        "--query-structure", str(query_structure),
        "--job-name", job_name,
        "--match-cutoff", str(match_cutoff),
    ]
    if template_atoms_csv:
        cmd += ["--template-atoms-csv", str(template_atoms_csv)]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        return True, proc.stdout.strip()
    detail = (proc.stderr or proc.stdout or "").strip()
    return False, detail


def ensure_matches_dir() -> Path:
    path = Path(os.getcwd()) / "matches"
    path.mkdir(exist_ok=True)
    return path


def read_alignment_metrics(matches_dir: Path, job_name: str) -> Tuple[Optional[float], Optional[int], Optional[int], Optional[float]]:
    """
    Pull RMSD from alignment_report.csv and network match stats from job summary.
    Returns (rmsd, matched, total, percent).
    """
    rmsd: Optional[float] = None
    matched: Optional[int] = None
    total: Optional[int] = None
    percent: Optional[float] = None

    job_dir = matches_dir / job_name
    report_path = job_dir / "alignment_report.csv"
    if report_path.is_file():
        with report_path.open() as f:
            for row in csv.DictReader(f):
                if (row.get("job_name") or "") == job_name:
                    try:
                        rmsd = float(row.get("rmsd", ""))
                    except Exception:
                        rmsd = None
                    break

    summary_path = matches_dir / job_name / f"{job_name}_results.md"
    if summary_path.is_file():
        with summary_path.open() as f:
            for line in f:
                line = line.strip()
                if line.startswith("- Matched"):
                    # Format: "- Matched X/Y network atoms (Z%)."
                    m = re.search(r"Matched\s+(\d+)\s*/\s*(\d+)", line)
                    if m:
                        matched = int(m.group(1))
                        total = int(m.group(2))
                        percent = (matched / total * 100.0) if total else None
                    break

    # Fallback: derive network stats from network_match.csv if summary parse failed
    if matched is None or total is None:
        net_path = matches_dir / job_name / f"{job_name}_network_match.csv"
        if net_path.is_file():
            try:
                with net_path.open() as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                total = len(rows)
                matched = sum(1 for r in rows if (r.get("matched_within_cutoff") or "").upper() == "Y")
                percent = (matched / total * 100.0) if total else None
            except Exception:
                pass

    return rmsd, matched, total, percent


def compute_job_metrics(matches_dir: Path, job_name: str) -> Dict[str, Optional[float]]:
    """Compute rmsd and network stats directly from job outputs (no markdown parsing)."""
    job_dir = matches_dir / job_name
    metrics: Dict[str, Optional[float]] = {
        "rmsd": None,
        "matched": None,
        "total": None,
        "percent": None,
    }

    # RMSD from alignment_report.csv
    report_path = job_dir / "alignment_report.csv"
    if report_path.is_file():
        try:
            with report_path.open() as f:
                for row in csv.DictReader(f):
                    if (row.get("job_name") or "") == job_name:
                        try:
                            metrics["rmsd"] = float(row.get("rmsd", ""))
                        except Exception:
                            pass
                        break
        except Exception:
            pass

    # Network stats from network_match.csv
    net_path = job_dir / f"{job_name}_network_match.csv"
    if net_path.is_file():
        try:
            with net_path.open() as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            total = len(rows)
            matched = sum(1 for r in rows if (r.get("matched_within_cutoff") or "").upper() == "Y")
            percent = (matched / total * 100.0) if total else None
            metrics["matched"] = matched
            metrics["total"] = total
            metrics["percent"] = percent
        except Exception:
            pass

    return metrics


def write_manifest(
    manifest_path: Path,
    rows: List[Dict[str, str]]
) -> None:
    header = [
        "run_dir",
        "best_seed",
        "best_sample",
        "ranking_score",
        "model_path",
        "job_name",
        "status",
        "message",
        "rmsd",
        "network_matched",
        "network_total",
        "network_percent",
        "iptm",
        "ptm",
        "af_ranking_score",
        "fraction_disordered",
        "has_clash",
    ]
    write_header = not manifest_path.is_file()
    with manifest_path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(summary_path: Path, rows: List[Dict[str, str]]) -> None:
    header = ["job_name", "run_dir", "model_path", "rmsd", "network", "iptm", "ptm", "af_ranking_score", "fraction_disordered", "has_clash"]
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def backfill_summary(summary_path: Path, matches_dir: Path) -> None:
    """
    Ensure network and rmsd columns are populated by recomputing from job outputs.
    This mirrors the manual fix users were running.
    """
    if not summary_path.is_file():
        return
    with summary_path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        header = reader.fieldnames or []
    updated_rows: List[Dict[str, str]] = []
    for row in rows:
        job_name = row.get("job_name", "")
        if not job_name:
            updated_rows.append(row)
            continue
        metrics = compute_job_metrics(matches_dir, job_name)
        rmsd = metrics.get("rmsd")
        matched = metrics.get("matched")
        total = metrics.get("total")
        percent = metrics.get("percent")
        pct = percent if percent is not None else ((matched / total * 100.0) if matched is not None and total else None)
        if matched is not None and total is not None:
            row["network"] = f"{matched}/{total} ({pct:.1f}%)" if pct is not None else f"{matched}/{total}"
        if rmsd is not None:
            row["rmsd"] = f"{rmsd:.4f}"
        updated_rows.append(row)

    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(updated_rows)


def main() -> None:
    args = parse_args()
    template_csv = Path(args.template_csv).resolve()
    template_atoms_csv = Path(args.template_atoms_csv).resolve() if args.template_atoms_csv else None
    af3_root = Path(args.af3_root).resolve()

    if not template_csv.is_file():
        raise FileNotFoundError(f"Template CSV not found: {template_csv}")
    if template_atoms_csv and not template_atoms_csv.is_file():
        raise FileNotFoundError(f"Template atoms CSV not found: {template_atoms_csv}")
    if not af3_root.is_dir():
        raise FileNotFoundError(f"AF3 root not found: {af3_root}")
    if not STSQC_SCRIPT.is_file():
        raise FileNotFoundError(f"STSQNC script not found: {STSQC_SCRIPT}")

    run_dirs = find_run_dirs(af3_root)
    if not run_dirs:
        raise RuntimeError(f"No ranking_scores.csv found under {af3_root}")

    matches_dir = ensure_matches_dir()
    manifest_path = matches_dir / f"{args.job_prefix}_batch_manifest.csv"
    manifest_rows: List[Dict[str, str]] = []
    summary_rows: List[Dict[str, str]] = []

    for run_dir in run_dirs:
        ranking_csv = run_dir / "ranking_scores.csv"
        best = pick_best_sample(ranking_csv)
        if best is None:
            manifest_rows.append({
                "run_dir": str(run_dir),
                "best_seed": "",
                "best_sample": "",
                "ranking_score": "",
                "model_path": "",
                "job_name": "",
                "status": "skipped",
                "message": "No valid rows in ranking_scores.csv",
            })
            continue

        seed, sample, score = best
        model_path = locate_model(run_dir, seed, sample)
        if model_path is None:
            manifest_rows.append({
                "run_dir": str(run_dir),
                "best_seed": seed,
                "best_sample": sample,
                "ranking_score": f"{score:.6f}",
                "model_path": "",
                "job_name": "",
                "status": "skipped",
                "message": "Model file not found for best sample",
                "rmsd": "",
                "network_matched": "",
                "network_total": "",
                "network_percent": "",
            })
            continue

        job_name = f"{args.job_prefix}__{run_dir.name}__seed{seed}_sample{sample}"
        ok, detail = run_stsqnc(
            template_csv=template_csv,
            template_atoms_csv=template_atoms_csv,
            query_structure=model_path,
            job_name=job_name,
            match_cutoff=float(args.match_cutoff),
        )
        metrics = compute_job_metrics(matches_dir, job_name) if ok else {"rmsd": None, "matched": None, "total": None, "percent": None}
        rmsd = metrics.get("rmsd")
        matched = metrics.get("matched")
        total = metrics.get("total")
        percent = metrics.get("percent")
        conf = load_confidences(run_dir, seed, sample)
        iptm = conf.get("iptm")
        ptm = conf.get("ptm")
        af_ranking = conf.get("ranking_score")
        frac_dis = conf.get("fraction_disordered")
        has_clash = conf.get("has_clash")
        pct = percent if percent is not None else ((matched / total * 100.0) if matched is not None and total else None)
        if ok:
            summary_rows.append({
                "job_name": job_name,
                "run_dir": str(run_dir),
                "model_path": str(model_path),
                "rmsd": f"{rmsd:.4f}" if rmsd is not None else "",
                "network": f"{matched}/{total} ({pct:.1f}%)" if matched is not None and total is not None else "",
                "iptm": f"{iptm:.4f}" if isinstance(iptm, (int, float)) else "",
                "ptm": f"{ptm:.4f}" if isinstance(ptm, (int, float)) else "",
                "af_ranking_score": f"{af_ranking:.4f}" if isinstance(af_ranking, (int, float)) else "",
                "fraction_disordered": f"{frac_dis:.4f}" if isinstance(frac_dis, (int, float)) else "",
                "has_clash": str(has_clash) if has_clash is not None else "",
            })
        manifest_rows.append({
            "run_dir": str(run_dir),
            "best_seed": seed,
            "best_sample": sample,
            "ranking_score": f"{score:.6f}",
            "model_path": str(model_path),
            "job_name": job_name,
            "status": "ok" if ok else "failed",
            "message": detail,
            "rmsd": f"{rmsd:.4f}" if rmsd is not None else "",
            "network_matched": str(matched) if matched is not None else "",
            "network_total": str(total) if total is not None else "",
            "network_percent": f"{pct:.1f}" if pct is not None else "",
            "iptm": f"{iptm:.4f}" if isinstance(iptm, (int, float)) else "",
            "ptm": f"{ptm:.4f}" if isinstance(ptm, (int, float)) else "",
            "af_ranking_score": f"{af_ranking:.4f}" if isinstance(af_ranking, (int, float)) else "",
            "fraction_disordered": f"{frac_dis:.4f}" if isinstance(frac_dis, (int, float)) else "",
            "has_clash": str(has_clash) if has_clash is not None else "",
        })

    write_manifest(manifest_path, manifest_rows)
    summary_path = matches_dir / f"{args.job_prefix}_batch_summary.csv"
    write_summary(summary_path, summary_rows)
    backfill_summary(summary_path, matches_dir)

    total = len(manifest_rows)
    ok_count = sum(1 for r in manifest_rows if r["status"] == "ok")
    skipped = sum(1 for r in manifest_rows if r["status"] == "skipped")
    failed = total - ok_count - skipped

    print(f"[DONE] Processed {total} runs under {af3_root}.")
    print(f"       ok={ok_count}, skipped={skipped}, failed={failed}")
    print(f"[INFO] Manifest: {manifest_path}")
    if summary_rows:
        print("[SUMMARY] Alignment metrics (best sample per run):")
        for row in summary_rows:
            line = f" - {row['job_name']}: RMSD={row['rmsd'] or 'NA'}"
            if row["network"]:
                line += f", Network={row['network']}"
            if row.get("iptm"):
                line += f", ipTM={row['iptm']}"
            if row.get("ptm"):
                line += f", pTM={row['ptm']}"
            if row.get("af_ranking_score"):
                line += f", AF rank score={row['af_ranking_score']}"
            print(line)
        print(f"[INFO] Batch summary: {summary_path}")


if __name__ == "__main__":
    main()
