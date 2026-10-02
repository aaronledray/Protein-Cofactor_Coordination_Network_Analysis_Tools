#!/usr/bin/env python3
"""Write label-sensitivity tables for the representative benchmark cases.

Water is not a variant: there is no water option; waters are ordinary shell
members (motif ``water``), so there is no water policy to vary. Output goes to the ignored ``benchmark_output/``.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from batch_smoke_benchmark import CASES, ROOT  # noqa: E402

from modules.sensitivity import sensitivity_table  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shells", type=int, default=3)
    parser.add_argument("--output-dir", default=str(ROOT / "benchmark_output"))
    args = parser.parse_args()

    frames = []
    for case in CASES:
        if not Path(case["path"]).is_file():
            continue
        options = {k: v for k, v in case["options"].items() if k != "site_mode"}
        table = sensitivity_table(
            case["path"], case["cofactor"], shells=args.shells, baseline_options=options
        )
        frames.append(table.assign(case=case["name"]))
    result = pd.concat(frames, ignore_index=True)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "label_sensitivity.csv", index=False)
    print(result.drop(columns="error").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
