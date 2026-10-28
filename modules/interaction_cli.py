"""Command line entry point for opt-in interaction typing."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .interaction_typing import run_arpeggio, run_arpeggio_batch


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--arpeggio", default="pdbe-arpeggio", help="Arpeggio 1.4.4 executable")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--vdw-compensation", type=float, default=0.1)
    parser.add_argument("--interacting-cutoff", type=float, default=5.0)
    parser.add_argument("--mapped-edges", type=Path, help="Mapped per-structure edge ledger CSV with sequence/group metadata")
    parser.add_argument("--selection", action="append", help="Arpeggio selection (repeatable), e.g. /A/272")
    parser.add_argument("--no-hydrogenated-output", action="store_true")
    args = parser.parse_args()
    common = dict(executable=args.arpeggio, ph=args.ph,
                  vdw_compensation=args.vdw_compensation,
                  interacting_cutoff_A=args.interacting_cutoff)
    if len(args.structure) == 1:
        manifest = run_arpeggio(
            args.structure[0], args.output, mapped_edges_csv=args.mapped_edges,
            selections=args.selection, write_hydrogenated=not args.no_hydrogenated_output, **common)
    else:
        manifest = run_arpeggio_batch(
            args.structure, args.output, mapped_edges_csv=args.mapped_edges, **common)
    print(json.dumps(manifest, indent=2))
    return 0
