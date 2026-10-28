#!/usr/bin/env python3
"""Run the optional Arpeggio-backed typed interaction analysis."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.interaction_cli import main


if __name__ == "__main__":
    raise SystemExit(main())
