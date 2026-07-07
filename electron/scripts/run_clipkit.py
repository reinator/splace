#!/usr/bin/env python3
"""Small wrapper used by the SPLACE desktop build to run ClipKIT safely.

Usage:
    python scripts/run_clipkit.py input.fasta output.fasta [extra ClipKIT args...]
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 3:
        print("Usage: run_clipkit.py input.fasta output.fasta [extra ClipKIT args...]", file=sys.stderr)
        return 2

    input_fasta = Path(sys.argv[1])
    output_fasta = Path(sys.argv[2])
    extra_args = sys.argv[3:]

    if not input_fasta.exists():
        print(f"Input FASTA not found: {input_fasta}", file=sys.stderr)
        return 2

    cmd = [
        sys.executable,
        "-m",
        "clipkit",
        str(input_fasta),
        "-o",
        str(output_fasta),
        *extra_args,
    ]

    completed = subprocess.run(cmd, text=True, capture_output=True)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="", file=sys.stderr)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
