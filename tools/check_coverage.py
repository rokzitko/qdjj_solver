"""Enforce Python line coverage; report branch coverage without a branch gate."""

import argparse
import json
from pathlib import Path


def check_report(report, minimum=95.):
    """Check coverage.py JSON counts, rather than its combined line/branch score."""
    totals = report["totals"]
    statements, covered = totals["num_statements"], totals["covered_lines"]
    if not report["files"] or statements <= 0 or not 0 <= covered <= statements:
        raise ValueError("coverage report is empty or has invalid line counts")
    if not report["meta"]["branch_coverage"]:
        raise ValueError("collect branch coverage as well as line coverage")
    branches, covered_branches = totals["num_branches"], totals["covered_branches"]
    if branches <= 0 or not 0 <= covered_branches <= branches:
        raise ValueError("coverage report has no branches or invalid branch counts")
    lines = 100.*covered/statements
    print(f"Python lines: {lines:.2f}% ({covered}/{statements}); required: {minimum:g}%")
    print(f"Python branches: {100.*covered_branches/branches:.2f}% "
          f"({covered_branches}/{branches}); diagnostic only")
    return lines >= minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    return 0 if check_report(json.loads(args.report.read_text(encoding="utf-8"))) else 1


if __name__ == "__main__":
    raise SystemExit(main())
