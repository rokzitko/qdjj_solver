"""Read-only retention diagnostics for test4's full-diagonalization SPU1 runs.

Conventions and pilot scope: README.md and README.model_and_parameters.md in
this directory. Format/semantics follow NRG Ljubljana's c++/step.hpp
(Step::infostring, one-based labels), truncation.hpp (highest_retained,
truncate_prepare and truncate_stats), core.hpp (do_diag retries), time_mem.hpp
(MemTime, Timing), and portabil.hpp (Linux memoryused reads VmPeak).

Cmax is the selection boundary in shell units, NOT necessarily the largest
retained energy: lastall overrides the selected count afterwards. Selection
uses count(E <= keepenergy) + 1, then keepmin/keep clamps, then safeguards;
ties can retain still more states. The log does not give the exact window
count, so no inverse clamp or n_window estimate is reported. Mechanisms label
the observed count boundary, not proof of a unique unclamped window count.
Silent logs do not report symmetry-block dimensions. nrall is their SUM.

Only stdlib imports; no solver, HDF5, NRG or Mathematica execution. This is not
a general NRG/DM/partial-diagonalization or checkpoint/retry log reader.
"""

from collections import Counter
import math
from pathlib import Path
import re

if __package__:
    from . import nrg
else:
    import nrg


NUMBER = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
HEADER = re.compile(rf"\*+\s*\[NRG\]\s+Iteration ([0-9]+)/([0-9]+) \(scale ({NUMBER})\)\s*\*+")
COUNTS = re.compile(rf"Cmax=({NUMBER})\s+nrkept=([0-9]+)\s+nrkeptmult=([0-9]+)\s+"
                    r"nrall=([0-9]+)\s+nrallmult=([0-9]+)")
KEPT = re.compile(rf"Kept:\s+([0-9]+) out of ([0-9]+), ratio=({NUMBER})")
MEMORY = re.compile(r"(?:Memory used|Peak usage):\s+([0-9]+(?:'[0-9]{3})*)\s+(?:MiB|MB)(?:\s+|$)")
WINDOW_REASON = ("Unknown: the log reports Cmax and post-selection counts, not count(E <= keepenergy); "
                 "the +1 boundary state, clamps, safeguards and ties cannot in general be inverted.")


def parse_log(text, numerical, gap=1.0) -> dict:
    """Parse a single full NRG pass, returning JSON-compatible diagnostics.

    numerical requires typed nmax, keep, keepmin, keepenergy. Optional lastall
    defaults to True (test4's deck); keepall is the actual comma-separated,
    ZERO-based param string, default empty, never inferred from untruncated.
    Scale and Cmax*scale are divided by the positive gap in deck units.

    Truncated logs return complete=False, with missing shell fields set to
    None. Malformed/duplicate records, invalid counters, inconsistent SPU1
    counts, DM passes and retries raise ValueError. Completion requires all
    nmax shells with Cmax and matching Kept records, the completion marker,
    and the timing footer's All/Other. diag/recalc can legitimately be absent
    (NRG folds components <=1% into Other); missing times remain None.

    elapsed_seconds is cumulative, never summed across shells. Memory reports
    labelled MB or MiB are Linux VmPeak in MiB, NOT RSS, converted by /1024.
    Neither this log nor complete=True certifies numerical convergence or
    the integrity of separate output artifacts such as raw.h5 and DONE.
    """
    for key, minimum in (("nmax", 1), ("keep", 1), ("keepmin", 0)):
        if type(numerical.get(key)) is not int or numerical[key] < minimum:
            raise ValueError(f"numerical.{key} must be an integer >= {minimum}")
    nmax, keep, keepmin = (numerical[key] for key in ("nmax", "keep", "keepmin"))
    if keepmin > keep:
        raise ValueError("keepmin must not exceed keep")
    if (isinstance(gap, bool) or not isinstance(gap, (int, float))
            or not math.isfinite(gap) or gap <= 0):
        raise ValueError("gap must be finite and positive")
    window = numerical.get("keepenergy")
    if isinstance(window, bool) or not isinstance(window, (int, float)) or not math.isfinite(window):
        raise ValueError("numerical.keepenergy must be finite and real")
    lastall = numerical.get("lastall", True)
    keepall = numerical.get("keepall", "")
    if type(lastall) is not bool:
        raise ValueError("lastall must be boolean")
    if numerical.get("untruncated") and "keepall" not in numerical:
        raise ValueError("untruncated requires the explicit zero-based keepall param list")
    if not isinstance(keepall, str) or not re.fullmatch(r"(?:[0-9]+(?:,[0-9]+)*)?", keepall):
        raise ValueError("keepall must be an empty or comma-separated zero-based iteration list")
    keepall_steps = [int(value) for value in keepall.split(",") if value]
    if len(set(keepall_steps)) != len(keepall_steps) or any(i >= nmax for i in keepall_steps):
        raise ValueError("keepall contains duplicate or out-of-range zero-based iterations")

    def finite(value, label):
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{label} must be finite and nonnegative")
        return value

    shells, records = [], []
    row, seen = None, set()
    completed = timing_report = False
    timings = dict.fromkeys(("All", "diag", "recalc", "Other"))
    vmpeak, elapsed = None, None
    # Validate relevant optional parameter echoes too, rather than trusting a
    # full-diagonalization deck when the supplied log belongs to another run.
    expected = {"Ninit": "0", "symtype": "SPU1", "diagratio": "1", "dm": "0", "fdm": "0",
                "resume": "0", "substeps": "0", "floquet": "0", "keepall": keepall,
                "lastall": str(int(lastall)), "keep": str(keep), "keepmin": str(keepmin),
                "keepenergy": str(window), "Nmax": str(nmax)}
    for line_number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if any(token in line for token in ("[DM]", "Restarting this iteration", "Insufficient number of states computed",
                                            "Partial diagonalisation retry", "Loading from checkpoint")):
            raise ValueError(f"line {line_number}: DM, partial diagonalization and retries are unsupported")
        echo = re.fullmatch(r"([A-Za-z0-9_]+)=(.*?)\s*(?:\*)?", line)
        if echo and echo[1] in expected:
            key, value = echo[1], echo[2].strip()
            target = expected[key]
            if key not in ("symtype", "keepall"):
                try:
                    matches = float(value) == float(target)
                except ValueError:
                    matches = False
            else:
                matches = value == target
            if not matches:
                raise ValueError(f"line {line_number}: unsupported or mismatched {key}={value}")
        if "[NRG]" in line or (line.startswith("*") and "Iteration" in line and "/" in line):
            match = HEADER.fullmatch(line)
            if not match:
                raise ValueError(f"line {line_number}: malformed NRG iteration header")
            iteration, total = int(match[1]), int(match[2])
            if completed or timing_report or total != nmax or iteration != len(shells) + 1 or iteration > nmax:
                raise ValueError(f"line {line_number}: invalid, duplicate or out-of-order iteration counter")
            scale = finite(match[3], "scale") / gap
            if not math.isfinite(scale) or scale <= 0:
                raise ValueError("scale_over_gap must be finite and positive")
            row = dict.fromkeys(("kept_reduced_states", "kept_physical_states", "pretruncation_reduced_dimension",
                                 "pretruncation_physical_dimension", "selection_boundary_shell",
                                 "selection_boundary_over_gap", "elapsed_seconds"))
            row.update(iteration=iteration, scale_over_gap=scale,
                       window_over_gap=finite(window * scale, "window_over_gap") if window > 0 else None,
                       safeguard_added_reduced_states=0, mechanism="ambiguous",
                       exact_window_count=None, exact_window_count_reason=WINDOW_REASON)
            shells.append(row)
            seen = set()
            records.append(seen)
            continue
        if line.startswith(("Cmax", "Kept:", "Safeguard:")):
            if row is None or completed or timing_report:
                raise ValueError(f"line {line_number}: retention record outside an active shell")
            kind = line.split(":" if ":" in line else "=", 1)[0]
            if kind in seen:
                raise ValueError(f"iteration {row['iteration']}: duplicate {kind}; retries are unsupported")
            seen.add(kind)
            if kind == "Safeguard":
                match = re.fullmatch(r"Safeguard: keep additional ([0-9]+) states", line)
                if not match or "Cmax" in seen or int(match[1]) == 0:
                    raise ValueError(f"line {line_number}: malformed or misplaced Safeguard record")
                row["safeguard_added_reduced_states"] = int(match[1])
            elif kind == "Cmax":
                match = COUNTS.fullmatch(line)
                if not match:
                    raise ValueError(f"line {line_number}: malformed Cmax/count record")
                boundary = finite(match[1], "Cmax")
                kept, physical, total, total_physical = map(int, match.groups()[1:])
                if not 0 < kept <= total or physical != kept or total_physical != total:
                    raise ValueError(f"iteration {row['iteration']}: inconsistent SPU1 reduced/physical counts")
                if row["safeguard_added_reduced_states"] >= kept:
                    raise ValueError(f"iteration {row['iteration']}: safeguard addition must be smaller than kept count")
                row.update(kept_reduced_states=kept, kept_physical_states=physical,
                           pretruncation_reduced_dimension=total, pretruncation_physical_dimension=total_physical,
                           selection_boundary_shell=boundary,
                           selection_boundary_over_gap=finite(boundary * row["scale_over_gap"], "boundary_over_gap"))
            else:
                match = KEPT.fullmatch(line)
                if not match or "Cmax" not in seen:
                    raise ValueError(f"line {line_number}: malformed Kept record or missing preceding Cmax")
                if (int(match[1]), int(match[2])) != (row["kept_reduced_states"], row["pretruncation_reduced_dimension"]):
                    raise ValueError(f"iteration {row['iteration']}: Kept crosscheck disagrees with Cmax counts")
                if finite(match[3], "Kept ratio") > 1:
                    raise ValueError("Kept ratio must be <= 1")
            continue
        if line == "** Iteration completed.":
            if completed:
                raise ValueError("duplicate completion marker; multiple passes are unsupported")
            completed = True
        if line.startswith(("Memory used:", "Peak usage:")):
            match = MEMORY.match(line)
            if not match:
                raise ValueError(f"line {line_number}: malformed memory report")
            value = int(match[1].replace("'", "")) / 1024
            vmpeak = max(vmpeak or 0, value)
            line = line[match.end():].strip()
            if line and not line.startswith("Time elapsed:"):
                raise ValueError(f"line {line_number}: malformed memory/time report")
        if line.startswith("Time elapsed:"):
            match = re.fullmatch(rf"Time elapsed:\s+({NUMBER})\s+s", line)
            if not match:
                raise ValueError(f"line {line_number}: malformed elapsed time")
            value = finite(match[1], "elapsed time")
            if elapsed is not None and value < elapsed:
                raise ValueError("cumulative elapsed time decreased; retries are unsupported")
            elapsed = value
            if row is not None:
                row["elapsed_seconds"] = value
        if line == "Timing report":
            if timing_report or not completed:
                raise ValueError("duplicate or misplaced Timing report")
            timing_report = True
        if re.match(r"(?:All|diag|recalc|Other):", line):
            match = re.fullmatch(rf"(All|diag|recalc|Other):\s+({NUMBER})\s+s", line)
            if not match or not timing_report:
                raise ValueError(f"line {line_number}: malformed or misplaced timing component")
            if timings[match[1]] is not None:
                raise ValueError(f"duplicate timing component {match[1]}")
            timings[match[1]] = finite(match[2], "timing component")

    for row, seen in zip(shells, records, strict=True):
        if not {"Cmax", "Kept"} <= seen:
            continue
        kept, total = row["kept_reduced_states"], row["pretruncation_reduced_dimension"]
        added = row["safeguard_added_reduced_states"]
        base = kept - added
        if row["iteration"] == nmax and lastall:
            mechanism = "lastall_override"
        elif row["iteration"] - 1 in keepall_steps:
            mechanism = "keepall_override"
        elif kept == total and total <= keep and not added:
            mechanism = "available_space_exhausted"
        elif window <= 0 and base == keep:
            mechanism = "fixed_count"
        elif window > 0 and base == keep and keepmin < keep:
            mechanism = "keep_cap"
        elif window > 0 and base == keepmin and keepmin < keep:
            mechanism = "keepmin_floor_or_boundary"
        elif window > 0 and keepmin < base < keep and row["selection_boundary_shell"] >= window:
            mechanism = "energy_window"
        else:
            mechanism = "ambiguous"
        if mechanism in ("lastall_override", "keepall_override") and kept != total:
            raise ValueError(f"iteration {row['iteration']}: {mechanism} must retain all states")
        row["mechanism"] = mechanism

    return {"complete": (len(shells) == nmax and all({"Cmax", "Kept"} <= seen for seen in records)
                         and completed and timings["All"] is not None and timings["Other"] is not None),
            "shells": shells, "mechanism_counts": dict(Counter(row["mechanism"] for row in shells)),
            "solver_vmpeak_gib": vmpeak, "elapsed_seconds": elapsed, "timings_seconds": timings,
            "max_pretruncation_dimension": max((row["pretruncation_reduced_dimension"] for row in shells
                                                if row["pretruncation_reduced_dimension"] is not None), default=None),
            "largest_block_dimension": None,
            "largest_block_dimension_reason": "Unknown: silent test4 logs do not report symmetry-block dimensions; nrall is a sum."}


def read_diagnostics(directory) -> dict:
    """Strictly read param and nrg.log; raise unless all requested shells succeeded.

    nrg.read_param validates test4's fixed full-diagonalization, SPU1, no-DM
    deck. The actual keepall string is used, not the lossy untruncated flag.
    Missing files raise OSError; incomplete or inconsistent logs raise ValueError.
    """
    directory = Path(directory)
    deck = nrg.read_param(directory / "param")
    numerical = dict(deck["numerical"], lastall=deck["param"]["lastall"] == "true",
                     keepall=deck["param"]["keepall"])
    result = parse_log((directory / "nrg.log").read_text(), numerical, gap=deck["physical"]["gap"])
    if not result["complete"]:
        raise ValueError(f"incomplete nrg.log: expected {numerical['nmax']} successful shells and completion/timing footer")
    return result
