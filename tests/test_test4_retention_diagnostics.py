"""Bounded text-only tests; no archive, solver, NRG, Mathematica or HDF5 needed."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from NRG_comparisons.test4 import nrg, retention_diagnostics as diagnostics


NUMERICAL = dict(nmax=20, keep=1000, keepmin=200, keepenergy=14.)
FOOTER = """
Total energy: -317.180107157130976
** Iteration completed.
Peak usage: 3'340 MB
Timing report
         All: 130.927 s
        diag: 126.966 s
      recalc: 2.440 s
       Other: 1.522 s
"""


def shell(iteration=1, nmax=2, *, kept=1000, total=16000, boundary=4., scale=2., added=0):
    return (f" ***** [NRG] Iteration {iteration}/{nmax} (scale {scale}) ***** \n"
            + (f"Safeguard: keep additional {added} states\n" if added else "")
            + f"Cmax={boundary} nrkept={kept} nrkeptmult={kept} nrall={total} nrallmult={total} \n"
            + f"Kept: {kept} out of {total}, ratio={kept / total:.3g}\n")


def pilot_text():
    # Hand-literal diagnostic columns from pilot-v2/nrg-k1000/nrg.log, not a
    # fixture dependency on ignored runs. Format provenance is in the parser.
    boundaries = [6.1313553228803315, 4.1003166226705554, 4.9428836617780636,
                  4.3079666149401739, 5.0825466927774325, 5.2102278665395882,
                  7.2164834493536238, 10.177129370579411, 14.181992760372491,
                  20.093151963066752, 38.78456243989703, 76.642174480940781,
                  151.75080923642801, 302.9162765532339, 605.7092679750599,
                  1211.3706798877304, 2422.7269624637543, 4391.8094113685866,
                  8570.46148069632, 17133.511584453365]
    counts = [(1000, 1024, 0), (1002, 16000, 2), (1000, 16032, 0),
              *[(1000, 16000, 0)] * 5, (444, 16000, 0), (200, 7104, 0),
              (200, 3200, 0), (200, 3200, 0), (202, 3200, 2), (200, 3232, 0),
              (200, 3200, 0), (200, 3200, 0), (210, 3200, 10), (200, 3360, 0),
              (200, 3200, 0), (3200, 3200, 3)]
    text = "NRG Ljubljana\ndiagratio=1 *\ndm=0 *\nsymtype=SPU1 *\nkeepall= *\nlastall=1 *\n"
    for i, (boundary, (kept, total, added)) in enumerate(zip(boundaries, counts, strict=True), 1):
        text += shell(i, 20, kept=kept, total=total, boundary=boundary,
                      scale=54.10106403333613 / 2**(i - 1), added=added)
    # Cumulative time need not be present at every shell.
    return text + "Memory used: 3'340 MB Time elapsed: 130.873 s\n" + FOOTER


def test_twenty_shell_pilot_columns_and_semantics():
    result = diagnostics.parse_log(pilot_text(), NUMERICAL)
    assert result["complete"] is True
    assert [row["iteration"] for row in result["shells"]] == list(range(1, 21))
    assert result["mechanism_counts"] == {"keep_cap": 8, "energy_window": 1,
                                          "keepmin_floor_or_boundary": 10, "lastall_override": 1}
    assert result["max_pretruncation_dimension"] == 16032
    assert result["largest_block_dimension"] is None
    assert "sum" in result["largest_block_dimension_reason"]
    assert result["solver_vmpeak_gib"] == 3340 / 1024
    assert result["elapsed_seconds"] == 130.873
    assert result["timings_seconds"] == {"All": 130.927, "diag": 126.966, "recalc": 2.440, "Other": 1.522}
    first, second, *_, final = result["shells"]
    assert first["pretruncation_reduced_dimension"] == first["pretruncation_physical_dimension"] == 1024
    assert first["kept_reduced_states"] == first["kept_physical_states"] == 1000
    assert first["scale_over_gap"] == 54.10106403333613
    assert first["selection_boundary_shell"] == 6.1313553228803315
    assert first["selection_boundary_over_gap"] == pytest.approx(331.712848529371)
    assert first["window_over_gap"] == 14 * first["scale_over_gap"]
    assert first["elapsed_seconds"] is None
    assert second["safeguard_added_reduced_states"] == 2
    assert second["mechanism"] == "keep_cap"
    assert final["safeguard_added_reduced_states"] == 3
    assert final["mechanism"] == "lastall_override"
    assert final["selection_boundary_shell"] == 17133.511584453365
    assert all(row["exact_window_count"] is None for row in result["shells"])
    assert all("Unknown" in row["exact_window_count_reason"] for row in result["shells"])
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("keep,keepmin", [(1500, 300), (2000, 300), (2000, 400)])
def test_first_shell_available_space_and_final_override(keep, keepmin):
    numerical = dict(nmax=2, keep=keep, keepmin=keepmin, keepenergy=14.)
    text = shell(kept=1024, total=1024, boundary=7.) + shell(2, kept=16384, total=16384, added=3) + FOOTER
    result = diagnostics.parse_log(text, numerical)
    assert result["complete"]
    assert [row["mechanism"] for row in result["shells"]] == ["available_space_exhausted", "lastall_override"]


@pytest.mark.parametrize("kept,total,added,boundary,settings,mechanism", [
    (1002, 16000, 2, 4., {}, "keep_cap"),
    (210, 3200, 10, 2422., {}, "keepmin_floor_or_boundary"),
    (444, 16000, 0, 14.18, {}, "energy_window"),
    (1002, 16000, 2, 4., {"keepenergy": -1.}, "fixed_count"),
    (1000, 16000, 0, 4., {"keepenergy": 0.}, "fixed_count"),
    (1005, 16000, 2, 4., {}, "ambiguous"),  # Unreported boundary ties, not invented window counts.
    (300, 16000, 0, 4., {}, "ambiguous"),
    (1000, 16000, 0, 4., {"keepmin": 1000}, "ambiguous"),
])
def test_mechanism_boundaries(kept, total, added, boundary, settings, mechanism):
    numerical = dict(NUMERICAL, nmax=1, lastall=False, **settings)
    row = diagnostics.parse_log(shell(1, 1, kept=kept, total=total, added=added, boundary=boundary) + FOOTER,
                                numerical)["shells"][0]
    assert row["mechanism"] == mechanism
    assert row["exact_window_count"] is None
    if numerical["keepenergy"] <= 0:
        assert row["window_over_gap"] is None


def test_zero_based_keepall_list_and_lastall_priority():
    text = shell(1, 3, kept=1024, total=1024) + shell(2, 3) + shell(3, 3, kept=16000, total=16000) + FOOTER
    numerical = dict(NUMERICAL, nmax=3, keepall="0,2")
    result = diagnostics.parse_log(text, numerical)
    assert [row["mechanism"] for row in result["shells"]] == ["keepall_override", "keep_cap", "lastall_override"]
    result = diagnostics.parse_log(text, dict(numerical, lastall=False))
    assert result["shells"][-1]["mechanism"] == "keepall_override"


@pytest.mark.parametrize("memory", ["3'340 MB", "3'340 MiB", "3340 MiB", "1'003'340 MB"])
def test_memory_units_apostrophes_and_cumulative_not_summed(memory):
    text = (shell() + "Memory used: 1'678 MB Time elapsed: 0.053 s\n"
            + shell(2, kept=16000) + f"Memory used: {memory}\nTime elapsed: 130.873 s\n" + FOOTER)
    result = diagnostics.parse_log(text, dict(NUMERICAL, nmax=2))
    assert result["solver_vmpeak_gib"] == int(memory.split()[0].replace("'", "")) / 1024
    assert result["elapsed_seconds"] == 130.873
    assert [row["elapsed_seconds"] for row in result["shells"]] == [.053, 130.873]
    assert not any("rss" in key.lower() for key in result)


def test_omitted_timing_components_and_memory_remain_unknown():
    text = shell(1, 1, kept=1024, total=1024) + "** Iteration completed.\nTiming report\nAll: 1 s\nOther: 1 s\n"
    result = diagnostics.parse_log(text, dict(NUMERICAL, nmax=1), gap=2.)
    assert result["complete"]
    assert result["timings_seconds"] == {"All": 1., "diag": None, "recalc": None, "Other": 1.}
    assert result["solver_vmpeak_gib"] is None
    assert result["elapsed_seconds"] is None
    row = result["shells"][0]
    assert row["scale_over_gap"] == 1.
    assert row["selection_boundary_over_gap"] == 4.
    assert row["window_over_gap"] == 14.


@pytest.mark.parametrize("text", [
    "", " ***** [NRG] Iteration 1/20 (scale 2) *****\n",
    shell(1, 20), pilot_text().replace("** Iteration completed.", "").split("Timing report")[0],
    pilot_text().split("Timing report")[0], pilot_text().replace("         All: 130.927 s", ""),
    pilot_text().replace("       Other: 1.522 s", ""),
    pilot_text().replace("Kept: 1002 out of 16000, ratio=0.0626\n", ""),
])
def test_incomplete_logs(text):
    result = diagnostics.parse_log(text, NUMERICAL)
    assert result["complete"] is False
    if result["shells"] and result["shells"][0]["kept_reduced_states"] is None:
        assert result["shells"][0]["selection_boundary_over_gap"] is None
        assert result["shells"][0]["mechanism"] == "ambiguous"


@pytest.mark.parametrize("old,new,match", [
    ("Iteration 1/20", "Iteration 0/20", "iteration counter"),
    ("Iteration 2/20", "Iteration 1/20", "iteration counter"),
    ("Iteration 2/20", "Iteration 3/20", "iteration counter"),
    ("Iteration 1/20", "Iteration 1/19", "iteration counter"),
    ("Iteration 1/20", "Iteration 1/", "header"),
    ("scale 54.10106403333613", "scale nan", "header"),
    ("scale 54.10106403333613", "scale 0", "positive"),
    ("nrkeptmult=1000", "nrkeptmult=2000", "SPU1"),
    ("nrallmult=1024", "nrallmult=2048", "SPU1"),
    ("nrall=1024", "nrall=999", "SPU1"),
    ("nrkept=1000", "nrkept=0", "SPU1"),
    ("nrallmult=1024", "", "Cmax/count"),
    ("Kept: 1000 out of 1024", "Kept: 999 out of 1024", "crosscheck"),
    ("Kept: 1000 out of 1024", "Kept: 1000 out of 1025", "crosscheck"),
    ("Safeguard: keep additional 2 states", "Safeguard: keep additional 1002 states", "safeguard"),
    ("Safeguard: keep additional 2 states", "Safeguard: keep additional 0 states", "Safeguard"),
    ("[NRG]", "[DM]", "unsupported"),
    ("diagratio=1 *", "diagratio=0.5 *", "diagratio"),
    ("dm=0 *", "dm=1 *", "dm="),
    ("symtype=SPU1 *", "symtype=SU2 *", "symtype"),
    ("keepall= *", "keepall=0 *", "keepall"),
    ("lastall=1 *", "lastall=0 *", "lastall"),
    ("3'340 MB", "3'34 MB", "memory"),
    ("Time elapsed: 130.873 s", "Time elapsed: nan s", "elapsed"),
    ("All: 130.927 s", "All: -1 s", "nonnegative"),
])
def test_reject_inconsistent_or_malformed_records(old, new, match):
    with pytest.raises(ValueError, match=match):
        diagnostics.parse_log(pilot_text().replace(old, new, 1), NUMERICAL)


@pytest.mark.parametrize("extra", [
    "Restarting this iteration step. diagratio=1\n",
    "Insufficient number of states computed.\n",
    "Cmax=4 nrkept=1000 nrkeptmult=1000 nrall=16000 nrallmult=16000\n",
    "Kept: 1000 out of 16000, ratio=0.0625\n",
])
def test_reject_retry_or_duplicate_count_instead_of_selecting_an_attempt(extra):
    with pytest.raises(ValueError, match="retries|duplicate"):
        diagnostics.parse_log(shell() + extra, dict(NUMERICAL, nmax=2))


@pytest.mark.parametrize("text,match", [
    ("Cmax=4 nrkept=1 nrkeptmult=1 nrall=16 nrallmult=16", "outside"),
    (" ***** [NRG] Iteration 1/2 (scale 2) *****\nKept: 1 out of 16, ratio=0.0625", "missing preceding Cmax"),
    (shell() + "Time elapsed: 2 s\nTime elapsed: 1 s\n", "decreased"),
    (shell() + FOOTER + "** Iteration completed.\n", "duplicate completion"),
    (shell() + FOOTER + "All: 131 s\n", "duplicate timing"),
    (shell() + FOOTER + "Timing report\n", "Timing report"),
    (shell() + FOOTER + shell(2), "iteration counter"),
    (shell() + "Timing report\n", "Timing report"),
    (shell(2), "iteration counter"),
])
def test_reject_record_order_and_multiple_passes(text, match):
    with pytest.raises(ValueError, match=match):
        diagnostics.parse_log(text, dict(NUMERICAL, nmax=2))


@pytest.mark.parametrize("settings", [dict(lastall=True), dict(lastall=False, keepall="0")])
def test_overrides_require_all_states(settings):
    with pytest.raises(ValueError, match="must retain all"):
        diagnostics.parse_log(shell(1, 1) + FOOTER, dict(NUMERICAL, nmax=1, **settings))


@pytest.mark.parametrize("settings", [
    {"nmax": 0}, {"keep": True}, {"keepmin": -1}, {"keepmin": 1001}, {"keepenergy": float("nan")},
    {"lastall": "false"}, {"keepall": "even"}, {"keepall": "1,1"}, {"keepall": "20"},
    {"keepall": True}, {"untruncated": True},
])
def test_invalid_settings(settings):
    with pytest.raises(ValueError):
        diagnostics.parse_log("", dict(NUMERICAL, **settings))


@pytest.mark.parametrize("gap", [0, -1, True, "1", float("nan"), float("inf")])
def test_invalid_gap(gap):
    with pytest.raises(ValueError, match="gap"):
        diagnostics.parse_log("", NUMERICAL, gap=gap)


def write_run(directory, *, untruncated=False, units="gap"):
    physical = dict(geometry="junction", gap=1., bandwidth=100., u=2., gamma=.4,
                    nu=.75, field=.03, phi=1.8849555921538759, wn=0., ws=0.)
    numerical = dict(NUMERICAL, **{"lambda": 4., "z": 1., "untruncated": untruncated})
    (directory / "param").write_text(nrg.render_param(physical, numerical, units=units))
    (directory / "nrg.log").write_text(pilot_text())


def test_strict_reader_uses_validated_deck_and_log(tmp_path):
    write_run(tmp_path)
    result = diagnostics.read_diagnostics(tmp_path)
    assert diagnostics.nrg is nrg
    assert result == diagnostics.parse_log(pilot_text(), NUMERICAL)
    (tmp_path / "nrg.log").write_text(shell(1, 20))
    with pytest.raises(ValueError, match="incomplete nrg.log.*20"):
        diagnostics.read_diagnostics(tmp_path)
    (tmp_path / "param").write_text((tmp_path / "param").read_text().replace("diagratio=1", "diagratio=0.5"))
    with pytest.raises(ValueError, match="diagratio"):
        diagnostics.read_diagnostics(tmp_path)


def test_reader_uses_deck_gap_and_actual_keepall(tmp_path):
    write_run(tmp_path, untruncated=True, units="bandwidth")
    text = "".join(shell(i, 20, kept=1024, total=1024, scale=.02) for i in range(1, 21)) + FOOTER
    (tmp_path / "nrg.log").write_text(text)
    result = diagnostics.read_diagnostics(tmp_path)
    assert result["mechanism_counts"] == {"keepall_override": 19, "lastall_override": 1}
    assert result["shells"][0]["scale_over_gap"] == 2.
    assert result["shells"][0]["selection_boundary_over_gap"] == 8.


def test_missing_input_files(tmp_path):
    with pytest.raises(FileNotFoundError):
        diagnostics.read_diagnostics(tmp_path)
    write_run(tmp_path)
    (tmp_path / "nrg.log").unlink()
    with pytest.raises(FileNotFoundError):
        diagnostics.read_diagnostics(tmp_path)


def test_standalone_import_is_stdlib_only(tmp_path):
    # No site packages or checkout imports available; direct-import fallback
    # must work with only these two source files, as in an exported case.
    for module in (nrg, diagnostics):
        source = Path(module.__file__)
        (tmp_path / source.name).write_bytes(source.read_bytes())
    check = subprocess.run([sys.executable, "-B", "-S", "-c",
                            "import retention_diagnostics as d; "
                            "assert d.nrg.__name__ == 'nrg'; "
                            "assert not d.parse_log('', dict(nmax=1, keep=2, keepmin=0, keepenergy=14))['complete']"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=15)
    assert check.returncode == 0, check.stderr
