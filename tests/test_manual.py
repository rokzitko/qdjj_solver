"""The PDF checker enforces diagnostics, without imposing a document length limit."""

import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.docs

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_manual", ROOT/"tools"/"check_manual.py")
check_manual = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_manual)


@pytest.mark.parametrize("pages", [1, 100, 101, 10000])
def test_clean_manual_has_no_page_ceiling(tmp_path, pages):
    (tmp_path/"manual.log").write_text(f"Output written on manual.pdf ({pages} pages, 12345 bytes).\n")
    assert check_manual.problems(tmp_path) == ([], pages)


@pytest.mark.parametrize("diagnostic", [
    "LaTeX Warning: Reference `eq:missing' undefined.",
    "LaTeX Warning: Citation `missing' undefined.",
    "Package biblatex Warning: Please rerun Biber.",
    "Class scrbook Warning: unsupported option.",
    r"Overfull \hbox (4.0pt too wide) in paragraph at lines 3--4",
])
def test_long_manual_still_rejects_typesetting_problems(tmp_path, diagnostic):
    (tmp_path/"manual.log").write_text(diagnostic+"\nOutput written on manual.pdf (1000 pages, 12345 bytes).\n")
    assert check_manual.problems(tmp_path) == ([diagnostic], 1000)


def test_successful_output_and_bibliography_are_checked(tmp_path):
    (tmp_path/"manual.log").write_text("LaTeX stopped before writing output.\n")
    (tmp_path/"manual.blg").write_text("INFO - Reading sources\nWARN - Missing entry\nERROR - Bad input\n")
    assert check_manual.problems(tmp_path) == (["no PDF written", "WARN - Missing entry", "ERROR - Bad input"], 0)


def test_cli_reports_pages_and_propagates_diagnostics(tmp_path, monkeypatch, capsys):
    inspect = check_manual.problems
    monkeypatch.setattr(check_manual, "problems", lambda *args: inspect(tmp_path))
    monkeypatch.setattr(check_manual, "artifact_problems", lambda *args, **kwargs: [])
    (tmp_path/"manual.log").write_text("Output written on manual.pdf (1000 pages, 12345 bytes).\n")
    assert check_manual.main([]) == 0
    assert "1000 pages; 0 problem(s)" in capsys.readouterr().out
    (tmp_path/"manual.blg").write_text("WARN - Missing entry\n")
    assert check_manual.main([]) == 1
    assert "Missing entry" in capsys.readouterr().out
    with pytest.raises(SystemExit) as caught:
        check_manual.main(["--max-pages", "100"])
    assert caught.value.code == 2
