"""Tutorials in examples/: every quick profile reproduces its reference (marker ``tutorial``)."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_tutorials", ROOT/"tools"/"check_tutorials.py")
check_tutorials = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_tutorials)


@pytest.mark.tutorial
@pytest.mark.timeout(600)
@pytest.mark.parametrize("directory", [entry["dir"] for entry in check_tutorials.catalog()])
def test_quick_profile_matches_reference(directory, tmp_path):
    assert check_tutorials.check(directory, "quick", tmp_path/directory) == []


@pytest.mark.docs
def test_reference_comparison():
    reference = dict(values=dict(energy=dict(value=1., atol=1e-6), residual=dict(max=1e-9),
                                 weight=dict(value=.5, rtol=.1)))
    assert check_tutorials.compare(reference, dict(energy=1. + 1e-7, residual=1e-10, weight=.52)) == []
    problems = check_tutorials.compare(reference, dict(energy=1.1, residual=1e-8, weight=float("nan")))
    assert len(problems) == 3
    assert "missing" in check_tutorials.compare(reference, dict(energy=1., residual=0.))[0]
