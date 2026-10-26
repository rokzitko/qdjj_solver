"""Small checks for the coverage gate and distribution boundary."""

import importlib.util
import io
from pathlib import Path
import re
import stat
import tarfile
import tomllib
from urllib.parse import urlsplit
import zipfile

import pytest


NRG_SOURCE_ASSETS = (
    "NRG_comparisons/README.md",
    "NRG_comparisons/test1/README.md",
    "NRG_comparisons/test1/run.py",
    "NRG_comparisons/test1/nrg.py",
    "NRG_comparisons/test1/control.py",
    "NRG_comparisons/test1/analysis.py",
    "NRG_comparisons/test1/full_analysis.py",
    "NRG_comparisons/test1/targeted_nrg.py",
    "NRG_comparisons/test1/nrg_dataset.py",
    "NRG_comparisons/test1/twist_analysis.py",
    "NRG_comparisons/test1/focused_internal.py",
    "NRG_comparisons/test1/generate.py",
    "NRG_comparisons/test1/input/physical.json",
    "NRG_comparisons/test1/input/profiles.json",
    "NRG_comparisons/test1/input/nrg-targeted.json",
    "NRG_comparisons/test1/input/qp-dmrg-focused.json",
    "NRG_comparisons/test1/input/final-baths.json",
    "NRG_comparisons/test1/output/results.json",
    "NRG_comparisons/test1/output/observables.csv",
    "NRG_comparisons/test1/output/comparisons.csv",
    "NRG_comparisons/test1/output/comparison.svg",
    "NRG_comparisons/test1/output/convergence.svg",
    "NRG_comparisons/test2/README.md",
    "NRG_comparisons/test2/generate.py",
    "NRG_comparisons/test2/input/physical.json",
    "NRG_comparisons/test2/input/baths.json",
    "NRG_comparisons/test2/reference/provenance.json",
    "NRG_comparisons/test2/output/results.json",
    "NRG_comparisons/test2/output/observables.csv",
    "NRG_comparisons/test2/output/comparisons.csv",
    "NRG_comparisons/test2/output/comparison.svg",
    "NRG_comparisons/test2/output/convergence.svg",
    "NRG_comparisons/test2/manual_run/README.md",
    "NRG_comparisons/test2/manual_run/.gitignore",
    "NRG_comparisons/test2/manual_run/param.loop",
    "NRG_comparisons/test2/manual_run/1_run",
    "NRG_comparisons/test2/manual_run/run_one",
    "NRG_comparisons/test2/manual_run/2_proc",
    "NRG_comparisons/test2/manual_run/process.py",
    "NRG_comparisons/test2/manual_run/extract.py",
    "NRG_comparisons/test2/manual_run/reference.json",
    "NRG_comparisons/test3/README.md",
    "NRG_comparisons/test3/run.py",
    "NRG_comparisons/test3/finite.py",
    "NRG_comparisons/test3/analyze.py",
    "NRG_comparisons/test3/generate.py",
    "NRG_comparisons/test3/input/physical.json",
    "NRG_comparisons/test3/input/baths.json",
    "NRG_comparisons/test3/reference/provenance.json",
    "NRG_comparisons/test3/output/report.md",
    "NRG_comparisons/test3/output/results.json",
    "NRG_comparisons/test3/output/observables.csv",
    "NRG_comparisons/test3/output/comparisons.csv",
    "NRG_comparisons/test3/manual_run/README.md",
    "NRG_comparisons/test3/manual_run/.gitignore",
    "NRG_comparisons/test3/manual_run/nrg.py",
    "NRG_comparisons/test3/manual_run/physics.py",
    "NRG_comparisons/test3/manual_run/process.py",
    "NRG_comparisons/test3/manual_run/param.loop",
    "NRG_comparisons/test3/manual_run/1_run",
    "NRG_comparisons/test3/manual_run/2_proc",
    "NRG_comparisons/test3/manual_run/run_one",
    "NRG_comparisons/test4/model.m",
    "NRG_comparisons/test4/modeloperators.m",
    "NRG_comparisons/test4/nrg.py",
    "NRG_comparisons/test4/physics.py",
    "NRG_comparisons/test4/finite.py",
    "NRG_comparisons/test4/run.py",
    "NRG_comparisons/test4/manual.py",
    "NRG_comparisons/test4/generate.py",
    "NRG_comparisons/test4/retention_diagnostics.py",
    "NRG_comparisons/test4/dmrg_states.py",
    "NRG_comparisons/test4/continuum_finite.py",
    "NRG_comparisons/test4/README.md",
    "NRG_comparisons/test4/FINAL_REPORT.md",
    "NRG_comparisons/test4/README.model_and_parameters.md",
    "NRG_comparisons/test4/.gitignore",
    "NRG_comparisons/test4/input/physical.json",
    "NRG_comparisons/test4/input/baths.json",
    "NRG_comparisons/test4/output/results.json",
    "NRG_comparisons/test4/output/report.md",
    "NRG_comparisons/test4/output/observables.csv",
    "NRG_comparisons/test4/output/comparisons.csv",
    "NRG_comparisons/test4/output/z_averages.csv",
)
TEST4_HISTORY_DIRS = (
    "retention", "high-retention", "followup", "continuum", "continuum-parallel",
    "continuum-continue", "continuum-lambda", "final",
)
NRG_PRIVATE_DIRS = (
    "private", "private_archive", "private-archive", "private_archives", "private-archives",
    "worktree", "worktrees/research", ".venv", "venv", "shared_venv", "shared-venv", ".git",
)
NRG_REJECTED_ASSETS = tuple(f"NRG_comparisons/test1/{name}" for name in (
    "residual_study.py", "input/unreviewed.json", "FULLTEST1.md", "output/pilot-v2.json",
    "nested/runs/results.json", "checkpoints/state.h5", "output/state.hdf5",
    "output/state.npz", "output/unreviewed.json",
)) + tuple(f"NRG_comparisons/test2/manual_run/{name}" for name in (
    "1/param", "8/raw.h5", "sweep.json", "processed/summary.json", "failed-3/data",
)) + tuple(f"NRG_comparisons/test3/{name}" for name in (
    "runs/pilot-v1/campaign.json", "runs/pilot-v1/chain-nrg/result.json",
    "runs/pilot-v1/chain-nrg/raw.h5", "checkpoints/state.h5", "output/state.npz",
    "output/unreviewed.json", "campaign.json", "output/campaign.json",
)) + tuple(f"NRG_comparisons/test3/manual_run/{name}" for name in (
    "1/param", "8/raw.h5", "sweep.json", "processed/per_z.csv",
    "processed/averages.csv", "processed/summary.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/pilot-v1/campaign.json", "runs/pilot-v1/junction-nrg/model.m",
    "runs/pilot-v1/junction-nrg/modeloperators.m", "runs/pilot-v1/junction-nrg/result.json",
    "runs/pilot-v1/junction-nrg/raw.h5", "checkpoints/state.h5", "output/state.npz",
    "output/unreviewed.json", "campaign.json", "output/campaign.json",
    "runs/retention-v1/campaign.json", "runs/retention-v1/controls.json", "runs/retention-v1/decisions.json",
    "runs/retention-v1/nrg-k1000/request.json", "runs/retention-v1/nrg-k1000/result.json",
    "runs/retention-v1/nrg-k1000/retention.json", "runs/retention-v1/nrg-k1500-m300/model.m",
    "runs/retention-v1/nrg-k1500-m300/modeloperators.m", "runs/retention-v1/nrg-k1500-m300/raw.h5",
    "runs/retention-v1/nrg-k1500-m300/nrg.log", "runs/retention-v1/nrg-k1500-m300/nrginit.log",
    "runs/retention-v1/nrg-k1500-m300/request.json", "runs/retention-v1/nrg-k1500-m300/result.json",
    "runs/retention-v1/nrg-k1500-m300/retention.json", "runs/retention-v1/output/results.json",
    "runs/retention-1500-2000-v1.session.json", "retention-1500-2000-v1.session.json",
    "output/retention-v1/results.json", "output/retention/snapshots/results.json",
    "output/retention/raw.h5", "output/retention/nrg.log", "output/retention/nrginit.log",
    "output/retention/campaign.json", "output/retention/retention.json",
    "output/retention/retention-1500-2000-v1.session.json", "output/retention/source/model.m",
    "output/retention/reused/nrg-k1000/result.json", "output/retention/controls/qp/result.json",
    "runs/high-retention-v1/campaign.json", "runs/high-retention-v1/controls.json",
    "runs/high-retention-v1/decisions.json", "runs/high-retention-v1/nrg-k2000-m400/request.json",
    "runs/high-retention-v1/nrg-k2000-m400/result.json", "runs/high-retention-v1/nrg-k2000-m400/retention.json",
    "runs/high-retention-v1/output/results.json", "runs/retention-3000-4000-v1.session.json",
    "retention-3000-4000-v1.session.json", "output/high-retention-v1/results.json",
    "output/high-retention/snapshots/results.json", "output/high-retention/raw.h5",
    "output/high-retention/raw.hdf5", "output/high-retention/nrg.log", "output/high-retention/nrginit.log",
    "output/high-retention/campaign.json", "output/high-retention/controls.json",
    "output/high-retention/decisions.json", "output/high-retention/request.json",
    "output/high-retention/result.json", "output/high-retention/retention.json", "output/high-retention/launch.json",
    "output/high-retention/retention-3000-4000-v1.session.json",
    "output/high-retention/model.m", "output/high-retention/modeloperators.m",
    "output/high-retention/source/model.m", "output/high-retention/source/modeloperators.m",
    "output/high-retention/reused/nrg-k2000-m400/result.json", "output/high-retention/controls/qp/result.json",
    "runs/next-priorities-v1.session.json", "next-priorities-v1.session.json",
    "runs/next-priorities-v1/output/results.json", "output/next-priorities-v1/results.json",
    "runs/next-priorities-v1/cosh6-qp-q4/request.json", "runs/next-priorities-v1/cosh6-qp-q4/result.json",
)) + tuple(f"NRG_comparisons/test4/runs/high-retention-v1/{task}/{name}"
           for task in ("nrg-k3000-m400", "nrg-k4000-m400") for name in (
    "model.m", "modeloperators.m", "param", "data", "raw.h5", "DONE", "nrg.log", "nrginit.log",
    "request.json", "result.json", "retention.json", "launch.json", "worker.log",
)) + tuple(f"NRG_comparisons/test4/{directory}/{name}"
           for directory in ("runs/next-priorities-v1", "output/followup") for name in (
    "campaign.json", "controls.json", "budgets.json", "decisions.json", "task_gates.json",
    "request.json", "result.json", "launch.json", "worker.log", "next-priorities-v1.session.json",
    "eigenstates.json", "eigenstates.mps.h5", "eigenstates.a1b2c3d4.mps.h5",
    "raw.h5", "raw.hdf5", "state.npz", "snapshots/results.json", "source/model.m",
    "controls/qp/request.json", "controls/dmrg/request.json", "reused/nrg-k4000-m400/result.json",
    "archive/pilot-v2/campaign.json", "archive/high-retention-v1/campaign.json",
)) + tuple(f"NRG_comparisons/test4/runs/next-priorities-v1/cosh4-dmrg-even-chi64/{name}" for name in (
    "request.json", "result.json", "launch.json", "worker.log",
    "eigenstates.json", "eigenstates.mps.h5", "eigenstates.a1b2c3d4.mps.h5",
)) + tuple(f"NRG_comparisons/test4/runs/next-priorities-v1/{task}/{name}"
           for task in ("nrg-length22-k4000", "nrg-window16-k4000") for name in (
    "model.m", "modeloperators.m", "param", "data", "raw.h5", "DONE", "nrg.log", "nrginit.log",
    "request.json", "result.json", "retention.json", "launch.json", "worker.log",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/continuum-v1.session.json", "continuum-v1.session.json", "runs/continuum-v1.controller.log",
    "runs/continuum-v1/campaign.json", "runs/continuum-v1/status.json", "runs/continuum-v1/clock.json",
    "runs/continuum-v1/decisions.json", "runs/continuum-v1/work/cosh6-dmrg-even-chi128/eigenstates.json",
    "runs/continuum-v1/work/cosh6-dmrg-even-chi128/eigenstates.a1b2c3d4.mps.h5",
    "runs/continuum-v1/work/cosh6-qp-q10/qp-sz+0.json",
    "runs/continuum-v1/work/nrg-z4-k4000/raw.h5", "runs/continuum-v1-presentation/results.json",
    "output/continuum/campaign.json", "output/continuum/clock.json", "output/continuum/eigenstates.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/continuum-parallel-v1.amendment.json", "runs/continuum-parallel-v1.controller.log",
    "runs/.continuum-cluster-mathematica.lock", "runs/continuum-parallel-v1/campaign.json",
    "runs/continuum-parallel-v1/clock.json", "runs/continuum-parallel-v1/status.json",
    "runs/continuum-parallel-v1/scheduler_preflight.json", "runs/continuum-parallel-v1/decisions.json",
    "runs/continuum-parallel-v1/jobs/nrg-z2-k4000.json",
    "runs/continuum-parallel-v1/jobs/nrg-z2-k4000.submission.json",
    "runs/continuum-parallel-v1/jobs/nrg-z2-k4000.out", "runs/continuum-parallel-v1/jobs/nrg-z2-k4000.err",
    "runs/continuum-parallel-v1/work/nrg-z2-k4000/raw.h5",
    "runs/continuum-parallel-v1/work/cosh8-dmrg-even-chi256/eigenstates.json",
    "runs/continuum-parallel-presentation/results.json",
    "output/continuum-parallel/campaign.json", "output/continuum-parallel/clock.json",
    "output/continuum-parallel/jobs/nrg-z2-k4000.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/continuum-continue-v1.amendment.json", "runs/continuum-continue-v1.controller.log",
    "runs/continuum-continue-v1/campaign.json", "runs/continuum-continue-v1/clock.json",
    "runs/continuum-continue-v1/jobs/nrg-z1-k5000.json",
    "runs/continuum-continue-v1/jobs/nrg-z1-k5000.submission.json",
    "runs/continuum-continue-v1/work/nrg-z1-k5000/raw.h5",
    "runs/continuum-continue-v1/work/cosh12-dmrg-even-chi256/eigenstates.json",
    "runs/continuum-continue-presentation/results.json",
    "output/continuum-continue/campaign.json", "output/continuum-continue/jobs/nrg-z1-k5000.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/continuum-lambda-v1.amendment.json", "runs/continuum-lambda-v1.controller.log",
    "runs/continuum-lambda-v1/campaign.json", "runs/continuum-lambda-v1/clock.json",
    "runs/continuum-lambda-v1/jobs/nrg-l8-control.json",
    "runs/continuum-lambda-v1/jobs/nrg-l8-control.submission.json",
    "runs/continuum-lambda-v1/work/nrg-l8-control/raw.h5",
    "runs/continuum-lambda-v1/work/nrg-l4-z1-k7000/raw.h5",
    "runs/continuum-lambda-presentation/results.json",
    "output/continuum-lambda/campaign.json", "output/continuum-lambda/jobs/nrg-l8-control.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "runs/final-endpoint-v1.amendment.json", "runs/final-endpoint-v1.controller.log",
    "runs/final-endpoint-v1/campaign.json", "runs/final-endpoint-v1/clock.json",
    "runs/final-endpoint-v1/jobs/nrg-l8-z1-k5000.json",
    "runs/final-endpoint-v1/jobs/nrg-l8-z1-k5000.submission.json",
    "runs/final-endpoint-v1/work/nrg-l8-z1-k5000/raw.h5",
    "runs/final-endpoint-v1/work/nrg-l8-z1-k5000/launch.json",
    "output/final/campaign.json", "output/final/jobs/nrg-l8-z1-k5000.json",
)) + tuple(f"NRG_comparisons/test4/manual_run/{name}" for name in (
    "README.md", "model.m", "modeloperators.m", "nrg.py", "physics.py", "manual.py",
    "param.loop", "1_run", "run_one", "2_proc", "bundle.json", "sweep.json",
    "1/param", "1/model.m", "1/modeloperators.m", "1/data", "1/raw.h5", "8/raw.h5",
    "1/started.json", "1/completed.json", "processed/per_z.csv", "processed/means.csv",
    "processed/summary.json",
)) + tuple(f"NRG_comparisons/test4/{name}" for name in (
    "analyze.py", "retention.py", "retention_analysis.py", "high_retention.py",
    "followup.py", "followup_analysis.py", "continuum.py", "continuum_analysis.py",
    "continuum_parallel.py", "continuum_cluster_worker.py", "continuum_parallel_analysis.py",
    "continuum_parallel_assessment.py", "continuum_continue.py", "continuum_continue_analysis.py",
    "continuum_lambda.py", "continuum_lambda_control.py", "continuum_lambda_analysis.py",
    "final_endpoint.py", "final_endpoint_analysis.py",
    "input/nrg-retention.json", "input/nrg-high-retention.json", "input/next-priorities.json",
    "input/followup-baths.json", "input/continuum.json", "input/continuum-baths.json",
    "input/continuum-parallel.json", "input/continuum-continue.json", "input/continuum-lambda.json",
    "input/final-endpoint.json", "output/retention/shells.csv", "output/high-retention/shells.csv",
    "output/continuum-parallel/recovery.json",
)) + tuple(f"NRG_comparisons/test4/output/{directory}/{name}"
           for directory in TEST4_HISTORY_DIRS for name in (
    "results.json", "observables.csv", "comparisons.csv", "report.md", "z_averages.csv",
)) + tuple(f"NRG_comparisons/{case}{name}"
           for case in ("", "test2/", "test3/", "test4/") for name in (
    "private/results.json", "private_archive/results.json", "private-archive/results.json",
    "private_archives/results.json", "private-archives/results.json",
    "worktree/.git", "worktrees/research/.git", "research.bundle",
    ".venv/pyvenv.cfg", "venv/pyvenv.cfg", "shared_venv/pyvenv.cfg", "shared-venv/pyvenv.cfg",
    ".git/config", ".env", "credentials.json", "secrets.json",
))


def load_tool(name):
    path = Path(__file__).resolve().parents[1]/"tools"/f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("covered,passed", [(94, False), (95, True), (100, True)])
def test_line_gate_is_independent_of_branch_percentage(covered, passed):
    report = dict(meta={"branch_coverage": True}, files={"module.py": {}},
                  totals=dict(num_statements=100, covered_lines=covered,
                              num_branches=100, covered_branches=20))
    assert load_tool("check_coverage").check_report(report) is passed


@pytest.mark.parametrize("change", ["empty", "no_branches", "unmeasured_branches", "invalid_lines"])
def test_incomplete_coverage_reports_fail(change):
    report = dict(meta={"branch_coverage": True}, files={"module.py": {}},
                  totals=dict(num_statements=100, covered_lines=100,
                              num_branches=10, covered_branches=10))
    if change == "empty":
        report["files"] = {}
    elif change == "no_branches":
        report["totals"].update(num_branches=0, covered_branches=0)
    elif change == "unmeasured_branches":
        report["meta"]["branch_coverage"] = False
    else:
        report["totals"]["covered_lines"] = 101
    with pytest.raises(ValueError):
        load_tool("check_coverage").check_report(report)


@pytest.mark.parametrize("member,message", [
    ("reports/results.xml", "unexpected wheel member"),
    *[(name, "unexpected wheel member") for name in NRG_SOURCE_ASSETS + NRG_REJECTED_ASSETS],
    ("qdjj_solver/__init__.py", "no native QP core"),
])
def test_distribution_checker_rejects_incomplete_or_leaking_wheels(tmp_path, member, message):
    wheel = tmp_path/"broken.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(member, "placeholder")
    with pytest.raises(ValueError, match=message):
        load_tool("check_distribution").inspect(wheel)


@pytest.mark.parametrize("extra,missing,valid", [
    (None, None, True),
    ("applications/case/output/paper/current.csv", None, True),
    ("APPLICATIONS.md", None, False),
    ("session-example.md", None, False),
    ("retention-1500-2000-v1.session.json", None, False),
    ("docs/retention-1500-2000-v1.session.json", None, False),
    ("src/qdjj_solver/retention-1500-2000-v1.session.json", None, False),
    ("runs/retention-v1/campaign.json", None, False),
    ("retention-3000-4000-v1.session.json", None, False),
    ("docs/retention-3000-4000-v1.session.json", None, False),
    ("src/qdjj_solver/retention-3000-4000-v1.session.json", None, False),
    ("runs/high-retention-v1/campaign.json", None, False),
    ("next-priorities-v1.session.json", None, False),
    ("docs/next-priorities-v1.session.json", None, False),
    ("src/qdjj_solver/next-priorities-v1.session.json", None, False),
    ("runs/next-priorities-v1/campaign.json", None, False),
    ("applications/case/runs/checkpoint.json", None, False),
    *[(name, None, False) for name in NRG_REJECTED_ASSETS],
    *[(None, name, False) for name in NRG_SOURCE_ASSETS],
    (None, "tools/check_publication.py", False),
    ("NRG_comparisons/test5/README.md", None, False),
    ("large_Gamma/output/observables.csv", None, True),
    ("large_Gamma/runs/checkpoint.json", None, False),
    ("docs/manual/build/manual.pdf", None, False),
    ("docs/manual/figures/primer_devices.pdf", None, False),
    (None, "docs/handbook/index.md", False),
    (None, "docs/formats.toml", False),
    (None, "examples/01_subgap_states/input/reference.json", False),
])
def test_research_source_distribution_boundary(tmp_path, extra, missing, valid):
    tool = load_tool("check_distribution")
    path = tmp_path/"source.tar.gz"
    names = tool.REQUIRED.copy()
    if extra:
        names.add(extra)
    if missing:
        names.discard(missing)
    with tarfile.open(path, "w:gz") as archive:
        for name in names:
            info = tarfile.TarInfo("qdjj_solver/"+name)
            data = b"fixture\n"
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    if valid:
        tool.inspect(path)
    else:
        with pytest.raises(ValueError):
            tool.inspect(path)


@pytest.mark.parametrize("suffix", [".tar.gz", ".whl"])
@pytest.mark.parametrize("member", [
    "/README.md", "../README.md", "./README.md", "qdjj_solver/../README.md",
    "qdjj_solver/./README.md", "qdjj_solver//README.md", "qdjj_solver/../../",
    "qdjj_solver\\README.md", "C:/README.md", "qdjj_solver/C:README.md",
])
def test_distribution_checker_rejects_unsafe_raw_paths(tmp_path, suffix, member):
    path = tmp_path/("unsafe" + suffix)
    if suffix == ".tar.gz":
        with tarfile.open(path, "w:gz") as archive:
            info = tarfile.TarInfo(member)
            if member.endswith("/"):
                info.type = tarfile.DIRTYPE
            archive.addfile(info)
    else:
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(member, "")
    with pytest.raises(ValueError, match="unsafe archive member path"):
        load_tool("check_distribution").inspect(path)


@pytest.mark.parametrize("extra,kind,valid", [
    ("qdjj_solver/src/", tarfile.DIRTYPE, True),
    ("other/src/", tarfile.DIRTYPE, False),
    ("other/README.md", tarfile.REGTYPE, False),
    ("qdjj_solver", tarfile.REGTYPE, False),
])
def test_source_distribution_requires_one_directory_root(tmp_path, extra, kind, valid):
    tool = load_tool("check_distribution")
    path = tmp_path/"roots.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        root = tarfile.TarInfo("qdjj_solver/")
        root.type = tarfile.DIRTYPE
        archive.addfile(root)
        for name in tool.REQUIRED:
            archive.addfile(tarfile.TarInfo("qdjj_solver/" + name))
        info = tarfile.TarInfo(extra)
        info.type = kind
        archive.addfile(info)
    if valid:
        tool.inspect(path)
    else:
        with pytest.raises(ValueError, match="top-level directory"):
            tool.inspect(path)


@pytest.mark.parametrize("kind", [
    tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE, tarfile.CHRTYPE,
    tarfile.BLKTYPE, tarfile.CONTTYPE, tarfile.GNUTYPE_SPARSE, b"Z",
])
def test_source_distribution_rejects_unsupported_entries(tmp_path, kind):
    path = tmp_path/"unsupported.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        info = tarfile.TarInfo("qdjj_solver/src/alias")
        info.type = kind
        if kind in {tarfile.SYMTYPE, tarfile.LNKTYPE}:
            info.linkname = "qdjj_solver/README.md"
        archive.addfile(info)
    with pytest.raises(ValueError, match="unsupported tar member type"):
        load_tool("check_distribution").inspect(path)


@pytest.mark.parametrize("directory", [
    "private", "NRG_comparisons/test1/runs", "src/__pycache__",
    "NRG_comparisons/test4/runs", "NRG_comparisons/test4/manual_run",
    "NRG_comparisons/test4/manual_run/1", "NRG_comparisons/test4/manual_run/processed",
    "runs", "NRG_comparisons/test4/runs/retention-v1", "NRG_comparisons/test4/output/retention-v1",
    "NRG_comparisons/test4/output/retention/snapshots", "NRG_comparisons/test4/output/retention/source",
    "NRG_comparisons/test4/output/retention/reused", "NRG_comparisons/test4/output/retention/controls",
    "NRG_comparisons/test4/runs/high-retention-v1", "NRG_comparisons/test4/output/high-retention-v1",
    "NRG_comparisons/test4/runs/high-retention-v1/nrg-k3000-m400",
    "NRG_comparisons/test4/runs/high-retention-v1/nrg-k4000-m400",
    "NRG_comparisons/test4/output/high-retention/snapshots", "NRG_comparisons/test4/output/high-retention/source",
    "NRG_comparisons/test4/output/high-retention/reused", "NRG_comparisons/test4/output/high-retention/controls",
    "NRG_comparisons/test4/runs/next-priorities-v1", "NRG_comparisons/test4/output/next-priorities-v1",
    "NRG_comparisons/test4/runs/next-priorities-v1/cosh4-dmrg-even-chi64",
    "NRG_comparisons/test4/runs/next-priorities-v1/nrg-length22-k4000",
    "NRG_comparisons/test4/runs/next-priorities-v1/nrg-window16-k4000",
    "NRG_comparisons/test4/output/followup/snapshots", "NRG_comparisons/test4/output/followup/source",
    "NRG_comparisons/test4/output/followup/reused", "NRG_comparisons/test4/output/followup/controls",
    "NRG_comparisons/test4/output/followup/checkpoints", "NRG_comparisons/test4/output/followup/archive",
    *[f"NRG_comparisons/test4/output/{name}" for name in TEST4_HISTORY_DIRS],
    *[f"NRG_comparisons/{case}{name}"
      for case in ("", "test2/", "test3/", "test4/") for name in NRG_PRIVATE_DIRS],
])
def test_source_distribution_rejects_unreviewed_empty_directories(tmp_path, directory):
    tool = load_tool("check_distribution")
    path = tmp_path/"directories.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name in tool.REQUIRED:
            archive.addfile(tarfile.TarInfo("qdjj_solver/" + name))
        info = tarfile.TarInfo("qdjj_solver/" + directory)
        info.type = tarfile.DIRTYPE
        archive.addfile(info)
    with pytest.raises(ValueError, match="source boundary mismatch"):
        tool.inspect(path)


@pytest.mark.parametrize("kind", [stat.S_IFLNK, stat.S_IFIFO])
def test_wheel_rejects_unsupported_entries(tmp_path, kind):
    path = tmp_path/"unsupported.whl"
    with zipfile.ZipFile(path, "w") as archive:
        info = zipfile.ZipInfo("qdjj_solver/alias")
        info.create_system = 3
        info.external_attr = (kind | 0o644) << 16
        archive.writestr(info, "qdjj_solver/__init__.py")
    with pytest.raises(ValueError, match="unsupported wheel member type"):
        load_tool("check_distribution").inspect(path)


@pytest.mark.parametrize("session", [
    None, "qdjj_solver/retention.session.json", "qdjj_solver-0.2.0.dist-info/retention.session.json",
    "qdjj_solver/retention-3000-4000-v1.session.json", "qpsolver/retention-3000-4000-v1.session.json",
    "qdjj_solver-0.2.0.dist-info/retention-3000-4000-v1.session.json",
    "qdjj_solver/next-priorities-v1.session.json", "qpsolver/next-priorities-v1.session.json",
    "qdjj_solver-0.2.0.dist-info/next-priorities-v1.session.json",
])
def test_wheel_allows_regular_files_and_directories_but_not_sessions(tmp_path, session):
    path = tmp_path/"valid.whl"
    with zipfile.ZipFile(path, "w") as archive:
        for name in ("qdjj_solver/", "qdjj_solver/qp_solver/_core.so",
                     "qdjj_solver/common/problem.py", "qdjj_solver/dmrg_solver/solver.py",
                     "qpsolver/__init__.py", "qdjj_solver-0.2.0.dist-info/METADATA"):
            info = zipfile.ZipInfo(name)
            info.create_system = 3
            info.external_attr = (stat.S_IFDIR if name.endswith("/") else stat.S_IFREG) << 16
            archive.writestr(info, "")
        if session:
            archive.writestr(session, "{}")
    if session:
        with pytest.raises(ValueError, match="local runtime session files leaked"):
            load_tool("check_distribution").inspect(path)
    else:
        load_tool("check_distribution").inspect(path)


def test_wheel_directory_cannot_satisfy_required_file(tmp_path):
    path = tmp_path/"directory.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("qdjj_solver/qp_solver/_core.so", "")
        info = zipfile.ZipInfo("qdjj_solver/common/problem.py")
        info.create_system = 3
        info.external_attr = stat.S_IFDIR << 16
        archive.writestr(info, "")
    with pytest.raises(ValueError, match="wheel is missing qdjj_solver/common/problem.py"):
        load_tool("check_distribution").inspect(path)


def test_nrg_source_assets_match_build_configuration():
    root = Path(__file__).resolve().parents[1]
    with (root / "pyproject.toml").open("rb") as handle:
        settings = tomllib.load(handle)["tool"]["scikit-build"]
    tool = load_tool("check_distribution")
    assets = set(NRG_SOURCE_ASSETS)
    assert len(NRG_SOURCE_ASSETS) == len(assets)
    assert len(NRG_REJECTED_ASSETS) == len(set(NRG_REJECTED_ASSETS))
    assert assets.isdisjoint(NRG_REJECTED_ASSETS)
    assert tool.NRG_SOURCE_FILES == assets
    assert {name for name in tool.REQUIRED if name.startswith("NRG_comparisons/")} == assets
    assert set(settings["sdist"]["include"]) == assets
    assert len(settings["sdist"]["include"]) == len(assets)
    assert {"NRG_comparisons", "runs", "APPLICATIONS.md", "session-*.md", "*.session.json"} <= set(settings["sdist"]["exclude"])
    assert settings["wheel"]["packages"] == ["src/qdjj_solver", "src/qpsolver"]
    assert "tools/check_publication.py" in tool.REQUIRED


def test_nrg_source_assets_exist():
    root = Path(__file__).resolve().parents[1]
    assert all((root / name).is_file() for name in NRG_SOURCE_ASSETS)


@pytest.mark.parametrize("workflow,job,archives", [
    ("tests.yml", "linux", ["dist/*.tar.gz dist/*.whl"]),
    ("tests.yml", "compatibility", ["dist/*.whl"]),
    ("tests.yml", "minimum-dependencies", []),
    ("tests.yml", "sanitizers", []),
    ("tests.yml", "native-coverage", []),
    ("tests.yml", "docs", []),
    ("tests.yml", "docs-extended", []),
    ("release.yml", "sdist", ["dist/*.tar.gz"]),
    ("release.yml", "wheels", ["dist/*.tar.gz", "wheelhouse/*.whl"]),
    ("release.yml", "documentation", ["dist/*.tar.gz"]),
])
def test_publication_ci_gates_precede_execution_and_artifact_upload(workflow, job, archives):
    root = Path(__file__).resolve().parents[1]
    text = (root / ".github/workflows" / workflow).read_text()
    body = re.split(r"\n  [\w-]+:\n", text.split(f"\n  {job}:\n", 1)[1], maxsplit=1)[0]
    steps = re.split(r"^      - ", body, flags=re.MULTILINE)[1:]
    source = next(step for step in steps if "id: publication-source\n" in step)
    executing = [step for step in steps if re.search(r"^(?:run:|        run:)", step, re.MULTILINE)]
    assert executing[0] == source
    assert "if:" not in source
    for step in steps:
        if job == "docs-extended" and step.startswith("name: External links (report only)\n"):
            assert "run: python -m sphinx -b linkcheck" in step
        else:
            assert "continue-on-error:" not in step
    command = "python -B tools/check_publication.py"
    assert f"run: {command}\n" in source
    assert body.index("uses: actions/setup-python@") < body.index(command)
    assert re.findall(re.escape(command) + r"[^\n]*", body) == [
        command, *(f"{command} {target}" for target in archives),
    ]
    for target in archives:
        boundary = body.index(f"python tools/check_distribution.py {target}")
        privacy = body.index(f"{command} {target}")
        assert boundary < privacy < body.index("uses: actions/upload-artifact@")
        if job in ("wheels", "documentation") and target == "dist/*.tar.gz":
            consume = "uses: pypa/cibuildwheel@" if job == "wheels" else "name: Install from the source archive"
            assert privacy < body.index("tar -xzf") < body.index(consume)
        else:
            build = "uses: pypa/cibuildwheel@" if job == "wheels" else "python -m build"
            assert body.index(build) < boundary
    for step in steps:
        if "uses: actions/upload-artifact@" in step:
            if "path: reports/" in step or "manual-logs/" in step:
                assert "if: ${{ always() && steps.publication-source.outcome == 'success' }}" in step
            else:
                assert "if:" not in step


@pytest.mark.parametrize("name", ["README.md", "NRG_comparisons/README.md",
                                  "NRG_comparisons/test1/README.md", "NRG_comparisons/test2/README.md",
                                  "NRG_comparisons/test2/manual_run/README.md",
                                  "NRG_comparisons/test3/README.md",
                                  "NRG_comparisons/test3/manual_run/README.md",
                                  "NRG_comparisons/test4/README.md",
                                  "NRG_comparisons/test4/README.model_and_parameters.md"])
def test_nrg_documentation_local_links(name):
    root = Path(__file__).resolve().parents[1]
    document = root / name
    for link in re.findall(r"\[[^\]\n]+\]\(([^\s)]+)\)", document.read_text()):
        url = urlsplit(link)
        if url.scheme or url.netloc:
            continue
        target = (document.parent / url.path).resolve() if url.path else document
        assert target.is_relative_to(root) and target.is_file(), f"{document}: {link}"
        if url.fragment:
            headings = re.findall(r"^#+\s+(.+)$", target.read_text(), re.MULTILINE)
            anchors = {re.sub(r"[^\w -]", "", heading.lower()).replace(" ", "-") for heading in headings}
            assert url.fragment in anchors, f"{document}: {link}"
