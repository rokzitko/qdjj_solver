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
)
NRG_REJECTED_ASSETS = tuple(f"NRG_comparisons/test1/{name}" for name in (
    "residual_study.py", "input/unreviewed.json", "FULLTEST1.md", "output/pilot-v2.json",
    "nested/runs/results.json", "checkpoints/state.h5", "output/state.hdf5",
    "output/state.npz", "output/unreviewed.json",
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
    ("session-ses_f317.md", None, False),
    ("applications/case/runs/checkpoint.json", None, False),
    *[(name, None, False) for name in NRG_REJECTED_ASSETS],
    *[(None, name, False) for name in NRG_SOURCE_ASSETS],
    ("NRG_comparisons/test2/README.md", None, False),
    ("large_Gamma/output/observables.csv", None, True),
    ("large_Gamma/runs/checkpoint.json", None, False),
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


@pytest.mark.parametrize("directory", ["private", "NRG_comparisons/test1/runs", "src/__pycache__"])
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


def test_wheel_allows_regular_files_and_directories(tmp_path):
    path = tmp_path/"valid.whl"
    with zipfile.ZipFile(path, "w") as archive:
        for name in ("qdjj_solver/", "qdjj_solver/qp_solver/_core.so",
                     "qdjj_solver/common/problem.py", "qdjj_solver/dmrg_solver/solver.py",
                     "qpsolver/__init__.py", "qdjj_solver-0.2.0.dist-info/METADATA"):
            info = zipfile.ZipInfo(name)
            info.create_system = 3
            info.external_attr = (stat.S_IFDIR if name.endswith("/") else stat.S_IFREG) << 16
            archive.writestr(info, "")
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
    assert tool.NRG_SOURCE_FILES == assets
    assert {name for name in tool.REQUIRED if name.startswith("NRG_comparisons/")} == assets
    assert set(settings["sdist"]["include"]) == assets
    assert len(settings["sdist"]["include"]) == len(assets)
    assert all((root / name).is_file() for name in assets)
    assert {"NRG_comparisons", "APPLICATIONS.md", "session-*.md"} <= set(settings["sdist"]["exclude"])
    assert settings["wheel"]["packages"] == ["src/qdjj_solver", "src/qpsolver"]


@pytest.mark.parametrize("name", ["README.md", "NRG_comparisons/README.md",
                                  "NRG_comparisons/test1/README.md"])
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
