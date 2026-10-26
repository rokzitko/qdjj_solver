"""Review source/wheel members against the standalone software packaging boundary."""

import argparse
from pathlib import Path, PurePosixPath
import stat
import tarfile
import tomllib
import zipfile


SOURCE_ROOTS = {"src", "cpp", "docs", "examples", "applications",
                "large_Gamma", "NRG_comparisons", "tests", "tools", ".github"}
SOURCE_FILES = {"CMakeLists.txt", "pyproject.toml", "README.md", "LICENSE", "CITATION.cff",
                "CONTRIBUTING.md", "AGENTS.md", ".gitattributes", ".gitignore", "PKG-INFO"}
NRG_SOURCE_FILES = {"NRG_comparisons/README.md"} | {
    f"NRG_comparisons/test1/{asset}" for asset in (
        "README.md", "run.py", "nrg.py", "control.py", "analysis.py", "full_analysis.py",
        "targeted_nrg.py", "nrg_dataset.py", "twist_analysis.py", "focused_internal.py", "generate.py",
        "input/physical.json", "input/profiles.json", "input/nrg-targeted.json",
        "input/qp-dmrg-focused.json", "input/final-baths.json",
        "output/results.json", "output/observables.csv", "output/comparisons.csv",
         "output/comparison.svg", "output/convergence.svg")}
NRG_SOURCE_FILES.update(f"NRG_comparisons/test2/{asset}" for asset in (
    "README.md", "generate.py", "input/physical.json", "input/baths.json", "reference/provenance.json",
    "output/results.json", "output/observables.csv", "output/comparisons.csv",
    "output/comparison.svg", "output/convergence.svg"))
NRG_SOURCE_FILES.update(f"NRG_comparisons/test2/manual_run/{asset}" for asset in (
    "README.md", ".gitignore", "param.loop", "1_run", "run_one", "2_proc",
    "process.py", "extract.py", "reference.json"))
NRG_SOURCE_FILES.update(f"NRG_comparisons/test3/{asset}" for asset in (
    "README.md", "run.py", "finite.py", "analyze.py", "generate.py",
    "input/physical.json", "input/baths.json", "reference/provenance.json",
    "output/report.md", "output/results.json", "output/observables.csv", "output/comparisons.csv"))
NRG_SOURCE_FILES.update(f"NRG_comparisons/test3/manual_run/{asset}" for asset in (
    "README.md", ".gitignore", "nrg.py", "physics.py", "process.py", "param.loop",
    "1_run", "2_proc", "run_one"))
NRG_SOURCE_FILES.update(f"NRG_comparisons/test4/{asset}" for asset in (
    "model.m", "modeloperators.m", "nrg.py", "physics.py", "finite.py", "run.py", "manual.py", "generate.py",
    "retention_diagnostics.py", "dmrg_states.py", "continuum_finite.py",
    "README.md", "FINAL_REPORT.md", "README.model_and_parameters.md", ".gitignore",
    "input/physical.json", "input/baths.json",
    "output/results.json", "output/report.md", "output/observables.csv", "output/comparisons.csv",
    "output/z_averages.csv"))
REQUIRED = {".gitattributes", "CMakeLists.txt", "cpp/core.cpp", "pyproject.toml", "LICENSE", "README.md",
            "src/qdjj_solver/common/problem.py", "src/qdjj_solver/qp_solver/solver.py",
            "src/qdjj_solver/dmrg_solver/solver.py", "src/qpsolver/__init__.py",
            "docs/numerics.md", "tests/test_physics.py", "applications/README.md",
            "docs/bath_representations.md", "tools/benchmark_baths.py", "tools/check_publication.py",
            "tools/plot_bath_benchmarks.py", "tools/refine_bath_chain.py",
            "docs/benchmarks/baths/input.json", "docs/benchmarks/baths/refinement-input.json",
            "docs/benchmarks/baths/chain-input.json",
            "docs/benchmarks/baths/manifest.json", "docs/benchmarks/baths/baths.json",
            "docs/benchmarks/baths/measurements.json", "docs/benchmarks/baths/accuracy.svg"}
REQUIRED.update(NRG_SOURCE_FILES)
REQUIRED.update(f"examples/00_first_calculation/{asset}" for asset in (
    "input/qp.json", "input/dmrg.json", "input/exact.json",
    "reference/qp.json", "reference/dmrg.json", "reference/exact.json", "reference/comparison.json"))
# Every catalogued tutorial ships its README, runner and references.
REQUIRED.update({"examples/catalog.toml", "examples/README.md", "tools/check_tutorials.py"})
REQUIRED.update(f"examples/{entry['dir']}/{asset}"
                for entry in tomllib.loads((Path(__file__).resolve().parents[1]/"examples"/"catalog.toml")
                                           .read_text())["tutorial"]
                for asset in ("README.md", "run.py", "reference/quick.json", "reference/full.json"))
REQUIRED.update({"docs/conf.py", "docs/index.md", "docs/_ext/github_markdown.py", "docs/conventions.toml",
                 "docs/concordance.toml", "docs/troubleshooting.toml", "tools/gen_docs.py", "tools/check_docs.py"})
REQUIRED.update(f"docs/benchmarks/baths/{asset}" for asset in (
    "README.md", "finite-controls-input.json", "edge-refinement-input.json",
    "kernel.csv", "quadratic.csv", "many_body.csv", "layouts.csv", "failures.json",
    "summary.json", "tables.md", "refinement.csv", "time_to_accuracy.csv",
    "interacting_comparison.csv", "coordinate_comparison.csv", "layout_comparison.csv",
    "fit_sensitivity.csv", "construction_and_window.svg", "interacting_cost.svg",
    "cost_breakdown.svg", "odd_even.svg"))
for section in ("refinement", "chain", "controls", "edge"):
    REQUIRED.update(f"docs/benchmarks/baths/{section}/{asset}" for asset in (
        "input.json", "manifest.json", "baths.json", "measurements.json", "many_body.csv"))
for case in ("zitko_2023_knight_shift", "zonda_2023_double_dot", "zalom_2024_multiterminal",
             "paaske_2023_surrogate_spectrum", "bargerbos_2022_parity_diagram", "choi_2004_kondo_josephson",
             "zonda_2016_unequal_gaps", "kadlecova_2017_asymmetry", "hecht_2008_gap_edge",
             "bobok_2025_chain_expansion"):
    REQUIRED.update(f"applications/{case}/{asset}" for asset in (
        "README.md", "run.py", "input/parameters.json", "reference/provenance.json",
        "output/quick/manifest.json", "output/paper/comparison.svg"))
REQUIRED.update({"applications/zonda_2016_unequal_gaps/reference/figure8.csv",
                 "applications/kadlecova_2017_asymmetry/reference/figure1.csv",
                 "applications/hecht_2008_gap_edge/reference/spectra.csv",
                 "applications/hecht_2008_gap_edge/spectral.py",
                 "tests/test_applications_third_set.py"})
REQUIRED.update({"docs/chain_expansion.md", "tests/test_chain_expansion.py",
                 "tests/test_application_chain_expansion.py"})
REQUIRED.update(f"applications/bobok_2025_chain_expansion/{asset}" for asset in (
    "models.py", "chain_check.py", "input/dmrg_check.json", "reference/figure8.csv",
    "reference/figure17.csv", "output/paper/shared_lead.svg", "output/validation.json"))
REQUIRED.update({"docs/qp_convergence.md", "examples/qp_convergence.py",
                 "tools/benchmark_qp_convergence.py", "tools/plot_qp_convergence.py",
                 "tests/test_qp_convergence.py"})
REQUIRED.update(f"docs/benchmarks/qp_convergence/{asset}" for asset in (
    "README.md", "input.json", "manifest.json", "baths.json", "measurements.json",
    "summary.json", "convergence.csv", "table.md", "observables.svg", "errors.svg", "weights.svg"))
REQUIRED.update(f"docs/benchmarks/qp_convergence/current/{asset}" for asset in (
    "input.json", "manifest.json", "baths.json", "measurements.json"))
REQUIRED.update(f"large_Gamma/{asset}" for asset in (
    "README.md", "run.py", "study.py", "analyze.py", "campaign.py", "input/study.json", "input/pilot.json",
    "deadline.py", "deadline_report.py", "input/deadline.json", "output/report.md",
    "output/summary.json", "output/manifest.json", "output/baths.json", "output/observables.csv",
    "output/convergence.csv", "output/comparisons.csv", "output/breakdown.csv"))
REQUIRED.add("tests/test_large_gamma.py")
REQUIRED.update({"docs/floquet.md", "examples/floquet.py", "tests/test_floquet.py",
                 "src/qdjj_solver/common/floquet.py", "src/qdjj_solver/qp_solver/floquet.py",
                 "tests/test_floquet_transitions.py", "src/qdjj_solver/common/probes.py",
                 "src/qdjj_solver/common/spectra.py", "src/qdjj_solver/qp_solver/transitions.py"})
REQUIRED.update(f"docs/manual/{asset}" for asset in (
    "README.md", "Makefile", "latexmkrc", "manual.tex", "preamble.tex", "references.bib",
    "chapters/04_floquet_theory.tex", "chapters/07_probes.tex", "scripts/manual_common.py", "scripts/fig_phase.py",
    "scripts/quoted.py", "data/phase.json", "generated/quoted.tex", "generated/concordance.tex",
    "chapters/01a_background.tex", "chapters/02a_correlations.tex", "chapters/02b_experiments.tex",
    "scripts/fig_primer.py", "data/primer.json"))
REQUIRED.update(f"docs/{asset}" for asset in (
    "handbook/index.md", "_ext/handbook.py", "parameters.toml", "formats.toml", "user_manual.md",
    "getting_started.md", "static_calculations.md", "reservoir_workflow.md", "qp_workflow.md",
    "dmrg.md", "workflows.md", "custom_models.md", "spectroscopy.md",
    "reference/models.md", "reference/configuration.md", "reference/records.md",
    "manual/chapters/03a_practical_calculations.tex", "manual/generated/practical_code.tex",
    "manual/scripts/edition.py"))
REQUIRED.update({"tools/build_manuals.py", "tools/check_manual.py", "tools/run_doc_examples.py",
                 "tools/docs-constraints.txt", "examples/01_subgap_states/input/reference.json",
                 "examples/11_dmrg_convergence/expected.png"})
REQUIRED.update(f"examples/handbook/{asset}" for asset in (
    "static.py", "scan.py", "custom.py", "custom.json", "reference.json"))


def member_parts(name):
    # Validate before pathlib can normalize away empty or dot components.
    parts = name.removesuffix("/").split("/")
    if "\\" in name or any(part in {"", ".", ".."} or ":" in part for part in parts):
        raise ValueError(f"unsafe archive member path: {name!r}")
    return parts


def inspect(path):
    path = Path(path)
    if path.name.endswith(".tar.gz"):
        names, roots, directories = set(), set(), set()
        with tarfile.open(path) as archive:
            for member in archive:
                parts = member_parts(member.name)
                roots.add(parts[0])
                if member.isdir():
                    if len(parts) > 1:
                        directories.add("/".join(parts[1:]))
                    continue
                if member.type not in {tarfile.REGTYPE, tarfile.AREGTYPE}:
                    raise ValueError(f"unsupported tar member type: {member.name!r} ({member.type!r})")
                if len(parts) < 2:
                    raise ValueError("source members must be inside a top-level directory")
                names.add("/".join(parts[1:]))
        if len(roots) != 1:
            raise ValueError("source archive must have a single top-level directory")
        if any(not any(name.startswith(directory + "/") for name in names) for directory in directories):
            raise ValueError("source boundary mismatch: directory has no reviewed files")
        unexpected = {name for name in names if name not in SOURCE_FILES and
                      PurePosixPath(name).parts[0] not in SOURCE_ROOTS}
        unexpected.update(name for name in names if name.startswith("NRG_comparisons/")
                           and name not in NRG_SOURCE_FILES)
        unexpected.update(name for name in names if name.startswith(("docs/manual/build/", "docs/manual/figures/"))
                           or name.startswith("docs/") and name.endswith((".aux", ".log", ".bcf", ".fdb_latexmk", ".fls", ".dvi")))
        if unexpected or REQUIRED-names:
            raise ValueError(f"source boundary mismatch: unexpected={unexpected}, missing={REQUIRED-names}")
        names.update(directories)
    elif path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = set()
            for member in archive.infolist():
                top = member_parts(member.filename)[0]
                if stat.S_IFMT(member.external_attr >> 16) not in {0, stat.S_IFREG, stat.S_IFDIR}:
                    raise ValueError(f"unsupported wheel member type: {member.filename!r}")
                if top not in {"qdjj_solver", "qpsolver"} and not top.endswith(".dist-info"):
                    raise ValueError(f"unexpected wheel member: {member.filename}")
                if not member.is_dir() and stat.S_IFMT(member.external_attr >> 16) != stat.S_IFDIR:
                    names.add(member.filename)
        if not any(name.startswith("qdjj_solver/qp_solver/_core") and name.endswith((".so", ".pyd")) for name in names):
            raise ValueError("wheel has no native QP core")
        for name in ("qdjj_solver/common/problem.py", "qdjj_solver/dmrg_solver/solver.py", "qpsolver/__init__.py"):
            if name not in names:
                raise ValueError(f"wheel is missing {name}")
    else:
        raise ValueError("expected an sdist .tar.gz or wheel .whl")
    if any(name.endswith(".session.json") for name in names):
        raise ValueError("local runtime session files leaked into distribution")
    if any("__pycache__" in PurePosixPath(name).parts for name in names):
        raise ValueError("unexpected generated or non-relative member")
    if any("APPLICATIONS.md" in PurePosixPath(name).parts or
           (name.startswith(("applications/", "large_Gamma/")) and "runs" in PurePosixPath(name).parts)
           for name in names):
        raise ValueError("local application research files leaked into distribution")
    print(f"{path.name}: checked {len(names)} software members")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", nargs="+", type=Path)
    args = parser.parse_args()
    for archive in args.archives:
        inspect(archive)


if __name__ == "__main__":
    main()
