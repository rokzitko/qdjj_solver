"""Public dispatch, installed command lines, and isolated optional-dependency failures."""

import importlib
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import textwrap

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import Hamiltonian, Sector, __version__, number, solve
from qdjj_solver import cli
from qdjj_solver.common.io import write_json
from qdjj_solver.qp_solver import SolverOptions
from qdjj_solver.qp_solver import cli as qp_cli


def _command(tmp_path, *args):
    return subprocess.run(args, cwd=tmp_path, capture_output=True, text=True, encoding="utf-8",
                          env=dict(os.environ, PYTHONIOENCODING="utf-8"), timeout=90)


def _module(tmp_path, name, *args):
    return _command(tmp_path, sys.executable, "-B", "-m", name, *map(str, args))


def _reference_config():
    return dict(model=dict(geometry="single", u=1.6, gamma=.2), bath=dict(pairs=1, bandwidth=3.))


@pytest.mark.parametrize("options", [None, {"method": "dense"}, SolverOptions(method="dense")])
def test_api_qp_options_give_the_same_physical_result(options):
    h = Hamiltonian(.7*number(0)+1.1*number(1), 2, (1, -1))
    result = solve(h, options=options)
    assert_allclose(result.energies, [.7], rtol=0, atol=1e-15)
    assert result.hamiltonian is h


def test_api_rejects_backend_and_invalid_option_dictionary():
    h = Hamiltonian(number(0), 2, (1, -1))
    with pytest.raises(ValueError, match="backend must be"):
        solve(h, backend="missing")
    with pytest.raises(TypeError, match="unexpected keyword"):
        solve(h, options={"chi_max": 8})


@pytest.mark.dmrg
def test_api_rejects_cross_backend_options():
    pytest.importorskip("tenpy")
    from qdjj_solver.dmrg_solver import SolverOptions as DMRGOptions

    h = Hamiltonian(number(0), 2, (1, -1))
    with pytest.raises(ValueError, match="options do not belong"):
        solve(h, options=SolverOptions(), backend="dmrg")
    with pytest.raises(ValueError, match="options do not belong"):
        solve(h, options=DMRGOptions(), backend="qp")


def test_module_version_and_parser_error_outside_checkout(tmp_path):
    completed = _module(tmp_path, "qdjj_solver.cli", "--version")
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == __version__
    missing = _module(tmp_path, "qdjj_solver.cli")
    assert missing.returncode == 2 and "required" in missing.stderr


@pytest.mark.parametrize("program, expected", [("qdjj-solver", "qdjj-eigenstates"),
                                               ("qpsolver", "bcs-qp-eigenstates")])
def test_installed_cli_solve_saves_states(tmp_path, program, expected):
    # The wheel's interpreter may be invoked without activating its virtualenv.
    executable = shutil.which(program, path=sysconfig.get_path("scripts"))
    assert executable is not None
    workdir = tmp_path/"working directory Δ"
    workdir.mkdir()
    config, output = workdir/"model input.json", workdir/"state result.json"
    write_json(config, dict(_reference_config(), solver={"eigenpairs": 2}))
    flag = "--save-states" if program == "qdjj-solver" else "--save-vectors"
    completed = _command(workdir, executable, "solve", str(config), "--output", str(output), flag)
    assert completed.returncode == 0, completed.stderr
    assert "State 0:" in completed.stdout and "State 1:" in completed.stdout
    assert "Saved" in completed.stdout and "energy reference:" in completed.stdout.lower()
    record = json.loads(output.read_text())
    assert record["format"] == expected and len(record["energies"]) == 2
    assert (output.parent/record["state_amplitudes"]["file"]).is_file()
    if program == "qpsolver":
        assert "QP probabilities:" in completed.stdout
        assert "impurity_charge=" in completed.stdout


@pytest.mark.parametrize("cutoff", [None, 2])
def test_legacy_dimension_cli_matches_independent_bit_count(tmp_path, capsys, cutoff):
    config = tmp_path/"model.json"
    write_json(config, dict(_reference_config(), cutoff=cutoff))
    assert qp_cli.main(["dimension", str(config)]) == 0
    output = capsys.readouterr().out
    q = 4 if cutoff is None else cutoff
    count = 0
    for bits in range(1 << 6):
        occupied = [i for i in range(6) if bits & (1 << i)]
        if (len(occupied) % 2 == 1 and sum(1-2*(i % 2) for i in occupied) == 1
                and sum(i >= 2 for i in occupied) <= q):
            count += 1
    assert "Impurity spin orbitals: 2; bath QP modes: 4" in output
    assert f"QP cutoff: {q}" in output
    assert f"Hilbert-space dimension: {count}" in output


def test_fit_command_output_matches_independent_one_pole_least_squares(tmp_path, capsys):
    output = tmp_path/"surrogate.json"
    assert qp_cli.main(["fit", "--levels", "1", "--delta", "1", "--bandwidth", "3",
                        "--frequency-cutoff", "5", "--starts", "1", "--output", str(output)]) == 0
    record = json.loads(output.read_text())
    assert record["xi"] == [0.]
    frequency = np.geomspace(1e-3, 5., 1000)
    pole = 1/(frequency**2+1)
    target = 2/np.pi*np.arctan(3/np.hypot(frequency, 1))/np.hypot(frequency, 1)
    residue = np.dot(pole, target)/np.dot(pole, pole)
    assert_allclose(record["weights"], [np.pi/6*residue], rtol=2e-9, atol=0)
    text = capsys.readouterr().out
    assert "Saved 1-level surrogate" in text and "Maximum relative g-function error" in text
    assert qp_cli.main(["fit", "--levels", "0", "--output", str(output)]) == 2
    assert "levels must be a positive integer" in capsys.readouterr().err


@pytest.mark.parametrize("main, case, message", [
    # Match the reported filename rather than platform/localized OS error text.
    (cli.main, "missing", "model.json"),
    (cli.main, "json", "Expecting"),
    (qp_cli.main, "unknown", "unknown configuration keys"),
    (qp_cli.main, "memory", "permitted maximum"),
])
def test_cli_reports_configuration_and_calculation_errors(tmp_path, capsys, main, case, message):
    config, output = tmp_path/"model.json", tmp_path/"result.json"
    data = _reference_config()
    if case == "unknown":
        data["typo"] = True
    elif case == "memory":
        data["solver"] = {"max_dimension": 1}
    if case == "json":
        config.write_text("{broken")
    elif case != "missing":
        write_json(config, data)
    assert main(["solve", str(config), "--output", str(output)]) == 2
    captured = capsys.readouterr()
    assert "Calculation failed:" in captured.err and message in captured.err
    assert "Saved" not in captured.out and not output.exists()


@pytest.mark.parametrize("main", [cli.main, qp_cli.main])
@pytest.mark.parametrize("data,message", [
    ([], "configuration must be an object"),
    (None, "configuration must be an object"),
    ("model", "configuration must be an object"),
    ({"model": []}, "model must be an object"),
    ({"bath": None}, "bath must be an object"),
    ({"sector": []}, "sector must be an object"),
    ({"solver": []}, "solver must be an object"),
    ({"model": {"compress": "false"}}, "model.compress must be a boolean"),
    ({"solver": {"mixer": "false"}}, "solver.mixer must be a boolean"),
    ({"sector": {"parity": True}}, "sector.parity must be an integer"),
    ({"cutoff": True}, "cutoff must be a nonnegative integer"),
    ({"format_version": True}, "unsupported configuration version"),
    ({"sector": {"twice_szz": 1}}, "unknown sector keys"),
    ({"solver": {"tollerance": 1e-9}}, "unknown solver keys"),
])
def test_wrong_shaped_and_mistyped_configs_fail_before_solving(tmp_path, capsys, monkeypatch, main, data, message):
    def forbidden(*args, **kwargs):
        pytest.fail("invalid configuration reached the numerical solver")
    monkeypatch.setattr(cli, "solve", forbidden)
    monkeypatch.setattr(qp_cli, "solve", forbidden)
    config, output = tmp_path/"bad.json", tmp_path/"result.json"
    write_json(config, data)
    assert main(["solve", str(config), "--output", str(output)]) == 2
    captured = capsys.readouterr()
    assert message in captured.err and "Traceback" not in captured.err
    assert not output.exists()


@pytest.mark.parametrize("config, extra, message", [
    ({"backend": "other"}, [], "backend must be qp or dmrg"),
    ({"backend": "dmrg", "cutoff": 1}, [], "DMRG requires cutoff=null"),
    ({"backend": "qp"}, ["--resume", "unused.json"], "resume currently accepts DMRG"),
])
def test_unified_backend_and_resume_errors(tmp_path, capsys, config, extra, message):
    path = tmp_path/"model.json"
    write_json(path, dict(_reference_config(), **config))
    assert cli.main(["solve", str(path), "--output", str(tmp_path/"result.json"), *extra]) == 2
    assert message in capsys.readouterr().err


def test_cli_backend_override_for_user_hamiltonian(tmp_path, capsys):
    h = Hamiltonian(.7*number(0)+1.1*number(1), 2, (1, -1))
    config, output = tmp_path/"model.json", tmp_path/"result.json"
    write_json(config, dict(backend="dmrg", model=dict(kind="hamiltonian", **h.record()), cutoff=None))
    assert cli.main(["solve", str(config), "--backend", "qp", "--output", str(output)]) == 0
    assert "Backend: qp; energy reference: as specified" in capsys.readouterr().out
    assert json.loads(output.read_text())["hamiltonian_sha256"] == h.fingerprint()


@pytest.mark.dmrg
@pytest.mark.parametrize("probability_tolerance", [None, 1e-22])
def test_unified_real_dmrg_cli_checkpoint_and_resume(tmp_path, capsys, probability_tolerance):
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    from qdjj_solver.dmrg_solver import SolverOptions as DMRGOptions
    from qdjj_solver.dmrg_solver.io import load_result

    config, first, second = tmp_path/"model.json", tmp_path/"first.json", tmp_path/"second.json"
    data = dict(_reference_config(), backend="dmrg",
                solver=dict(seed_trials=1, chi_max=16))
    if probability_tolerance is not None:
        data["solver"]["lanczos_probability_tolerance"] = probability_tolerance
    write_json(config, data)
    assert cli.main(["solve", str(config), "--output", str(first), "--save-states"]) == 0
    text = capsys.readouterr().out
    assert "Backend: dmrg" in text and "Finite-problem checks passed:" in text
    assert "residual=" in text and "not measured" not in text
    assert cli.main(["solve", str(config), "--output", str(second), "--save-states", "--resume", str(first)]) == 0
    a, b = load_result(first), load_result(second)
    assert_allclose(a.energies, b.energies, rtol=0, atol=1e-12)
    assert_allclose(abs(a.overlaps(b)), [[1.]], rtol=0, atol=1e-12)
    assert b.metadata["roots"][0]["seed_trials"][0]["seed"] == "continuation"
    expected = 1e-14 if probability_tolerance is None else probability_tolerance
    for path, result in ((first, a), (second, b)):
        saved_options = json.loads(path.read_text())["calculation"]["options"]
        assert saved_options["lanczos_probability_tolerance"] == expected
        assert result.metadata["options"]["lanczos_probability_tolerance"] == expected
        assert DMRGOptions(**saved_options).lanczos_probability_tolerance == expected


def _blocked_import_script(module, body):
    return textwrap.dedent(f"""
        import importlib.abc
        import sys
        class Block(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == {module!r} or fullname.startswith({module!r} + '.'):
                    raise ModuleNotFoundError('blocked ' + fullname, name={module!r})
        sys.meta_path.insert(0, Block())
    """)+textwrap.dedent(body)


def test_missing_tenpy_is_lazy_and_actionable_in_isolated_process(tmp_path):
    config, output = tmp_path/"model.json", tmp_path/"result.json"
    write_json(config, dict(_reference_config(), backend="dmrg"))
    script = _blocked_import_script("tenpy", f"""
        import qdjj_solver, qpsolver
        import qdjj_solver.dmrg_solver as backend
        from qpsolver.reference import chain_coefficients
        assert 'tenpy' not in sys.modules
        assert chain_coefficients(qdjj_solver.cosh_grid(1))[0].size == 2
        try:
            backend.solve
        except ImportError as error:
            assert 'Install qdjj-solver[dmrg]' in str(error)
            assert isinstance(error.__cause__, ModuleNotFoundError)
        else:
            raise AssertionError('the unavailable backend was silently accepted')
        from qdjj_solver.cli import main
        assert main(['solve', {str(config)!r}, '--output', {str(output)!r}]) == 2
    """)
    completed = _command(tmp_path, sys.executable, "-B", "-c", script)
    assert completed.returncode == 0, completed.stderr
    assert "Calculation failed:" in completed.stderr and "qdjj-solver[dmrg]" in completed.stderr
    assert not output.exists()


@pytest.mark.dmrg
def test_missing_h5py_allows_scalar_results_and_reports_checkpoint_error(tmp_path):
    pytest.importorskip("tenpy")
    config, output = tmp_path/"model.json", tmp_path/"result.json"
    write_json(config, dict(_reference_config(), backend="dmrg", solver={"seed_trials": 1}))
    script = _blocked_import_script("h5py", f"""
        from pathlib import Path
        from qdjj_solver.cli import main
        args = ['solve', {str(config)!r}, '--output', {str(output)!r}]
        assert main(args) == 0
        before = Path({str(output)!r}).read_bytes()
        assert main(args + ['--save-states']) == 2
        assert Path({str(output)!r}).read_bytes() == before
        assert not list(Path({str(tmp_path)!r}).glob('*.h5'))
    """)
    completed = _command(tmp_path, sys.executable, "-B", "-c", script)
    assert completed.returncode == 0, completed.stderr
    assert "Calculation failed:" in completed.stderr and "h5py" in completed.stderr


def test_compatibility_shims_export_identical_modules_and_functions():
    import qpsolver
    from qpsolver import cli as legacy_cli, nrg as legacy_nrg, reference
    from qdjj_solver.qp_solver import nrg
    from qdjj_solver.dmrg_solver import legacy

    assert qpsolver.Sector is Sector and qpsolver.Hamiltonian is Hamiltonian
    assert legacy_cli.main is qp_cli.main and legacy_nrg.compare is nrg.compare
    for name in ("algebra", "baths", "models"):
        assert importlib.import_module(f"qpsolver.{name}") is importlib.import_module(f"qdjj_solver.common.{name}")
    for name in ("basis", "solver", "io", "reference", "_core"):
        assert importlib.import_module(f"qpsolver.{name}") is importlib.import_module(f"qdjj_solver.qp_solver.{name}")
    assert reference.chain_coefficients is legacy.chain_coefficients
    assert reference.dmrg_reference is legacy.dmrg_reference
    assert reference.dmrg_chain_reference is legacy.dmrg_chain_reference
    assert qpsolver.SolverOptions is SolverOptions
