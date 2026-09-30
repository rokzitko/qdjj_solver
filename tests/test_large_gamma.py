"""Independent quadratic physics and bounded large-Gamma workflow regressions."""

from copy import deepcopy
from datetime import datetime, timedelta, timezone
import os
import subprocess
import sys
from time import monotonic

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import Impurity, cosh_grid
from qdjj_solver.common.io import fingerprint, write_json
from electron_oracle import electron_hamiltonian, electron_operators, sector_indices
from large_Gamma import analyze, campaign, deadline, deadline_report, run, study


@pytest.mark.parametrize("gamma,phi", [(.1, 0.), (10., 0.), (10., 3*np.pi/5)])
def test_gaussian_spin_charge_energy_against_independent_electron_ed(gamma, phi):
    bath = cosh_grid(1, bandwidth=100.)
    point = dict(u=0., gamma=gamma, phi=phi)
    gaussian = study.quadratic_finite(bath, point)
    t = np.eye(2)*np.sqrt(gamma/(2*np.pi*bath.rho))
    matrix = electron_hamiltonian(Impurity.anderson(0.), [bath, bath], [t, t], [-phi/2, phi/2])
    cs = electron_operators(10)
    ns = [c.getH() @ c for c in cs[:2]]
    for parity in (0, 1):
        selected = sector_indices(10, parity, parity)
        energies, vectors = np.linalg.eigh(matrix[selected][:, selected].toarray())
        vector = vectors[:, 0]

        def expectation(op, selected=selected, vector=vector):
            return float(np.vdot(vector, op[selected][:, selected] @ vector).real)

        assert_allclose(gaussian[f"energy_{'even' if parity == 0 else 'odd'}"], energies[0], atol=3e-10, rtol=0)
        assert_allclose(gaussian["P_2" if parity else "even_P_2"], expectation(ns[0] @ ns[1]), atol=3e-11, rtol=0)
        if parity:
            assert_allclose(gaussian["q_d"], expectation(ns[0]-ns[1]), atol=3e-11, rtol=0)


@pytest.mark.parametrize("gamma,phi", [(.01, 0.), (.4, 3*np.pi/5), (10., 0.), (10., 3*np.pi/5)])
def test_continuum_oracle_against_refined_bdg_and_phase_derivatives(gamma, phi):
    point = dict(u=0., gamma=gamma, phi=phi)
    exact = study.quadratic_continuum(point)["values"]
    bath = cosh_grid(48, bandwidth=100.)
    finite = study.quadratic_finite(bath, point)
    for name in study.VALUES:
        assert_allclose(finite[name], exact[name], atol=2e-10, rtol=0)
    if phi:
        step = 1e-4
        shifted = [study.quadratic_continuum(dict(point, phi=phi+sign*step))["values"] for sign in (-1, 1)]
        for branch in ("even", "odd"):
            derivative = (shifted[1][f"energy_{branch}"]-shifted[0][f"energy_{branch}"])/(2*step)
            assert_allclose(derivative, exact[f"current_{branch}"], atol=2e-9, rtol=0)


def test_runner_exact_bath_resume_and_report_does_not_invent_convergence(tmp_path):
    config = study.small_config()
    identity = run.prepare_archive(tmp_path, config)
    bath = run.prepare_bath(tmp_path, config["baths"][0], 2, config)
    point = dict(u=0., gamma=.4, phi=3*np.pi/5)
    task = run.task_record(point, bath, "qp")
    path = run.result_path(tmp_path, task)
    request = dict(config=config, task=task, identity=identity, bath=bath["bath"])
    run.execute_task(request, path)
    result = run.load_task(path, identity, task)
    exact = next(r for r in result["records"] if r["method"] == "BdG")
    ed = next(r for r in result["records"] if r["method"] == "ED")
    full = next(r for r in result["records"] if r["method"] == "6QP")
    for name in study.VALUES:
        assert_allclose(ed["values"][name], exact["values"][name], atol=2e-10, rtol=0)
        assert_allclose(full["values"][name], exact["values"][name], atol=2e-10, rtol=0)
    energies = [next(r for r in result["records"] if r["method"] == f"{q}QP")["values"]["energy_odd"]
                for q in (2, 4, 6)]
    assert np.all(np.diff(energies) < 2e-10)
    run.execute_task(request, path)
    assert run.load_task(path)["records"] == result["records"]
    summary = analyze.export(tmp_path, tmp_path/"export")
    assert not summary["complete"]
    assert summary["coverage"]["6QP"]["measured_grid_points"] == 1
    assert summary["coverage"]["6QP"]["all_observables_bath_converged_grid_points"] == 0
    assert analyze.verify_archive(tmp_path/"export")
    (tmp_path/"export/observables.csv").write_text("corrupted\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        analyze.verify_archive(tmp_path/"export")
    changed = deepcopy(config)
    changed["qp_options"]["seed"] += 1
    with pytest.raises(ValueError, match="resume requires"):
        run.prepare_archive(tmp_path, changed)
    result["records"][0]["values"]["q_d"] += .01
    write_json(path, result)
    with pytest.raises(ValueError, match="checksum"):
        run.load_task(path)


def test_reference_resolution_and_nonmonotone_threshold_brackets():
    config = study.small_config()
    values = {name: .5 for name in study.VALUES}
    rows = [dict(values={k: v+offset for k, v in values.items()}, accepted=True)
            for offset in (1e-10, 5e-11, 0.)]
    assert max(study.empirical_resolution(rows, config).values()) == 1e-8
    rows[-1]["accepted"] = False
    assert study.empirical_resolution(rows, config) is None
    comparisons = [dict(u=2., phi=0., family="cosh-grid", bandwidth=100., method="2QP",
                        observable="q_d", scope="matched_bath", levels=8, gamma=g,
                        absolute_error=e, reference_resolution=1e-8)
                   for g, e in ((.1, 1e-5), (1., 1e-2), (10., 1e-5))]
    brackets = analyze.threshold_brackets(comparisons, [1e-3])
    assert [r["status"] for r in brackets] == ["below_to_above", "above_to_below"]
    proposed = analyze.refinement_points(comparisons, config)
    assert {p["gamma"] for p in proposed} == {.55, 5.5}
    unresolved = [dict(r, scope="continuum", qp_bath_resolved=False) for r in comparisons]
    assert analyze.threshold_brackets(unresolved, [1e-3]) == []


def test_process_worker_and_resource_preflight(tmp_path):
    config = study.small_config()
    config["qp_options"]["max_dimension"] = 1
    config["resources"]["poll_seconds"] = .05
    identity = run.prepare_archive(tmp_path, config)
    bath = run.prepare_bath(tmp_path, config["baths"][0], 2, config)
    task = run.task_record(dict(u=2., gamma=.4, phi=3*np.pi/5), bath, "qp")
    run.execute_queue([(task, bath, 1.)], tmp_path, identity, config, workers=1)
    result = run.load_task(run.result_path(tmp_path, task), identity, task)
    assert result["complete"]
    assert all(row["status"] == "resource_limited" for row in result["records"])
    # Repeated submission preserves the completed file byte for byte.
    before = run.result_path(tmp_path, task).read_bytes()
    run.execute_queue([(task, bath, 1.)], tmp_path, identity, config, workers=1)
    assert run.result_path(tmp_path, task).read_bytes() == before


def test_control_derivative_and_dark_mode_comparison(tmp_path):
    config = study.small_config()
    identity = run.prepare_archive(tmp_path, config)
    bath = run.prepare_bath(tmp_path, config["baths"][0], 2, config)
    point = dict(u=2., gamma=.4, phi=3*np.pi/5)
    tasks = [run.task_record(point, bath, "qp", variant) for variant in ("standard", "uncompressed")]
    step = config["controls"]["phase_step"]
    tasks += [dict(run.task_record(dict(point, phi=point["phi"]+sign*step), bath, "qp", variant), anchor=point)
              for sign, variant in ((-1, "phase_minus"), (1, "phase_plus"))]
    for task in tasks:
        run.execute_task(dict(config=config, identity=identity, task=task, bath=bath["bath"]),
                         run.result_path(tmp_path, task))
    _, _, rows, failures, _ = analyze.collect(tmp_path)
    assert not failures
    checks = analyze.control_checks(rows, config)
    exact = [row for row in checks if row["method"] == "ED"]
    assert {row["check"] for row in exact} == {"phase_derivative", "uncompressed"}
    assert max(max(row["absolute_changes"].values()) for row in exact) < 1e-8


@pytest.mark.dmrg
def test_dmrg_bond_ladder_against_quadratic_control(tmp_path):
    pytest.importorskip("tenpy")
    config = study.small_config()
    identity = run.prepare_archive(tmp_path, config)
    bath = run.prepare_bath(tmp_path, config["baths"][0], 2, config)
    point = dict(u=0., gamma=10., phi=3*np.pi/5)
    task = run.task_record(point, bath, "dmrg")
    path = run.result_path(tmp_path, task)
    run.execute_task(dict(config=config, task=task, identity=identity, bath=bath["bath"]), path)
    result = run.load_task(path)
    rows = [r for r in result["records"] if r.get("method") == "DMRG"]
    assert len(rows) == 3 and rows[-1]["bond_converged"]
    exact = study.quadratic_finite(cosh_grid(1, bandwidth=100.), point)
    for name in study.VALUES:
        assert_allclose(rows[-1]["values"][name], exact[name], atol=2e-9, rtol=0)
    assert max(s["residual"] for s in rows[-1]["sectors"]) < 1e-7
    assert fingerprint(result["records"]) == result["records_sha256"]


def test_shipped_archive_scope_integrity_and_small_file_budget():
    directory = study.ROOT/"large_Gamma/output"
    assert analyze.verify_archive(directory)
    summary = study.read_json(directory/"summary.json")
    config = study.read_json(study.DEFAULT_INPUT)
    assert summary["expected_base_grid_points"] == 288
    rows = analyze.archive_rows(directory, "observables.csv")
    for method in ("2QP", "4QP", "6QP"):
        measured = {(float(row["u"]), float(row["gamma"]), float(row["phi"])) for row in rows
                    if row["method"] == method and row["variant"] == "standard"}
        expected = {(p["u"], p["gamma"], p["phi"]) for p in study.points(config)}
        assert expected <= measured
    sizes = [path.stat().st_size for path in directory.rglob("*") if path.is_file()]
    assert max(sizes) < 5*1024**2
    assert sum(sizes) < 20*1024**2
    assert not any(path.suffix in (".npz", ".h5", ".hdf5") for path in directory.rglob("*"))
    if summary["complete"]:
        assert all(c["all_observables_bath_converged_grid_points"] == 288 for c in summary["coverage"].values())


def test_condensed_csv_partition_and_reassembly(tmp_path):
    rows = [dict(index=i, value=.1*i) for i in range(100)]
    parts = analyze.write_csv(tmp_path/"measurements.csv", rows, max_bytes=120)
    assert len(parts) > 1
    assert all((tmp_path/name).stat().st_size <= 120 for name in parts)
    write_json(tmp_path/"manifest.json", dict(datasets={"measurements.csv": parts}))
    restored = analyze.archive_rows(tmp_path, "measurements.csv")
    assert [(int(r["index"]), float(r["value"])) for r in restored] == [(r["index"], r["value"]) for r in rows]


def test_campaign_records_publication_failure(tmp_path, monkeypatch):
    class CompletedController:
        def wait(self):
            return 0

    def fail_export(*args, **kwargs):
        raise subprocess.CalledProcessError(1, "archive-export")

    monkeypatch.setattr(sys, "argv", ["campaign", "--output", str(tmp_path), "--publish", str(tmp_path/"published")])
    monkeypatch.setattr(campaign.subprocess, "Popen", lambda *a, **k: CompletedController())
    monkeypatch.setattr(campaign.subprocess, "run", fail_export)
    monkeypatch.setattr(campaign.signal, "signal", lambda *a: None)
    campaign.main()
    status = study.read_json(tmp_path/"campaign.json")
    assert status["state"] == "campaign_failed"
    assert status["exception"] == "CalledProcessError"


def deadline_fixture(tmp_path):
    config = study.small_config()
    config["u_values"] = [0.]
    config["gamma_values"] = [.4]
    original_input = tmp_path/"base.json"
    write_json(original_input, config)
    profile = study.read_json(deadline.DEFAULT_PROFILE)
    profile.update(base_input=str(original_input), source_run=str(tmp_path/"original"),
                   output=str(tmp_path/"deadline"), publish=str(tmp_path/"published"), chi_values=[8, 16, 32])
    source = tmp_path/"original"
    identity = run.prepare_archive(source, config)
    bath = run.prepare_bath(source, config["baths"][0], 2, config)
    point = dict(u=0., gamma=.4, phi=3*np.pi/5)
    task = run.task_record(point, bath, "qp")
    run.execute_task(dict(config=config, identity=identity, task=task, bath=bath["bath"]), run.result_path(source, task))
    new_config, new_identity = deadline.prepare(tmp_path/"deadline", profile)
    return profile, new_config, new_identity, bath, point


def test_deadline_import_preserves_values_and_exact_reference_resolution(tmp_path):
    profile, config, identity, bath, point = deadline_fixture(tmp_path)
    task = run.task_record(point, bath, "qp")
    original = run.load_task(run.result_path(tmp_path/"original", task))
    copied = run.load_task(run.result_path(tmp_path/"deadline", task), identity)
    assert original["records"] == copied["records"]
    assert copied["imported_from"]["records_sha256"] == original["records_sha256"]
    assert copied["identity"] != original["identity"]
    _, _, rows, _, _ = analyze.collect(tmp_path/"deadline")
    comparisons, _ = analyze.compare(rows, config)
    assert all(r["reference_resolution"] == 1e-8 for r in comparisons if r["scope"] == "matched_bath")
    assert config["convergence"]["target"] == .001
    assert config["convergence"]["finite_target"] == 1e-5
    assert profile["report_deadline"].endswith("+02:00")


@pytest.mark.dmrg
def test_deadline_exact_reference_validates_dmrg_without_repeated_full_bond_ladder(tmp_path):
    pytest.importorskip("tenpy")
    profile, config, identity, bath, point = deadline_fixture(tmp_path)
    _, _, rows, _, _ = analyze.collect(tmp_path/"deadline")
    oracle = deadline.references(rows)[(study.point_key(point), bath["bath_key"], "standard")]
    assert len(deadline.selected_sectors(oracle, 1e-8)) == 2
    task = run.task_record(point, bath, "dmrg")
    request = dict(config=config, profile=profile, task=task, identity=identity,
                   bath=bath["bath"], oracle=oracle,
                   numerical_scripts={name: study.file_hash(study.ROOT/"large_Gamma"/name)
                                      for name in ("deadline.py", "run.py", "study.py")})
    destination = run.result_path(tmp_path/"deadline", task)
    deadline.worker(request, destination)
    result = run.load_task(destination, identity)
    measurements = [r for r in result["records"] if r.get("method") == "DMRG"]
    assert len(measurements) == 1
    assert measurements[0]["observable_converged"]
    assert measurements[0]["sector_selection"] == "exact reference minima"
    assert max(measurements[0]["observable_errors"].values()) <= profile["finite_observable_target"]


def test_deadline_observable_acceptance_does_not_claim_wavefunction_convergence():
    profile = study.read_json(deadline.DEFAULT_PROFILE)
    values = {name: .5 for name in study.VALUES}
    row = dict(values=values, accepted=False, sectors=[dict(residual=.2, finite_problem_converged=False)])
    oracle = dict(values={k: v+2e-6 for k, v in values.items()}, method="ED")
    assessed = deadline.assess(row, oracle, profile)
    assert assessed["observable_converged"]
    assert not assessed["accepted"]
    assert assessed["sectors"][0]["residual"] == .2
    assert assessed["accuracy_evidence"] == "exact same-bath observable comparison"
    assert max(assessed["method_resolution"].values()) < 1e-5
    assessed = deadline.assess(row, dict(oracle, values={k: v+.01 for k, v in values.items()}), profile)
    assert not assessed["observable_converged"]
    with pytest.raises(ValueError, match="timezone"):
        deadline.timestamp("2026-09-24T07:00:00")


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group termination check")
def test_hard_deadline_kills_even_a_term_ignoring_child(tmp_path, monkeypatch):
    profile, config, identity, bath, point = deadline_fixture(tmp_path)
    directory = tmp_path/"deadline"
    profile.update(compute_stop=(datetime.now(timezone.utc)+timedelta(seconds=.8)).isoformat(),
                   termination_grace_seconds=.1, poll_seconds=.02)
    command = [sys.executable, "-B", "-c", "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)"]
    monkeypatch.setattr(deadline, "worker_command", lambda *args: command)
    task = run.task_record(point, bath, "dmrg")
    started = monotonic()
    result = deadline.schedule([dict(task=task, bath=bath, memory_gib=1., timeout=60)],
                               directory, config, identity, profile)
    assert monotonic()-started < 5
    assert result["stopped_tasks"] == 1
    saved = run.load_task(run.result_path(directory, task), identity)
    assert saved["interruption"] == "hard computation deadline"
    assert saved["complete"]


def test_deadline_report_finality_keeps_unresolved_physics_explicit(tmp_path, monkeypatch):
    deadline_fixture(tmp_path)
    output = tmp_path/"published"
    summary = analyze.export(tmp_path/"deadline", output)
    readme = tmp_path/"large_Gamma/README.md"
    readme.parent.mkdir()
    readme.write_text("# Preserved report\n"+deadline_report.BEGIN+"\nold\n"+deadline_report.END+"\nUser text\n", encoding="utf-8")
    monkeypatch.setattr(study, "ROOT", tmp_path)
    deadline_report.write_report(tmp_path/"deadline", output, summary, final=True)
    saved = study.read_json(output/"summary.json")
    assert saved["deadline_final"]
    assert not saved["complete"]
    assert "remaining bath-resolution limits" in saved["status"]
    assert "$10^{-5}$" in (output/"report.md").read_text(encoding="utf-8")
    assert "User text" in readme.read_text(encoding="utf-8")
    analyze.verify_archive(output)
