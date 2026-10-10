"""Offline synthetic privacy checks; never scan the current research datasets."""

import base64
import importlib.util
import io
import json
from pathlib import Path
import stat
import subprocess
import tarfile
import zipfile

import pytest


SPEC = importlib.util.spec_from_file_location(
    "check_publication", Path(__file__).resolve().parents[1] / "tools/check_publication.py"
)
publication = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publication)


def github_token():
    # Construct inert patterns at runtime, never publish a token-shaped literal.
    return "gh" + "p_" + "Ab7cD9" * 6


def categories(diagnostics):
    return {item.category for item in diagnostics}


def json_bytes(value):
    return json.dumps(value).encode()


def write(path, content=b"{}"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def make_archive(path, members):
    if path.suffix == ".whl":
        with zipfile.ZipFile(path, "w") as archive:
            for name, data in members:
                archive.writestr(name, data)
    else:
        with tarfile.open(path, "w:gz") as archive:
            for name, data in members:
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
    return path


def test_nested_metadata_locations_and_values_redacted():
    value = {"records": [{"parent_campaign": "private-campaign-label",
                          "nested": {"SLURM_JOB_ID": 987654321}}]}
    issues = publication.check_text("results.json", json_bytes(value))
    assert categories(issues) == {"operational-key"}
    output = "\n".join(map(str, issues))
    assert "$[key#0][0][key#0:parent_campaign]" in output
    assert "$[key#0][0][key#1][key#0:SLURM_*]" in output
    assert "private-campaign-label" not in output
    assert "987654321" not in output
    assert "SLURM_JOB_ID" not in output


@pytest.mark.parametrize("key", sorted(publication.FORBIDDEN_KEYS))
def test_all_operational_keys(key):
    issues = publication.check_text("results.json", json_bytes({"nested": {key: None}}))
    assert "operational-key" in categories(issues)


def test_json_keys_duplicate_keys_freeform_and_encoded_secrets():
    token = github_token()
    data = ('{"nested": [{"' + token + '": "ordinary"}], "value": "' + token +
            '", "value": "overwritten", "stderr_hex": "' + token.encode().hex() + '"}').encode()
    issues = publication.check_text("results.json", data)
    assert sum(item.category == "credential/github" for item in issues) == 3
    assert any(":key:" in item.location for item in issues)
    assert any(":decoded:" in item.location for item in issues)
    assert token not in str(issues)
    assert token.encode().hex() not in str(issues)


def test_credentials():
    header = base64.urlsafe_b64encode(b'{"alg":"HS256","typ":"JWT"}').decode().rstrip("=")
    claims = base64.urlsafe_b64encode(b'{"sub":"synthetic-subject"}').decode().rstrip("=")
    jwt = ".".join((header, claims, "A9bc" * 10))
    pem = ("-----BEGIN " + "PRIVATE KEY-----\n" +
           base64.b64encode(b"synthetic-not-a-real-key" * 5).decode() +
           "\n-----END " + "PRIVATE KEY-----")
    samples = [
        (github_token(), "credential/github"),
        ("github_" + "pat_" + "Ab3" * 8 + "_" + "X9y" * 20, "credential/github"),
        ("AK" + "IA" + "A1B2" * 4, "credential/aws"),
        ("xo" + "xb-" + "1234567890-" * 2 + "a1b2" * 8, "credential/slack"),
        ("https://hooks.slack.com/" + "services/T" + "A1" * 5 + "/B" + "B2" * 5 +
         "/" + "a3" * 12, "credential/slack"),
        (jwt, "credential/jwt"), (pem, "credential/pem"),
        ("https://synthetic-user:" + "not-a-real-password-83" + "@service.example.org/repo",
         "credential/url"),
    ]
    for value, expected in samples:
        issues = publication.check_text("notes.txt", value.encode())
        assert expected in categories(issues)
        assert value not in str(issues)


def test_science_public_provenance_and_documentation_are_allowed():
    document = {
        "field": 0.01, "detuning": 0.0, "methods": {"cpu_threads": 4, "mem_gib": 3.2,
        "wall_seconds": 12, "software_versions": {"python": "3.11"}},
        "source_fingerprint": "a1b2c3d4" * 8,
        "source_file": "NRG_comparisons/test4/physics.py",
        "authors": ["Ada Lovelace", "Emmy Noether"], "doi": "10.1234/example.5678",
        "repository": "https://github.com/example/scientific-project",
        "license": ["MIT", "CC-BY"], "error": "No valid password found",
        "instructions": "/absolute/path/to/nrg /path/to/nrg /tmp/work /home/<user>/nrg "
                        "/project/${USER}/nrg /scratch/{account}/job ssh user@host "
                        "ssh example.org --account=<account>",
        "pem_documentation": "-----BEGIN " + "PRIVATE KEY-----",
    }
    assert publication.check_text("results.json", json_bytes(document)) == []
    assert publication.check_text("README.md", json_bytes(document)) == []
    assert publication.check_text("script.py", b'import subprocess\nsubprocess.run(["ssh", host])') == []
    assert publication.check_text("README.md", b"https://user:password@example.org/repo") == []
    assert publication.check_text("README.md", b"/projected/field /homework/solution /scratchpad") == []


@pytest.mark.parametrize("prefix", ["/home", "/project", "/scratch", "/Users", "/share/apps"])
def test_private_paths_in_text_keys_and_decoded_logs(prefix):
    value = prefix + "/example-user/software/nrg"
    issues = publication.check_text("output.json", json_bytes({value: [value],
                                     "stdout_hex": value.encode().hex()}))
    assert sum(item.category == "private-path" for item in issues) == 3
    assert value not in str(issues)
    assert "example-user" not in str(issues)


@pytest.mark.parametrize("prefix", ["/home", "/project", "/scratch", "/Users", "/share/apps"])
def test_bare_roots_are_not_identifying_paths(prefix):
    for value in (prefix, prefix + "/"):
        assert publication.check_text("notes.txt", value.encode()) == []
        assert publication.check_text("data.json", json_bytes({value: value})) == []
        assert publication.check_text("data.json", json_bytes(value)) == []
        assert publication.check_text("runner.py", f"root = {value!r}\n".encode()) == []
        issues = publication.check_text("data.json", json_bytes({"parent_campaign": value}))
        assert categories(issues) == {"operational-key"}


@pytest.mark.parametrize("pattern", [
    r"sha256|task_id|budget|attempt-|/project/|\b[0-9a-f]{64}\b",
    r"(?:/project/|/home/|/tmp/|/share/|SLURM_|SBATCH_|Traceback)",
    r"^/project/(?:[^/]+/)?$", r"/project/([a-z-]+)", r"/project/[a-z-]+",
    r"/project/\w+", r"/project/.*", r"/project/.+", r"/project/$",
])
def test_denial_regexes_are_not_concrete_private_paths(pattern):
    for literal in (repr(pattern), 'r"' + pattern + '"'):
        source = f"import re\nDENY = re.compile(\n    {literal}\n)\n"
        assert publication.check_text("generate.py", source.encode()) == []
    assert publication.check_text("data.json", json_bytes({"pattern": pattern})) == []


@pytest.mark.parametrize("prefix", ["/home", "/project", "/scratch", "/Users", "/share/apps"])
def test_concrete_paths_still_fail_in_python_regexes_and_delimited_data(prefix):
    value = prefix + "/example-user/software/nrg"
    for source in (f"path = {value!r}\n", f'path = "{value}"\n',
                   f'import re\nDENY = re.compile(r"{value}/.*")\n'):
        issues = publication.check_text("runner.py", source.encode())
        assert "private-path" in categories(issues)
        assert value not in str(issues)
    for ending in ("\n", '"', "'", r"|\bother", r"\s+", "(capture)"):
        text = "before\n" + value + ending
        issues = publication.check_text("notes.txt", text.encode())
        assert categories(issues) == {"private-path"}
        assert issues[0].location == "text:line 2"
        issues = publication.check_text("data.json", json_bytes({value + ending: value + ending}))
        assert sum(item.category == "private-path" for item in issues) == 2
        assert value not in str(issues)


def test_python_literals_site_account_and_escaped_path():
    code = ('path = "\\x2fhome/example-user/software"\n'
            'command = ["ssh", "researcher@compute.private-cluster"]\n'
            'account = ["sbatch", "--account", "synthetic-allocation"]\n'
            'SSH_HOST = "compute.private-cluster"\n'
            'SLURM_ACCOUNT = "synthetic-allocation"\n'
            'metadata = {"SLURM_JOB_ID": job_id}\n')
    issues = publication.check_text("runner.py", code.encode())
    assert {"private-path", "ssh-site", "scheduler-account", "site-literal"} <= categories(issues)
    assert "operational-key" not in categories(issues)
    assert "compute.private-cluster" not in str(issues)
    assert "synthetic-allocation" not in str(issues)


def test_environment_strings_and_site_defaults_not_generic_slurm_code():
    value = "USER=example-user SLURM_JOB_ID=12345 HOSTNAME=compute.private-cluster"
    issues = publication.check_text("results.json", json_bytes({"notes": value,
                                     "stdout_hex": value.encode().hex()}))
    assert "environment-string" in categories(issues)
    assert "example-user" not in str(issues)
    assert "12345" not in str(issues)
    assert "site-literal" in categories(publication.check_text(
        "runner.py", b'parser.add_argument("--account", default="synthetic-allocation")'))
    assert publication.check_text("runner.py", b'SLURM_ACCOUNT = os.environ["ACCOUNT"]') == []
    assert publication.check_text("runner.py", b'SSH_HOST = "cluster.example.org"') == []
    assert publication.check_text("runner.py", b'command = ["ssh", "-o", "BatchMode=yes", host, "echo"]') == []
    assert publication.check_text("README.md", b'Use `ssh -p 22 -i "key file" user@host`.') == []
    assert "ssh-site" in categories(publication.check_text(
        "runner.py", b'command = ["ssh", "-p", "22", "researcher@compute.private-cluster"]'))


def test_malformed_jwt_and_pem_candidates_do_not_crash_or_false_positive():
    payload = base64.urlsafe_b64encode(b'{"sub":"synthetic-subject"}').decode().rstrip("=")
    for header in (b'{"alg":[]}', b'{"alg":"none"}', b'not a JSON header'):
        value = ".".join((base64.urlsafe_b64encode(header).decode().rstrip("="), payload, "A9bc" * 10))
        assert publication.check_text("notes.txt", value.encode()) == []
    empty_pem = "-----BEGIN " + "PRIVATE KEY-----\n" + " " * 100 + "-----END " + "PRIVATE KEY-----"
    assert publication.check_text("notes.txt", empty_pem.encode()) == []


def test_invalid_inputs_and_binary():
    assert categories(publication.check_text("a.json", b"{broken")) == {"invalid-json"}
    assert categories(publication.check_text("a.json", b"\0\xff")) == {"invalid-json"}
    assert categories(publication.check_text("a.py", b"def (")) == {"invalid-python"}
    assert "invalid-hex" in categories(publication.check_text("a.json", b'{"stdout_hex":"not hex"}'))
    assert publication.check_text("a.bin", b"\0\xff") == []
    assert "credential/github" in categories(publication.check_text(
        "a.txt", github_token().encode("utf-16")))


def test_test1_exemption_is_semantic_only():
    name = "package/NRG_comparisons/test1/output/results.json"
    data = json_bytes({"parent_campaign": "/home/example-user/private"})
    assert publication.check_text(name, data, max_json_bytes=1) == []
    issues = publication.check_text(name, json_bytes({"value": github_token()}))
    assert categories(issues) == {"credential/github"}
    assert "operational-key" in categories(publication.check_text(
        "test1/NRG_comparisons/test4/results.json", b'{"uid":1}'))


def test_default_git_selection_includes_manual_code_not_ignored_runs(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True, capture_output=True)
    write(tmp_path / ".gitignore", b"**/runs/\nNRG_comparisons/test4/manual_run/\n")
    tracked = write(tmp_path / "NRG_comparisons/test2/manual_run/run_one", b"/home/example-user/nrg")
    subprocess.run(["git", "add", "--", str(tracked)], cwd=tmp_path, check=True, capture_output=True)
    write(tmp_path / "NRG_comparisons/test3/manual_run/nrg.py", b'path = "/scratch/example-user/nrg"')
    write(tmp_path / "NRG_comparisons/test4/output/results.json")
    for name in ("test2/runs/private.json", "test3/runs/private.json", "test4/manual_run/raw.json",
                 "test1/output/results.json"):
        write(tmp_path / "NRG_comparisons" / name, github_token().encode())
    report = publication.check_paths(root=tmp_path)
    assert report.files == 3
    assert categories(report.diagnostics) == {"private-path"}
    assert {"test2", "test3"} == {item.path.split("/")[1] for item in report.diagnostics}
    # Forced tracked private material is rejected rather than hidden by ignores.
    raw = tmp_path / "NRG_comparisons/test2/runs/private.json"
    subprocess.run(["git", "add", "-f", "--", str(raw)], cwd=tmp_path, check=True, capture_output=True)
    assert "private-member" in categories(publication.check_paths(root=tmp_path).diagnostics)


@pytest.mark.parametrize("read_error", [None, FileNotFoundError, PermissionError])
def test_default_git_selection_excludes_only_snapshotted_deletions(tmp_path, monkeypatch, read_error):
    deleted = "NRG_comparisons/test4/analyze.py"
    tracked = "NRG_comparisons/test2/output/results.json"
    untracked = "NRG_comparisons/test4/new.json"
    write(tmp_path / tracked)
    write(tmp_path / untracked, b'{"field":0.1}')
    unknown = write(tmp_path / "unknown-untracked.txt", github_token().encode())
    calls = []

    def fake_git(command, **kwargs):
        assert kwargs == {"cwd": tmp_path, "capture_output": True, "check": True}
        if not calls:
            assert command == ["git", "ls-files", "--deleted", "-z", "--", *publication.DEFAULT_PATHS]
            stdout = deleted.encode() + b"\0"
        else:
            assert command == ["git", "ls-files", "--cached", "--others", "--exclude-standard",
                               "-z", "--", *publication.DEFAULT_PATHS]
            stdout = b"\0".join(name.encode() for name in (tracked, deleted, untracked)) + b"\0"
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, stdout=stdout)

    original_open = publication.os.open

    def open_file(path, flags, *args, **kwargs):
        if read_error is not None and path == tmp_path / tracked:
            assert len(calls) == 2  # The file becomes unreadable after discovery.
            raise read_error("synthetic private error details")
        return original_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(publication.subprocess, "run", fake_git)
    monkeypatch.setattr(publication.os, "open", open_file)
    report = publication.check_paths(root=tmp_path)
    assert len(calls) == 2
    assert report.excluded == 1
    assert report.files == (2 if read_error is None else 1)
    assert categories(report.diagnostics) == (set() if read_error is None else {"read-error"})
    assert all(item.path == tracked for item in report.diagnostics)
    assert "synthetic private error details" not in str(report.diagnostics)
    assert not (tmp_path / deleted).exists()
    assert (tmp_path / untracked).read_bytes() == b'{"field":0.1}'
    assert unknown.read_bytes() == github_token().encode()


def test_explicit_offline_walk_and_symlinks(tmp_path):
    root = tmp_path / "export"
    write(root / "NRG_comparisons/test2/manual_run/run_one", b"generic nrg command")
    write(root / "NRG_comparisons/test3/manual_run/nrg.py", b"field = 0.1\n")
    for name in ("test2/runs/raw.json", "test2/manual_run/1/param",
                 "test3/manual_run/processed/summary.json", "test4/manual_run/raw.json"):
        write(root / "NRG_comparisons" / name, github_token().encode())
    write(root / "NRG_comparisons/test4/output/plot.svg", b'<svg><text>Delta</text></svg>')
    write(root / "NRG_comparisons/test4/output/plot.png", b"\x89PNG\0synthetic")
    report = publication.check_paths([root], root=tmp_path)
    assert report.diagnostics == []
    assert (report.files, report.text_files, report.binary_skipped, report.excluded) == (4, 3, 1, 4)
    outside = write(tmp_path / "outside.txt", github_token().encode())
    (root / "linked.txt").symlink_to(outside)
    (root / "linked-dir").symlink_to(tmp_path, target_is_directory=True)
    report = publication.check_paths([root], root=tmp_path)
    assert categories(report.diagnostics) == {"symlink"}
    assert len(report.diagnostics) == 2
    assert outside.read_bytes() == github_token().encode()
    report = publication.check_paths([root / "linked-dir/outside.txt"], root=tmp_path)
    assert categories(report.diagnostics) == {"symlink"}


def test_size_limits_before_read_and_case_totals(tmp_path):
    paths = [write(tmp_path / "NRG_comparisons/test4/output" / name, b" " * 20)
             for name in ("a.json", "b.json")]
    report = publication.check_paths(paths, root=tmp_path, max_json_bytes=10, max_case_bytes=30)
    assert [item.category for item in report.diagnostics].count("json-size") == 2
    assert [item.category for item in report.diagnostics].count("case-size") == 1
    assert report.case_bytes == {"test4": 40}
    assert report.text_files == 0
    report = publication.Report()

    def must_not_read(_size):
        pytest.fail("oversized JSON must not be read")

    report.content("test4/results.json", 166 * publication.MIB, must_not_read)
    assert categories(report.diagnostics) == {"json-size", "case-size"}
    # Standalone exports outside --root retain their case identity.
    report = publication.check_paths(paths, root=tmp_path / "unrelated", max_case_bytes=30)
    assert "case-size" in categories(report.diagnostics)


@pytest.mark.parametrize("suffix", [".tar.gz", ".whl"])
def test_archives_content_semantics_global_credentials_and_no_extraction(tmp_path, monkeypatch, suffix):
    def no_extract(*args, **kwargs):
        pytest.fail("archives must never be extracted")

    monkeypatch.setattr(tarfile.TarFile, "extractall", no_extract)
    monkeypatch.setattr(tarfile.TarFile, "extract", no_extract)
    monkeypatch.setattr(zipfile.ZipFile, "extractall", no_extract)
    monkeypatch.setattr(zipfile.ZipFile, "extract", no_extract)
    members = [
        ("pkg/NRG_comparisons/test4/results.json", json_bytes({"parent_campaign": "hidden"})),
        ("pkg/NRG_comparisons/test1/results.json", json_bytes({"parent_campaign": "untouched"})),
        ("pkg/elsewhere/notes.txt", github_token().encode()),
        ("pkg/NRG_comparisons/test4/plot.svg", b'<svg><text>/home/example-user/nrg</text></svg>'),
        ("pkg/lib/core.so", b"\x00native"),
    ]
    path = make_archive(tmp_path / ("public" + suffix), members)
    before = sorted(tmp_path.iterdir())
    report = publication.check_archive(path)
    assert categories(report.diagnostics) == {"operational-key", "credential/github", "private-path"}
    assert report.files == 5 and report.binary_skipped == 1
    assert all("test1" not in item.path for item in report.diagnostics)
    assert sorted(tmp_path.iterdir()) == before
    assert github_token() not in str(report.diagnostics)


@pytest.mark.parametrize("suffix", [".tar.gz", ".whl"])
def test_archive_private_members_at_any_root_and_traversal(tmp_path, suffix):
    names = ["APPLICATIONS.md", "arbitrary/runs/result.json", "nested/checkpoints/state.dat",
             "root/shared_venv/lib/module.py", "root/private_archive/data.txt", "root/raw.h5",
             "root/private.zip", "root/NRG_comparisons/test4/manual_run/nrg.py",
             "root/NRG_comparisons/test2/manual_run/1/param", "../outside.txt",
             "/absolute.txt", "root/./file.txt", "root\\file.txt", "C:/file.txt"]
    path = make_archive(tmp_path / ("unsafe" + suffix), [(name, b"x") for name in names])
    report = publication.check_archive(path)
    assert len(report.diagnostics) == len(names)
    assert categories(report.diagnostics) == {"private-member", "unsafe-member"}
    assert report.files == 0


@pytest.mark.parametrize("kind", ["symlink", "hardlink", "fifo", "zip-symlink"])
def test_archive_special_members(tmp_path, kind):
    if kind == "zip-symlink":
        path = tmp_path / "links.whl"
        with zipfile.ZipFile(path, "w") as archive:
            info = zipfile.ZipInfo("pkg/link")
            info.create_system = 3
            info.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive.writestr(info, "/home/example-user/private")
    else:
        path = tmp_path / "links.tar.gz"
        with tarfile.open(path, "w:gz") as archive:
            info = tarfile.TarInfo("pkg/link")
            info.type = {"symlink": tarfile.SYMTYPE, "hardlink": tarfile.LNKTYPE,
                         "fifo": tarfile.FIFOTYPE}[kind]
            info.linkname = "/home/example-user/private"
            archive.addfile(info)
    assert categories(publication.check_archive(path).diagnostics) == {"special-member"}


def test_archive_budget_duplicates_and_errors(tmp_path):
    path = make_archive(tmp_path / "budget.tar.gz", [
        ("pkg/NRG_comparisons/test3/a.json", b"{}"),
        ("pkg/NRG_comparisons/test3/b.json", b"{}"),
        ("pkg/NRG_comparisons/test3/a.json", b"{}"),
    ])
    report = publication.check_archive(path, max_json_bytes=1, max_case_bytes=3)
    assert categories(report.diagnostics) == {"json-size", "case-size", "duplicate-member"}
    broken = write(tmp_path / "broken.whl", b"not a zip")
    assert categories(publication.check_archive(broken).diagnostics) == {"archive-read-error"}
    assert categories(publication.check_archive(tmp_path / "bad.txt").diagnostics) == {"unsupported-archive"}
    empty = make_archive(tmp_path / "empty.whl", [])
    assert categories(publication.check_archive(empty).diagnostics) == {"empty-archive"}


def test_diagnostic_labels_are_redacted_too():
    token = github_token()
    issues = publication.check_text("notes/" + token + ".txt", token.encode())
    assert token not in str(issues)
    assert "<redacted>" in issues[0].path
    issues = publication.check_text("/home/example-user/private/results.json", b'{"uid":1}')
    assert "example-user" not in str(issues)


def test_source_names_missing_files_and_empty_selection(tmp_path):
    path = write(tmp_path / (github_token() + ".txt"), b"ordinary notes")
    report = publication.check_paths([path], root=tmp_path)
    assert categories(report.diagnostics) == {"credential/github"}
    assert github_token() not in str(report.diagnostics)
    report = publication.check_paths([tmp_path / "missing.txt"], root=tmp_path)
    assert categories(report.diagnostics) == {"read-error"}
    empty = tmp_path / "empty"
    empty.mkdir()
    assert categories(publication.check_paths([empty], root=tmp_path).diagnostics) == {"empty-selection"}


def test_cli_only_explicit_targets_counts_exit_status_and_no_git(tmp_path, monkeypatch, capsys):
    def no_git(*args, **kwargs):
        pytest.fail("explicit targets must work without Git")

    monkeypatch.setattr(publication.subprocess, "run", no_git)
    path = write(tmp_path / "test4/results.json", b'{"field":0.1}')
    monkeypatch.chdir(tmp_path)
    assert publication.main(["test4/results.json"]) == 0
    assert "files=1, text=1, binary_skipped=0, excluded=0, failures=0" in capsys.readouterr().out
    path.write_bytes(json_bytes({"parent_campaign": github_token()}))
    assert publication.main([str(path)]) == 1
    output = capsys.readouterr().out
    assert "operational-key" in output and "credential/github" in output
    assert github_token() not in output
    path = make_archive(tmp_path / "public.whl", [("pkg/notes.txt", b"public notes")])
    assert publication.main([str(path)]) == 0
    assert "failures=0" in capsys.readouterr().out
    with pytest.raises(SystemExit) as error:
        publication.main(["--max-json-mib", "nan", str(path)])
    assert error.value.code == 2


def test_default_selection_fails_closed_without_git(tmp_path, monkeypatch):
    def missing_git(*args, **kwargs):
        raise FileNotFoundError("sensitive local details")

    monkeypatch.setattr(publication.subprocess, "run", missing_git)
    report = publication.check_paths(root=tmp_path)
    assert categories(report.diagnostics) == {"git-unavailable"}
    assert "sensitive local details" not in str(report.diagnostics)


def test_checker_and_synthetic_tests_ship_without_literal_credentials():
    for path in (Path(publication.__file__), Path(__file__)):
        assert publication.check_text(path.name, path.read_bytes(), semantic=False) == []
