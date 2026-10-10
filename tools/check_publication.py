"""Read-only, offline publication privacy gate (not a dataset schema validator).

With no targets, inspect Git-visible files in NRG_comparisons/test{2,3,4}:
    python -B tools/check_publication.py
Explicit files/directories work without Git; raw working directories are pruned:
    python -B tools/check_publication.py export/NRG_comparisons/test4
Check fresh artifacts before uploading them (no extraction or source-tree scan):
    python -B tools/check_publication.py dist/*.tar.gz dist/*.whl

Archives have no ignore exclusions: private members are failures. Credentials
are checked in all archive text; semantic metadata and dataset size policies
apply only to Test2/3/4. Test1 is exempt from these policies. Reads of other text
are capped at 32 MiB per file as a resource safeguard. Binary files are
counted as skipped, not certified free of secrets or image/EXIF metadata. No
network authentication, writes, imports of inspected code, or secret values in
diagnostics. This complements, rather than replaces, check_distribution.py.
"""

import argparse
import ast
import base64
from dataclasses import dataclass, field
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import stat
import subprocess
import tarfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATHS = tuple(f"NRG_comparisons/test{i}" for i in (2, 3, 4))
MIB = 1024 * 1024
MAX_JSON_BYTES = 8 * MIB
MAX_CASE_BYTES = 16 * MIB
MAX_TEXT_BYTES = 32 * MIB

FORBIDDEN_KEYS = frozenset({
    "parent_campaign", "execution_metadata", "job_id", "job_file", "scheduler",
    "scheduler_preflight", "source_jobmeta_path", "authorization", "source_archive_paths",
    "raw_command", "rawcommand", "stdout_hex", "stderr_hex", "stdout", "stderr",
    "command_stdout", "command_stderr", "hostname", "host_name", "user", "username",
    "uid", "gid", "environment", "environ", "env", "full_environment",
    "environment_variables", "ssh_host", "remote_host", "scheduler_account",
})
PRIVATE_DIRS = frozenset({
    "runs", "raw_runs", "raw-runs", "nrg_runs", "nrg-runs", "raw", "checkpoints",
    "nrg_raw", "nrg-raw", "raw_nrg", "raw-nrg",
    "checkpoint", ".venv", "venv", "shared_venv", "shared-venv", ".git",
    "__pycache__", ".pytest_cache", ".ruff_cache", "private", "private_archive",
    "private-archive", "private_archives", "private-archives",
})
MANUAL_FILES = {
    "test2": {"README.md", ".gitignore", "param.loop", "1_run", "run_one", "2_proc",
              "process.py", "extract.py", "reference.json"},
    "test3": {"README.md", ".gitignore", "param.loop", "1_run", "run_one", "2_proc",
              "nrg.py", "physics.py", "process.py"},
}
BINARY_SUFFIXES = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".webp", ".ico", ".pdf", ".so", ".pyd",
    ".dll", ".dylib", ".a", ".o", ".ttf", ".woff", ".woff2",
})
PRIVATE_PATH = re.compile(
    # Bare roots and regex-only components do not identify an installation.
    r"(?<![\w./])/(?:home|project|scratch|Users|share/apps)/"
    r"(?P<component>[\w.-]*\w[\w.-]*)(?:/[^\s\"'`;,()\[\]{}|\\^$*?+<>]*)?"
)
PLACEHOLDER = re.compile(r"^(?:<[^>]+>|\$\{[^}]+\}|\{[^}]+\})(?:/|$)")
ACCOUNT = re.compile(r"(?:--account(?:=|\s+)|#SBATCH\s+-A\s+)([\w.${}<>-]+)")
SSH = re.compile(r"\bssh[ \t]+([^\r\n`]+)")
ENVIRONMENT = re.compile(r"\b(?:SLURM_[A-Z0-9_]+|HOSTNAME|USER|LOGNAME|UID|GID)="
                         r"([\w./@${}<>-]+)")
SITE_NAMES = frozenset({"ssh_host", "remote_host", "cluster_host", "slurm_account",
                        "scheduler_account", "cluster_account", "account"})
CREDENTIAL_PATTERNS = (
    ("credential/github", re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{36,255}|"
                                     r"github_pat_[A-Za-z0-9]{20,}_[A-Za-z0-9]{40,})\b")),
    ("credential/aws", re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b")),
    ("credential/slack", re.compile(r"\bxox[baprs]-(?:[0-9]{8,}-){2,3}[A-Za-z0-9]{16,}\b")),
    ("credential/slack", re.compile(r"https://hooks\.slack\.com/services/"
                                    r"T[A-Z0-9]{8,}/B[A-Z0-9]{8,}/[A-Za-z0-9]{20,}")),
    ("credential/pem", re.compile(
        r"-----BEGIN ((?:(?:RSA|EC|DSA|OPENSSH|ENCRYPTED) )?PRIVATE KEY)-----"
        r"([A-Za-z0-9+/=\r\n\t ]{64,16384})-----END \1-----")),
)
JWT = re.compile(r"\b[A-Za-z0-9_-]{12,}\.[A-Za-z0-9_-]{8,}\.[A-Za-z0-9_-]{16,}\b")
CREDENTIAL_URL = re.compile(
    r"\b[a-z][a-z0-9+.-]*://([^\s/:@\"'<>`]+):([^\s/@\"'<>`]+)@[^\s/\"'<>`]+", re.I
)


def _placeholder(value):
    return (bool(PLACEHOLDER.match(value)) or value.startswith("$") or
            value in {"USER", "USERNAME", "HOST", "SSH_HOST", "ACCOUNT", "SLURM_ACCOUNT",
                      "user", "username", "host", "hostname", "account", "password",
                      "changeme", "example", "..."})


def _site_placeholder(value):
    host = value.rsplit("@", 1)[-1]
    return (_placeholder(value) or _placeholder(host) or
            host in {"localhost", "127.0.0.1", "::1"} or
            bool(re.fullmatch(r"(?:[\w-]+\.)*example\.(?:com|org|net)", host)))


def _matches(text, semantic):
    """Yield category and spans only; matched bytes never enter a diagnostic."""
    for category, pattern in CREDENTIAL_PATTERNS:
        for match in pattern.finditer(text):
            if category == "credential/pem":
                body = re.sub(r"\s", "", match[2])
                try:
                    if len(body) < 64 or not base64.b64decode(body, validate=True):
                        continue
                except ValueError:
                    continue
            yield category, match.start(), match.end()
    for match in JWT.finditer(text):
        try:
            header, payload, _ = match.group().split(".")
            header = json.loads(base64.urlsafe_b64decode(header + "=" * (-len(header) % 4)))
            payload = json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))
        except (ValueError, UnicodeError, RecursionError):
            continue
        if (isinstance(header, dict) and isinstance(header.get("alg"), str) and
                header["alg"].lower() != "none" and
                isinstance(payload, dict)):
            yield "credential/jwt", match.start(), match.end()
    for match in CREDENTIAL_URL.finditer(text):
        if not _placeholder(match[2]):
            yield "credential/url", match.start(), match.end()
    if not semantic:
        return
    for match in PRIVATE_PATH.finditer(text):
        if not _placeholder(match["component"]):
            yield "private-path", match.start(), match.end()
    for match in ACCOUNT.finditer(text):
        if not _site_placeholder(match[1]):
            yield "scheduler-account", match.start(), match.end()
    for match in SSH.finditer(text):
        try:
            arguments = iter(shlex.split(match[1], comments=True))
        except ValueError:
            continue
        for argument in arguments:
            if argument in {"-B", "-b", "-c", "-D", "-E", "-e", "-F", "-I", "-i", "-J",
                            "-L", "-l", "-m", "-O", "-o", "-p", "-Q", "-R", "-S", "-W", "-w"}:
                next(arguments, None)
                continue
            if argument.startswith("-"):
                continue
            if not _site_placeholder(argument):
                yield "ssh-site", match.start(), match.end()
            break


def _label(name):
    # Absolute caller paths are not public identifiers. Archive/member labels
    # also need redaction: credentials and private paths can occur in names.
    name = os.fspath(name)
    if name.startswith("/"):
        marker = name.find("/NRG_comparisons/")
        name = name[marker + 1:] if marker >= 0 else "<absolute>/" + Path(name).name
    spans = sorted((start, end) for _, start, end in _matches(name, True))
    result, end = "", 0
    for start, stop in spans:
        if start >= end:
            result += name[end:start] + "<redacted>"
        end = max(end, stop)
    return (result + name[end:]).encode("unicode_escape").decode("ascii")


@dataclass(frozen=True)
class Diagnostic:
    path: str
    category: str
    location: str

    def __str__(self):
        return f"{self.path}: {self.category} at {self.location} (value redacted)"


def _case(name):
    parts = PurePosixPath(os.fspath(name)).parts
    for parent, part in zip(parts, parts[1:], strict=False):
        if parent == "NRG_comparisons" and part in {"test1", "test2", "test3", "test4"}:
            return part
    for part in parts:
        if part in {"test1", "test2", "test3", "test4"}:
            return part
    return None


class _Object(list):
    """Preserve duplicate JSON keys so an overwritten secret cannot disappear."""


def _decode(data):
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16")
    if b"\0" in data:
        raise UnicodeError
    return data.decode("utf-8-sig")


def check_text(name, data: bytes, *, semantic=True,
               max_json_bytes=MAX_JSON_BYTES) -> list[Diagnostic]:
    """Check text/JSON/Python bytes; return redacted path/category/location records.

    JSON locations use zero-based key ordinals, not arbitrary key contents.
    Known prohibited keys are included in the location; values never are.
    Binary input returns no diagnostics; the file/archive APIs count skips.
    """
    semantic = semantic and _case(name) != "test1"
    label = _label(name)
    diagnostics = []

    def add(category, location):
        diagnostics.append(Diagnostic(label, category, location))

    def scan(text, location, *, environment=False):
        for category, start, _ in _matches(text, semantic):
            add(category, f"{location}:line {text.count(chr(10), 0, start) + 1}")
        if semantic and environment:
            for match in ENVIRONMENT.finditer(text):
                if not _placeholder(match[1]):
                    add("environment-string", location)

    is_json = PurePosixPath(name).suffix.lower() == ".json"
    if semantic and is_json and len(data) > max_json_bytes:
        add("json-size", "file")
        return diagnostics
    try:
        text = _decode(data)
    except UnicodeError:
        if semantic and is_json:
            add("invalid-json", "file")
        return diagnostics
    if not is_json:
        scan(text, "text")
    else:
        try:
            document = json.loads(text, object_pairs_hook=_Object)
        except (ValueError, RecursionError):
            scan(text, "text")
            if semantic:
                add("invalid-json", "file")
            return diagnostics
        stack = [(document, "$", None)]
        while stack:
            value, location, key = stack.pop()
            if isinstance(value, _Object):
                for index, (child_key, child) in enumerate(value):
                    normalized = child_key.lower().replace("-", "_")
                    forbidden = normalized in FORBIDDEN_KEYS or normalized.startswith("slurm_")
                    known = normalized if normalized in FORBIDDEN_KEYS else "SLURM_*"
                    child_location = f"{location}[key#{index}{':' + known if forbidden else ''}]"
                    scan(child_key, child_location + ":key")
                    if semantic and forbidden:
                        add("operational-key", child_location)
                    stack.append((child, child_location, normalized))
            elif isinstance(value, list):
                stack.extend((child, f"{location}[{i}]", None) for i, child in enumerate(value))
            elif isinstance(value, str):
                scan(value, location, environment=True)
                if key in {"stdout_hex", "stderr_hex"}:
                    try:
                        decoded = bytes.fromhex(value).decode("utf-8", errors="replace")
                    except ValueError:
                        if semantic:
                            add("invalid-hex", location)
                    else:
                        scan(decoded, location + ":decoded", environment=True)
    if PurePosixPath(name).suffix.lower() == ".py":
        try:
            tree = ast.parse(text)
        except (SyntaxError, ValueError, RecursionError):
            if semantic:
                add("invalid-python", "file")
        else:
            for node in ast.walk(tree):
                location = f"python:line {getattr(node, 'lineno', 1)}"
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    scan(node.value, location)
                if isinstance(node, ast.Constant) and isinstance(node.value, bytes):
                    scan(node.value.decode("utf-8", errors="replace"), location)
                if not semantic:
                    continue
                if isinstance(node, (ast.List, ast.Tuple)):
                    strings = []
                    for item in node.elts:
                        if not isinstance(item, ast.Constant) or not isinstance(item.value, str):
                            break
                        strings.append(item.value)
                    scan(shlex.join(strings), location)
                if isinstance(node, (ast.Assign, ast.AnnAssign, ast.keyword)):
                    targets = (node.targets if isinstance(node, ast.Assign) else
                               [node.target] if isinstance(node, ast.AnnAssign) else [])
                    names = [target.id.lower() for target in targets if isinstance(target, ast.Name)]
                    if isinstance(node, ast.keyword):
                        names.append((node.arg or "").lower())
                    if (set(names) & SITE_NAMES
                            and isinstance(node.value, ast.Constant)
                            and isinstance(node.value.value, str)
                            and node.value.value and not _site_placeholder(node.value.value)):
                        add("site-literal", location)
                if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                        and node.func.attr == "add_argument"
                        and any(arg.value in
                                {"--account", "--ssh-host", "--remote-host", "--host"}
                                for arg in node.args if isinstance(arg, ast.Constant)
                                and isinstance(arg.value, str))):
                    for keyword in node.keywords:
                        if (keyword.arg == "default" and isinstance(keyword.value, ast.Constant)
                                and isinstance(keyword.value.value, str) and keyword.value.value
                                and not _site_placeholder(keyword.value.value)):
                            add("site-literal", location)
    return list(dict.fromkeys(diagnostics))


def _private_member(name):
    parts = PurePosixPath(name).parts
    lower = [part.lower() for part in parts]
    if any(part in PRIVATE_DIRS or part.startswith("runs-") for part in lower):
        return True
    if any(part == "applications.md" or part.startswith(".env") or
           part.endswith((".session.json", ".h5", ".hdf5", ".npz", ".checkpoint",
                          ".gz", ".bz2", ".xz", ".tar", ".tgz", ".zip", ".7z", ".whl")) or
           part in {"nrg.log", "nrginit.log", "jobmeta.json", "pyvenv.cfg"} for part in lower):
        return True
    if "manual_run" in parts:
        index = parts.index("manual_run")
        case = _case(name)
        if case not in MANUAL_FILES:
            return True
        rest = parts[index + 1:]
        return bool(rest) and (len(rest) != 1 or rest[0] not in MANUAL_FILES[case])
    return False


@dataclass
class Report:
    diagnostics: list[Diagnostic] = field(default_factory=list)
    files: int = 0
    text_files: int = 0
    binary_skipped: int = 0
    excluded: int = 0
    case_bytes: dict[str, int] = field(default_factory=dict)
    max_json_bytes: int = MAX_JSON_BYTES
    max_case_bytes: int = MAX_CASE_BYTES

    def issue(self, name, category, location="file"):
        self.diagnostics.append(Diagnostic(_label(name), category, location))

    def content(self, name, size, read, *, semantic=True):
        """Check declared size before calling the bounded, regular-file reader."""
        self.files += 1
        case = _case(name)
        semantic = semantic and case != "test1"
        if case in {"test2", "test3", "test4"}:
            previous = self.case_bytes.get(case, 0)
            self.case_bytes[case] = previous + size
            if previous <= self.max_case_bytes < previous + size:
                self.issue(case, "case-size", "case total")
        if semantic and PurePosixPath(name).suffix.lower() == ".json" and size > self.max_json_bytes:
            self.issue(name, "json-size")
            return
        if PurePosixPath(name).suffix.lower() in BINARY_SUFFIXES:
            self.binary_skipped += 1
            return
        if size > MAX_TEXT_BYTES:
            self.issue(name, "text-size")
            return
        data = read(MAX_TEXT_BYTES + 1)
        if len(data) != size:
            self.issue(name, "size-mismatch")
            return
        try:
            _decode(data)
        except UnicodeError:
            self.binary_skipped += 1
        else:
            self.text_files += 1
        self.diagnostics.extend(check_text(name, data, semantic=semantic,
                                           max_json_bytes=self.max_json_bytes))


def check_paths(paths=None, *, root=ROOT, max_json_bytes=MAX_JSON_BYTES,
                max_case_bytes=MAX_CASE_BYTES) -> Report:
    """Default: cached + unignored untracked Git files, minus worktree deletions.

    Never descend through symlinks or raw directories. A private Git-visible
    candidate is a failure, not silently ignored. Explicit directory walks
    prune known private working material; explicit private files still fail.
    """
    report = Report(max_json_bytes=max_json_bytes, max_case_bytes=max_case_bytes)
    root = Path(root).absolute()
    candidates = []
    if paths is None:
        try:
            # Snapshot pending deletions first; later disappearances must fail.
            deleted_result = subprocess.run(
                ["git", "ls-files", "--deleted", "-z", "--", *DEFAULT_PATHS],
                cwd=root, capture_output=True, check=True,
            )
            result = subprocess.run(
                ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z",
                 "--", *DEFAULT_PATHS], cwd=root, capture_output=True, check=True,
            )
        except (OSError, subprocess.CalledProcessError):
            report.issue("<selection>", "git-unavailable", "use explicit paths without Git")
            return report
        deleted = {os.fsdecode(name) for name in deleted_result.stdout.split(b"\0") if name}
        names = {os.fsdecode(name) for name in result.stdout.split(b"\0") if name}
        report.excluded += len(names & deleted)
        candidates = [(root / name, name) for name in sorted(names - deleted)]
        if not candidates:
            report.issue("<selection>", "empty-selection")
    else:
        for target in paths:
            target = Path(target)
            target = target if target.is_absolute() else root / target
            try:
                prefix = target.relative_to(root).as_posix()
            except ValueError:
                parts = target.parts
                case = _case(target.as_posix())
                prefix = (PurePosixPath(*parts[parts.index(case):]).as_posix()
                          if case else target.name)
            if any(parent.is_symlink() for parent in (target, *target.parents)):
                report.issue(prefix, "symlink")
                continue
            if not target.is_dir():
                candidates.append((target, prefix))
                continue
            if _private_member(prefix):
                report.issue(prefix, "private-member")
                continue

            def walk_error(_error, prefix=prefix):
                report.issue(prefix, "read-error")

            for directory, dirs, files in os.walk(target, followlinks=False, onerror=walk_error):
                relative = Path(directory).relative_to(target)
                for entry in sorted(dirs[:]):
                    path = Path(directory) / entry
                    name = (PurePosixPath(prefix) / relative / entry).as_posix()
                    if path.is_symlink():
                        report.issue(name, "symlink")
                        dirs.remove(entry)
                    elif _private_member(name):
                        report.excluded += 1
                        dirs.remove(entry)
                for entry in sorted(files):
                    path = Path(directory) / entry
                    name = (PurePosixPath(prefix) / relative / entry).as_posix()
                    if _private_member(name):
                        report.excluded += 1
                    else:
                        candidates.append((path, name))
        if not candidates and not report.diagnostics:
            report.issue("<selection>", "empty-selection")
    seen = set()
    for path, name in candidates:
        if path in seen:
            continue
        seen.add(path)
        if _private_member(name):
            report.issue(name, "private-member")
            continue
        for category, _, _ in _matches(name, True):
            report.issue(name, category, "file name")
        try:
            if any(parent.is_symlink() for parent in (path, *path.parents)):
                report.issue(name, "symlink")
                continue
            info = path.stat()
            if not stat.S_ISREG(info.st_mode):
                report.issue(name, "special-file")
                continue
            # O_NOFOLLOW also prevents following a final-component link swapped
            # after lstat. Like Git/build tools, this assumes stable parent dirs.
            fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
            with os.fdopen(fd, "rb") as stream:
                report.content(name, info.st_size, stream.read)
        except OSError:
            report.issue(name, "read-error")
    return report


def check_archive(path, *, max_json_bytes=MAX_JSON_BYTES,
                  max_case_bytes=MAX_CASE_BYTES) -> Report:
    """Inspect .tar.gz/.whl member names and bounded contents, never extract."""
    report = Report(max_json_bytes=max_json_bytes, max_case_bytes=max_case_bytes)
    path = Path(path)
    seen = set()

    def accept(name, regular, directory):
        parts = name.removesuffix("/").split("/")
        if (any(part in {"", ".", ".."} or ":" in part for part in parts) or
                "\\" in name or any(ord(char) < 32 for char in name)):
            report.issue(name, "unsafe-member")
            return False
        if name in seen:
            report.issue(name, "duplicate-member")
            return False
        seen.add(name)
        if not regular and not directory:
            report.issue(name, "special-member")
            return False
        if _private_member(name):
            report.issue(name, "private-member")
            return False
        for category, _, _ in _matches(name, True):
            report.issue(name, category, "member name")
        return not directory

    try:
        if any(parent.is_symlink() for parent in (path, *path.parents)):
            report.issue(path.name, "symlink")
        elif path.name.endswith(".tar.gz"):
            with tarfile.open(path, "r:gz") as archive:
                for member in archive:
                    if not accept(member.name, member.isreg(), member.isdir()):
                        continue
                    with archive.extractfile(member) as stream:
                        report.content(member.name, member.size, stream.read,
                                       semantic=_case(member.name) in {"test2", "test3", "test4"})
        elif path.suffix == ".whl":
            with zipfile.ZipFile(path) as archive:
                for member in archive.infolist():
                    mode = stat.S_IFMT(member.external_attr >> 16)
                    directory = member.is_dir() and mode in {0, stat.S_IFDIR}
                    regular = not member.is_dir() and mode in {0, stat.S_IFREG}
                    if not accept(member.filename, regular, directory):
                        continue
                    with archive.open(member) as stream:
                        report.content(member.filename, member.file_size, stream.read,
                                       semantic=_case(member.filename) in {"test2", "test3", "test4"})
        else:
            report.issue(path.name, "unsupported-archive")
    except (OSError, ValueError, EOFError, RuntimeError, tarfile.TarError, zipfile.BadZipFile):
        # Exception messages may contain raw member names, passwords or paths.
        report.issue(path.name, "archive-read-error")
    if not seen and not report.diagnostics:
        report.issue(path.name, "empty-archive")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("targets", nargs="*", type=Path, help="source paths, .tar.gz or .whl")
    parser.add_argument("--root", type=Path, default=ROOT, help="root for default Git selection")
    parser.add_argument("--max-json-mib", type=float, default=8)
    parser.add_argument("--max-case-mib", type=float, default=16)
    args = parser.parse_args(argv)
    if any(not math.isfinite(value) or value <= 0
           for value in (args.max_json_mib, args.max_case_mib)):
        parser.error("size limits must be finite positive MiB values")
    limits = {"max_json_bytes": int(args.max_json_mib * MIB),
              "max_case_bytes": int(args.max_case_mib * MIB)}
    archives = [path for path in args.targets if path.name.endswith((".tar.gz", ".whl"))]
    sources = [path for path in args.targets if path not in archives]
    reports = []
    if sources or not args.targets:
        reports.append(check_paths([path.absolute() for path in sources] or None,
                                   root=args.root, **limits))
    reports.extend(check_archive(path, **limits) for path in archives)
    for report in reports:
        for diagnostic in report.diagnostics:
            print(diagnostic)
    failures = sum(len(report.diagnostics) for report in reports)
    print(f"Publication check: files={sum(r.files for r in reports)}, "
          f"text={sum(r.text_files for r in reports)}, "
          f"binary_skipped={sum(r.binary_skipped for r in reports)}, "
          f"excluded={sum(r.excluded for r in reports)}, failures={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
