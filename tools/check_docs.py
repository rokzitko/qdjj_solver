"""Fast consistency checks for the documentation (no TeX, no solver runs).

    python tools/check_docs.py              # run every check on the repository
    python tools/check_docs.py links labels # run selected checks

The user manual is built by Sphinx from the repository's Markdown pages (see
``docs/conf.py``), the physics manual lives in ``docs/manual`` and the
tutorials in ``examples``. These checks keep the three consistent: links and
anchors resolve; headings contain no mathematics, so that a page renders
identically on GitHub and on the site; physics-manual labels and quoted section
titles exist; generated files are current; public objects and option and result
fields have docstrings; tutorials match their catalogue; quoted numbers are
current; and the code named in the physics manual exists in the package. Each
check takes the repository root and returns a list of ``path:line: message``
strings.
"""

import argparse
import ast
import dataclasses
import importlib
import importlib.util
import inspect
import json
from fnmatch import fnmatch
from pathlib import Path
import re
import sys
import tomllib
from urllib.parse import unquote, urlsplit

from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[1]
SELF_URL = "https://github.com/rokzitko/qdjj_solver/blob/main/"
LABEL = re.compile(r"^(?:ch|sec|app|fig|eq|tab):[A-Za-z0-9-]+$")
LATEX_LABEL = re.compile(r"\\label\{([^}]+)\}")
EXTRA_LINKED = ("AGENTS.md",)


def _conf_lists(root):
    """``include_patterns`` and ``exclude_patterns`` read from docs/conf.py without importing it."""
    values = {}
    for node in ast.parse((root/"docs"/"conf.py").read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ("include_patterns", "exclude_patterns"):
                values[node.targets[0].id] = ast.literal_eval(node.value)
    return values["include_patterns"], values["exclude_patterns"]


def _excluded(relative, patterns):
    parts = relative.split("/")
    prefixes = ["/".join(parts[:i]) for i in range(1, len(parts) + 1)]
    return any(fnmatch(prefix, pattern) for prefix in prefixes for pattern in patterns)


def site_pages(root=ROOT):
    """The Markdown pages of the Sphinx site, as repository-relative POSIX paths."""
    include, exclude = _conf_lists(root)
    pages = set()
    for pattern in include:
        for path in root.glob(pattern):
            relative = path.relative_to(root).as_posix()
            if path.is_file() and not _excluded(relative, exclude):
                pages.add(relative)
    return sorted(pages)


def _parser():
    return MarkdownIt("commonmark").enable("table").enable("strikethrough")


def slugify(text):
    """GitHub's (and MyST's) heading anchor."""
    return re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")


class Page:
    """A parsed Markdown page: headings, anchors, links and inline code with line numbers."""

    def __init__(self, root, relative):
        self.relative = relative
        self.path = root/relative
        tokens = _parser().parse(self.path.read_text())
        self.anchors, self.headings, self.links, self.code = set(), [], [], []
        counts = {}
        for index, token in enumerate(tokens):
            if token.type == "heading_open":
                inline = tokens[index + 1]
                text = "".join(child.content for child in inline.children
                               if child.type in ("text", "code_inline"))
                line = token.map[0] + 1
                self.headings.append((line, inline.content))
                slug = slugify(text)
                count = counts.get(slug, 0)
                counts[slug] = count + 1
                self.anchors.add(slug if count == 0 else f"{slug}-{count}")
            if token.type == "inline":
                line = token.map[0] + 1 if token.map else 0
                text = ""
                for child in token.children or ():
                    if child.type == "link_open":
                        self.links.append((line, child.attrs.get("href", ""), "link"))
                    elif child.type == "image":
                        self.links.append((line, child.attrs.get("src", ""), "image"))
                    elif child.type == "code_inline":
                        self.code.append((line, child.content, text))
                        text += f"`{child.content}`"
                    elif child.type in ("softbreak", "hardbreak"):
                        text += " "
                    else:
                        text += child.content


def _pages(root):
    names = site_pages(root)
    linked = [name for name in EXTRA_LINKED if (root/name).is_file()]
    return {name: Page(root, name) for name in names + linked}, set(names)


def check_links(root=ROOT):
    """Relative links, images and anchors resolve; no directory links; site pages link to site pages."""
    pages, site = _pages(root)
    problems = []
    cache = dict(pages)
    for name, page in pages.items():
        for line, href, kind in page.links:
            where = f"{name}:{line}"
            url = urlsplit(href)
            if href.startswith(SELF_URL):
                target_name = unquote(href[len(SELF_URL):].split("#")[0])
                target = root/target_name
                fragment = url.fragment
            elif url.scheme or url.netloc:
                continue
            elif not url.path:
                target, target_name, fragment = page.path, name, url.fragment
            else:
                target = (page.path.parent/unquote(url.path)).resolve()
                if not target.is_relative_to(root.resolve()):
                    problems.append(f"{where}: {kind} leaves the repository: {href}")
                    continue
                target_name = target.relative_to(root.resolve()).as_posix()
                fragment = url.fragment
            if target.is_dir():
                problems.append(f"{where}: directory {kind}s do not work on the site; link a file: {href}")
                continue
            if not target.is_file():
                problems.append(f"{where}: missing {kind} target: {href}")
                continue
            if target.suffix == ".md":
                if name in site and target_name not in site:
                    problems.append(f"{where}: links a Markdown page that is not part of the site: {href}")
                    continue
                if fragment:
                    if target_name not in cache:
                        cache[target_name] = Page(root, target_name)
                    if unquote(fragment) not in cache[target_name].anchors:
                        problems.append(f"{where}: missing anchor #{fragment} in {target_name}")
            elif fragment and kind == "link" and target != page.path:
                problems.append(f"{where}: anchors only work on Markdown pages: {href}")
    return problems


def check_headings(root=ROOT):
    """Headings contain no mathematics (GitHub and the site would build different anchors)."""
    pages, _ = _pages(root)
    return [f"{name}:{line}: mathematics in a heading: {text}"
            for name, page in pages.items() for line, text in page.headings if "$" in text]


SECTIONING = re.compile(r"\\(?:chapter|section|subsection)\*?(?:\[[^\]]*\])?\{")
QUOTED_TITLE = re.compile(r"“([^”]+)”\s*\($")


def _balanced(text, start):
    """The brace group that opens at ``text[start-1]``; returns (content, end)."""
    depth, index = 1, start
    while depth and index < len(text):
        depth += {"{": 1, "}": -1}.get(text[index], 0)
        index += 1
    return text[start:index - 1], index


def _plain_title(title):
    """A LaTeX heading as plain text, or None when it holds mathematics or other macros."""
    title = re.sub(r"\\texorpdfstring\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}\{([^{}]*)\}", r"\1", title)
    title = title.replace("---", "—").replace("--", "–").replace("~", " ")
    return None if "\\" in title or "$" in title else " ".join(title.split())


def latex_titles(root=ROOT):
    """``\\label`` keys of the physics manual mapped to their heading (None for non-headings)."""
    titles = {}
    for path in sorted((root/"docs"/"manual").rglob("*.tex")):
        if "build" in path.relative_to(root/"docs"/"manual").parts:
            continue
        text = path.read_text()
        for label in LATEX_LABEL.findall(text):
            titles.setdefault(label, None)
        for match in SECTIONING.finditer(text):
            title, end = _balanced(text, match.end())
            label = re.match(r"\s*\\label\{([^}]+)\}", text[end:])
            if label:
                titles[label.group(1)] = _plain_title(title)
    return titles


def latex_labels(root=ROOT):
    """All ``\\label`` keys of the physics manual."""
    return set(latex_titles(root))


def check_labels(root=ROOT):
    """Physics-manual labels quoted in Markdown as inline code exist; a quoted title “...” matches."""
    pages, _ = _pages(root)
    titles = latex_titles(root)
    problems = []
    for name, page in pages.items():
        for line, content, before in page.code:
            if not LABEL.match(content):
                continue
            if content not in titles:
                problems.append(f"{name}:{line}: unknown physics-manual label `{content}`")
                continue
            quoted = QUOTED_TITLE.search(before)
            expected = titles[content]
            if quoted and expected is not None and " ".join(quoted.group(1).split()) != expected:
                problems.append(f"{name}:{line}: `{content}` is titled “{expected}”, not “{quoted.group(1)}”")
    return problems


def _tool(root, name):
    spec = importlib.util.spec_from_file_location(name, root/"tools"/f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_generated(root=ROOT):
    """Generated reference pages and LaTeX tables equal what tools/gen_docs.py produces."""
    generator = _tool(root, "gen_docs")
    problems = []
    for relative, content in generator.outputs().items():
        path = root/relative
        if not path.is_file() or path.read_text() != content:
            problems.append(f"{relative}: stale generated file; run python tools/gen_docs.py")
    return problems


PUBLIC_MODULES = ("qdjj_solver", "qdjj_solver.qp_solver", "qdjj_solver.dmrg_solver")
FIELD_STYLE = re.compile(r"\||:\w+:`")


def check_docstrings(root=ROOT):
    """Public objects have docstrings; field docstrings are plain text with ``literals`` only."""
    problems = []
    for module_name in PUBLIC_MODULES:
        module = importlib.import_module(module_name)
        for name in module.__all__:
            obj = getattr(module, name)
            if name.startswith("__") or not (inspect.isclass(obj) or inspect.isfunction(obj)):
                continue
            doc = obj.__doc__ or ""
            generated = dataclasses.is_dataclass(obj) and doc.startswith(f"{obj.__name__}(")
            if not doc.strip() or generated:
                problems.append(f"{module_name}.{name}: missing docstring")
    generator = _tool(root, "gen_docs")
    for _, _, relative, cls in generator.OPTION_CLASSES + generator.RESULT_CLASSES + generator.MODEL_CLASSES:
        for field in generator.dataclass_fields(root/"src"/"qdjj_solver"/relative, cls):
            where = f"src/qdjj_solver/{relative}:{field.line}"
            if not field.doc:
                problems.append(f"{where}: {cls}.{field.name} has no field docstring")
            elif FIELD_STYLE.search(field.doc):
                problems.append(f"{where}: {cls}.{field.name}: use plain text and ``literals`` (no | or roles)")
    for module_name in PUBLIC_MODULES:
        module = importlib.import_module(module_name)
        for name in module.__all__:
            obj = getattr(module, name)
            if not isinstance(obj, type):
                continue
            for member, value in vars(obj).items():
                if member.startswith("_"):
                    continue
                if (isinstance(value, (property, classmethod, staticmethod)) or inspect.isfunction(value)) \
                        and not inspect.getdoc(value):
                    problems.append(f"{module_name}.{name}.{member}: missing public member docstring")
    try:
        generator.parameter_entries()
    except ValueError as error:
        problems.append(str(error))
    return problems


TUTORIAL_HEADINGS = ("Goal", "Physics in brief", "Physics manual", "Run", "Expected results", "What to look at",
                     "Try next", "Files", "See also")
TUTORIAL_POINTER = re.compile(r"\\tutorial\{([^}]+)\}")


def _label_files(root):
    """Physics-manual label -> chapter file that defines it."""
    files = {}
    for path in sorted((root/"docs"/"manual"/"chapters").glob("*.tex")):
        for label in LATEX_LABEL.findall(path.read_text()):
            files[label] = path
    return files


def check_tutorials(root=ROOT):
    """Tutorials match the catalogue; READMEs follow the template; chapters point to their tutorial."""
    catalog_path = root/"examples"/"catalog.toml"
    if not catalog_path.is_file():
        return ["examples/catalog.toml: missing"]
    entries = tomllib.loads(catalog_path.read_text())["tutorial"]
    problems = []
    listed = [entry["dir"] for entry in entries]
    present = sorted(p.name for p in (root/"examples").iterdir() if p.is_dir() and re.match(r"\d\d_", p.name))
    if listed != present:
        problems.append(f"examples/catalog.toml lists {listed}, but the tutorial directories are {present}")
    files = _label_files(root)
    pointed = {}
    for path in (root/"docs"/"manual"/"chapters").glob("*.tex"):
        for name in TUTORIAL_POINTER.findall(path.read_text()):
            pointed.setdefault(name, set()).add(path)
    quoted_path = root/"docs"/"manual"/"generated"/"quoted.json"
    quoted = json.loads(quoted_path.read_text()) if quoted_path.is_file() else {}
    for entry in entries:
        name = entry["dir"]
        base = root/"examples"/name
        where = f"examples/{name}"
        for required in ("README.md", "run.py", "reference/quick.json", "reference/full.json"):
            if not (base/required).is_file():
                problems.append(f"{where}/{required}: missing")
        if any(not (base/required).is_file() for required in ("README.md", "reference/quick.json",
                                                              "reference/full.json")):
            continue
        headings = [text for line, text in Page(root, f"{where}/README.md").headings if not text.startswith("Tutorial")]
        order = [h for h in headings if h in TUTORIAL_HEADINGS]
        if order != list(TUTORIAL_HEADINGS):
            problems.append(f"{where}/README.md: sections must include, in order: {', '.join(TUTORIAL_HEADINGS)}")
        quick = json.loads((base/"reference"/"quick.json").read_text())["values"]
        full = json.loads((base/"reference"/"full.json").read_text())["values"]
        if not set(quick) <= set(full):
            problems.append(f"{where}/reference: quick keys {sorted(set(quick) - set(full))} are not in full.json")
        for page in entry.get("user", []):
            if not (root/page).is_file():
                problems.append(f"examples/catalog.toml, {name}: unknown user page {page}")
        primary = entry["physics"][0]
        if primary not in files:
            problems.append(f"examples/catalog.toml, {name}: unknown physics-manual label {primary}")
        elif files[primary] not in pointed.get(name, set()):
            problems.append(f"{files[primary].relative_to(root)}: add \\tutorial{{{name}}} (catalogue label {primary})")
        for shared in entry.get("shared", []):
            reference, manual = full.get(shared["reference"]), quoted.get(shared["manual"])
            if reference is None or manual is None or manual.get("value") is None:
                problems.append(f"examples/catalog.toml, {name}: shared quantity {shared} not found")
                continue
            tolerance = shared.get("rtol", 1e-6)*abs(manual["value"]) + shared.get("atol", 0.)
            if abs(reference["value"] - manual["value"]) > tolerance:
                problems.append(f"{where}: {shared['reference']}={reference['value']!r} disagrees with the physics "
                                f"manual's {shared['manual']}={manual['value']!r}")
    for name, paths in pointed.items():
        if name not in listed:
            problems.append(f"{sorted(paths)[0].relative_to(root)}: \\tutorial{{{name}}} is not in the catalogue")
    return problems


HARD_REFERENCE = re.compile(r"\b(?:Fig\.|Figure|Chapter|Section|Table|Eq\.)~?\s?\d+(?:\.\d+)+")
CODE = re.compile(r"\\code\{([^{}]*)\}")
CODE_ALLOWED = {"make", "python", "pip", "latexmk", "pdflatex", "biber", "True", "False", "None", "dict", "np"}


def source_identifiers(root=ROOT):
    """Every name, attribute, argument and identifier-like string constant of the package source."""
    names = set()
    for path in (root/"src").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Name):
                names.add(node.id)
            elif isinstance(node, ast.Attribute):
                names.add(node.attr)
            elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                names.add(node.name)
            elif isinstance(node, ast.arg):
                names.add(node.arg)
            elif isinstance(node, ast.keyword) and node.arg:
                names.add(node.arg)
            elif isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.isidentifier():
                names.add(node.value)
    return names


def check_latex(root=ROOT):
    """Physics-manual code matches the package; no hard-coded cross-reference numbers."""
    import numpy
    spec = importlib.util.spec_from_file_location("run_manual_listings", root/"tools"/"run_manual_listings.py")
    listings = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(listings)
    known = source_identifiers(root) | set(dir(numpy)) | set(dir(str)) | set(dir(list)) | set(dir(dict))
    found = listings.listings(root/"docs"/"manual")
    for _, _, code in found:   # names that the listings define themselves, such as model keywords
        try:
            known |= {node.arg for node in ast.walk(ast.parse(code)) if isinstance(node, ast.arg)}
        except SyntaxError:
            pass
    problems = []
    for path in sorted((root/"docs"/"manual"/"chapters").glob("*.tex")):
        relative = path.relative_to(root)
        for number, line in enumerate(path.read_text().split("\n"), 1):
            if line.lstrip().startswith("%"):
                continue
            for match in HARD_REFERENCE.finditer(line):
                problems.append(f"{relative}:{number}: hard-coded reference '{match.group(0)}'; use \\cref")
            for content in CODE.findall(line):
                if "/" in content or content.startswith("."):
                    if "/" in content and not ((root/content).exists() or (root/"docs"/"manual"/content).exists()):
                        problems.append(f"{relative}:{number}: \\code{{{content}}} is not a path of the repository")
                    continue
                leading = re.match(r"[A-Za-z_][\w.]*", content)
                if leading and not all(part in known or part in CODE_ALLOWED
                                       for part in leading.group(0).split(".") if part):
                    problems.append(f"{relative}:{number}: \\code{{{content}}} names nothing in the package")
    for path, line, code in found:
        where = f"{path.relative_to(root)}:{line}"
        try:
            tree = ast.parse(code)
        except SyntaxError as error:
            problems.append(f"{where}: listing does not parse: {error.msg}")
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("qdjj_solver"):
                module = importlib.import_module(node.module)
                problems += [f"{where}: {node.module} has no {alias.name}" for alias in node.names
                             if not hasattr(module, alias.name)]
            elif isinstance(node, ast.Attribute) and node.attr not in known:
                problems.append(f"{where}: attribute .{node.attr} does not occur in the package")
            elif isinstance(node, ast.keyword) and node.arg and node.arg not in known:
                problems.append(f"{where}: keyword {node.arg}= does not occur in the package")
    return problems


QUOTE = re.compile(r"\\q\{([^}]+)\}")


def check_quoted(root=ROOT):
    """The manual's quoted numbers are current, and every \\q key is defined and used."""
    scripts = root/"docs"/"manual"/"scripts"
    problems = []
    sys.path.insert(0, str(scripts))
    try:
        spec = importlib.util.spec_from_file_location("manual_quoted", scripts/"quoted.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        values = module.collect()
        tex, record = module.render(values)
    finally:
        sys.path.remove(str(scripts))
    generated = root/"docs"/"manual"/"generated"
    for name, content in (("quoted.tex", tex), ("quoted.json", record)):
        if not (generated/name).is_file() or (generated/name).read_text() != content:
            problems.append(f"docs/manual/generated/{name}: stale; run make quoted in docs/manual")
    used = {}
    for path in sorted((root/"docs"/"manual"/"chapters").glob("*.tex")):
        for number, line in enumerate(path.read_text().split("\n"), 1):
            for key in QUOTE.findall(line):
                used.setdefault(key, f"{path.relative_to(root)}:{number}")
    problems += [f"{where}: undefined quoted value \\q{{{key}}}" for key, where in used.items() if key not in values]
    problems += [f"docs/manual/scripts: quoted value {key} is never used" for key in sorted(set(values) - set(used))]
    return problems


CHECKS = {"links": check_links, "headings": check_headings, "labels": check_labels,
          "generated": check_generated, "docstrings": check_docstrings, "tutorials": check_tutorials,
          "quoted": check_quoted, "latex": check_latex}


def run(names=None, root=ROOT):
    problems = []
    for name in names or CHECKS:
        problems.extend(f"[{name}] {problem}" for problem in CHECKS[name](root))
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("checks", nargs="*", help=f"checks to run (default: all of {', '.join(CHECKS)})")
    arguments = parser.parse_args(argv)
    unknown = sorted(set(arguments.checks) - set(CHECKS))
    if unknown:
        parser.error(f"unknown checks: {', '.join(unknown)}")
    problems = run(arguments.checks or None)
    for problem in problems:
        print(problem)
    print(f"{len(problems)} documentation problem(s)" if problems else "documentation checks passed")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
