"""Publication regressions that do not require TeX, a network, or numerical scans."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("sphinx")
pypdf = pytest.importorskip("pypdf")
pytestmark = pytest.mark.docs
ROOT = Path(__file__).resolve().parents[1]


def module(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT/relative)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def small_pdf(path, text="Software version 0.2.0 Required chapter", title="User's Guide"):
    from pypdf.generic import DecodedStreamObject, DictionaryObject, NameObject
    writer = pypdf.PdfWriter()
    page = writer.add_blank_page(width=600, height=800)
    font = DictionaryObject({NameObject("/Type"): NameObject("/Font"),
                             NameObject("/Subtype"): NameObject("/Type1"),
                             NameObject("/BaseFont"): NameObject("/Helvetica")})
    page[NameObject("/Resources")] = DictionaryObject({NameObject("/Font"): DictionaryObject({NameObject("/F1"): font})})
    stream = DecodedStreamObject()
    stream.set_data(f"BT /F1 12 Tf 40 700 Td ({text}) Tj ET".encode("ascii"))
    page[NameObject("/Contents")] = stream
    writer.add_metadata({"/Title": title})
    writer.write(path)


def test_actual_pdf_identity_version_and_contents(tmp_path):
    checker = module("tools/check_manual.py")
    path = tmp_path/"book.pdf"
    assert checker.artifact_problems(path)
    small_pdf(path)
    assert checker.artifact_problems(path, title="User's Guide", version="0.2.0", pages=1,
                                     required=("Required chapter",)) == []
    for kwargs in ({"title": "physics manual"}, {"version": "99.0"}, {"pages": 2},
                   {"required": ("absent chapter",)}):
        assert checker.artifact_problems(path, **kwargs)
    path.write_bytes(b"not a PDF")
    assert "unreadable PDF" in checker.artifact_problems(path)[0]


def test_success_log_does_not_substitute_for_pdf(tmp_path):
    checker = module("tools/check_manual.py")
    (tmp_path/"book.log").write_text("Output written on book.pdf (500 pages).\n")
    assert checker.problems(tmp_path, "book") == ([], 500)
    assert checker.main(["--build-dir", str(tmp_path), "--basename", "book"]) == 1


def test_missing_glyph_is_a_publication_error(tmp_path):
    checker = module("tools/check_manual.py")
    (tmp_path/"book.log").write_text("Missing character: There is no x in font\nOutput written on book.pdf (1 page).\n")
    assert checker.problems(tmp_path, "book")[0] == ["Missing character: There is no x in font"]


def test_external_document_link_retains_revision_and_anchor():
    from docutils import nodes
    from docutils.utils import new_document
    from sphinx import addnodes
    extension = module("docs/_ext/handbook.py")
    config = SimpleNamespace(handbook_source_ref="v0.2.0")
    env = SimpleNamespace(config=config, docname="docs/models", found_docs={"docs/models", "docs/validation"},
                          toctree_includes={"docs/handbook/index": ["docs/models"]})
    document = new_document("fixture")
    document.settings.env = env
    external = addnodes.pending_xref(refdomain="doc", reftype="myst", refdoc="docs/models",
                                     reftarget="docs/validation", reftargetid="independent-physics-checks")
    external += nodes.literal("", "validation")
    document += external
    internal = addnodes.pending_xref(refdomain="doc", reftype="myst", refdoc="docs/models", reftarget="docs/models")
    document += internal
    extension.ExternalBookLinks(document).run()
    assert document[0]["refuri"].endswith("/v0.2.0/docs/validation.md#independent-physics-checks")
    assert document[1] is internal


def test_reference_table_print_adaptation_preserves_cells():
    from docutils import nodes
    from docutils.utils import new_document
    extension = module("docs/_ext/handbook.py")
    tree = new_document("fixture")
    table = nodes.table()
    table.source = str(ROOT/"docs/reference/options.md")
    group = nodes.tgroup(cols=2)
    head, body = nodes.thead(), nodes.tbody()
    for container, values in ((head, ("Option", "Meaning")), (body, ("chi_max", "Maximum MPS bond dimension"))):
        row = nodes.row()
        for text in values:
            row += nodes.entry("", nodes.paragraph("", text))
        container += row
    group += head
    group += body
    table += group
    tree += table
    app = SimpleNamespace(builder=SimpleNamespace(format="latex"), srcdir=ROOT,
                          env=SimpleNamespace(toctree_includes={}), config=SimpleNamespace(handbook_source_ref="v0.2.0"))
    extension.prepare(app, tree, "docs/reference/options")
    assert not list(tree.findall(nodes.table))
    assert "chi_max" in tree.astext() and "Maximum MPS bond dimension" in tree.astext()


def test_release_manuals_require_complete_inputs(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT/"tools"))
    builder = module("tools/build_manuals.py")
    with pytest.raises(SystemExit) as error:
        builder.main(["--manual", "user", "--release-tag", "v0.2.0"])
    assert error.value.code == 2
    with pytest.raises(ValueError, match="release tag"):
        builder.verify_versions("v999.0")


def test_autodoc_is_parsed_as_python_domain_objects(tmp_path):
    """A warning-free build can still print '.. py:function::' as plain text."""
    import io
    from sphinx.application import Sphinx
    from sphinx import addnodes
    generator = module("tools/gen_docs.py")
    source = tmp_path/"source"
    source.mkdir()
    (source/"conf.py").write_text("extensions=['myst_parser','sphinx.ext.autodoc']\n"
                                 "root_doc='index'\nautodoc_preserve_defaults=True\n")
    (source/"index.md").write_text(generator.api_page())
    for name in ("options", "results"):
        (source/f"{name}.md").write_text(f"---\norphan: true\n---\n\n# {name.title()}\n")
    warning = io.StringIO()
    app = Sphinx(str(source), str(source), str(tmp_path/"html"), str(tmp_path/"doctrees"),
                 "html", status=io.StringIO(), warning=warning, freshenv=True)
    app.build()
    assert not warning.getvalue(), warning.getvalue()
    tree = app.env.get_doctree("index")
    signatures = list(tree.findall(addnodes.desc_signature))
    assert any("qdjj_solver.solve" in node.get("ids", []) for node in signatures)
    assert any("qdjj_solver.DiscreteBath.record" in node.get("ids", []) for node in signatures)
    assert ".. py:function::" not in tree.astext()
