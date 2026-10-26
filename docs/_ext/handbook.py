"""Print presentation of canonical Sphinx nodes; no duplicated technical descriptions."""

from pathlib import Path, PurePosixPath
import posixpath
import re
from urllib.parse import quote, urlsplit, unquote

from docutils import nodes
from sphinx import addnodes
from sphinx.transforms.post_transforms import SphinxPostTransform


def source_url(ref, path, anchor=""):
    base = "https://github.com/rokzitko/qdjj_solver/blob/"
    return base+quote(ref, safe="")+"/"+quote(path, safe="/")+("#"+quote(anchor, safe="-") if anchor else "")


def book_documents(env, start="docs/handbook/index"):
    """Explicit assembled document set, independent of the website's broader tree."""
    found, pending = set(), [start]
    while pending:
        name = pending.pop()
        if name in found:
            continue
        found.add(name)
        pending.extend(env.toctree_includes.get(name, []))
    return found


def node_document(node, fallback, root):
    while node is not None:
        if node.source:
            path = Path(node.source)
            if path.suffix == ".md" and path.is_relative_to(root):
                return path.relative_to(root).with_suffix("").as_posix()
        if isinstance(node, addnodes.start_of_file):
            return node.get("docname", fallback)
        node = node.parent
    return fallback


class ExternalBookLinks(SphinxPostTransform):
    """Resolve links to documents outside this book before the domain resolver."""
    default_priority = 1
    formats = ("latex",)

    def run(self, **kwargs):
        included = book_documents(self.env)
        for node in list(self.document.findall(addnodes.pending_xref)):
            target = node.get("reftarget", "")
            if node.get("refdomain") == "doc" and target in self.env.found_docs and target not in included:
                node.replace_self(nodes.reference("", "", *node.children,
                    refuri=source_url(self.config.handbook_source_ref, target+".md", node.get("reftargetid") or "")))
                continue
            url = urlsplit(node.get("reftarget", ""))
            if url.scheme or url.netloc or not url.path.endswith(".md"):
                continue
            origin = node.get("refdoc", self.env.docname)
            path = unquote(url.path)
            path = path.lstrip("/") if path.startswith("/") else posixpath.normpath(
                posixpath.join(posixpath.dirname(origin), path))
            if path[:-3] in self.env.found_docs and path[:-3] not in included:
                node.replace_self(nodes.reference("", "", *node.children,
                    refuri=source_url(self.config.handbook_source_ref, path, url.fragment)))


def prepare(app, doctree, docname):
    if app.builder.format != "latex":
        return
    included = book_documents(app.env)
    # External sources and data remain real versioned hyperlinks in the offline
    # book, not HTML-only _downloads paths or nonexistent PDF destinations.
    for node in list(doctree.findall(addnodes.download_reference)):
        target = node.get("reftarget", "")
        origin = node.get("refdoc", docname)
        path = (PurePosixPath(origin).parent/target).as_posix()
        path = posixpath.normpath(path).lstrip("/")
        node.replace_self(nodes.reference("", "", *node.children,
                                         refuri=source_url(app.config.handbook_source_ref, path)))
    # The complete two-orbital JSON example fits on one page. Keep its frame
    # together: splitting it after the input-field table can overfill the first
    # fragment when framed.sty adds the continuation label.
    for block in list(doctree.findall(nodes.literal_block)):
        if node_document(block, docname, Path(app.srcdir)) == "docs/interfaces" and block.get("language") == "json":
            parent, index = block.parent, block.parent.index(block)
            parent.insert(index, nodes.raw("", r"\begingroup\sphinxverbatimwithminipagetrue", format="latex"))
            parent.insert(index + 2, nodes.raw("", r"\endgroup", format="latex"))
    # Long pipe tables are readable as repeated, labelled entries in print. The
    # HTML/GitHub tables and these entries have precisely the same cell content.
    for table in list(doctree.findall(nodes.table)):
        origin = node_document(table, docname, Path(app.srcdir))
        if origin.startswith("docs/reference/") or origin == "docs/troubleshooting":
            group = next(table.findall(nodes.tgroup), None)
            if group is None:
                continue
            head = next(group.findall(nodes.thead), None)
            body = next(group.findall(nodes.tbody), None)
            if head is None or body is None:
                continue
            headings = [cell.astext() for cell in head.findall(nodes.entry)]
            listing = nodes.definition_list()
            for row in body.children:
                cells = list(row.children)
                if not cells:
                    continue
                item = nodes.definition_list_item()
                term = nodes.term("", cells[0].astext())
                definition = nodes.definition()
                for label, cell in zip(headings[1:], cells[1:], strict=True):
                    para = nodes.paragraph()
                    para += nodes.strong("", label+": ")
                    for child in cell.children:
                        if isinstance(child, nodes.paragraph):
                            para.extend(child.children)
                        else:
                            para += child
                    definition += para
                item += term
                item += definition
                listing += item
            table.replace_self(listing)
        else:
            table["classes"].append("longtable")
    # Resolved out-of-book references are supplied by the missing-reference hook;
    # the assertion here catches an accidental inclusion of a research chapter.
    if docname == "docs/handbook/index":
        app.env.handbook_documents = sorted(included)
        section = next(doctree.findall(nodes.section))
        section.insert(1, nodes.paragraph("", f"Software version {app.config.release}. "
                                          f"{app.config.handbook_build_label}."))
        for start in doctree.findall(addnodes.start_of_file):
            if start.get("docname", "").startswith("docs/reference/"):
                start.parent.insert(start.parent.index(start), nodes.raw("", r"\appendix", format="latex"))
                break
        for section in doctree.findall(nodes.section):
            origin = node_document(section, docname, Path(app.srcdir))
            if not origin.startswith("docs/reference/") or not section.children:
                continue
            title = section.children[0].astext()
            if re.fullmatch(r"[a-z][a-z0-9_]*", title):
                identity = f"handbook-parameter-{app.env.new_serialno('handbook-index')}"
                section.insert(1, nodes.target("", "", ids=[identity]))
                section.insert(2, addnodes.index(entries=[("single", title, identity, "", None)]))


def missing_reference(app, env, node, contnode):
    if app.builder.format != "latex":
        return None
    target = node.get("reftarget", "")
    if node.get("refdomain") == "std" and node.get("reftype") == "doc":
        name = posixpath.normpath(posixpath.join(posixpath.dirname(node.get("refdoc", "")), target))
        if target.startswith("/"):
            name = target.lstrip("/")
        if name not in book_documents(env) and name in env.found_docs:
            return nodes.reference("", "", contnode,
                                   refuri=source_url(app.config.handbook_source_ref, name+".md"))
    return None


def setup(app):
    app.add_config_value("handbook_source_ref", "main", "env")
    app.add_config_value("handbook_build_label", "Development snapshot", "env")
    app.connect("doctree-resolved", prepare)
    app.connect("missing-reference", missing_reference)
    app.add_post_transform(ExternalBookLinks)
    return {"parallel_read_safe": True, "parallel_write_safe": True, "version": "1"}
