"""Render the repository's GitHub Markdown conventions in Sphinx.

GitHub writes inline mathematics that contains escaped punctuation as
``$`...`$``; MyST's dollarmath would keep the backticks inside the formula.
Absolute links to files of this repository on GitHub (used where a page is
also shown on PyPI) are turned into relative links, so the site stays
self-contained. Both rewrites skip code fences and inline code spans and keep
line numbers unchanged.
"""

import posixpath
import re

REPOSITORY = "https://github.com/rokzitko/qdjj_solver/blob/main/"
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})")
INLINE = re.compile(r"(?P<code>`+).+?(?P=code)"
                    r"|\$`(?P<math>[^`\n]+)`\$"
                    r"|\]\(" + re.escape(REPOSITORY) + r"(?P<path>[^)\s#]+)(?P<anchor>#[^)\s]*)?\)")


def rewrite_line(line, directory):
    """Rewrite GitHub math spans and self-links in one line outside code fences."""
    def replace(match):
        if match.group("code"):
            return match.group(0)
        if match.group("math") is not None:
            return "$" + match.group("math") + "$"
        target = posixpath.relpath(match.group("path"), directory or ".")
        return "](" + target + (match.group("anchor") or "") + ")"
    return INLINE.sub(replace, line)


def rewrite(text, docname):
    """Apply :func:`rewrite_line` to every line that is not inside a fenced block."""
    directory = posixpath.dirname(docname)
    lines, fence = [], None
    for line in text.split("\n"):
        opening = FENCE.match(line)
        if fence is None:
            if opening:
                fence = opening.group(1)
                lines.append(line)
            else:
                lines.append(rewrite_line(line, directory))
        else:
            if opening and opening.group(1)[0] == fence[0] and len(opening.group(1)) >= len(fence) \
                    and not line.strip().strip(fence[0]):
                fence = None
            lines.append(line)
    return "\n".join(lines)


def _source_read(app, docname, source):
    source[0] = rewrite(source[0], docname)


def setup(app):
    app.connect("source-read", _source_read)
    return {"parallel_read_safe": True, "parallel_write_safe": True, "version": "1"}
