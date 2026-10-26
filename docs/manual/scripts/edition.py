"""Write transient edition metadata for the physics PDF; no numerical calculation."""

import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[3]
label = os.environ.get("QDJJ_DOCS_LABEL")
if not label:
    if (ROOT/".git").exists():
        commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "--short=12", "HEAD"], text=True).strip()
        label = f"Development snapshot {commit}"
    else:
        label = "Development source edition"
label = re.sub(r"([\\{}_$%&#])", r"\\\1", label)
target = ROOT/"docs"/"manual"/"build"/"edition.tex"
target.parent.mkdir(parents=True, exist_ok=True)
content = "% Generated build identity; do not commit.\n\\newcommand{\\manualedition}{"+label+"}\n"
if not target.is_file() or target.read_text() != content:
    target.write_text(content)
