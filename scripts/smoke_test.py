#!/usr/bin/env python
"""Smoke tests for the scVelo + CellRank2 fate-mapping skill.

These tests do not require a real AnnData object or heavy dependencies. They
verify skill packaging, script syntax, core placeholders, and failure guards.
"""

from __future__ import annotations

import argparse
import json
import py_compile
import tempfile
import re
import sys
from pathlib import Path


def assert_true(cond, msg):
    if not cond:
        raise AssertionError(msg)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skill-root", default=".")
    args = parser.parse_args()

    root = Path(args.skill_root).resolve()
    skill = root / "SKILL.md"
    assert_true(skill.exists(), "SKILL.md missing")
    text = skill.read_text(encoding="utf-8")
    assert_true(text.startswith("---\n"), "SKILL.md missing YAML frontmatter")
    front = text.split("---", 2)[1]
    assert_true("name:" in front and "description:" in front, "frontmatter must include name and description")
    assert_true((root / "references").exists(), "references/ directory missing")
    assert_true(not (root / "README.md").exists(), "README.md should not be present in final skill package")

    scripts = sorted((root / "scripts").glob("*.py"))
    assert_true(scripts, "scripts missing")
    for script in scripts:
        # Compile into a temporary .pyc so the smoke test also works when the
        # skill directory is read-only or mounted without __pycache__ writes.
        with tempfile.TemporaryDirectory() as td:
            cfile = Path(td) / f"{script.stem}.pyc"
            py_compile.compile(str(script), cfile=str(cfile), doraise=True)
        stext = script.read_text(encoding="utf-8")
        if script.name not in {"smoke_test.py"}:
            assert_true("NotImplementedError" not in stext, f"{script.name} still contains NotImplementedError")

    template = (root / "templates" / "report_template.md").read_text(encoding="utf-8")
    placeholders = re.findall(r"\{[^{}\n]+\}", template)
    # Template may contain placeholders, but report generator must remove them.
    report_script = (root / "scripts" / "06_generate_report.py").read_text(encoding="utf-8")
    assert_true("re.sub" in report_script and "unfilled_placeholders" in report_script, "report generator lacks placeholder guard")

    result = {
        "status": "PASS",
        "checked_scripts": [str(p.relative_to(root)) for p in scripts],
        "template_placeholders_detected": placeholders,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
