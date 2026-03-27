import subprocess
import tempfile
import os

from sybil import Sybil, Example
from sybil.parsers.rest import ClearNamespaceParser, PythonCodeBlockParser, CodeBlockParser


def ruff_check(example: Example) -> str | None:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        f.write(example.parsed)
        tmpfile = f.name
    try:
        result = subprocess.run(
            ["ruff", "check", tmpfile],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return result.stdout.replace(tmpfile, "<example>")
    finally:
        os.unlink(tmpfile)

linting = Sybil(
    name="linting",
    parsers=[
        CodeBlockParser(language="python", evaluator=ruff_check),
    ],
    patterns=["*.py", "*.rst"],
)


tests = Sybil(
    name="tests",
    parsers=[
        PythonCodeBlockParser(),
        ClearNamespaceParser(),
    ],
    patterns=["*.py", "*.rst"],
)

pytest_collect_file = (linting + tests).pytest()
