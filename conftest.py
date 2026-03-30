import os
import subprocess
import tempfile

from sybil import Example, Sybil
from sybil.parsers.rest import ClearNamespaceParser, CodeBlockParser, PythonCodeBlockParser


tests = Sybil(
    name="tests",
    parsers=[
        PythonCodeBlockParser(),
        ClearNamespaceParser(),
    ],
    patterns=["*.py", "*.rst"],
)

pytest_collect_file = (linting + tests).pytest()
