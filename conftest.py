from sybil import Sybil
from sybil.parsers.rest import ClearNamespaceParser, PythonCodeBlockParser

pytest_collect_file = Sybil(
    name="tests",
    parsers=[
        PythonCodeBlockParser(),
        ClearNamespaceParser(),
    ],
    patterns=["*.py", "*.rst"],
).pytest()
