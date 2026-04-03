from sybil import Sybil
from sybil.parsers.rest import PythonCodeBlockParser

pytest_collect_file = Sybil(
    name="tests",
    parsers=[
        PythonCodeBlockParser(),
    ],
    patterns=["*.py", "*.rst"],
).pytest()
