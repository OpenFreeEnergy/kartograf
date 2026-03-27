import os
import subprocess
import tempfile

from sybil import Example, Sybil
from sybil.parsers.rest import ClearNamespaceParser, CodeBlockParser, PythonCodeBlockParser


def ruff_fix(example: Example) -> str | None:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        f.write(example.parsed)
        tmpfile = f.name
    try:
        # Fix what can be fixed
        subprocess.run(["ruff", "check", "--fix", tmpfile], capture_output=True, text=True)

        with open(tmpfile) as f:
            fixed = f.read()

        # Write fixed code back to the rst file with correct indentation
        if fixed != example.parsed:
            source_path = str(example.path)
            with open(source_path) as f:
                source_lines = f.readlines()

            # Detect indentation from the first non-empty line of the block
            block_start = example.line  # 1-indexed
            indent = ""
            for line in source_lines[block_start:]:
                if line.strip():
                    indent = line[: len(line) - len(line.lstrip())]
                    break

            # Re-indent the fixed code
            fixed_indented = "".join(
                indent + line if line.strip() else "\n" for line in fixed.splitlines(keepends=True)
            )

            # Replace the original block in the file
            original_indented = "".join(
                indent + line if line.strip() else "\n" for line in example.parsed.splitlines(keepends=True)
            )
            with open(source_path, "w") as f:
                f.write("".join(source_lines).replace(original_indented, fixed_indented, 1))

        # Report any remaining unfixable errors, ignoring F821 (undefined names from namespace)
        result = subprocess.run(
            ["ruff", "check", "--ignore", "F821", tmpfile],
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
        CodeBlockParser(language="python", evaluator=ruff_fix),
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
