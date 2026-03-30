# ruff: noqa
import os
import re
import subprocess
import sys
import tempfile

CODEBLOCK_RE = re.compile(
    r"(?P<directive>[ \t]*\.\. code-block:: python\n(?:[ \t]*:[^:]+:[^\n]*\n)*\n)"
    r"(?P<code>(?:[ \t]+[^\n]*\n|\n)+)"
)


def fix_block(code: str, indent: str) -> str:
    dedented = "".join(line.removeprefix(indent) for line in code.splitlines(keepends=True))
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        f.write(dedented)
        tmpfile = f.name
    try:
        subprocess.run(["ruff", "check", "--fix", "--ignore", "F821", tmpfile], capture_output=True)
        subprocess.run(["ruff", "format", tmpfile], capture_output=True)
        with open(tmpfile) as f:
            fixed = f.read()
    finally:
        os.unlink(tmpfile)
    return "".join(indent + line if line.strip() else "\n" for line in fixed.splitlines(keepends=True))


def process_file(path: str) -> bool:
    with open(path) as f:
        original = f.read()

    def replacer(m):
        code = m.group("code")
        # detect indentation from first non-empty line
        indent = ""
        for line in code.splitlines():
            if line.strip():
                indent = line[: len(line) - len(line.lstrip())]
                break
        fixed = fix_block(code, indent)
        return m.group("directive") + fixed

    fixed_content = CODEBLOCK_RE.sub(replacer, original)

    if fixed_content != original:
        with open(path, "w") as f:
            f.write(fixed_content)
        return True
    return False


if __name__ == "__main__":
    changed = [p for p in sys.argv[1:] if process_file(p)]
    if changed:
        print(f"Reformatted code blocks in: {', '.join(changed)}")
        sys.exit(1)  # pre-commit expects exit 1 when files are modified
