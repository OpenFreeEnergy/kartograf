# ruff: noqa
import os
import re
import subprocess
import sys
import tempfile

CODEBLOCK_RE = re.compile(
    r"(?P<directive>[ \t]*\.\.[ \t]+(?:code-block|code|sourcecode)::[ \t]*python3?\n"
    r"(?:[ \t]*:[^:]+:[^\n]*\n)*\n)"
    r"(?P<code>(?:[ \t]+[^\n]*\n|\n)+)"
)


def check_ruff_available() -> None:
    result = subprocess.run(["ruff", "--version"], capture_output=True)
    if result.returncode != 0:
        print("Error: ruff not found on PATH", file=sys.stderr)
        sys.exit(2)


def fix_block(code: str, indent: str) -> str:
    dedented = "".join(line.removeprefix(indent) for line in code.splitlines(keepends=True))
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfile = os.path.join(tmpdir, "block.py")
        with open(tmpfile, "w", encoding="utf-8") as f:
            f.write(dedented)

        result = subprocess.run(
            ["ruff", "check", "--fix", "--ignore", "F821", tmpfile],
            capture_output=True,
        )
        if result.returncode not in (0, 1):  # ruff exits 1 on unfixable lint errors, which is ok
            raise RuntimeError(f"ruff check failed:\n{result.stderr.decode()}")

        result = subprocess.run(["ruff", "format", tmpfile], capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(f"ruff format failed:\n{result.stderr.decode()}")

        with open(tmpfile, encoding="utf-8") as f:
            fixed = f.read()

    reindented = "".join(indent + line if line != "\n" else "\n" for line in fixed.splitlines(keepends=True))

    # We need to add a SINGLE blank line after the code block, but ruff trims
    # the empty new line so we remove the one ruff adds, then we add one to the
    # line and then another to make a blank line
    return reindented.rstrip("\n") + "\n\n"


def process_file(path: str, check_only: bool = False) -> bool:
    with open(path, encoding="utf-8") as f:
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
        if not check_only:
            with open(path, "w", encoding="utf-8") as f:
                f.write(fixed_content)
        return True
    return False


if __name__ == "__main__":
    check_only = "--check" in sys.argv
    paths = [p for p in sys.argv[1:] if p != "--check"]

    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        print(f"Error: files not found: {', '.join(missing)}", file=sys.stderr)
        sys.exit(2)

    check_ruff_available()

    changed = [p for p in paths if process_file(p, check_only=check_only)]
    if changed:
        verb = "Would reformat" if check_only else "Reformatted"
        print(f"{verb} code blocks in: {', '.join(changed)}")
        sys.exit(1)  # pre-commit expects exit 1 when files are modified
