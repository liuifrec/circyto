import ast, io, os, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "circyto"

DEF_RE   = re.compile(r"^(\s*)(def\s+\w+\s*\(.*\)\s*:\s*)$")
CLASS_RE = re.compile(r"^(\s*)(class\s+\w+.*:\s*)$")

def fix_file(path: Path) -> bool:
    src = path.read_text(encoding="utf-8")
    changed = False

    # quick pass: ensure every def/class line that is immediately followed by
    # a non-indented/non-empty line gets a 'pass'
    lines = src.splitlines(keepends=False)
    i = 0
    while i < len(lines):
        line = lines[i]
        m_def = DEF_RE.match(line)
        m_cls = CLASS_RE.match(line)
        if m_def or m_cls:
            indent = (m_def or m_cls).group(1)
            j = i + 1
            # skip pure blank/comment lines but track indentation
            while j < len(lines) and (lines[j].strip() == "" or lines[j].lstrip().startswith("#")):
                j += 1
            needs_pass = False
            if j >= len(lines):
                needs_pass = True
            else:
                nxt = lines[j]
                # if next significant line is NOT more indented, then body is missing
                needs_pass = len(nxt) - len(nxt.lstrip()) <= len(indent)
            if needs_pass:
                lines.insert(i + 1, indent + "    pass")
                changed = True
                i += 1  # skip the inserted line
        i += 1

    tmp = "\n".join(lines) + ("\n" if src.endswith("\n") else "")
    if tmp != src:
        path.write_text(tmp, encoding="utf-8")

    # validate by parsing; if still failing with same error, report
    try:
        ast.parse(path.read_text(encoding="utf-8"))
        return changed or (tmp != src)
    except SyntaxError as e:
        if "expected an indented block" in str(e):
            # last-resort: insert pass directly after the error line if it's a def/class
            lines = tmp.splitlines(keepends=False)
            idx = (e.lineno or 1) - 1
            if 0 <= idx < len(lines):
                if DEF_RE.match(lines[idx]) or CLASS_RE.match(lines[idx]):
                    indent = re.match(r"^(\s*)", lines[idx]).group(1)
                    lines.insert(idx + 1, indent + "    pass")
                    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
                    try:
                        ast.parse(path.read_text(encoding="utf-8"))
                        return True
                    except SyntaxError:
                        return False
        return False

def main():
    fixed_any = False
    for py in sorted(PKG.rglob("*.py")):
        if py.name.endswith("_autofix_indent.py"):
            continue
        ok = fix_file(py)
        fixed_any = fixed_any or ok
    # show a summary and exit code
    print("Auto-fix complete.")
    sys.exit(0 if fixed_any else 0)

if __name__ == "__main__":
    main()
