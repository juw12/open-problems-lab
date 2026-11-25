"""Check every Lean file of this folder with Lean itself and write one verdict.

    python lean_check.py <lake project dir> <lean dir> <out dir> [file.lean ...]

Each file is copied, `#print axioms` is added for every theorem in it, and the copy is run with
`lake env lean`.  A file counts as checked when Lean exits with 0, prints no error, uses no
`sorry`, and every theorem depends only on the three standard axioms of Lean with Mathlib
(propext, Classical.choice, Quot.sound).  Anything else is written down as it is.
"""
import json
import re
import subprocess
import sys
import time
from pathlib import Path

STANDARD = {"propext", "Classical.choice", "Quot.sound"}
TIMEOUT_S = 150 * 60


def audit_copy(text):
    namespace = re.search(r"^namespace\s+(\S+)", text, re.M)
    prefix = namespace.group(1) + "." if namespace else ""
    names = re.findall(r"^theorem\s+([^\s({\[:]+)", text, re.M)
    lines = [text.rstrip("\n"), "", "-- added by lean_check.py"]
    lines += [f"#print axioms {prefix}{name}" for name in names]
    return "\n".join(lines) + "\n", [prefix + name for name in names]


def check(project, lean_file, out_dir):
    text = lean_file.read_text(encoding="utf-8")
    audit_text, names = audit_copy(text)
    audit_file = out_dir / (lean_file.stem + "__audit.lean")
    audit_file.write_text(audit_text, encoding="utf-8")
    started = time.time()
    try:
        run = subprocess.run(
            ["lake", "env", "lean", str(audit_file.resolve())],
            cwd=project, capture_output=True, text=True, timeout=TIMEOUT_S,
        )
        output, exit_code, timed_out = run.stdout + run.stderr, run.returncode, False
    except subprocess.TimeoutExpired as stop:
        output = (stop.stdout or b"").decode("utf-8", "replace") if isinstance(stop.stdout, bytes) else (stop.stdout or "")
        exit_code, timed_out = 124, True
    seconds = round(time.time() - started, 1)
    (out_dir / (lean_file.stem + ".out.txt")).write_text(output, encoding="utf-8")

    axioms = {}
    # a theorem name may itself end in a prime (adj_flip'), so the name runs to the last quote
    for name, rest in re.findall(r"^'(.+)' depends on axioms: \[([^\]]*)\]", output, re.M):
        axioms[name] = [a.strip() for a in rest.replace("\n", " ").split(",") if a.strip()]
    for name in re.findall(r"^'(.+)' does not depend on any axioms", output, re.M):
        axioms[name] = []
    errors = [line for line in output.splitlines() if re.search(r":\d+:\d+: error", line)]
    uses_sorry = "sorryAx" in output or "declaration uses 'sorry'" in output
    extra = sorted({a for used in axioms.values() for a in used} - STANDARD)
    missing = [name for name in names if name not in axioms]
    checked = (exit_code == 0 and not errors and not uses_sorry and not extra and not missing)
    return {
        "file": lean_file.name,
        "verdict": "checked" if checked else "failed",
        "exit": exit_code,
        "timed_out": timed_out,
        "seconds": seconds,
        "theorems": len(names),
        "theorems_with_axiom_report": len(names) - len(missing),
        "axioms_used": sorted({a for used in axioms.values() for a in used}),
        "axioms_beyond_the_standard_three": extra,
        "uses_sorry": uses_sorry,
        "errors": errors[:40],
        "error_count": len(errors),
    }


def main():
    project, lean_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = sys.argv[4:]
    files = [lean_dir / name for name in wanted] if wanted else sorted(lean_dir.glob("*.lean"))
    results = [check(project, lean_file, out_dir) for lean_file in files]
    for result in results:
        (out_dir / (Path(result["file"]).stem + ".json")).write_text(
            json.dumps(result, indent=1), encoding="utf-8")
        print(json.dumps(result, indent=1))
    return 0 if all(r["verdict"] == "checked" for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
