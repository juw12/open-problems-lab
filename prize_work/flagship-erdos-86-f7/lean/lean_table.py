"""One table from the per-file results of lean_check.py:  python lean_table.py <dir> <out.md>"""
import json
import sys
from pathlib import Path

results = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(Path(sys.argv[1]).glob("*.json"))
           if p.name != "lean_check.json"]
rows = ["| File | Lean verdict | Theorems | Axioms used | Seconds |", "|---|---|---|---|---|"]
for r in results:
    axioms = ", ".join(r["axioms_used"]) or "none"
    verdict = r["verdict"] + (" (timed out)" if r["timed_out"] else "")
    if r["error_count"]:
        verdict += f", {r['error_count']} error(s)"
    rows.append(f"| `{r['file']}` | {verdict} | {r['theorems_with_axiom_report']} of {r['theorems']} | {axioms} | {r['seconds']} |")
all_checked = bool(results) and all(r["verdict"] == "checked" for r in results)
head = "# Lean check of this folder\n\n" + (
    "Every file below was accepted by Lean: no error, no `sorry`, and no axiom beyond the three "
    "standard ones (propext, Classical.choice, Quot.sound).\n"
    if all_checked else "Not every file is checked - see the table.\n")
Path(sys.argv[2]).write_text(head + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
Path(sys.argv[1], "lean_check.json").write_text(
    json.dumps({"all_checked": all_checked, "files": results}, indent=1), encoding="utf-8")
print(head + "\n" + "\n".join(rows))
