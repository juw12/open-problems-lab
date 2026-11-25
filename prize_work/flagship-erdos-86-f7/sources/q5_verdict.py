"""Reads the JSON files written by q5_inequality_cases.py and states one verdict on (*).

(*): every C4-free subgraph H of Q5 has sum_v g(d(v)) <= 800, g = (0, 12, 19, 21, 29, 45).

  FALSE        some run found a graph with sum >= 801 (the graph is re-checked here from its edge list).
  HOLDS        for at least one solver, the cases with target 801 that came back INFEASIBLE cover every
               C4-free subgraph of Q5, and the control (target 800 must be FEASIBLE in d4 and in d5) passed.
               This is a solver verdict: it is as good as the solver, it is not a formal proof.
  NOT DECIDED  anything else.

A solver covers everything when `all` is INFEASIBLE, or when d4 and d5 are both covered; d4 is
covered by the run `d4` or by the four runs d4 top=1..4, d5 by the run `d5` or by d5 top=1..5.
(Largest degree <= 3 needs no run: the sum is at most 32 * 21 = 672.)

Usage: python q5_verdict.py <folder with the json files> <out.json> [out.md]
"""
import json
import sys
from itertools import combinations
from pathlib import Path

G = [0, 12, 19, 21, 29, 45]
N = 5
TOPS = {"d4": [1, 2, 3, 4], "d5": [1, 2, 3, 4, 5]}


def graph_sum(edge_list):
    """C4-free subgraph of Q5?  Returns the sum of g over its degrees, or raises."""
    chosen = {(min(a, b), max(a, b)) for a, b in edge_list}
    assert all(0 <= a < b < 32 and bin(a ^ b).count("1") == 1 for a, b in chosen), "not edges of Q5"
    for i, j in combinations(range(N), 2):
        for v in range(32):
            if not v & (1 << i) and not v & (1 << j):
                square = [(v, v | (1 << i)), (v | (1 << i), v | (1 << i) | (1 << j)),
                          (v | (1 << j), v | (1 << i) | (1 << j)), (v, v | (1 << j))]
                assert sum(edge in chosen for edge in square) <= 3, "4-cycle"
    degrees = [sum((min(v, v ^ (1 << i)), max(v, v ^ (1 << i))) in chosen for i in range(N)) for v in range(32)]
    return sum(G[d] for d in degrees)


def main():
    folder, out = Path(sys.argv[1]), Path(sys.argv[2])
    runs = []
    for path in sorted(folder.rglob("*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except ValueError:
            continue
        if "verdict" in report and "case" in report:
            report["file"] = path.name
            runs.append(report)

    counterexamples = []
    for run in runs:
        if run["verdict"] == "FEASIBLE" and "graph" in run:
            total = graph_sum(run["graph"]["edge_list"])
            assert total >= run["target"], f"{run['file']}: graph sum {total} below target {run['target']}"
            if total >= 801:
                counterexamples.append({"file": run["file"], "sum": total, "edge_list": run["graph"]["edge_list"]})

    def infeasible(solver, case, top=None):
        return any(run["solver"] == solver and run["case"] == case and run.get("top_degree") == top
                   and run["target"] == 801 and run["verdict"] == "INFEASIBLE" for run in runs)

    def covered(solver, case):
        return infeasible(solver, case) or all(infeasible(solver, case, top) for top in TOPS[case])

    solvers = sorted({run["solver"] for run in runs})
    holds_by = [solver for solver in solvers
                if infeasible(solver, "all") or (covered(solver, "d4") and covered(solver, "d5"))]
    control = {case: sorted({run["solver"] for run in runs if run["case"] == case and run["target"] == 800
                             and run.get("top_degree") is None and run["verdict"] == "FEASIBLE"})
               for case in ("d4", "d5")}
    control_passed = all(control.values())

    if counterexamples:
        verdict = "FALSE"
    elif holds_by and control_passed:
        verdict = "HOLDS"
    else:
        verdict = "NOT DECIDED"

    table = [{"file": run["file"], "solver": run["solver"], "case": run["case"], "top_degree": run.get("top_degree"),
              "target": run["target"], "plain": run.get("plain", False), "verdict": run["verdict"],
              "seconds": run["seconds"], "time_limit": run["time_limit"]} for run in runs]
    summary = {"claim": "sum_v g(d(v)) <= 800 for every C4-free subgraph of Q5, g = (0, 12, 19, 21, 29, 45)",
               "verdict": verdict, "holds_by_solver": holds_by, "control_feasible_at_800": control,
               "control_passed": control_passed, "counterexamples": counterexamples, "runs": table}
    out.write_text(json.dumps(summary, indent=1), encoding="utf-8")

    lines = [f"# Q5 inequality: {verdict}", "",
             f"Solvers that cover every case at target 801: {', '.join(holds_by) or 'none'}.  "
             f"Control (800 is reachable): d4 by {', '.join(control['d4']) or 'nobody'}, "
             f"d5 by {', '.join(control['d5']) or 'nobody'}.", "",
             "| solver | case | top | target | verdict | seconds | limit |", "|---|---|---|---|---|---|---|"]
    for row in table:
        lines.append(f"| {row['solver']} | {row['case']}{' plain' if row['plain'] else ''} | {row['top_degree'] or ''} | "
                     f"{row['target']} | {row['verdict']} | {row['seconds']} | {row['time_limit']:.0f} |")
    text = "\n".join(lines) + "\n"
    if len(sys.argv) > 3:
        Path(sys.argv[3]).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
