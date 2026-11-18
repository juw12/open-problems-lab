#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Erdos #86 -- job lab_erdos_261117_2: run the three legs and cross-check them.

  leg A  legA_gate_slack_parity.py  verbatim copy of job lab_erdos_261114_3's gate
                                    (its G1 also imports tools/gate_odsquare_cap.py,
                                     the committed gate tool -- a fourth code path)
  leg B  fresh_mod8_table.py        this job's CRT / closed-form descent
  leg C  audit_mod8_caps.py         this job's independent audit (cross-multiplication
                                    minimum, downward-scan S_max, complete support
                                    search, and a decode of the real n=12 field)

Each leg recomputes the two-point programme and Lemma 12 from its own code.  This
driver only compares the JSON they write; it contains no arithmetic of its own.
Writes crosscheck.json and prints the three tables side by side.

Exit 0 iff every leg exits 0 and all three agree on every compared number.
"""
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEGS = [("A", "legA_gate_slack_parity.py", "slack_table.json"),
        ("B", "fresh_mod8_table.py", "cap_table.json"),
        ("C", "audit_mod8_caps.py", "audit_caps.json")]
OPEN_N = [12, 14, 18, 20, 22, 24, 25]
PUB_EXACT = {3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328,
             11: 7296, 13: 33792, 15: 153600, 16: 327680, 17: 688128,
             19: 3047424, 21: 13369344, 23: 58195968, 26: 520093696}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rows_of(tag, data):
    """Normalise each leg's JSON to {n: (cap_plain, cap_quant, S_quant)}."""
    out = {}
    if tag == "A":
        for r in data["rows"]:
            out[int(r["n"])] = (int(r["cap_plain"]), int(r["cap_quant"]), int(r["S_quant"]))
    elif tag == "B":
        for r in data["rows"]:
            out[int(r["n"])] = (int(r["cap_plain"]), int(r["cap_quant"]), int(r["S_quant"]))
    else:
        for k, r in data["rows"].items():
            out[int(k)] = (int(r["cap_plain"]), int(r["cap_quant"]), int(r["S_quant"]))
    return out


def main():
    print("job lab_erdos_261117_2 -- three-leg cross-check of the mod-8 caps\n")
    rc, tables, meta = {}, {}, {}
    for tag, script, out in LEGS:
        path = os.path.join(HERE, script)
        print("=== leg %s: %s  (sha256 %s)" % (tag, script, sha256(path)[:16]))
        p = subprocess.run([sys.executable, script], cwd=HERE, capture_output=True, text=True)
        rc[tag] = p.returncode
        log = os.path.join(HERE, "leg%s_run.log" % tag)
        with open(log, "w", encoding="utf-8") as f:
            f.write(p.stdout + "\n--- stderr ---\n" + p.stderr)
        tail = [l for l in p.stdout.strip().splitlines() if l.startswith(("GATES", "FAILURES"))]
        print("    rc=%d  %s" % (p.returncode, " | ".join(tail)))
        if p.returncode != 0:
            print(p.stdout[-2000:])
        with open(os.path.join(HERE, out), encoding="utf-8") as f:
            tables[tag] = rows_of(tag, json.load(f))
        meta[tag] = dict(script=script, sha256=sha256(path), rc=p.returncode, log=log)

    ns = sorted(set().union(*[set(t) for t in tables.values()]))
    mismatch = {}
    for n in ns:
        vals = {tag: tables[tag].get(n) for tag in tables}
        if len(set(vals.values())) != 1:
            mismatch[n] = vals

    print("\n--- the three legs on the seven open dims (cap_plain / cap_quant / S_quant) ---")
    print("| n | A | B | C | agree |")
    print("|---|---|---|---|---|")
    for n in OPEN_N:
        v = [tables[t][n] for t in "ABC"]
        print("| %d | %s | %s | %s | %s |" % (n, v[0], v[1], v[2],
                                              "yes" if len(set(v)) == 1 else "NO"))

    print("\n--- published exact rows: all three legs give cap_quant >= g, equal at the attained dims ---")
    gate = {}
    for tag in "ABC":
        below = [n for n, g in PUB_EXACT.items() if tables[tag][n][1] < g]
        att = [n for n, g in PUB_EXACT.items() if tables[tag][n][1] == g]
        gate[tag] = dict(below=below, attained=sorted(att), n_attained=len(att))
        print("leg %s: below=%s  attained=%d dims %s" % (tag, below or "none",
                                                         len(att), sorted(att)))

    ok = all(v == 0 for v in rc.values()) and not mismatch \
        and all(gate[t] == gate["A"] for t in "ABC") \
        and len(gate["A"]["attained"]) == 16
    print("\nVERDICT: %s" % ("all three legs agree on every compared number"
                             if ok else "MISMATCH -- see crosscheck.json"))

    out = dict(job="lab_erdos_261117_2", legs=meta, rc=rc, mismatch=mismatch,
               open_caps={n: dict(zip("ABC", [tables[t][n] for t in "ABC"])) for n in OPEN_N},
               published_row_gate=gate, verdict=("agree" if ok else "mismatch"))
    with open(os.path.join(HERE, "crosscheck.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
