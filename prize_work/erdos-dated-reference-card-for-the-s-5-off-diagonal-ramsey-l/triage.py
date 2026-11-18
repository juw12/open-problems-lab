"""Triage helper: print title + abstract for a set of bare_ids from arxiv_hits.json."""
import json, sys
from pathlib import Path
HERE = Path(__file__).parent
hits = json.loads((HERE / "arxiv_hits.json").read_text(encoding="utf-8"))
byid = {d["bare_id"]: d for d in hits}
WANT = sys.argv[1:] or [
    "2610.06081", "2610.15293", "2610.17891",
    "2610.03358", "2610.14829", "2610.32644",
]
for bid in WANT:
    d = byid.get(bid)
    if not d:
        print(f"\n=== {bid}: NOT IN HITS ==="); continue
    print(f"\n=== {bid} ({d['v1_date']}) ===")
    print("TITLE:", d["title"])
    print("AUTHORS:", ", ".join(d["authors"][:4]))
    print("ABS:", d["summary"][:700])
    print()
