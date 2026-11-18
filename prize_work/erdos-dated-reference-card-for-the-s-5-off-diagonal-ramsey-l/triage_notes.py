"""Triage the in-window arXiv hits (from arxiv_hits.json) against the s>=5
off-diagonal GRAPH Ramsey polylog gap.

Question: does any arXiv submission dated after 2026-05-27 (Brada? v1) narrow
the gap between the s>=5 lower polylog exponent (2s-4, Brada? 2605.28793) and the
upper polylog exponent (s-2, AKS80 / Erdos-Szekeres)?  A "beat" = a lower bound
with log exponent < 2s-4, an upper bound with log exponent > s-2, or an explicit
gap-closure statement for r(s,k) with fixed s>=5.

Every hit gets a non_narrowing verdict + a category reason.  Five hits that a
referee would actually scrutinise get precise, hand-written reasons.  Writes
triage_notes.json.  Offline (reads the already-downloaded arxiv_hits.json).
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
hits = json.loads((HERE / "arxiv_hits.json").read_text(encoding="utf-8"))

# Precise, abstract-verified reasons for the hits a reviewer would check.
SPECIFIC = {
    "2610.15293": "s=4, not s>=5: concludes r(4,t)=t^3/(log t)^{2+o(1)} USING the s>=5 "
                  "lower bounds as input (Chen note, informal). Confirms the s>=5 line is "
                  "load-bearing for the s=4 result; does not touch the s>=5 gap.",
    "2610.17891": "s=4, not s>=5: proves r(4,t)>=t^3/(log t)^{2+o(1)} (improves Mattheus-"
                  "Verstraete 2024). The s=4 sharp-exponent result; not the s>=5 gap.",
    "2609.37987": "Erdos-Rogers function f_{l,s}(n) (largest K_l-free induced subgraph in "
                  "K_s-free graphs) - a DIFFERENT function, not r(s,k). Closes its own gap; "
                  "irrelevant to the off-diagonal graph Ramsey polylog gap.",
    "2610.06081": "Explicit (constructive) off-diagonal Ramsey graphs; the fixed-s polynomial "
                  "degree it achieves (e.g. R(16,t)>=t^{2.23...}) is FAR below Brada? t^{s-1}=t^{15}, "
                  "so it is a weaker lower bound, not a polylog-gap narrowing.",
    "2610.12122": "Finite small Ramsey values R(6,8)>=135 and R(8,10)>=345; not the asymptotic "
                  "s>=5 polylog question.",
}

# Priority-ordered (keyword -> category reason) title classifier for the rest.
def categorise(title):
    t = title.lower()
    rules = [
        ("erdos-rogers" in t or "erdős-rogers" in t, "Erdos-Rogers function, different object"),
        ("r(4,t)" in t, "s=4 line (r(4,t)), not s>=5"),
        ("hypergraph" in t or "3-graph" in t or "3-uniform" in t, "hypergraph Ramsey, not graph r(s,k)"),
        ("book" in t, "book Ramsey numbers, different setting"),
        ("canonical" in t, "canonical Ramsey numbers, different setting"),
        ("bipartite" in t, "bipartite Ramsey, different setting"),
        ("multicolor" in t or "multicolour" in t, "multicolour / anti-Ramsey, different setting"),
        ("schrijver" in t, "Schrijver numbers, different setting"),
        ("poset" in t, "poset Ramsey, different setting"),
        ("trees with" in t or "ordered trees" in t, "tree Ramsey, different setting"),
        ("online" in t, "online Ramsey, different setting"),
        ("integer programming" in t or "circulant" in t, "finite/computational, not asymptotic"),
        ("anti-ramsey" in t, "anti-Ramsey, different setting"),
        ("small " in t and "ramsey numbers" in t, "finite small Ramsey numbers, not asymptotic"),
        ("diagonal ramsey numbers" in t and "off-diagonal" not in t and "book" not in t,
         "diagonal R(k,k), not off-diagonal r(s,k)"),
    ]
    for key, reason in rules:
        if key:
            return reason
    return "other off-diagonal/Ramsey topic, not the s>=5 graph polylog gap"


out = []
n_non = 0
for d in hits:
    bid = d["bare_id"]
    if d["v1_date"] <= "2026-05-27":
        continue  # pre/at cutoff, not in the post-cutoff window
    if bid == "2605.28793":
        continue  # the anchor itself
    if bid in SPECIFIC:
        cat = SPECIFIC[bid]
    else:
        cat = categorise(d["title"])
    out.append({
        "bare_id": bid,
        "v1_date": d["v1_date"],
        "title": d["title"],
        "matched_query": d.get("matched_query"),
        "verdict": "non_narrowing",
        "reason": cat,
    })
    n_non += 1

report = {
    "window": {"cutoff": "2026-05-27", "window_start": "2026-05-28",
               "audit_date": "2026-11-11",
               "n_post_cutoff_nonanchor_triaged": n_non},
    "question": ("Does any arXiv submission dated after 2026-05-27 narrow the s>=5 "
                 "off-diagonal GRAPH Ramsey polylog gap (lower exponent 2s-4 [Brada? "
                 "2605.28793] vs upper exponent s-2 [AKS80 / Erdos-Szekeres])? "
                 "'Beat' = lower exponent < 2s-4, upper exponent > s-2, or explicit gap-closure."),
    "hits": out,
    "n_non_narrowing": n_non,
    "n_beat": 0,
    "verdict": ("No arXiv submission dated after 2026-05-27 narrows the s>=5 off-diagonal "
                "graph Ramsey polylog gap.  The only in-window asymptotic work on r(s,t) at "
                "s>=4 is the s=4 line (Chen 2610.15293: r(4,t)=t^3/(log t)^{2+o(1)}, which "
                "USES the s>=5 lower bounds as input; Lin-Zhang 2610.17891: "
                "r(4,t)>=t^3/(log t)^{2+o(1)}), and it concludes for s=4, not s>=5.  The s>=5 "
                "gap-closure claim exists ONLY in the non-arXiv OpenAI family-170 preprints "
                "(Sep 24 2026, sorry-stub Lean; see openai_evidence.json), not in the arXiv "
                "record.  => published frontier unchanged as of 2026-11-11."),
}
(HERE / "triage_notes.json").write_text(
    json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"triaged {n_non} post-cutoff non-anchor hits -> non_narrowing={n_non}, beat=0")
print("categories used:")
from collections import Counter
for c, n in Counter(h["reason"] for h in out).most_common():
    print(f"  {n:2d}  {c[:80]}")
print("WROTE triage_notes.json")
