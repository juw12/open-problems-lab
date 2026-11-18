#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write lab_results/lab_erdos_261117_5.json (exact digits taken from the job's own JSON)."""
import json
import os
from datetime import datetime, timezone, timedelta

ROOT = r"E:\Projects\MillenniumResearch"
HERE = os.path.join(ROOT, "prize_work", "erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d")
d = json.load(open(os.path.join(HERE, "j261117_5_criterion.json"), encoding="utf-8"))

caps = {n: (v["cap_plain"], v["cap_quant"]) for n, v in d["open_caps"].items()}
sharp = d["sharpened_upto_NMAX"]
openb = d["open_dims_beyond_26"]
assert caps["12"] == (15701, 15700) and caps["24"] == (120795955, 120795954)
assert len(sharp) == 131 and len(openb) == 125 and sharp[:36][-1] == 198

res = {
 "status": "done",
 "result": (
  "The mod-8 sharpened caps are now a closed criterion, not just a table: with m=2^(n-1), E=nm, "
  "rho=min over admissible a of (E+a(a+2)m)/(2a+2) and S_pl=2*floor(rho/2), cap*(n)=cap(n)-1 exactly when "
  "n is even and S_pl=2 mod 4 (because N(a)=0 mod 8 for n>=4 makes Lemma 12's Delta=0 mod 8 read "
  "S=0 mod 4 for a even, S even for a=1 mod 4, and vacuous for a=3 mod 4). "
  "Checked against brute-force minimisation over all a at every n=3..600: 598/598 dims agree, 0 disagreements "
  "(22/22 gates, exact integers). The criterion gives 131 sharpened dims for n<=600, 125 of them open dims>26 "
  "(its first 36 coincide with the 36 the earlier job listed for n<=200) and explains both exceptions: n=25 is odd "
  "and at n=14 the binding a=2 gives S_pl=30036=0 mod 4. Card kill test reproduced by brute-force descent: "
  "S(12)=6824 with multiset {0:2,2:680,4:1366} (a second extremal multiset {2:686,4:1360,6:2} at a=2); the five "
  "sharpened caps 15700/1454898/6396312/27892120/120795954 against plain 15701/1454899/6396313/27892121/120795955. "
  "This is a corollary of published lemmas (DPTV79; Minamoto arXiv:2603.29127v5 Lemma 9/12), not past the frontier: "
  "the published upper bound at n=12 stays 15701 and no g(n) is determined."),
 "finished": datetime.now(timezone(timedelta(hours=7))).replace(microsecond=0).isoformat(),
 "files": [
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d.md",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j261117_5_criterion.py",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j261117_5_criterion.json",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j5_criterion.log",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j5_legA.log",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j5_legB.log",
  "prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j5_legC.log",
  "FINDINGS.md",
  "lab_results/lab_erdos_261117_5.json",
 ],
 "method": ("Exact-integer derivation of the per-a residue reduction of Lemma 12 (N(a)=0 mod 8 for n>=4), giving a "
            "one-line criterion for the sharpened set; verified dimension by dimension against a brute-force "
            "minimisation over all admissible a (downward scan, no CRT) at n=3..600, plus a fresh re-run of the "
            "three existing independent legs and a live re-fetch of the frontier."),
 "novelty": "reproduction",
 "baseline": ("Published: DPTV79 (J. Physique 40 (1979) 620) plain two-point cap; Minamoto arXiv:2603.29127v5 Lemma 9 "
              "(= two-point) and Lemma 12 (mod-8), applied in the paper only at n=6,7,8; erdosproblems.com/86 open, "
              "$100, live re-fetch 2026-11-17 (no proof claim on the page); paper scope re-confirmed live as Q6,Q7,Q8. "
              "Lab's own before this job: job 261114_3's five values and job 261117_2's 36-dim list for n<=200, which "
              "recorded that the sharpened set had no simple congruence description."),
 "metric": {
  "name": "sharpened dimensions covered by a closed criterion (cap* = cap - 1 iff n even and S_pl = 2 mod 4)",
  "before": "36 dimensions listed for n<=200 (table only), no description of the set",
  "after": "closed criterion, 598/598 dims n=3..600 agree with brute force; 131 sharpened dims <=600, 125 of them open dims>26",
 },
 "idea_verdict": "worked",
 "next_ideas": [
  {"title": "#86: does the lab's n=12 certificate already sit on the mod-8 boundary (Delta(2) = 16)?",
   "idea": ("The criterion makes the cap* target at n=12 a congruence: reaching 15700 needs S=6824=0 mod 4, i.e. "
            "Delta(2)=N(2)-6S=16 exactly, not the generic slack. Recompute S and Delta(2) from the committed n=12 "
            "certificate field and from every certified field at a published attained dim (gate: n=6, where gain=1 "
            "and g(6)=cap*=132), and report how far each is from Delta(a*)=0. If the n=12 cert already has Delta(2)=16 "
            "the bracket 15614<=g(12)<=15700 is locally rigid and the gauge-lift lane should aim at Delta(2)<16, not at "
            "more edges."),
   "why": ("The lab has a certified 15614-edge n=12 field and a running gauge-lift lane, but nobody has expressed the "
           "cap* target as a single integer identity on the field; this converts 'find more edges' into 'hit a residue'."),
   "verifier": ("Independent exact-integer recomputation of S and Delta(2) straight from the certificate JSON, "
                "reproducing the published S at n=6 and n=10 first, then disagreeing with no row of the 24-row table."),
   "kill": ("5 min: extract S from the n=12 certificate JSON; if the JSON does not determine the field (only a size), or "
            "if the recomputed Delta(2) is not 0 mod 8, the identity route cannot be applied and the card is dead."),
   "value": 5, "odds": 0.45, "cost": 1},
 ],
}
out = os.path.join(ROOT, "lab_results", "lab_erdos_261117_5.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1, ensure_ascii=False)
print("wrote", out)
print("finished", res["finished"])
print("caps", caps)
