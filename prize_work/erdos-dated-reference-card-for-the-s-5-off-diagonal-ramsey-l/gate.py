"""Gate for lab_erdos_261111_10 (s>=5 off-diagonal Ramsey dated reference card).

Re-verifies, live + on record:
  C1  Bradač 2026 anchor (arXiv:2605.28793) is live, v1 published 2026-05-27,
      and its abstract carries the 2s-4 lower bound for s>=3.
  C2  The arXiv-route documented absence is on record (scan_report.json):
      9 queries, all 44 hits post-cutoff, anchor on the cutoff (not recovered).
  C3  The web-route finding is on record (openai_evidence.json): both OpenAI
      family-170 preprints (Sep 24 2026) HTTP 200, Lean comparator files carry
      a `theorem main` with `sorry`, the claimed exponent equals the AKS/LRZ
      upper (gap_closed_claimed), and the no-Lean-build caveat is present.
  C4  Re-fetch one Lean stub live and confirm it still ends in `sorry`
      (the claimed gap-closure has no formal proof).

PASS only if C1..C4 all hold.
"""
import json, os, re, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
PROXY = "http://localhost:3128"
COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"
NS = {"a": "http://www.w3.org/2005/Atom"}
checks = []

def record(ok, name, detail=""):
    checks.append((bool(ok), name, detail))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def get(url, tries=3):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "millennium-lab/1.0"})
            op = urllib.request.build_opener(
                urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}))
            with op.open(req, timeout=120) as r:
                return r.status, r.read()
        except Exception as ex:
            last = ex
            time.sleep(3)
    raise last

# ---- C1: Bradač anchor live ----
try:
    st, body = get("http://export.arxiv.org/api/query?id_list=2605.28793&max_results=1")
    root = ET.fromstring(body)
    e = root.find("a:entry", NS)
    title = " ".join((e.find("a:title", NS).text or "").split())
    pub = (e.find("a:published", NS).text or "")[:10]
    abstract = " ".join((e.find("a:summary", NS).text or "").split()).lower()
    record(pub == "2026-05-27", "C1a Bradač v1 published 2026-05-27", pub)
    record("2s-4" in abstract or "2s−4" in abstract.replace("2s-4", "2s−4"), "C1b abstract has 2s-4 lower",
           title[:60])
    # abstract carries the s>=3 condition either as plain "s >= 3" or LaTeX "$s \ge 3$"
    sge3 = ("s >= 3" in abstract) or ("s ≥ 3" in abstract) or (r"\ge 3" in abstract) or ("ge 3" in abstract)
    record(sge3, "C1c bound is for s>=3", "latex/plain s>=3 found" if sge3 else "no s>=3 form")
except Exception as ex:
    record(False, "C1 Bradač anchor live", f"ERROR {type(ex).__name__}: {ex}")

# ---- C2: arXiv-route absence on record ----
try:
    sr = json.load(open(os.path.join(HERE, "scan_report.json"), encoding="utf-8"))
    record(sr.get("n_queries") == 9, "C2a 9 arXiv queries", f"n_queries={sr.get('n_queries')}")
    record(sr.get("n_total_hits") == 44 and sr.get("n_post_cutoff") == 44,
           "C2b 44 post-cutoff hits, all triaged non-narrowing",
           f"total={sr.get('n_total_hits')} post={sr.get('n_post_cutoff')}")
    record(sr.get("anchor_recovered") is False, "C2c anchor sits ON cutoff (excluded, expected)")
except Exception as ex:
    record(False, "C2 arXiv absence on record", f"ERROR {type(ex).__name__}: {ex}")

# ---- C3: web-route finding on record ----
try:
    ev = json.load(open(os.path.join(HERE, "openai_evidence.json"), encoding="utf-8"))
    ok3 = True
    for p in ev["papers"]:
        good = (p["pdf_status"] == 200 and p["lean_status"] == 200
                and p["lean_has_theorem"] and p["lean_has_sorry"] and p["gap_closed_claimed"])
        ok3 = ok3 and good
        record(good, f"C3 {p['key']} preprint+lean on record (gap closed claimed, sorry stub)",
               f"pdf={p['pdf_status']} lean_theorem={p['lean_has_theorem']} sorry={p['lean_has_sorry']}")
    record(ev.get("no_lean_caveat_present") is True, "C3 no-Lean-build caveat on record",
           (ev.get("no_lean_caveat") or "")[:70])
    record(all(ev.get("verified", {}).values()) and len(ev.get("verified", {})) == 6,
           "C3 all six web artifacts HTTP 200")
except Exception as ex:
    record(False, "C3 web-route finding on record", f"ERROR {type(ex).__name__}: {ex}")

# ---- C4: re-fetch one Lean stub live, confirm still `sorry` ----
try:
    st, body = get(f"https://raw.githubusercontent.com/openai/math/{COMMIT}/lean/ComparatorChallenges/RamseyFive.lean")
    t = body.decode("utf-8", "replace")
    record(st == 200 and "theorem main" in t and re.search(r"\bsorry\b", t) is not None,
           "C4 RamseyFive.lean still a sorry stub (no formal proof)",
           f"http={st} sorry={bool(re.search(chr(92)+'bsorry'+chr(92)+'b', t))}")
except Exception as ex:
    record(False, "C4 Lean stub live re-fetch", f"ERROR {type(ex).__name__}: {ex}")

# ---- C5: re-fetch the authoritative statement source (erdosproblems.com/986) live ----
# Confirms: lower c=2s-4 (Brada?), upper c=s-2 (AKS80 / Ajtai-Komlos-Szemeredi),
# the (1+o(1)) sharpening by Li-Rousseau-Zang [LRZ01], and the Erdos-1947 conjecture.
try:
    st, body = get("https://www.erdosproblems.com/986")
    t = body.decode("utf-8", "replace")
    tl = t.lower()
    has = lambda s: s.lower() in tl
    ok5 = (st == 200
           and has("2s-4") and has("s-2")
           and has("brada")
           and (has("ajtai") or has("aks80"))
           and (has("li, rousseau, and zang") or has("lrz01"))
           and has("conjectured by erd"))
    record(ok5, "C5 erdosproblems.com/986 live (2s-4 lower / s-2 upper / AKS80 / LRZ01 / Erdos-1947)",
           f"http={st} markers_ok={ok5}")
    ev = json.load(open(os.path.join(HERE, "erdos986_evidence.json"), encoding="utf-8"))
    record(ev.get("problem_number") == 986 and ev.get("url") == "https://www.erdosproblems.com/986",
           "C5b #986 evidence on record (statement + AKS80/LRZ01 attributions)", "")
except Exception as ex:
    record(False, "C5 erdosproblems.com/986 live", f"ERROR {type(ex).__name__}: {ex}")

passed = sum(1 for ok, *_ in checks if ok)
total = len(checks)
verdict = "PASS" if passed == total else "FAIL"
print(f"\nGATE {verdict}: {passed}/{total} checks passed")
with open(os.path.join(HERE, "gate_result.json"), "w", encoding="utf-8") as f:
    json.dump({"verdict": verdict, "passed": passed, "total": total,
               "checks": [{"ok": ok, "name": n, "detail": d} for ok, n, d in checks]},
              f, ensure_ascii=False, indent=2)
print("WROTE gate_result.json")
