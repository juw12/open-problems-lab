"""Multi-route arXiv absence scan: windowed AFTER 2026-05-27 (Brada\u010d v1).

Question: does ANY submission dated after 2026-05-27 narrow the log-gap
between the s>=5 off-diagonal Ramsey lower exponent (2s-4, Brada\u010d) and the
upper exponent (s-2, AKS80 / LRZ01)?  A "beat" = a lower bound with a log
exponent < 2s-4, an upper bound with a log exponent > s-2, or an explicit
gap-closing statement.

Writes arxiv_hits.json (deduped entries) and scan_report.json (windowed + triage scaffold).
Single-threaded, well under 20 min.
"""
import json, re, sys, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from pathlib import Path

PROXY = "http://localhost:3128"
HERE = Path(__file__).parent
HITS = HERE / "arxiv_hits.json"
REPORT = HERE / "scan_report.json"
NS = {"a": "http://www.w3.org/2005/Atom", "arxiv": "http://arxiv.org/schemas/atom"}

CUTOFF = "2026-05-27"          # Brada\u010d v1 date (the frontier anchor)
WINDOW_START = "2026-05-28"    # post-cutoff window start
AUDIT = "2026-11-11"

# anchor we expect to (re)recover inside/near the window
ANCHOR_ID = "2605.28793"

# (label, search_query) -- each will be combined with the submittedDate window
QUERIES = [
    ("not-rsk",        'all:"r(s,k)"'),
    ("offdiag-term",   'all:"off-diagonal Ramsey"'),
    ("bradac-author",  'all:"Brada"'),
    ("ramsey-lb",      'ti:"Ramsey" AND abs:"lower bound"'),
    ("log-exponent",   'abs:Ramsey AND abs:"independent set" AND abs:"log"'),
    ("ks-offdiag",     'abs:"K_s" AND abs:Ramsey AND abs:"off-diagonal"'),
    ("aks-followup",   'all:"Alon" AND all:"Kahn" AND all:"Ramsey" AND all:"Szab"'),
    ("polylog",        'abs:Ramsey AND abs:"polylog" AND abs:"off-diagonal"'),
    ("erdos-szekeres", 'abs:"Erdos" AND abs:"Szekeres" AND abs:Ramsey AND abs:"lower bound"'),
]

def fetch(url, timeout=60, retries=3):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "millennium-lab/1.0 (absence scan)"})
            opener = urllib.request.build_opener(urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}))
            with opener.open(req, timeout=timeout) as r:
                return r.read()
        except Exception as ex:
            last = ex
            time.sleep(4 * (i + 1))
    raise last

def parse_entry(e):
    def t(tag):
        el = e.find("a:" + tag, NS)
        return (el.text or "").strip() if el is not None and el.text else ""
    authors = [a.find("a:name", NS).text.strip() for a in e.findall("a:author", NS)
               if a.find("a:name", NS) is not None and a.find("a:name", NS).text]
    prim = e.find("arxiv:primary_category", NS)
    return {
        "id": t("id"),
        "title": re.sub(r"\s+", " ", t("title")),
        "published": t("published"),
        "updated": t("updated"),
        "v1_date": t("published")[:10],
        "primary_category": prim.get("term") if prim is not None else "",
        "authors": authors,
        "summary": re.sub(r"\s+", " ", t("summary")),
    }

def run_query(label, sq, maxr=40):
    # windowed query: submittedDate after cutoff, before audit
    win = f"submittedDate:[{WINDOW_START.replace('-','')}0000 TO {AUDIT.replace('-','')}2359]"
    full = f"({sq}) AND {win}"
    q = urllib.parse.quote(full)
    url = f"http://export.arxiv.org/api/query?search_query={q}&start=0&max_results={maxr}&sortBy=submittedDate&sortOrder=ascending"
    out = {"label": label, "search_query": full, "n": 0, "entries": [], "error": None}
    try:
        raw = fetch(url).decode("utf-8", "replace")
        root = ET.fromstring(raw)
        for e in root.findall("a:entry", NS):
            d = parse_entry(e)
            # skip the arXiv "Error" pseudo-entry (no author / empty)
            if d["id"].startswith("http://arxiv.org/api/errors") or not d["title"]:
                continue
            d["matched_query"] = label
            out["entries"].append(d)
            out["n"] += 1
    except Exception as ex:
        out["error"] = f"{type(ex).__name__}: {ex}"
    return out

def main():
    results = []
    seen = {}
    for label, sq in QUERIES:
        r = run_query(label, sq)
        results.append({k: v for k, v in r.items() if k != "entries"})
        for d in r["entries"]:
            key = d["id"].split("/abs/")[-1].split("v")[0]  # bare id
            if key not in seen:
                d = dict(d); d["bare_id"] = key
                seen[key] = d
        time.sleep(3)  # arXiv API courtesy

    hits = list(seen.values())
    for d in hits:
        d["in_window_after_cutoff"] = d["v1_date"] > CUTOFF
        d["is_anchor"] = (d["bare_id"] == ANCHOR_ID)

    post = [d for d in hits if d["in_window_after_cutoff"]]
    post_nonanchor = [d for d in post if not d["is_anchor"]]

    HITS.write_text(json.dumps(hits, indent=2, ensure_ascii=False), encoding="utf-8")
    report = {
        "cutoff": CUTOFF,
        "window_start": WINDOW_START,
        "audit_date": AUDIT,
        "anchor_id": ANCHOR_ID,
        "n_queries": len(QUERIES),
        "queries": results,
        "n_total_hits": len(hits),
        "n_post_cutoff": len(post),
        "anchor_recovered": any(d["is_anchor"] for d in hits),
        "post_cutoff_nonanchor": [
            {"bare_id": d["bare_id"], "v1_date": d["v1_date"], "title": d["title"],
             "matched_query": d.get("matched_query"), "triage": ""}
            for d in post_nonanchor
        ],
        "verdict": "",
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({
        "n_queries": report["n_queries"],
        "n_total_hits": report["n_total_hits"],
        "n_post_cutoff": report["n_post_cutoff"],
        "anchor_recovered": report["anchor_recovered"],
        "post_cutoff_nonanchor": report["post_cutoff_nonanchor"],
    }, indent=2, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    main()
