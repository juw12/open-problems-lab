"""Fetch live arXiv metadata for the s>=5 off-diagonal Ramsey anchors.

Anchor (the card's load-bearing paper):
  Brada\u010d 2026, arXiv:2605.28793, v1 2026-05-27
    lower bound: r(s,k) >= k^{s-1} / (log k)^{2s-4}   (s>=5)
Classical anchors (no arXiv, cited by page):
  Alon-Kahn-Szab\u00f3 1980 (AKS): r(s,k) << k^{s-1} / (log k)^{s-2}
  Li-Rousseau-Zang 2001 (LRZ01): sharpen AKS to <= (1+o(1)) * k^{s-1} / (log k)^{s-2}

Writes anchors.json. Runs single-threaded, well under 20 min.
"""
import json, re, sys, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from pathlib import Path

PROXY = "http://localhost:3128"
OUT = Path(__file__).parent / "anchors.json"
NS = {"a": "http://www.w3.org/2005/Atom",
      "arxiv": "http://arxiv.org/schemas/atom"}

# arXiv IDs to verify live (the card's anchors that have an arXiv home)
ANCHORS = {
    "bradac2026": {
        "id": "2605.28793",
        "who": "Brada\u010d 2026 (s>=5 lower bound)",
        "expected_v1": "2026-05-27",
        "claim": "r(s,k) >= k^{s-1} / (log k)^{2s-4}  for s>=5",
    },
}

def fetch(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "millennium-lab/1.0 (frontier audit)"})
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}))
    with opener.open(req, timeout=timeout) as r:
        return r.read()

def parse_entry(entry):
    def t(tag):
        el = entry.find("a:" + tag, NS)
        return (el.text or "").strip() if el is not None and el.text else ""
    authors = [a.find("a:name", NS).text.strip()
               for a in entry.findall("a:author", NS)
               if a.find("a:name", NS) is not None and a.find("a:name", NS).text]
    prim = entry.find("arxiv:primary_category", NS)
    return {
        "title": re.sub(r"\s+", " ", t("title")),
        "published": t("published"),   # v1 submission
        "updated": t("updated"),
        "primary_category": prim.get("term") if prim is not None else "",
        "authors": authors,
        "summary": re.sub(r"\s+", " ", t("summary")),
    }

def main():
    out = {}
    ok = True
    ids = ",".join(a["id"] for a in ANCHORS.values())
    url = f"http://export.arxiv.org/api/query?id_list={ids}&max_results=50"
    try:
        raw = fetch(url).decode("utf-8", "replace")
        root = ET.fromstring(raw)
        entries = {e.find("a:id", NS).text: e for e in root.findall("a:entry", NS)}
    except Exception as ex:
        out["error"] = f"arXiv fetch failed: {type(ex).__name__}: {ex}"
        OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        print(json.dumps(out, indent=2, ensure_ascii=False))
        sys.exit(2)

    for key, a in ANCHORS.items():
        base = f"arxiv.org/abs/{a['id']}"
        # entry id may carry a version suffix
        match = None
        for e in entries.values():
            eid = e.find("a:id", NS).text or ""
            if a["id"] in eid:
                match = e
                break
        if match is None:
            # could be a "no results" entry with an error tag
            out[key] = {"id": a["id"], "who": a["who"], "status": "NOT_FOUND_IN_API",
                        "expected_v1": a["expected_v1"], "claim": a["claim"]}
            ok = False
            continue
        d = parse_entry(match)
        d.update({"id": a["id"], "who": a["who"], "expected_v1": a["expected_v1"], "claim": a["claim"]})
        d["v1_date"] = d["published"][:10]
        d["date_ok"] = (d["v1_date"] == a["expected_v1"])
        out[key] = d

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(out, indent=2, ensure_ascii=False))
    print("ALL_ANCHORS_OK" if ok else "SOME_ANCHOR_MISSING", file=sys.stderr)
    sys.exit(0 if ok else 2)

if __name__ == "__main__":
    main()
