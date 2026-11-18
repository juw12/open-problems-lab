"""Capture the OpenAI 'family 170' (Sep 24 2026) post-cutoff evidence, reproducible.

Re-fetches (via proxy) the two preprint PDFs' URLs (existence), the two Lean
comparator files (content -> theorem + sorry finding), the family scope page, and
the latent-seminar review caveat. Writes openai_evidence.json.

This is the WEB-ROUTE finding that the arXiv-only scan could not see.
"""
import json, os, re, time, urllib.request

PROXY = "http://localhost:3128"
HERE = os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(HERE, "openai_preprints")
COMMIT = "adc7f1241b42e322a6451854ab7e4b4c146bf78a"

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

raw = f"https://raw.githubusercontent.com/openai/math/{COMMIT}"
PAPERS = [
    dict(key="r5t",
         title="The sharp logarithmic exponent of r(5,t)",
         pdf_rel="preprints/The-Sharp-Logarithmic-Exponent-of-r-5-t-September-24-2026/paper.pdf",
         lean_rel="lean/ComparatorChallenges/RamseyFive.lean",
         claim="r(5,t)=t^4/(log t)^{3+o(1)};  lim (4 log t - log r(5,t))/log log t = 3",
         claimed_exponent=3, s="5", upper_exponent=3, bradac_exponent=6),
    dict(key="fixed_s",
         title="Sharp Logarithmic Exponents for Fixed Off-Diagonal Ramsey Numbers",
         pdf_rel="preprints/Sharp-Logarithmic-Exponents-for-Fixed-Off-Diagonal-Ramsey-Numbers-September-24-2026/paper.pdf",
         lean_rel="lean/ComparatorChallenges/SharpLogRamsey.lean",
         claim="for every fixed s>=6, r(s,t)=t^{s-1}/(log t)^{s-2+o(1)}; lim ((s-1) log t - log r(s,t))/log log t = s-2",
         claimed_exponent="s-2", s=">=6", upper_exponent="s-2", bradac_exponent="2s-4"),
]

out = {
    "audited_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "family": "OpenAI 'math' collection, family 170 (release dated Oct 6 2026; manuscripts dated Sep 24 2026)",
    "commit": COMMIT,
    "papers": [],
    "verified": {},
}

for p in PAPERS:
    rec = dict(p)
    # PDF existence
    st_pdf, pdf_bytes = get(raw + "/" + p["pdf_rel"])
    rec["pdf_url"] = raw + "/" + p["pdf_rel"]
    rec["pdf_status"] = st_pdf
    rec["pdf_bytes"] = len(pdf_bytes)
    out["verified"][f"{p['key']}_pdf_http200"] = (st_pdf == 200)
    # local copy already downloaded
    lp = os.path.join(PRE, p["key"] + ".pdf")
    rec["pdf_local"] = ("openai_preprints/" + p["key"] + ".pdf") if os.path.exists(lp) else None
    # Lean comparator file
    st_lean, lean_bytes = get(raw + "/" + p["lean_rel"])
    lean_text = lean_bytes.decode("utf-8", "replace")
    rec["lean_url"] = raw + "/" + p["lean_rel"]
    rec["lean_status"] = st_lean
    out["verified"][f"{p['key']}_lean_http200"] = (st_lean == 200)
    rec["lean_has_theorem"] = bool(re.search(r"\btheorem\s+main\b", lean_text))
    rec["lean_has_sorry"] = bool(re.search(r"\bsorry\b", lean_text))
    rec["lean_n_sorry"] = len(re.findall(r"\bsorry\b", lean_text))
    # extract the theorem line(s)
    thm = re.findall(r"theorem\s+main[^\n]*", lean_text)
    rec["lean_statement"] = " ".join(x.strip() for x in thm)
    rec["gap_closed_claimed"] = (str(p["claimed_exponent"]) == str(p["upper_exponent"]))
    out["papers"].append(rec)

# family scope page
st_sc, sc_bytes = get(raw + "/lean/docs/170.md")
out["verified"]["family_scope_http200"] = (st_sc == 200)
sc = sc_bytes.decode("utf-8", "replace")
out["family_scope_page"] = raw + "/lean/docs/170.md"
out["family_scope_notes_r5t"] = ("r(5,t)" in sc or "r(5,t)" in sc.replace("$", "").replace("{", "").replace("}", ""))
out["family_scope_notes_s6"] = ("s\\ge6" in sc or "s>=6" in sc)

# latent-seminar review caveat
LS = "https://lowrank.github.io/latent-seminar/OpenAI-Math/170-off-diagonal-ramsey"
st_ls, ls_bytes = get(LS)
ls = ls_bytes.decode("utf-8", "replace")
out["verified"]["latent_seminar_http200"] = (st_ls == 200)
out["latent_seminar_url"] = LS
m = re.search(r"No Lean build[^<]*", ls)
out["no_lean_caveat"] = m.group(0).strip() if m else None
out["no_lean_caveat_present"] = bool(m)
m2 = re.search(r"family 170 of the October 6, 2026 OpenAI collection", ls)
out["latent_seminar_confirms_family170"] = bool(m2)

with open(os.path.join(HERE, "openai_evidence.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print("papers:", [(r["key"], r["pdf_status"], "lean_theorem" if r["lean_has_theorem"] else "?",
                   "SORRY" if r["lean_has_sorry"] else "no-sorry", "gap_closed_claimed" if r["gap_closed_claimed"] else "?") for r in out["papers"]])
print("no_lean_caveat_present:", out["no_lean_caveat_present"])
print("verified:", out["verified"])
print("WROTE openai_evidence.json")
