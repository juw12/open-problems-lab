# -*- coding: utf-8 -*-
import json
o = json.load(open("sing_exp_N10000000.json"))
print("TOP KEYS:", sorted(o.keys()))
print("\n=== gates ===")
for k in ["gate_mean_phiN","gate_a_mass12","gate_b_growth12","gate_c_gamma12","all_gates_pass","n_divergent_job35","n_strong"]:
    print(" ", k, ":", o.get(k))
print("\n=== local-dim spectrum (60-core) ===")
ld = o.get("localdim_spectrum", {})
for k,v in ld.items():
    if k != "top10_lowest_p":
        print("  ", k, ":", v)
print("  top10_lowest_p:", ld.get("top10_lowest_p"))
print("\n=== p_local overall / by_absS ===")
print("  overall:", o.get("p_local_overall"))
print("  by_absS:", o.get("p_local_by_absS"))
print("\n=== mass spectrum ===")
ms = o.get("mass_spectrum", {})
for k,v in ms.items():
    if k != "top12":
        print("  ", k, ":", v)
print("  top12 (c, class, M_w0, p):", ms.get("top12"))
print("\n=== verifier_N2 (5e6) ===")
print(" ", o.get("verifier_N2"))
print("\n=== rows_N1e7 sample (first 3) ===")
for r in o.get("rows_N1e7", [])[:3]:
    print("  ", r)
print("  n rows N1e7:", len(o.get("rows_N1e7",[])), " n rows N5e6:", len(o.get("rows_N5e6",[])))

print("\n\n=== JOB-35 scan_N10000000.json structure ===")
j35 = json.load(open("scan_N10000000.json"))
print("TOP KEYS:", sorted(j35.keys())[:30])
dl = j35.get("divergent_list", [])
print("n divergent:", len(dl))
print("sample rows:")
for r in dl[:3]:
    print("  ", r)
# how many strong (growth>=20)
strong = [r for r in dl if r.get("growth",0) >= 20.0]
print("n strong (growth>=20):", len(strong))
# top-20 by growth
top20 = sorted(dl, key=lambda x: -x.get("growth",0))[:20]
print("\nTOP-20 by growth (c, class, |S|, growth):")
for r in top20:
    print("   c=%s class=%s |S|=%d growth=%.2f" % (r.get("c"), r.get("class"), len(r.get("S",[])), r.get("growth",0)))
