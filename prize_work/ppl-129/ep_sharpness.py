#!/usr/bin/env python3
"""EP lower-bound sharpness at the new anchors, BOTH definitions, script-verified.
  (a) card-literal  r_card(N) = f * loglog N / ((2/sqrt(e)) N sqrt(log N))   [card's '~1.2127 at 5000']
  (b) proper        r_prop(N) = f / EP_lower(N)
        where EP_lower(N) = (2/sqrt(e)) N sqrt(log N / loglog N)
        = f * sqrt(loglog N) / ((2/sqrt(e)) N sqrt(log N))          [-> 1  iff EP bound asymptotically tight]
"""
import json, os, math
D = os.path.dirname(os.path.abspath(__file__))
T = 2 / math.sqrt(math.e)
grid = {int(k): v for k, v in json.load(open(os.path.join(D, "fit_data_v2.json"))).items()}
push = json.load(open(os.path.join(D, "push_data.json")))
grid.update({int(k): v["f"] for k, v in push.items()})

rows = []
for n in [5000, 10000, 15000, 20000]:
    f = grid[n]
    L, LL = math.log(n), math.log(math.log(n))
    card = f * LL / (T * n * math.sqrt(L))
    prop = f * math.sqrt(LL) / (T * n * math.sqrt(L))
    rows.append(dict(n=n, f=f, ep_lower_asym=round(T * n * math.sqrt(L / LL), 1),
                     ratio_card=round(card, 6), ratio_proper=round(prop, 6)))
print(json.dumps(rows, indent=2))
json.dump(rows, open(os.path.join(D, "ep_sharpness.json"), "w"), indent=2)
