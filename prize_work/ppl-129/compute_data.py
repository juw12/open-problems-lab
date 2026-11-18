"""Resumable data generation for the fit. Reads/writes fit_data_v2.json, flushing after
each n so a server restart loses nothing. Recomputes all grid points (fresh file, no stale
off-by-one values)."""
import json, time, os
from fit_loglog import f_exact

grid = [100, 125, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1250, 1500, 1750,
        2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000]
F = "fit_data_v2.json"

data = {}
if os.path.exists(F):
    data = {int(k): v for k, v in json.load(open(F)).items()}

t0 = time.time()
for n in grid:
    if n in data:
        continue
    t1 = time.time()
    f = f_exact(n)
    data[n] = f
    with open(F, "w") as fh:
        json.dump({str(k): v for k, v in sorted(data.items())}, fh)
    done = sum(1 for x in grid if x in data)
    print(f"n={n:>5} f={f}  ({time.time()-t1:.1f}s)  [{done}/{len(grid)}]", flush=True)

if all(n in data for n in grid):
    print(f"ALL DONE in {time.time()-t0:.1f}s", flush=True)
else:
    print(f"INCOMPLETE: {sum(1 for x in grid if x in data)}/{len(grid)} (resume next run)", flush=True)
