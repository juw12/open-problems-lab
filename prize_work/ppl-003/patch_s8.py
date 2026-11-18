#!/usr/bin/env python3
# patch_s8.py — correct the |S|=8 class split in the N=5e6 artifact.
# The scan's underlying counts are exact; only the summary LABEL was wrong
# ("2pqrstuvw" vs the correct convention "2pqrstuv" = 2 + 7 letters, matching job1).
# has-2 = total - no-2 is an exact identity (every |S|=8 combo is 2-class or no-2),
# independently confirmed by diag_2class.py (2-class max growth 57.99, 2.18% ~ 989.7K).
import json
f = "E:/Projects/MillenniumResearch/prize_work/ppl-003/scan_s8_N5000000k8_8.json"
d = json.load(open(f))
no2 = d["new_s8_by_class"]["pqrstuvw"]
total = d["new_s8_total"]
has2 = total - no2
d["new_s8_by_class"] = {"pqrstuvw": no2, "2pqrstuv": has2}
d["new_s8_class_note"] = ("has-2 (2pqrstuv) = total - no-2, exact identity; "
                           "independently confirmed by diag_2class.py (2-class max growth 57.99, "
                           "2.181% of C(45,7) ~ 989.7K divergent). Label key corrected from the "
                           "buggy '2pqrstuvw' to '2pqrstuv' (2 + 7 letters).")
json.dump(d, open(f, "w"), indent=2)
print("CORRECTED |S|=8 class split @N=5e6:")
print("  total   =", total)
print("  pqrstuvw(no-2) =", no2)
print("  2pqrstuv(has-2) =", has2)
print("  sanity no2+has2 == total:", no2 + has2 == total)
print("  frac no-2  = %.5f  (of C(45,8)=215,553,195)" % (no2 / 215553195))
print("  frac has-2 = %.5f  (of C(45,7)=45,379,620)" % (has2 / 45379620))
