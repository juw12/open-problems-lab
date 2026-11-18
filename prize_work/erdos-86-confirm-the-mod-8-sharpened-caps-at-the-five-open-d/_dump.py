import json,io
t=json.load(io.open("cap_table.json",encoding="utf-8"))
w={r["n"]:r for r in t["wide_rows"]}
print("sharpened:",t["wide_sharpened"])
print("count:",len(t["wide_sharpened"]))
for n in t["wide_sharpened"]:
    if n>=27: print(n, w[n]["cap_plain"], w[n]["cap_quant"], w[n]["cap_plain"]-w[n]["cap_quant"])
for n in (100,150,200):
    print("anchor",n,w[n]["cap_plain"],w[n]["cap_quant"],w[n]["drop"])
