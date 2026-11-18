import sys
sys.path.insert(0,".")
from fresh_mod8_table import cap_plain, cap_quant
sharp=[]
rows=[]
for n in range(3,201):
    cp,ap,sp=cap_plain(n); cq,aq,SQ,W,N,D=cap_quant(n)
    rows.append((n,cp-cq,W,N-s_rem if False else D))
    if cp-cq: sharp.append(n)
print("sharpened up to 200:",sharp)
print("count",len(sharp))
