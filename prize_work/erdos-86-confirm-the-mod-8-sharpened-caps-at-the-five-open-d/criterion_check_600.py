# Criterion check extended: does  cap* = cap-1  iff (n even and S_pl = 2 mod 4)  hold for n=3..600?
# Independent of the main script's closed form: S_max by naive downward scan.
import sys
def NV(n): return 2**(n-1)
def NE(n): return n*2**(n-1)
def smax(n,a):
    W=2*a+2; N=NE(n)+a*(a+2)*NV(n); S=N//W
    while S>=0:
        if (N-W*S)%8==0: return S
        S-=1
def caps(n):
    E=NE(n); st=0 if n%2==0 else 1
    # plain: min over a of N/(2a+2) via cross multiplication
    best=None
    for a in range(st,n+1,2):
        p,q=E+a*(a+2)*NV(n),2*a+2
        if best is None or p*best[1]<best[0]*q: best=(p,q)
    Spl=2*((best[0]//best[1])//2)
    cap=(E+Spl)//2
    SQ=min(smax(n,a)-(smax(n,a)%2) for a in range(st,n+1,2))
    return cap,(E+SQ)//2,Spl
bad=[];sharp=[]
for n in range(3,601):
    cap,cq,Spl=caps(n)
    pred=(n%2==0) and (Spl%4==2)
    act=(cap-cq==1)
    if pred!=act: bad.append((n,cap,cq,Spl%4))
    if act: sharp.append(n)
print("rows 598, mismatches:",bad)
print("sharpened count:",len(sharp))
print("sharpened open (>26):",len([n for n in sharp if n>26]))
