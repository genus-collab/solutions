# Cor L1 by MV-free route: Xi_h(R)=Sum_k z^k Xi^(k), z=-1/sigma, from Mumford+DVV. Theorem L directly on Hodge side as polys in z.
# also degree<=h check, Cor L2, part 0, Cor L3.
from hodge import *
from t4_lists import multisets
import sys,time
hmax=int(sys.argv[1]); maxr=int(sys.argv[2])
fails=0;count=0;ctrl=0
for h in range(2,hmax+1):
    for D in range(1,2*h-3+3):
        for R in multisets(D,maxr=maxr):
            if h>=4 and len(R)>maxr-1: continue
            t0=time.time()
            p=Xi_poly(h,R)
            inr=D<=2*h-3
            deg=len(p)-1
            okdeg=deg<=h
            ok=(not p)==inr and okdeg
            count+=1
            if not inr: ctrl+=1
            print(f"h={h} R={R} D={D} {'IN ' if inr else 'OUT'} Xi^(k)={[str(c) for c in p] if len(str(p))<400 else 'nonzero(deg %d)'%deg} {'PASS' if ok else 'FAIL'} t={time.time()-t0:.1f}",flush=True)
            if not ok: fails+=1
print("TOTAL",count,"controls",ctrl,"fails",fails)
