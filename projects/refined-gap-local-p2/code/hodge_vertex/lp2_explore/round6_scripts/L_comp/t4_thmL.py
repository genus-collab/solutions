# Theorem L via finite formula s.8.2, Schur side only. Usage: t4_thmL.py exact|mod hmax
from schur_side import *
import sys,time
mode=sys.argv[1]; hmax=int(sys.argv[2]); maxr=int(sys.argv[3]) if len(sys.argv)>3 else 3
P61=(1<<61)-1
K=Fld(None) if mode=="exact" else Fld(P61)
taus=[Fr(2,7),Fr(-5,3),Fr(13,9)]
def multisets(D,minp=2,maxr=9):
    """multisets of parts>=2 with sum(part-1)=D, at most maxr parts"""
    def rec(D,maxpart,left):
        if D==0:
            yield ()
            return
        if left==0: return
        for k in range(min(D+1,maxpart),1,-1):
            for r in rec(D-(k-1),k,left-1):
                yield (k,)+r
    return list(rec(D,10**6,maxr))
res=[]
for tau in taus:
    for h in range(2,hmax+1):
        for D in range(1,2*h-3+2+1):  # in-range D<=2h-3, plus D=2h-2, 2h-1 as out-of-range controls
            for R in multisets(D,maxr=maxr):
                r=len(R); N=3*h-3+r; J=sum(R)
                if J+N>26: continue
                t0=time.time()
                P=2*h-2+r+2*N+J+1
                F,items=compute_F(R,N,P,tau,K)
                val=Xi_from_F(F,items,R,h,tau,K)
                inr=D<=2*h-3
                ok=(val==0)==inr
                print(f"tau={tau} h={h} R={R} D={D} {'IN ' if inr else 'OUT'} Xi={'0' if val==0 else ('nonzero: '+(str(val) if mode=='exact' and len(str(val))<60 else '...'))} {'PASS' if ok else 'FAIL'} t={time.time()-t0:.1f}",flush=True)
                res.append((tau,h,R,D,inr,val==0,ok))
print("TOTAL",len(res),"in-range",sum(1 for x in res if x[4]),"out-range",sum(1 for x in res if not x[4]),"fails",sum(1 for x in res if not x[6]))
