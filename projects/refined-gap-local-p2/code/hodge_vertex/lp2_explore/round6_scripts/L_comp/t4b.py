# exact s.8.2 at h=4,5,6 for the named multi-part R; three rational tau. usage: t4b.py num den
from schur_side import *
import sys,time
tau=Fr(int(sys.argv[1]),int(sys.argv[2]))
K=Fld(None)
cases=[(4,(2,)),(4,(2,2)),(4,(2,3)),(4,(2,2,2)),(4,(3,3,2)),(4,(2,2,2,2)),(4,(6,)),(4,(7,)),(4,(4,4)),(4,(3,3,3)),(4,(2,2,2,2,2)),(4,(2,2,2,2,2,2)),
       (5,(2,)),(5,(2,2)),(5,(2,3)),(5,(2,2,2)),(5,(3,3,3)),(5,(8,)),(5,(9,)),(5,(5,5)),(5,(4,4,3)),(5,(4,4,2)),(5,(3,3,3,2)),
       (6,(2,)),(6,(2,2)),(6,(2,3)),(6,(10,)),(6,(11,)),(6,(5,6)),(6,(6,6))]
bad=0
for h,R in cases:
    r=len(R);N=3*h-3+r;J=sum(R);D=J-r
    P=2*h-2+r+2*N+J+1
    t0=time.time()
    F,items=compute_F(R,N,P,tau,K)
    v=Xi_from_F(F,items,R,h,tau,K)
    inr=D<=2*h-3
    ok=(v==0)==inr
    bad+=not ok
    print(f"tau={tau} h={h} R={R} D={D} {'IN ' if inr else 'OUT'} Xi={v} {'PASS' if ok else 'FAIL'} t={time.time()-t0:.0f}",flush=True)
print("DONE fails",bad)
