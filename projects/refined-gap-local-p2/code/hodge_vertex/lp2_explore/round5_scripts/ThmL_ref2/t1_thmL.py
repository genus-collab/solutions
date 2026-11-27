# Theorem L directly on the Hodge side: all T in range, at genus h, several tau, plus outside-range controls.
import sys,time
from hodge import *
hmax=int(sys.argv[1]); rmax=int(sys.argv[2])
taus=[Fr(2,7),Fr(-5,3),Fr(-1,2),Fr(13,9)]
def Ts(D):
    # multisets of parts>=2 with Sum(t-1)=D
    for p in partitions(D):
        yield tuple(sorted(x+1 for x in p))
t0=time.time()
for h in range(2,hmax+1):
    nz=0;bad=[];ctrl=0;ctrln=0;skipped=[]
    for D in range(1,2*h-3+2+1):
        for T in Ts(D):
            if len(T)>rmax: skipped.append(T); continue
            for tau in taus:
                v=X(h,T,tau)
                if D<=2*h-3:
                    nz+=1
                    if v!=0: bad.append((T,tau,v))
                else:
                    ctrln+=1
                    if v!=0: ctrl+=1
            # with units
            if D<=2*h-3 and len(T)+1<=rmax:
                vv=X(h,T+(1,),taus[0])
                if vv!=0: bad.append((T+(1,),taus[0],vv))
                nz+=1
    print("h",h,"in-range tests",nz,"nonzero:",bad[:5],"outside-range nonzero",ctrl,"/",ctrln,"skipped r>rmax",len(skipped),"time",round(time.time()-t0,1),flush=True)
    # empty, part 0 identity, part-1 string
    t=taus[0]
    print("  Xi(empty)",X(h,(),t),"Xi(0)/( (2-2h) Xi(empty))",X(h,(0,),t)/((2-2*h)*X(h,(),t)),flush=True)
    print("  string test X(3,1)/X(3) (3 = J) :",X(h,(3,1),t)/X(h,(3,),t) if X(h,(3,),t)!=0 else ("X(3)=0, X(3,1)=",X(h,(3,1),t)),flush=True)
