from gsum import *
import time,sys
samples=[((Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7))),((Fr(2),Fr(-3),Fr(8)),(Fr(1),Fr(4))),((Fr(1,3),Fr(5,2),Fr(-3)),(Fr(11,5),Fr(-2,7)))]
mode=sys.argv[1]
if mode=='legs':
    Dmax=2
    for al,ep in samples:
        c=setup(al,ep,Dmax,6); gs=GS(c)
        for g in (1,2):
            for a in (0,1,2):
                gt=time.time()
                lhs=gs.total(g,leg=a)
                gw=[loc(al,ep,g,d,ins=(a,)) for d in range(Dmax+1)]
                rhs=q_to_series(c,gw)
                print("g",g,"H^%d"%a,"match",lhs==rhs,[str(x) for x in lhs],round(time.time()-gt,1),flush=True)
if mode=='g3':
    Dmax=int(sys.argv[2])
    for al,ep in samples[:int(sys.argv[3])]:
        t=time.time()
        c=setup(al,ep,Dmax,8); gs=GS(c)
        F=gs.total(3)
        print("graph sum",[str(x) for x in F],round(time.time()-t,1),flush=True)
        gw=[loc(al,ep,3,d) for d in range(Dmax+1)]
        R=q_to_series(c,gw)
        print("loc      ",[str(x) for x in R],"MATCH",F==R,round(time.time()-t,1),flush=True)
