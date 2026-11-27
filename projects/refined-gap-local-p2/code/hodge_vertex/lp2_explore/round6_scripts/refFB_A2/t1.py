from gsum import *
import time,sys
samples=[((Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7))),((Fr(2),Fr(-3),Fr(8)),(Fr(1),Fr(4))),((Fr(1,3),Fr(5,2),Fr(-4)),(Fr(11,5),Fr(-2,7)))]
al,ep=samples[0]
e12=ep[0]*ep[1]
print("g=0 d=1,2,3 *e1e2 (3,-45/8,244/9):",[loc(al,ep,0,d)*e12 for d in (1,2,3)])
print("g=2 d=0,1 loc:",loc(al,ep,2,0),loc(al,ep,2,1))
Dmax=int(sys.argv[1]) if len(sys.argv)>1 else 2
for al,ep in samples:
    t=time.time()
    c=setup(al,ep,Dmax,5)
    print(checks(c))
    gs=GS(c)
    F2=gs.total(2)
    gw=[loc(al,ep,2,d) if d<=Dmax else 0 for d in range(Dmax+1)]
    print("graph sum y-series",F2)
    print("loc -> y-series   ",q_to_series(c,gw))
    print("MATCH g=2:",F2==q_to_series(c,gw),time.time()-t)
    # control: c->-c
    c2=setup(al,ep,Dmax,5,cflip=True)
    print("control c->-c match?",GS(c2).total(2)==q_to_series(c,gw))
