# Part 1 at literal rational square samples, and at higher order mod p.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from part1 import *
from fractions import Fraction as Fr
found=[]
N=60
for a0 in list(range(-N,0))+list(range(1,N+1)):
    for a1 in range(-N,N+1):
        if a1==0 or a1==a0: continue
        for a2 in range(a1+1,N+1):
            if a2==0 or a2==a0: continue
            r1=Fr(-a1*(a2-a1),a0*(a2-a0))
            if r1<=0 or rat_sqrt(r1) is None: continue
            r2=Fr(a2*(a2-a1),a0*(a1-a0))
            if r2<=0 or rat_sqrt(r2) is None: continue
            b0=3*a0*(a1-a0)*(a2-a0)
            found.append(((Fr(a0),Fr(a1),Fr(a2)),(Fr(1),Fr(b0))))
print('rational samples with all W_i^0 squares (eps=(1,b0)):',len(found),[(tuple(map(str,a)),tuple(map(str,e))) for a,e in found[:6]],flush=True)
for al,ep in found[:3]:
    run(Fr,al,ep,3,5,'Q (literal square sample)')
    run(Fr,al,ep,4,7,'Q (literal square sample)')
plist=primes3mod4(2**127,400)
for al,ep in [((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7)))]:
    for p in plist:
        F=make_fp(p)
        P=P2Data(F,al,ep,1,1)
        if all(fsqrt(P.W0(i)) is not None for i in range(3)):
            run(F,al,ep,6,9,'F_p, high order, p=%d'%p)
            break
