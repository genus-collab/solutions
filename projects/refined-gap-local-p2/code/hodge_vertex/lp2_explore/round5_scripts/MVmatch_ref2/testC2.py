# (C) with R nonempty: literal Schur FF vs table, termination, values; controls on the substitution.
from fractions import Fraction as F
from testC import *
from hodge import Qk
from functools import reduce
import sys,time

def falling(J,r):
    v=F(1)
    # (J-1)(J-2)...(J-r+3)
    for i in range(1,r-2):
        v*=J-i
    return v

for tau in [F(5,3),F(-5,8)]:
    sig=tau*(1+tau)
    for R,M,G in [((2,),6,1),((3,),6,1),((4,),5,2),((2,),10,3),((5,),9,3),((2,2),8,2),((3,2),10,3),((2,2,2),7,1),((3,2,2),6,1),((3,3,2,2),5,0),((4,),14,4)]:
        t0=time.time()
        ff=FF(R,tau,M,G)
        r=len(R);J=sum(R)
        msg=[f"tau={tau} R={R} M={M} G={G} forb={ff.forbidden()}"]
        for g in range(G+1):
            q=ff.QQ(g)
            n=2*g-2+r
            if n<=0:
                if g==0 and r==1:
                    ex=[F(0)]*(M+1);ex[0]=F(1,J*J)
                    if M>=1: ex[1]=F(1,J*J)-F(1,R[0]*(R[0]+1))
                    msg.append(f"(0,(j)) {q==ex}")
                if g==0 and r==2:
                    ex=[F(0)]*(M+1);ex[0]=F(1,J)
                    msg.append(f"(0,(j,k)) {q==ex}")
                continue
            N=3*g-3+r
            if N>M:
                msg.append(f"g={g} N={N}>M skip");continue
            tail=all(v==0 for v in q[N+1:])
            xi=sum((-1)**m*q[m] for m in range(N+1))
            lit=None
            if g<=3 and N<=10:
                lit=all(q[m]==Qk(g,tau,R,m) for m in range(N+1))
            part=f"g={g} N={N} tail0={tail}({M-N}) lit={lit} Xi={'0' if xi==0 else xi} D={J-r} range={J-r<=2*g-3}"
            if g==0:
                part+=f" g0value={xi==falling(J,r)}"
            if g==1 and r==1:
                ex=[F(R[0]-1)/24-1/(24*sig),F(1,24)]
                part+=f" g1table={q[:2]==ex}"
            msg.append(part)
        print(" | ".join(msg),round(time.time()-t0,1),flush=True)
    # controls: must break: no x' ; no e^{-u|nu|}; sign of s ; no Pi
    for name,kw in [("x'=x",dict(xshift=False)),("no e^-u",dict(eu=False)),("s->-s",dict(sgn=-1)),("no Pi",dict(usePi=False))]:
        ff=FF((3,),tau,5,1,**kw)
        q=ff.QQ(1)
        ex=[F(3-1)/24-1/(24*sig),F(1,24)]+[F(0)]*4
        print("control",name,"(1,(3)) table holds:",q==ex,"forb",ff.forbidden(),flush=True)
