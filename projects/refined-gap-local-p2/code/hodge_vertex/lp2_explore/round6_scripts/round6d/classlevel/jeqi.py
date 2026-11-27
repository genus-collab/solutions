# eq:J=I checked against own genus-0 GP, through y^3 and 1/z^6, as series in s=1/z at each P_j.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from rows_p2 import P2Data
from gp import GP
from fractions import Fraction as Fr
for al,ep in [((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7)))]:
    D=3;Ns=6
    P=P2Data(Fr,al,ep,D,1);G=GP(al,ep)
    t=Ctx(Fr,D,Ns)   # z slot used as s=1/z
    Qser=None
    L={d:Fr(P.Ay[d])/d for d in range(1,D+1)}
    u=Bi(t,{(d,0):L[d] for d in L})
    Q=ymono(t,1,1)*bexp(u)
    Qp=[const(t,1)]
    for d in range(1,D+1): Qp.append(Qp[-1]*Q)
    etai=P.etai_f
    ok=True
    for j in range(3):
        a=P.al[j]
        I=const(t,1)
        for d in range(1,D+1):
            I=I+Bi(t,{(d,n):x for n,x in enumerate(P.c_at_inf(d,j,Ns)) if x!=0})
        J=const(t,1)
        for d in range(1,D+1):
            co={}
            for k in range(0,Ns-1):
                v=Fr(0)
                for b in range(3):
                    cor=G.value(0,[(k,b)],d,'desc')
                    for c in range(3):
                        v+=cor*etai[b][c]*(-a)**c
                co[(0,k+2)]=v
            J=J+Qp[d]*Bi(t,co).clean()
        # exp(-alpha u s)
        ex=Bi(t,{(d,1):-a*L[d] for d in L})
        J=bexp(ex)*J
        diff=J-I
        if not diff.is_zero():
            ok=False;print('J!=I at',j,sorted(diff.c.items())[:5])
    print('eq:J=I against genus-0 GP, weights',al,ep,'through y^3, s^6:',ok,flush=True)
