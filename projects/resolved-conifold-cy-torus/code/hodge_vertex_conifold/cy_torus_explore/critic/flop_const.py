# Critic: closed form on T_CY implies F_g(Q=oo) = F_g(Q=0). At Q = oo the two idempotents have L = a1, a2, i.e. the fixed points of the flop:
# wt'_1 = (a1-a2, a2, a0+a2, -e1, -e2), wt'_2 = (a2-a1, a1, a0+a1, -e1, -e2). Test the Hodge-integral identity, g = 2..6, CY and non-CY.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
from cylib import *
from gsum import Omega

def cm(wts, g):
    r = 0
    for wt in wts:
        wt = [fp(x) for x in wt]
        e = 1
        for x in wt: e = e * x % P
        r += pow(e, g - 1, P) * Omega(wt).integral(g, ())
    return r % P

for (a0, a1, a2, e1, e2) in [(Fr(7,5),Fr(-3,8),Fr(13,6),Fr(11,9),None),(Fr(-5,4),Fr(19,7),Fr(2,3),Fr(17,5),None),(Fr(1),Fr(2),Fr(4),Fr(3),Fr(5))]:
    if e2 is None: e2 = a0 + a1 + a2 - e1
    for g in range(2, 7):
        v0 = cm([(a0,a1,a2,-e1,-e2),(-a0,a0+a1,a0+a2,-e1,-e2)], g)
        voo = cm([(a1-a2,a2,a0+a2,-e1,-e2),(a2-a1,a1,a0+a1,-e1,-e2)], g)
        print((str(a0),str(a1),str(a2),str(e1),str(e2)), 'g', g, 'F_g(0) == F_g^flop(0):', v0 == voo)
