# -*- coding: utf-8 -*-
"""rf_ur.py -- REFEREE: read-only wrapper of the PROJECT module natural_form/theory/ur_gauge.py (the stored regular gauge), and the
canonical normalisation done here (not with gi_core).  The conifold frame is chosen by rebinding ur_gauge.Conifold IN MEMORY (no file
is modified).  Output: p[k][b] = list of Fractions (u^0 .. u^(prec-1)) of the entry (b + k, b) of [z^k] P, prec common to all."""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory"))
import fq_conifold
import ur_gauge


def stored_gauge(K, NU, al=0, be=0):
    orig = ur_gauge.Conifold
    ur_gauge.Conifold = lambda R, Nu=40: fq_conifold.Conifold(R, Nu=Nu, alpha=Fr(al), beta=Fr(be))
    try:
        res = ur_gauge.gauge(K, NU)
    finally:
        ur_gauge.Conifold = orig
    p = res['p']
    prec = min(x.prec for pk in p for x in pk)
    out = []
    for k in range(K + 1):
        row = []
        for b in range(5):
            s = p[k][b]
            assert (not s.c) or s.v >= 0, "stored gauge not regular at (k, b) = (%d, %d)" % (k, b)
            row.append([Fr(s.coef(n)) for n in range(prec)])
        out.append(row)
    return out, prec


def p11_zero(p):
    """{k: p_{k,1}(0)} for k = 0 mod 5 (the series P_11(0, z))"""
    return {k: p[k][1][0] for k in range(0, len(p), 5)}


def normalise(p, c):
    """P c(z)^-1 for a scalar series c = {5q: Fraction} in z^5 (c[0] != 0): pc[k][b] = Sum_j cinv_j p[k-j][b]"""
    K = len(p) - 1
    ci = {0: 1 / c[0]}
    for k in range(5, K + 1, 5):
        ci[k] = -sum(c.get(j, 0) * ci[k - j] for j in range(5, k + 1, 5)) / c[0]
    out = []
    for k in range(K + 1):
        row = []
        for b in range(5):
            acc = [Fr(0)] * len(p[k][b])
            for j, x in ci.items():
                if j <= k and x:
                    acc = [s + x * y for s, y in zip(acc, p[k - j][b])]
            row.append(acc)
        out.append(row)
    return out, ci
