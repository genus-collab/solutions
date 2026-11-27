# -*- coding: utf-8 -*-
"""rf_tables.py -- flat tables from the referee amplitudes, and exact values by CRT.

flat_table(model, h, n): {((a_1,k_1),..,(a_n,k_n)) ordered: <prod tau_k(H^a)>_h mod p}, by contracting the legs of the amputated
amplitude one at a time with L[a][j][c] = (R^-1_j Psi^-1 H^a)_c.  exact(h, ins): rational from both primes.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as rc
import rf_gsum as gsm
import rf_flat as rf

_T = {}


def flat_table(m, h, n):
    key = (m.p, h, n)
    if key in _T: return _T[key]
    p = m.p; A = gsm.amp(m, h, n)
    L = [[[sum(m.Rinv[j][c][b] * m.flat[a][b] for b in range(5)) % p for c in range(5)] for j in range(m.K + 1)] for a in range(5)]
    cur = {}
    for (cols, exps), val in A.items():
        cur[((), cols, exps)] = val
    for l in range(n):
        nxt = {}
        for (done, cols, exps), val in cur.items():
            c = cols[0]; d = exps[0]
            for a in range(5):
                for j in range(d + 1):
                    x = L[a][j][c]
                    if x:
                        k2 = (done + ((a, d - j),), cols[1:], exps[1:])
                        nxt[k2] = (nxt.get(k2, 0) + val * x) % p
        cur = nxt
    out = {k[0]: v for k, v in cur.items() if v}
    _T[key] = out
    return out


def exact_val(h, ins):
    m1, m2 = rf.models()
    t1 = flat_table(m1, h, len(ins)); t2 = flat_table(m2, h, len(ins))
    x1 = t1.get(tuple(ins), 0); x2 = t2.get(tuple(ins), 0)
    Mod = m1.p * m2.p
    x = (x1 + m1.p * ((x2 - x1) * pow(m1.p, m2.p - 2, m2.p) % m2.p)) % Mod
    return rf.ratrec(x, Mod)


def Fder(h, n):
    """d_v^n F_h at v = 1: F_0 = -(5/2) v^2 log v (n >= 3), F_1 = -(1/12) log v (n >= 1), F_h = c_h v^(2-2h)"""
    if h == 0:
        c = Fr(-5)
        for j in range(n - 3): c *= -(j + 1)
        return c
    if h == 1:
        c = Fr(-1, 12)
        for j in range(1, n): c *= -j
        return c
    B = {2: Fr(-1, 30), 3: Fr(1, 42), 4: Fr(-1, 30), 5: Fr(5, 66)}[h]
    c = (-1) ** (h - 1) * B / (2 * h * (2 * h - 2)) / 5 ** (h - 1); e = 2 - 2 * h
    for j in range(n): c *= (e - j)
    return c
