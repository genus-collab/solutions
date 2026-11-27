# -*- coding: utf-8 -*-
"""ms_corr.py -- helpers on top of ms_engine: admissible insertion lists, the gap constants, the vanishing covector PiC.

Conventions (ms_engine): an insertion is (a, k) = tau_k(H^a), H^a the flat class (H = d/dv), correlators are the ANCESTOR
correlators Int_{Mbar_{h,n}} V_{h,n}(H^a1, ..) prod psi^k of the frozen theory V = M_0.omega^ = the LG model at v = 1;
at a general v:  <..>_h(v) = <..>_h(1) v^e,  e = ms_engine.vexp(h, ins).
F_h(v) = c_h v^(2-2h) (h >= 2),  F_1 = -(1/12) log v;  c_h = (-1)^(h-1) B_2h/(2h(2h-2) 5^(h-1)).
PiC(tau_k(H^b)) = ms_model.pic(b, k) (coefficient at v = 1, power of v): the vanishing covector (pairing with S(-z)H^2/5).
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr
import ms_engine as me
import ms_model as mm


def cg(h):
    if h == 0 or h == 1: raise ValueError("c_h only for h >= 2")
    sys.path.insert(1, mr.PROJ)
    from fq_fp import bernoulli
    return (-1) ** (h - 1) * Fr(bernoulli(2 * h)) / (2 * h * (2 * h - 2)) / 5 ** (h - 1)


def Fder(h, m, v=Fr(1)):
    """m-th v-derivative of F_h at v (h >= 1; for h = 1 and m = 0 not defined)"""
    if h == 1:
        if m == 0: raise ValueError("F_1 = -(1/12) log v")
        x = Fr(-1, 12)                     # d/dv: -(1/12)/v ; d^m: -(1/12) (-1)^(m-1) (m-1)! v^-m
        f = 1
        for j in range(1, m): f *= j
        return x * (-1) ** (m - 1) * f / v ** m
    c = cg(h); e = 2 - 2 * h
    for j in range(m): c *= (e - j)
    return c * v ** (e - m)


def pic(ins):
    """PiC of one insertion (a, k) at v = 1"""
    return mm.pic(ins[0], ins[1])[0]


def charge(ins):
    return sum(a + k - 1 for a, k in ins) % 5


def all_insertions(h, n, kmax=None):
    """sorted tuples of n insertions (a, k), a in 0..4, Sum k <= 3h - 3 + n, total charge 0"""
    D = 3 * h - 3 + n
    if kmax is None: kmax = D
    single = [(a, k) for k in range(kmax + 1) for a in range(5)]
    out = []
    for combo in itertools.combinations_with_replacement(single, n):
        if sum(k for a, k in combo) > D: continue
        if charge(combo) != 0: continue
        out.append(tuple(sorted(combo)))
    return sorted(set(out))


def corr(h, ins):
    return me.corr(h, tuple(ins))


def vexp(h, ins):
    e = me.vexp(h, tuple(ins)); assert e.denominator == 1
    return int(e)


def fmt_ins(ins):
    return " ".join("t%d(H%d)" % (k, a) for a, k in ins)
