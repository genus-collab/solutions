# -*- coding: utf-8 -*-
"""cyc_ref_ops.py -- referee cyc_formal: the differential ring R = K[l] and operators over K (stdlib only).

  REl   element of R = K[l], K = Q((u)): dict k -> Ser (coefficient of l^k), with D l = dl (default DL = 1 - 1/u);
        shift_l(c) is the K-algebra map l -> l + c (c = 1: the automorphism tau of CYC s.5).
  Op    operator Sum_k a_k D^k, coefficients on the LEFT; composition via D^j o b = Sum_m C(j,m) D^m(b) D^(j-m);
        adjoint() = Sum_k (-D)^k o a_k.
Imports only cyc_ref_core (this referee's own module).
"""
from fractions import Fraction as Fr
from math import comb
from cyc_ref_core import Ser, ONE, U, DL, EX


class REl:
    def __init__(self, F, dl=DL):
        self.F = dict(F); self.dl = dl

    @staticmethod
    def of(s, dl=DL):
        return REl({0: s}, dl)

    def __add__(a, b):
        d = dict(a.F)
        for k, s in b.F.items():
            d[k] = d[k] + s if k in d else s
        return REl(d, a.dl)

    def __neg__(a):
        return REl({k: -s for k, s in a.F.items()}, a.dl)

    def __sub__(a, b):
        return a + (-b)

    def mulS(a, s):
        return REl({k: s * x for k, x in a.F.items()}, a.dl)

    def sc(a, x):
        return REl({k: s.sc(x) for k, s in a.F.items()}, a.dl)

    def __mul__(a, b):
        d = {}
        for k1, s1 in a.F.items():
            for k2, s2 in b.F.items():
                t = s1 * s2; k = k1 + k2
                d[k] = d[k] + t if k in d else t
        return REl(d, a.dl)

    def D(a):
        d = {}
        for k, s in a.F.items():
            t = s.D(); d[k] = d[k] + t if k in d else t
            if k:
                t = (s * a.dl).sc(Fr(k)); d[k - 1] = d[k - 1] + t if (k - 1) in d else t
        return REl(d, a.dl)

    def shift_l(a, c):                       # the K-algebra map l -> l + c (c a Ser; c = 1 is the automorphism tau)
        d = {}
        for k, s in a.F.items():
            for i in range(k + 1):           # (l + c)^k = Sum_i C(k,i) c^(k-i) l^i
                cp = ONE
                for _ in range(k - i):
                    cp = cp * c
                t = (s * cp).sc(Fr(comb(k, i))); d[i] = d[i] + t if i in d else t
        return REl(d, a.dl)

    def zero_below(a, M):
        return all(s.zero_below(M) for s in a.F.values())


class Op:
    def __init__(self, A):
        self.A = {k: s for k, s in A.items()}

    @staticmethod
    def fn(s):
        return Op({0: s})

    @staticmethod
    def Dk(k=1):
        return Op({k: ONE})

    def __add__(a, b):
        d = dict(a.A)
        for k, s in b.A.items():
            d[k] = d[k] + s if k in d else s
        return Op(d)

    def __neg__(a):
        return Op({k: -s for k, s in a.A.items()})

    def __sub__(a, b):
        return a + (-b)

    def scal(a, x):
        return Op({k: s.sc(Fr(x)) for k, s in a.A.items()})

    def __mul__(a, b):                       # composition a o b
        d = {}
        for j, x in a.A.items():
            for k, y in b.A.items():
                dy = y
                for m in range(j + 1):
                    t = (x * dy).sc(Fr(comb(j, m))); key = j - m + k
                    d[key] = d[key] + t if key in d else t
                    if m < j:
                        dy = dy.D()
        return Op(d)

    def __call__(a, f):                      # apply to an REl
        out, cur = None, f
        for k in range(max(a.A) + 1):
            if k in a.A:
                t = cur.mulS(a.A[k]); out = t if out is None else out + t
            cur = cur.D()
        return out

    def adjoint(a):                          # Sum_k (-D)^k o a_k
        out = Op({})
        for k, x in a.A.items():
            out = out + (Op.Dk(k) * Op.fn(x)).scal((-1) ** k)
        return out

    def diff(a, b, M):                       # a - b has every coefficient zero below u^M
        c = a - b
        return all(s.zero_below(M) for s in c.A.values())
