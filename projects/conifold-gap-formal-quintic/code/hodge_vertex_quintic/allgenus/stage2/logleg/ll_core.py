# -*- coding: utf-8 -*-
"""ll_core.py -- leg pairings of the TR correlators of the model curve x = W = p^4/4 + u/p, y = p (u = 1), stage2/logleg.

Engine: the allgenus/tr modules (imported READ-ONLY): tr_field (Q(zeta5), F_p), tr_core.ModelCurve, tr_rec.TR (recursion (2-28) of
arXiv:1509.06954 in the basis e_{j,k}(p) = dp/(p - c_j)^(k+1), k >= 1, c_j = zeta^-j), tr_v0.at_infinity (e-basis -> moments at p = oo).
Pairing of a leg Phi (a function) with one slot:  <Phi, w> := Sum_b Res_{p=c_b} Phi(p) w(p).  On the e-basis:
  L-leg  p^m (m any integer):      Res_{c_j} p^m e_{j,k} = binom(m, k) c_j^(m-k)          [= coefficient of dp/p^(m+1) at oo for m >= 1]
  log leg -log p:                  Res_{c_j} (-log p) e_{j,k} = (-1)^k / (k c_j^k)        [= Int_0^oo e_{j,k}; cut along a ray 0 -> oo]
  control -log(p - a):             (-1)^k / (k (c_j - a)^k)
Vector field xi = p/(p^5 - u) d/dp (d_v|_p = d_v|_x + Lie_xi, critic's convention) and the Euler derivative d_u|_p at u = 1:
  d_u|_p  [C prod e_{j_i,k_i}] = C [ (E/5) prod e + Sum_i ((k_i+1) c_{j_i}/5) e_{j_i,k_i+1} prod_{l != i} e ],  E = 5(2-2g-n) + Sum k_i
  Lie_xi e_{j,k} = d(g),  g = p/((p^5-1)(p-c_j)^(k+1)) = Sum_{l,r} a_{l,r}/(p-c_l)^r  =>  Lie_xi e_{j,k} = Sum -r a_{l,r} e_{l,r}.
All functions take a field F (tr_field) and the branch points c = curve.c.  py -I -B only (stdlib).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
TRDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'tr'))
if TRDIR not in sys.path: sys.path.insert(0, TRDIR)
from tr_field import QzField, make_fp, P89, P61          # noqa: E402  (read-only import of allgenus/tr)
from tr_core import ModelCurve, binom, smul, sinv        # noqa: E402
from tr_rec import TR                                    # noqa: E402
from tr_v0 import at_infinity                            # noqa: E402


def setup(F, N):
    cur = ModelCurve(F, N, 1); return cur, TR(cur)


def fpow(F, a, e):
    if e < 0: a = F.one / a; e = -e
    r = F.one
    for _ in range(e): r = r * a
    return r


class Legs(object):
    """cached one-slot functionals on the e-basis"""
    def __init__(self, F, c):
        self.F = F; self.c = c; self._c = {}

    def logleg(self, jk):
        key = ('log', jk); r = self._c.get(key)
        if r is None:
            j, k = jk; r = self.F.of(Fr((-1) ** k, k)) * fpow(self.F, self.c[j], -k); self._c[key] = r
        return r

    def shiftlog(self, a):
        F = self.F; a = F.of(a)
        def f(jk):
            j, k = jk
            return F.of(Fr((-1) ** k, k)) * fpow(F, self.c[j] - a, -k)
        return f

    def mono(self, m):
        F = self.F
        def f(jk):
            j, k = jk; key = ('m', m, jk); r = self._c.get(key)
            if r is None:
                r = F.of(binom(m, k)) * fpow(F, self.c[j], m - k); self._c[key] = r
            return r
        return f


def contract_last(W, funcs):
    """contract the LAST len(funcs) slots of the e-basis tensor W with the given one-slot functionals"""
    r = len(funcs); out = {}
    if r == 0: return dict(W)
    for key, v in W.items():
        n = len(key); x = v
        for i, f in enumerate(funcs):
            x = x * f(key[n - r + i])
            if x.iszero(): break
        if x.iszero(): continue
        kk = key[:n - r]
        out[kk] = out[kk] + x if kk in out else x
    return {k: v for k, v in out.items() if not v.iszero()}


def pair_table(F, c, W, n, klog, mtot, legs=None):
    """{(m_1..m_a): <p^m_1, .., p^m_a, (-log p)^klog>} for Sum m <= mtot (a = n - klog >= 1), or {(): value} if a = 0"""
    legs = legs or Legs(F, c)
    R = contract_last(W, [legs.logleg] * klog)
    a = n - klog
    if a == 0: return R
    return at_infinity(F, c, R, a, mtot)


def compositions(total, parts):
    """ordered tuples of `parts` positive integers with the given sum"""
    if parts == 1:
        if total >= 1: yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1): yield (first,) + rest


def bernoulli(n):
    from math import comb
    B = [Fr(1)]
    for m in range(1, n + 1): B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B[n]


def Fh_deriv(h, k):
    """k-th u-derivative at u = 1 of the TR free energy of x = W:  h >= 2: B_2h/(2h(2h-2)) u^(2-2h) (Gap Theorem, TR normalisation);
    h = 1: -(1/12) log u;  h = 0: (1/2) u^2 log u (k >= 3).  (frozen = (-5)^(1-h) x TR.)"""
    if h >= 2:
        r = bernoulli(2 * h) / (2 * h * (2 * h - 2)); e = 2 - 2 * h
        for i in range(k): r *= (e - i)
        return r
    if h == 1:
        assert k >= 1
        r = Fr(-1, 12); f = Fr(1)
        for i in range(1, k): f *= -i
        return r * f
    assert k >= 3
    f = Fr(1)
    for i in range(1, k - 2): f *= -i
    return f                                                    # d^3 [(1/2) u^2 log u] = 1/u


# ------------------------------------------------------------------ d_u at fixed p, Lie_xi, Int_0^oo in a slot
def dv_p(F, c, W, g, n):
    out = {}; fifth = F.of(Fr(1, 5))
    def add(k, v):
        if not v.iszero(): out[k] = out[k] + v if k in out else v
    for key, C in W.items():
        E = 5 * (2 - 2 * g - n) + sum(k for _, k in key)
        add(key, C * F.of(Fr(E, 5)))
        for i, (j, k) in enumerate(key):
            add(key[:i] + ((j, k + 1),) + key[i + 1:], C * c[j] * fifth * (k + 1))
    return {k: v for k, v in out.items() if not v.iszero()}


class Lie(object):
    """Lie_xi e_{j,k} in the e-basis, xi = p/(p^5 - 1) d/dp (u = 1)"""
    def __init__(self, F, c, sign=1):
        self.F = F; self.c = c; self.sign = sign; self._t = {}

    def table(self, j, k):
        key = (j, k); r = self._t.get(key)
        if r is not None: return r
        F = self.F; c = self.c; cj = c[j]; M = k + 2
        Q = [F.zero] * (M + 1)                                   # Q_j(c_j + t) = Sum_i c_j^(4-i) (c_j + t)^i
        for i in range(5):
            w = fpow(F, cj, 4 - i)
            for s in range(min(i, M) + 1): Q[s] = Q[s] + w * F.of(binom(i, s)) * fpow(F, cj, i - s)
        num = [cj, F.one] + [F.zero] * (M - 1)
        s = smul(num, sinv(Q, M, F), M, F)                      # (c_j + t)/Q_j(c_j + t)
        out = {}
        for r_ in range(1, k + 3):                              # a_{j,r} = s_{k+2-r}
            a = s[k + 2 - r_]
            if not a.iszero(): out[(j, r_)] = a * (-r_) * self.sign
        for l in range(5):
            if l == j: continue
            a = c[l] / (fpow(F, c[l], 4) * 5 * fpow(F, c[l] - cj, k + 1))
            out[(l, 1)] = out.get((l, 1), F.zero) + a * (-1) * self.sign
        self._t[key] = out
        return out

    def apply(self, W):
        out = {}
        for key, C in W.items():
            for i, jk in enumerate(key):
                for e, a in self.table(*jk).items():
                    kk = key[:i] + (e,) + key[i + 1:]; v = C * a
                    out[kk] = out[kk] + v if kk in out else v
        return {k: v for k, v in out.items() if not v.iszero()}


def sub(F, A, B):
    out = dict(A)
    for k, v in B.items(): out[k] = out[k] - v if k in out else -v
    return {k: v for k, v in out.items() if not v.iszero()}


def fmt(F, x):
    if hasattr(F, 'ratrec'): return F.fmt(x)
    return repr(x.rational()) if x.is_rational() else repr(x)


def rat(F, x):
    """Fraction (exact field: rational part; F_p: rational reconstruction, None if it fails)"""
    if hasattr(F, 'ratrec'): return F.ratrec(x)
    return x.rational() if x.is_rational() else None
