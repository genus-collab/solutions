# -*- coding: utf-8 -*-
"""cy_core.py -- referee cyc_adversary (Stage 5): exact algebra over Q for an independent test of (CYC).
Own code; nothing is imported from the project.

  Ser   truncated Laurent series in u over Q:  Sum_i c[i] u^(v+i) + O(u^prec); stored coefficients past the end are 0
        (prec = INF for exact polynomials).
  Rel   elements of R = K[l], K = Q((u)), l = log u, D l = 1 - 1/u: list of Ser (coefficient of l^k).
  operators: lists [e_0, .., e_n] of Ser, meaning Sum_k e_k D^k, D = q d/dq = -(1-u) d/du.
  PF4, its closed form on u^x and the Frobenius basis: cy_frob.py (split off to keep every file under 250 lines).
"""
from fractions import Fraction as Fr
from math import comb

INF = 10 ** 9


class Ser(object):
    __slots__ = ('v', 'c', 'prec')

    def __init__(self, v, c, prec=INF):
        n = prec - v
        c = [Fr(x) for x in c]
        if len(c) > n: c = c[:max(n, 0)]
        while c and c[-1] == 0: c.pop()
        self.v, self.c, self.prec = v, c, prec

    @staticmethod
    def poly(d, prec=INF):
        if not d: return Ser(0, [], prec)
        v = min(d); return Ser(v, [d.get(k, 0) for k in range(v, max(d) + 1)], prec)

    def coef(self, k):
        if k >= self.prec: raise ValueError('coefficient u^%d beyond precision %d' % (k, self.prec))
        i = k - self.v
        return self.c[i] if 0 <= i < len(self.c) else Fr(0)

    def add(a, b):
        v = min(a.v, b.v); p = min(a.prec, b.prec)
        top = min(p, max(a.v + len(a.c), b.v + len(b.c)))
        return Ser(v, [a.coef(k) + b.coef(k) for k in range(v, top)], p)

    def neg(a): return Ser(a.v, [-x for x in a.c], a.prec)
    def sub(a, b): return a.add(b.neg())
    def scale(a, x): x = Fr(x); return Ser(a.v, [x * y for y in a.c], a.prec)
    def shift(a, k): return Ser(a.v + k, a.c, a.prec + k)            # times u^k
    def trunc(a, p): return Ser(a.v, a.c, min(a.prec, p))

    def val(a):
        for i, x in enumerate(a.c):
            if x: return a.v + i
        return None

    def mul(a, b):
        v = a.v + b.v; p = min(a.v + b.prec, b.v + a.prec)
        n = min(p - v, len(a.c) + len(b.c) - 1)
        out = [Fr(0)] * max(n, 0)
        for i, x in enumerate(a.c):
            if not x or i >= n: continue
            for j, y in enumerate(b.c):
                if i + j >= n: break
                if y: out[i + j] += x * y
        return Ser(v, out, p)

    def inv(a, cap=None):
        w = a.val()
        if w is None: raise ZeroDivisionError('inverse of a series with no known non-zero coefficient')
        r = a.prec - w
        if r > 10 ** 6: r = cap
        c = a.c[w - a.v:]; c0 = c[0]
        out = [Fr(0)] * r; out[0] = 1 / c0
        for n in range(1, r):
            s = Fr(0)
            for k in range(1, min(n, len(c) - 1) + 1): s += c[k] * out[n - k]
            out[n] = -s / c0
        return Ser(-w, out, -w + r)

    def deriv(a): return Ser(a.v - 1, [(a.v + i) * x for i, x in enumerate(a.c)], a.prec - 1)

    def D(a):
        d = a.deriv(); return d.shift(1).sub(d)                       # -(1-u) d/du

    def integ0(a):
        """Int_0^u; requires no coefficient of u^k, k < 0, within the stored list"""
        assert all(a.coef(k) == 0 for k in range(a.v, min(0, a.prec))), 'integrand not a power series'
        top = min(a.prec, a.v + len(a.c))
        cs = [a.coef(k) / (k + 1) for k in range(max(a.v, 0), top)]
        return Ser(max(a.v, 0) + 1, cs, a.prec + 1)

    def is_zero(a): return all(x == 0 for x in a.c)
    def head(a, n=4): return ' + '.join('(%s)u^%d' % (x, a.v + i) for i, x in enumerate(a.c[:n]) if x) or '0'


ONE_M = Ser.poly({0: 1, 1: -1})          # 1 - u
DL = Ser.poly({-1: -1, 0: 1})            # D log u = 1 - 1/u
U1 = Ser.poly({1: 1})


class Rel(object):
    """Sum_k cs[k] l^k, l = log u"""
    def __init__(self, cs): self.cs = list(cs)

    @staticmethod
    def of(f): return Rel([f])

    def add(a, b):
        n = max(len(a.cs), len(b.cs)); z = Ser(0, [])
        return Rel([(a.cs[k] if k < len(a.cs) else z).add(b.cs[k] if k < len(b.cs) else z) for k in range(n)])

    def neg(a): return Rel([c.neg() for c in a.cs])
    def scale(a, x): return Rel([c.scale(x) for c in a.cs])
    def mulser(a, f): return Rel([c.mul(f) for c in a.cs])

    def mul(a, b):
        out = {}
        for i, x in enumerate(a.cs):
            for j, y in enumerate(b.cs):
                t = x.mul(y); out[i + j] = t if i + j not in out else out[i + j].add(t)
        return Rel([out[k] for k in range(max(out) + 1)])

    def D(a):
        out = [c.D() for c in a.cs]
        for k in range(1, len(a.cs)): out[k - 1] = out[k - 1].add(a.cs[k].mul(DL).scale(k))
        return Rel(out)

    def is_zero(a): return all(c.is_zero() for c in a.cs)


# ------------------------------------------------------------------ operators
def op_add(A, B):
    n = max(len(A), len(B)); z = Ser(0, [])
    return [(A[k] if k < len(A) else z).add(B[k] if k < len(B) else z) for k in range(n)]


def op_D_left(E):
    """D o E"""
    out = [e.D() for e in E] + [Ser(0, [])]
    for k, e in enumerate(E): out[k + 1] = out[k + 1].add(e)
    return out


def op_compose(A, B):
    """A o B by Leibniz:  a D^i o b D^j = a Sum_m C(i,m) (D^m b) D^(i-m+j)"""
    out = {}
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            Db = b
            for m in range(i + 1):
                t = a.mul(Db).scale(comb(i, m)); k = i - m + j
                out[k] = t if k not in out else out[k].add(t)
                Db = Db.D()
    return [out[k] for k in range(max(out) + 1)]


def op_adjoint(A):
    """Sum_k (-D)^k o A_k"""
    out = [Ser(0, [])]
    for k, a in enumerate(A):
        T = [a]
        for _ in range(k): T = op_D_left(T)
        if k % 2: T = [t.neg() for t in T]
        out = op_add(out, T)
    return out


def op_apply(E, f):
    """E f for a Ser or a Rel f"""
    rel = isinstance(f, Rel); g = f; acc = None
    for k, e in enumerate(E):
        t = g.mulser(e) if rel else g.mul(e)
        acc = t if acc is None else acc.add(t)
        g = g.D()
    return acc


def op_diff_zero(E, F):
    """(all coefficients of E - F vanish within precision, minimal precision of the difference)"""
    d = op_add(E, [f.neg() for f in F])
    return all(x.is_zero() for x in d), min(x.prec for x in d)
