# -*- coding: utf-8 -*-
"""g2_core.py -- self-contained exact arithmetic for the G2 gap checks (Theorem U, the gauge, colours).

Stdlib only.  No project module is imported.  Truncated Laurent series in u over Fraction:
  Ser(v, c, N): value Sum_i c[i] u^(v+i), known modulo u^N (absolute precision N).
Systems (notation of the note, NOTATION.md): the quintic one-pole system b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2),
t_c', -varpi_0/(1-u), -varpi_0/(1-u)) in the programme frame varpi_0 = f_0 + alpha varpi_1 + beta f_2, and the pulled-back
model bt = (v', -v'/v, v', -v', -v'), v = u tau^(-5/4), 1/tau = Sum 4u^n/(5n+4).  The Frobenius basis f_0, f_1, f_2 of PF4 and
the log obstruction (gamma) are computed here from PF4 = u D^4 - (1-u)(2 D^3 + (7/5) D^2 + (2/5) D + 24/625), D = -(1-u) d/du.
"""
from fractions import Fraction as Fr


class Ser(object):
    __slots__ = ("v", "c", "N")

    def __init__(self, v, c, N):
        c = [Fr(x) for x in c]
        while c and c[0] == 0:
            c.pop(0); v += 1
        if v >= N or not c:
            v, c = N, []          # the zero series carries v = N, so it never lowers the valuation of a sum
        c = c[:N - v]
        self.v, self.c, self.N = v, c, N

    @staticmethod
    def const(a, N):
        return Ser(0, [a], N)

    @staticmethod
    def mono(k, N, a=1):
        return Ser(k, [a], N)

    def iszero(self):
        return not self.c

    def coef(self, k):
        if k >= self.N:
            raise ValueError("precision: asked u^%d, known mod u^%d" % (k, self.N))
        return self.c[k - self.v] if self.v <= k < self.v + len(self.c) else Fr(0)

    def __add__(self, o):
        if not isinstance(o, Ser):
            o = Ser.const(o, self.N)
        N = min(self.N, o.N); v = min(self.v, o.v)
        return Ser(v, [self.coef(k) + o.coef(k) for k in range(v, N)], N) if v < N else Ser(N, [], N)
    __radd__ = __add__

    def __neg__(self):
        return Ser(self.v, [-x for x in self.c], self.N)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, a):
        return Ser(self.v, [a * x for x in self.c], self.N)

    def __mul__(self, o):
        if not isinstance(o, Ser):
            return self.scale(Fr(o))
        if not self.c or not o.c:
            return Ser(0, [], min(self.N + o.v, o.N + self.v))
        N = min(self.N + o.v, o.N + self.v); v = self.v + o.v
        if v >= N:
            return Ser(N, [], N)
        out = [Fr(0)] * (N - v)
        for i, x in enumerate(self.c):
            if x == 0:
                continue
            for j, y in enumerate(o.c):
                if i + j >= len(out):
                    break
                out[i + j] += x * y
        return Ser(v, out, N)
    __rmul__ = __mul__

    def inv(self):
        assert self.c, "inverse of zero"
        n = self.N - self.v; a0 = self.c[0]; out = [1 / a0]
        for k in range(1, n):
            s = Fr(0)
            for j in range(1, min(k, len(self.c) - 1) + 1):
                s += self.c[j] * out[k - j]
            out.append(-s / a0)
        return Ser(-self.v, out, n - self.v)

    def __truediv__(self, o):
        return self * o.inv() if isinstance(o, Ser) else self.scale(1 / Fr(o))

    def d(self):
        """d/du"""
        return Ser(self.v - 1, [(self.v + i) * x for i, x in enumerate(self.c)], self.N - 1)

    def shift(self, k):
        return Ser(self.v + k, self.c, self.N + k)

    def log1(self):
        """log of a power series with constant term 1"""
        assert self.v == 0 and self.c[0] == 1
        f = self.d() * self.inv()
        return Ser(0, [Fr(0)] + [f.coef(k) / (k + 1) for k in range(0, f.N)], f.N + 1)

    def exp0(self):
        """exp of a power series with zero constant term"""
        assert self.v >= 1 or not self.c
        N = self.N; out = [Fr(1)]
        for k in range(1, N):
            out.append(sum(j * self.coef(j) * out[k - j] for j in range(1, k + 1)) / k)
        return Ser(0, out, N)

    def powq(self, q):
        """f^q for f = 1 + O(u), q rational"""
        return self.log1().scale(Fr(q)).exp0()

    def resid(self):
        """coefficient of u^-1"""
        return self.coef(-1) if self.v <= -1 else Fr(0)


def D_op(f):
    """D = q d/dq = -(1-u) d/du on Laurent polynomials given as dict {exponent: coefficient}"""
    out = {}
    for n, a in f.items():
        if n != 0:
            out[n - 1] = out.get(n - 1, 0) - n * a
            out[n] = out.get(n, 0) + n * a
    return {k: x for k, x in out.items() if x != 0}


def mul1mu(f):
    out = {}
    for n, a in f.items():
        out[n] = out.get(n, 0) + a
        out[n + 1] = out.get(n + 1, 0) - a
    return {k: x for k, x in out.items() if x != 0}


def add(*fs):
    out = {}
    for f in fs:
        for k, x in f.items():
            out[k] = out.get(k, 0) + x
    return {k: x for k, x in out.items() if x != 0}


def sc(a, f):
    return {k: a * x for k, x in f.items() if a * x != 0}


def PF4_mono(n):
    """PF4 u^n as a dict"""
    f1 = D_op({n: Fr(1)}); f2 = D_op(f1); f3 = D_op(f2); f4 = D_op(f3)
    top = {k + 1: x for k, x in f4.items()}
    rest = add(sc(Fr(2), f3), sc(Fr(7, 5), f2), sc(Fr(2, 5), f1), sc(Fr(24, 625), {n: Fr(1)}))
    return add(top, sc(Fr(-1), mul1mu(rest)))


def frobenius(N):
    """f_0, f_1, f_2 (f_k = u^k + O(u^3)) as coefficient lists of length N; returns also the check that u^-3..u^-1 vanish"""
    pm = [PF4_mono(n) for n in range(N + 1)]
    basis = []
    for k in range(3):
        a = [Fr(0)] * N; a[k] = Fr(1)
        for n in range(3, N):
            s = sum(a[m] * pm[m].get(n - 3, 0) for m in range(max(0, n - 4), n))
            piv = pm[n].get(n - 3, 0)
            assert piv == n * (n - 1) ** 2 * (n - 2), "pivot"
            a[n] = -s / piv
        low = [sum(a[m] * pm[m].get(e, 0) for m in range(0, 4)) for e in (-3, -2, -1)]
        assert all(x == 0 for x in low), "power series obstruction"
        basis.append(a)
    return basis


def log_obstruction(basis):
    """[u^-3..u^-1] of PF4(f l) - l PF4(f) for f in basis (D l = 1 - 1/u): the matrix of S2(iii)"""
    def Dpair(A, B):  # D(A + B l) = DA + B(1 - 1/u) + (DB) l
        return add(D_op(A), B, sc(Fr(-1), {k - 1: x for k, x in B.items()})), D_op(B)
    cols = []
    for a in basis:
        B = {n: x for n, x in enumerate(a[:8]) if x != 0}; A = {}
        der = [(A, B)]
        for _ in range(4):
            der.append(Dpair(*der[-1]))
        top = {k + 1: x for k, x in der[4][0].items()}
        rest = add(sc(Fr(2), der[3][0]), sc(Fr(7, 5), der[2][0]), sc(Fr(2, 5), der[1][0]), sc(Fr(24, 625), der[0][0]))
        r = add(top, sc(Fr(-1), mul1mu(rest)))
        cols.append([r.get(e, Fr(0)) for e in (-3, -2, -1)])
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def systems(N, alpha=0, beta=0, gamma=Fr(7, 10)):
    """quintic b and pulled-back model bt, as lists of Ser known mod u^(N-2) or so; plus periods"""
    f0, f1, f2 = [Ser(0, a, N) for a in frobenius(N)]
    vp1 = f1 + f2.scale(gamma)
    var0 = f0 + vp1.scale(Fr(alpha)) + f2.scale(Fr(beta))
    tc = vp1 * var0.inv()
    one_m = Ser(0, [1, -1], N)
    u = Ser.mono(1, N)
    tcp = tc.d()
    b1 = -(u * var0 * var0 * one_m * one_m * one_m * tcp * tcp).inv()
    b3 = -(var0 * one_m.inv())
    b = [tcp, b1, tcp, b3, b3]
    tau = Ser(0, [Fr(4, 5 * n + 4) for n in range(N)], N).inv()
    v = tau.powq(Fr(-5, 4)).shift(1)
    vp = v.d()
    bt = [vp, -(vp * v.inv()), vp, -vp, -vp]
    return dict(b=b, bt=bt, tau=tau, v=v, var0=var0, tc=tc, vp1=vp1)


def prod(lst):
    out = lst[0]
    for x in lst[1:]:
        out = out * x
    return out
