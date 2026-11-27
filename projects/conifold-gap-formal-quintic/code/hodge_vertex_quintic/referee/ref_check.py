# -*- coding: utf-8 -*-
"""Referee's independent recomputation (fresh code, Fractions + Q(zeta_5)) of the formal-quintic genus-2 graph sum in the
point frame, as Laurent series in s (s^5 = u = 1 - 3125 q, L = 1/s), via the conifold Frobenius basis.
Nothing is imported from the project."""
import sys, math, itertools
from fractions import Fraction as Fr

# ---------------------------------------------------------------- Q(zeta_5)
class Z5(object):
    __slots__ = ('c',)
    def __init__(self, c):
        a4 = c[4]
        self.c = (c[0] - a4, c[1] - a4, c[2] - a4, c[3] - a4, Fr(0))
    @staticmethod
    def of(x): return Z5((Fr(x), Fr(0), Fr(0), Fr(0), Fr(0)))
    @staticmethod
    def zeta(k):
        c = [Fr(0)] * 5; c[k % 5] = Fr(1); return Z5(c)
    def __add__(self, o): return Z5(tuple(a + b for a, b in zip(self.c, o.c)))
    def __sub__(self, o): return Z5(tuple(a - b for a, b in zip(self.c, o.c)))
    def __neg__(self): return Z5(tuple(-a for a in self.c))
    def __mul__(self, o):
        if isinstance(o, (int, Fr)): return Z5(tuple(a * o for a in self.c))
        out = [Fr(0)] * 5
        for i, a in enumerate(self.c):
            if a == 0: continue
            for j, b in enumerate(o.c):
                if b == 0: continue
                out[(i + j) % 5] += a * b
        return Z5(out)
    def iszero(self): return all(a == 0 for a in self.c)
    def __eq__(self, o): return (self - o).iszero()
    def conj(self, k):
        out = [Fr(0)] * 5
        for j, a in enumerate(self.c): out[(j * k) % 5] += a
        return Z5(out)
    def inv(self):
        pr = self.conj(2) * self.conj(3) * self.conj(4)
        n = self * pr
        assert all(n.c[j] == 0 for j in range(1, 5)), "norm not rational"
        return pr * (1 / n.c[0])
    def rational(self):
        if all(self.c[j] == 0 for j in range(1, 5)): return self.c[0]
        return None
    def __repr__(self):
        r = self.rational()
        if r is not None: return str(r)
        return "(" + " + ".join("%s*z^%d" % (a, j) for j, a in enumerate(self.c) if a != 0) + ")"

ONE = Z5.of(1); ZERO = Z5.of(0)
def zeta(k): return Z5.zeta(k)

# ---------------------------------------------------------------- Laurent series in one variable, coefficients Z5
BIG = 10 ** 6
class S(object):
    __slots__ = ('v', 'c', 'prec')
    def __init__(self, v, c, prec):
        # strip leading zeros, cap at prec
        c = list(c)
        while c and c[0].iszero(): c.pop(0); v += 1
        if prec < BIG:
            n = prec - v
            if n < 0: c = []; v = prec
            else: c = c[:n]
        while c and c[-1].iszero(): c.pop()
        if not c: v = prec
        self.v = v; self.c = c; self.prec = prec
    @staticmethod
    def const(x, prec=BIG): return S(0, [x if isinstance(x, Z5) else Z5.of(x)], prec)
    @staticmethod
    def mono(k, x=1, prec=BIG): return S(k, [x if isinstance(x, Z5) else Z5.of(x)], prec)
    @staticmethod
    def from_dict(d, prec=BIG):
        if not d: return S(0, [], prec)
        lo = min(d); hi = max(d)
        c = [ZERO] * (hi - lo + 1)
        for k, x in d.items(): c[k - lo] = x if isinstance(x, Z5) else Z5.of(x)
        return S(lo, c, prec)
    def coef(self, k):
        if k >= self.prec: raise ValueError("coefficient %d beyond precision %d" % (k, self.prec))
        i = k - self.v
        return self.c[i] if 0 <= i < len(self.c) else ZERO
    def iszero(self): return not self.c
    def __add__(self, o):
        prec = min(self.prec, o.prec)
        if not self.c: return S(o.v, o.c, prec)
        if not o.c: return S(self.v, self.c, prec)
        lo = min(self.v, o.v); hi = max(self.v + len(self.c), o.v + len(o.c))
        c = [ZERO] * (hi - lo)
        for i, x in enumerate(self.c): c[self.v + i - lo] = c[self.v + i - lo] + x
        for i, x in enumerate(o.c): c[o.v + i - lo] = c[o.v + i - lo] + x
        return S(lo, c, prec)
    def __neg__(self): return S(self.v, [-x for x in self.c], self.prec)
    def __sub__(self, o): return self + (-o)
    def scale(self, x):
        if not isinstance(x, Z5): x = Z5.of(x)
        return S(self.v, [a * x for a in self.c], self.prec)
    def __mul__(self, o):
        if isinstance(o, (int, Fr, Z5)): return self.scale(o)
        if not self.c or not o.c:
            prec = min(self.v + o.prec, o.v + self.prec) if (self.c or o.c) else min(self.prec, o.prec)
            return S(0, [], prec)
        prec = min(self.v + o.prec, o.v + self.prec)
        n = prec - (self.v + o.v)
        c = [ZERO] * max(n, 0)
        for i, a in enumerate(self.c):
            if a.iszero(): continue
            for j, b in enumerate(o.c):
                if i + j >= n: break
                c[i + j] = c[i + j] + a * b
        return S(self.v + o.v, c, prec)
    def shift(self, k): return S(self.v + k, self.c, self.prec + k)
    def deriv(self):
        c = [x * (self.v + i) for i, x in enumerate(self.c)]
        return S(self.v - 1, c, self.prec - 1)
    def inv(self):
        assert self.c, "inverse of zero"
        a0 = self.c[0].inv(); n = self.prec - self.v if self.prec < BIG else len(self.c) + 60
        if self.prec >= BIG: n = max(n, 60)
        b = [ZERO] * n; b[0] = a0
        for k in range(1, n):
            acc = ZERO
            for j in range(1, min(k, len(self.c) - 1) + 1):
                acc = acc + self.c[j] * b[k - j]
            b[k] = -(acc * a0)
        return S(-self.v, b, -self.v + n)
    def pow(self, k):
        if k < 0: return self.inv().pow(-k)
        r = S.const(1, BIG)
        for _ in range(k): r = r * self
        return r
    def trunc(self, prec): return S(self.v, self.c, min(prec, self.prec))
    def fmt(self, kmax=None, n=6):
        t = []
        for i, x in enumerate(self.c):
            k = self.v + i
            if x.iszero(): continue
            if kmax is not None and k > kmax: break
            t.append("(%s) x^%d" % (x, k));
            if len(t) >= n: break
        return (" + ".join(t) if t else "0") + " + O(x^%d)" % self.prec

def D_s(f):
    """D = q d/dq = -((1 - s^5)/5) s^-4 d/ds on an s-series"""
    one_m = S.from_dict({0: 1, 5: -1})
    return (one_m * f.deriv().shift(-4)).scale(Fr(-1, 5))

def D_u(f):
    """D = -(1 - u) d/du on a u-series"""
    return -(S.from_dict({0: 1, 1: -1}) * f.deriv())

def u_to_s(F):
    """u-series -> s-series, u = s^5"""
    c = []
    for x in F.c: c.extend([x, ZERO, ZERO, ZERO, ZERO])
    return S(5 * F.v, c, 5 * F.prec)

def s_to_u(F):
    assert all(x.iszero() or (F.v + i) % 5 == 0 for i, x in enumerate(F.c)), "not a u-series: " + F.fmt()
    lo = -((-F.v) // 5)
    c = [F.coef(5 * m) for m in range(lo, (F.prec + 4) // 5) if 5 * m < F.prec]
    return S(lo, c, (F.prec + 4) // 5)

# ---------------------------------------------------------------- Laurent polynomials in L (Fractions), WKB at colour 0
def lp_add(a, b):
    out = dict(a)
    for k, v in b.items():
        x = out.get(k, Fr(0)) + v
        if x: out[k] = x
        elif k in out: del out[k]
    return out
def lp_scale(a, c): return {k: v * c for k, v in a.items() if v * c}
def lp_mul(a, b):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            out[k1 + k2] = out.get(k1 + k2, Fr(0)) + v1 * v2
    return {k: v for k, v in out.items() if v}
def lp_D(a):
    """D L^k = (k/5)(L^(k+5) - L^k)"""
    out = {}
    for k, v in a.items():
        if k == 0: continue
        c = v * Fr(k, 5)
        out[k + 5] = out.get(k + 5, Fr(0)) + c
        out[k] = out.get(k, Fr(0)) - c
    return {k: v for k, v in out.items() if v}

def zs_apply_MpL(f):
    """(M + L) f for a z-series f = {n: LP}, M = zD"""
    out = {}
    for n, a in f.items():
        out[n + 1] = lp_add(out.get(n + 1, {}), lp_D(a))
        out[n] = lp_add(out.get(n, {}), lp_mul(a, {1: Fr(1)}))
    return {n: a for n, a in out.items() if a}

def zs_add(f, g, c=Fr(1), shift=0):
    out = dict(f)
    for n, a in g.items():
        out[n + shift] = lp_add(out.get(n + shift, {}), lp_scale(a, c))
    return {n: a for n, a in out.items() if a}

Q_LP = {0: Fr(1, 3125), -5: Fr(-1, 3125)}        # q = (1 - L^-5)/3125

def Pcal(f):
    """Pcal f = (M+L)^5 f - f - q[3125 (M+L)^5 + 9375 z (M+L)^4 + 10625 z^2 (M+L)^3 + 5625 z^3 (M+L)^2 + 1370 z^4 (M+L) + 120 z^5] f"""
    A = [f]
    for k in range(5): A.append(zs_apply_MpL(A[-1]))
    out = zs_add(A[5], f, Fr(-1))
    for k, c, zp in ((5, 3125, 0), (4, 9375, 1), (3, 10625, 2), (2, 5625, 3), (1, 1370, 4), (0, 120, 5)):
        t = {n + zp: lp_mul(lp_scale(a, Fr(-c)), Q_LP) for n, a in A[k].items()}
        out = zs_add(out, t)
    return out

def integrate_D(r):
    """phi with D phi = r (Laurent polynomial in L), phi(L^0 coefficient) = 0; residue checks.
    r = (L^5 - 1) s + rem; phi = 5 Sum_{k != 0} s_k L^k / k.  Returns (phi, rem, s_0)."""
    r = dict(r); s = {}
    if r:
        kmin = min(r); kmax = max(r)
        for k in range(kmax, kmin + 4, -1):
            v = r.get(k, Fr(0))
            if not v: continue
            s[k - 5] = v; r[k] = Fr(0); r[k - 5] = r.get(k - 5, Fr(0)) + v
    rem = {k: v for k, v in r.items() if v}
    phi = {k: 5 * v / k for k, v in s.items() if k != 0 and v}
    return phi, rem, s.get(0, Fr(0))

def power_sums_wt0(kmax):
    """p_k = Sum_a w_a^-k for wt_0 = {1 - zeta^j (j = 1..4), -5}, computed in Q(zeta) (must be rational)"""
    winv = [(ONE - zeta(j)).inv() for j in range(1, 5)] + [Z5.of(Fr(-1, 5))]
    ps = {}
    for k in range(1, kmax + 1):
        tot = ZERO
        for x in winv:
            y = ONE
            for _ in range(k): y = y * x
            tot = tot + y
        r = tot.rational(); assert r is not None
        ps[k] = r
    return ps

def bernoulli(n):
    A = [Fr(0)] * (n + 1)
    for m in range(n + 1):
        A[m] = Fr(1, m + 1)
        for j in range(m, 0, -1): A[j - 1] = j * (A[j - 1] - A[j])
    return A[0]

def h_series(ps, zmax, sign=1):
    """h(z)^sign = exp(-sign Sum_k B_2k/(2k(2k-1)) z^(2k-1) p_(2k-1))"""
    ex = [Fr(0)] * (zmax + 1)
    k = 1
    while 2 * k - 1 <= zmax:
        ex[2 * k - 1] = -sign * bernoulli(2 * k) / (2 * k * (2 * k - 1)) * ps[2 * k - 1]; k += 1
    out = [Fr(0)] * (zmax + 1); out[0] = Fr(1)
    for n in range(1, zmax + 1):
        out[n] = sum(k * ex[k] * out[n - k] for k in range(1, n + 1)) / n
    return out

def wkb_colour0(Nz, hcoef):
    """psi_n(L), n = 0..Nz, with psi_0 = L and phi_n(L=1) = h_n"""
    psi = [{1: Fr(1)}]
    # check Pcal_1 (L phi) = 5 L D phi on a test phi, and Pcal_1 psi_0 = 0
    test = {3: Fr(2), -2: Fr(7)}
    res = Pcal({0: lp_mul({1: Fr(1)}, test)})
    assert not res.get(0), "symbol equation fails"
    assert res.get(1) == lp_scale(lp_D(test), Fr(5)), "Pcal_1 (L phi) != 5 D phi: %s" % res.get(1)
    assert not Pcal({0: psi[0]}).get(1), "Pcal_1 psi_0 != 0"
    for n in range(1, Nz + 1):
        f = {m: psi[m] for m in range(n)}
        r = Pcal(f).get(n + 1, {})                      # residual at z^(n+1) without psi_n
        # 5 D phi_n + r = 0  ->  D phi_n = -r/5
        rhs = lp_scale(r, Fr(-1, 5))
        phi, rem, s0 = integrate_D(rhs)
        assert not rem and s0 == 0, "WKB residue at order %d: rem=%s s0=%s" % (n, rem, s0)
        val1 = sum(phi.values(), Fr(0))
        phi[0] = hcoef[n] - val1
        phi = {k: v for k, v in phi.items() if v}
        psi.append(lp_mul(phi, {1: Fr(1)}))
    # full check: Pcal(psi) = O(z^(Nz+1))
    chk = Pcal({m: psi[m] for m in range(Nz + 1)})
    assert all(not chk.get(k) for k in range(Nz + 1)), "WKB verification failed"
    return psi

# ---------------------------------------------------------------- conifold Frobenius basis (u-series)
# P = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u),  D = -(1-u) d/du
def PF_u(f):
    one_m = S.from_dict({0: 1, 1: -1})
    j = [f]
    for k in range(4): j.append(D_u(j[-1]))
    return (S.mono(1) * j[4]) - (one_m * j[3]).scale(2) - (one_m * j[2]).scale(Fr(7, 5)) - (one_m * j[1]).scale(Fr(2, 5)) - (one_m * j[0]).scale(Fr(24, 625))

class LogPair(object):
    """A log u + B"""
    def __init__(self, A, B): self.A = A; self.B = B
    def D(self):
        one_m = S.from_dict({0: 1, 1: -1})
        dA = -(one_m * self.A.deriv())
        dB = -(one_m * self.B.deriv()) - (one_m * self.A).shift(-1)        # D log u = -(1-u)/u
        return LogPair(dA, dB)
def PF_log(lp):
    one_m = S.from_dict({0: 1, 1: -1})
    j = [lp]
    for k in range(4): j.append(j[-1].D())
    pk = {4: S.mono(1), 3: one_m.scale(-2), 2: one_m.scale(Fr(-7, 5)), 1: one_m.scale(Fr(-2, 5)), 0: one_m.scale(Fr(-24, 625))}
    A = None; B = None
    for k in range(5):
        tA = pk[k] * j[k].A; tB = pk[k] * j[k].B
        A = tA if A is None else A + tA; B = tB if B is None else B + tB
    return A, B

def indicial(n):
    return Fr(n) * (n - 1) ** 2 * (n - 2)

def solve_analytic(k, N):
    b = [Fr(0)] * N; b[k] = Fr(1)
    for n in range(3, N):
        f = S(0, [Z5.of(x) for x in b[:n]], n + 10)
        r = PF_u(f).coef(n - 3).rational()
        b[n] = -r / indicial(n)
    f = S(0, [Z5.of(x) for x in b], N)
    chk = PF_u(f)
    assert all(chk.coef(m).iszero() for m in range(chk.v, N - 4)), "analytic solution check"
    return f

def solve_log(g, N):
    """h with PF(g log u + h) = 0, h = O(u^3); returns (h, obstruction coefficients u^-3..u^-1 of PF(g log u) B-part)"""
    A, K = PF_log(LogPair(g, S(0, [], N)))
    obst = [K.coef(m).rational() for m in range(-3, 0)]
    b = [Fr(0)] * N
    for n in range(3, N):
        f = S(0, [Z5.of(x) for x in b[:n]], n + 10)
        r = PF_u(f).coef(n - 3).rational() + K.coef(n - 3).rational()
        b[n] = -r / indicial(n)
    return S(0, [Z5.of(x) for x in b], N), obst

def polar_in_t(Fu, tu):
    """Fu Laurent series in u, tu = t(u) = u + ...; returns Fu as Laurent series in t (compositional inverse by iteration)"""
    N = min(Fu.prec, tu.prec)
    # reversion: u(t) = t - a2 u^2 - ... iterate
    t = S.mono(1, 1, N)
    u = t
    for _ in range(N + 2):
        # u = t - (tu(u) - u)
        comp = compose(tu, u, N)
        u = t - (comp - u)
        u = u.trunc(N)
    w = u.shift(-1)            # u/t
    winv = w.inv()
    res = S(0, [], N)
    for i, x in enumerate(Fu.c):
        k = Fu.v + i
        if x.iszero(): continue
        term = (w.pow(k) if k >= 0 else winv.pow(-k)).shift(k).scale(x)
        res = res + term
    return res

def compose(f, g, N):
    """f(g) for power series f (v >= 0), g with v >= 1"""
    assert f.v >= 0 and g.v >= 1
    res = S(0, [], N)
    for i, x in enumerate(f.c):
        k = f.v + i
        if x.iszero(): continue
        res = res + (g.pow(k).trunc(N) if k > 0 else S.const(1, N)).scale(x)
    return res.trunc(N)

# ---------------------------------------------------------------- WK numbers used (hand values, checked against the string/dilaton eqs)
WK = {
    (2, (4,)): Fr(1, 1152), (2, (2, 3)): Fr(29, 5760), (2, (2, 2, 2)): Fr(7, 240),
    (1, (1,)): Fr(1, 24), (1, (0, 2)): Fr(1, 24), (1, (1, 1)): Fr(1, 24), (1, (0, 1, 2)): Fr(1, 12), (1, (0, 0, 3)): Fr(1, 24),
    (1, (0, 0, 2, 2)): Fr(1, 6), (0, (0, 0, 0)): Fr(1), (0, (0, 0, 0, 1)): Fr(1), (0, (0, 0, 0, 0, 2)): Fr(1),
}

def main():
    out = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    # ---- QRR constants h_{wt_0}
    ps = power_sums_wt0(9)
    log("power sums p_k of 1/wt_0:", {k: str(v) for k, v in ps.items()})
    h = h_series(ps, 7)
    log("h_{wt_0}(z) coefficients:", [str(x) for x in h])
    assert h[:5] == [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(-269, 80000), Fr(2823, 6400000)], "h_wt0 differs from the log"
    log("PASS h_{wt_0} = 1, -3/20, 9/800, -269/80000, 2823/6400000, ... (matches rows_check.log)")
    # ---- constant-map genus-2 term by hand: Int_{M_2} Hodge_wt = -e_3(x)/2880 - (p1 p2 - p3)/5760, x = 1/w
    e1 = ps[1]; e2 = (ps[1] ** 2 - ps[2]) / 2; e3 = (ps[1] ** 3 - 3 * ps[1] * ps[2] + 2 * ps[3]) / 6
    F20 = 5 * (-e3 / 2880 - (ps[1] * ps[2] - ps[3]) / 5760)
    log("F_2(0) = 5 * Int_{M_2} Hodge_{wt_0} =", F20, "(log: -19/9600)")
    assert F20 == Fr(-19, 9600)
    # genus 1: Sum_i zeta^i Int_{M_{1,1}} Hodge_{wt_i} = Sum_i zeta^i (-(1/24) p_1^{(i)}) with p_1^{(i)} = zeta^-i p_1 -> 5 * (-p1/24)
    H10 = 5 * (-ps[1] / 24)
    log("[Q^0] <H>_1 =", H10, "(log: -3/8)")
    assert H10 == Fr(-3, 8)
    # ---- WKB
    Nz = 4
    psi = wkb_colour0(Nz, h)
    for n in range(Nz + 1):
        log("psi_%d / L = %s" % (n, " + ".join("(%s) L^%d" % (v, k - 1) for k, v in sorted(psi[n].items()))))
    assert psi[1] == {5: Fr(-3, 20)} and psi[2] == {9: Fr(9, 800)}
    assert psi[3] == {3: Fr(67, 1250), 8: Fr(-126, 625), 13: Fr(11571, 80000)}
    log("PASS WKB psi_1..psi_3 match rows_check.log; top L-degree of psi_n/L = 4n:", [max(psi[n]) - 1 == 4 * n for n in range(Nz + 1)])
    # ---- conifold Frobenius basis
    NU = 12
    for n in range(8):
        Pu = PF_u(S.mono(n, 1, n + 12))
        lead = Pu.coef(n - 3).rational() if n >= 3 else None
        if n >= 3: assert lead == indicial(n), (n, lead)
        else: assert all(Pu.coef(k).iszero() for k in range(n - 3, 0))
    log("PASS indicial polynomial n (n-1)^2 (n-2): exponents {0, 1, 1, 2}")
    f = [solve_analytic(k, NU) for k in range(3)]
    for k in range(3): log("f_%d =" % k, f[k].fmt(5, 6))
    h1, ob1 = solve_log(f[1], NU); h2, ob2 = solve_log(f[2], NU)
    log("log obstructions (u^-3, u^-2, u^-1): f1:", ob1, " f2:", ob2)
    assert ob1[0] == ob1[1] == ob2[0] == ob2[1] == 0
    gamma = -ob1[2] / ob2[2]
    log("gamma =", gamma, "(claimed 7/10)"); assert gamma == Fr(7, 10)
    g = f[1] + f[2].scale(gamma)
    hlog, ob = solve_log(g, NU); assert all(x == 0 for x in ob)
    log("g = f1 + gamma f2 =", g.fmt(5, 6))
    tc = g * f[0].inv()
    log("t_c = g/f0 =", tc.fmt(4, 5))
    assert tc.coef(1).rational() == 1 and tc.coef(2).rational() == Fr(7, 10) and tc.coef(3).rational() == Fr(41, 75)
    log("PASS t_c = u + 7/10 u^2 + 41/75 u^3 + ...")
    # D on s-series vs D on u-series
    test_u = f[0];
    assert (u_to_s(D_u(test_u)) - D_s(u_to_s(test_u))).iszero(), "D_s != D_u"
    log("PASS D = -((1-s^5)/5) s^-4 d/ds agrees with D = -(1-u) d/du on f0")

    # ---- Cf images and rows, generic (alpha, beta)
    def run_frame(alpha, beta, tag):
        varpi0 = f[0] + g.scale(alpha) + f[2].scale(beta)
        tcu = g * varpi0.inv()
        C0 = u_to_s(varpi0); C1 = u_to_s(D_u(tcu)); C0inv = C0.inv(); C1inv = C1.inv()
        log("[%s] C0 = %s" % (tag, C0.fmt(10, 3)), " C1 = %s" % C1.fmt(10, 3))
        Ls = S.mono(-1)
        kinv = [C1inv, (C0 * C0 * C1 * C1 * C1).shift(5) * C1inv, C1inv, C0inv]       # 1/(C1 kappa_a), a = 0..3
        # X_0 = psi/C0 ; X_{a+1} = kinv[a] (D X_{a,n-1} + L X_{a,n}) ; P_a = (C0/L) X_a
        X = [S.from_dict({-k: v for k, v in psi[n].items()}) * C0inv for n in range(Nz + 1)]
        P = []
        for a in range(5):
            P.append([(C0 * S.mono(1)) * x for x in X])      # P_a = (C0/L) X_a = C0 s X_a
            if a == 4: break
            Xn = []
            for n in range(Nz + 1):
                t = Ls * X[n]
                if n >= 1: t = t + D_s(X[n - 1])
                Xn.append(kinv[a] * t)
            X = Xn
        log("[%s] valuations of P_{a,n}^{(0)}:" % tag, [[P[a][n].v if P[a][n].c else None for n in range(Nz + 1)] for a in range(5)])
        log("[%s] -4n + c'_a bound (c' = 0,-1,3,2,1) satisfied:" % tag,
            all(P[a][n].v >= -4 * n + (0, -1, 3, 2, 1)[a] for a in range(5) for n in range(Nz + 1) if P[a][n].c))
        # colour i rows: P^{(i)}_{a,n} = zeta^{i(a-n)} P^{(0)}_{a,n}
        def Prow(i, a, n): return P[a][n].scale(zeta(i * (a - n)))
        Delta = {i: (C0 * C0).shift(2).scale(zeta(3 * i)) for i in range(5)}
        Dinv = {i: Delta[i].inv() for i in range(5)}
        eta_inv = {(0, 3): Fr(1, 5), (1, 2): Fr(1, 5), (2, 1): Fr(1, 5), (3, 0): Fr(1, 5), (4, 4): Fr(1, 5)}
        # symplecticity check Sum_ab P_a^{(i)}(z) eta^{ab} P_b^{(j)}(-z) = delta_ij Delta_i to z^Nz
        for (i, j) in ((0, 0), (0, 1), (1, 1), (2, 4)):
            for n in range(Nz + 1):
                acc = S(0, [], BIG)
                for (a, b), e in eta_inv.items():
                    for m in range(n + 1):
                        acc = acc + (Prow(i, a, m) * Prow(j, b, n - m)).scale(e * (1 if (n - m) % 2 == 0 else -1))
                tgt = Delta[i] if (i == j and n == 0) else S(0, [], BIG)
                d = acc - tgt
                assert d.iszero() or d.v >= min(acc.prec, tgt.prec) - 1, ("symplecticity", i, j, n, d.fmt())
        log("[%s] PASS symplecticity of the rows (pairs (0,0),(0,1),(1,1),(2,4), z^<=%d) as s-series" % (tag, Nz))
        def G(i, j, m, n):
            acc = S(0, [], BIG)
            for (a, b), e in eta_inv.items():
                acc = acc + (Prow(i, a, m) * Prow(j, b, n)).scale(e)
            acc = acc.scale(-((-1) ** (m + n)))
            if i == j and m == 0 and n == 0: acc = acc + Delta[i]
            return acc
        Ecache = {}
        def E(i, j, m, n):
            key = (i, j, m, n)
            if key not in Ecache:
                acc = S(0, [], BIG)
                for k in range(n + 1):
                    t = G(i, j, m + 1 + k, n - k)
                    acc = acc + (t if k % 2 == 0 else -t)
                Ecache[key] = acc
            return Ecache[key]
        def Dil(i, d): return P[0][d - 1].scale(zeta(-i * (d - 1)) * ((-1) ** d))
        log("[%s] valuations of E^{(00)}_{mn}, E^{(01)}_{mn} (m+n <= 2):" % tag,
            {(m, n): (E(0, 0, m, n).v, E(0, 1, m, n).v) for m in range(3) for n in range(3) if m + n <= 2},
            " claimed -4(m+n)-2 bound holds:", all(E(i, j, m, n).v >= -4 * (m + n) - 2 for i in range(5) for j in range(5) for m in range(3) for n in range(3) if m + n <= 2))
        log("[%s] valuations of Dil_0(d), d = 2..5:" % tag, [Dil(0, d).v for d in range(2, Nz + 2)])
        if tag == "a=0,b=0":
            lead = lambda f: f.c[0]
            log("  leading coefficients (coefficient of s^val): E^{(00)}:", {(m, n): str(lead(E(0, 0, m, n))) for m in range(3) for n in range(3) if m + n <= 2})
            log("  E^{(01)}:", {(m, n): str(lead(E(0, 1, m, n))) for m in range(3) for n in range(3) if m + n <= 2})
            log("  Dil_0(2..4):", [str(lead(Dil(0, d))) for d in (2, 3, 4)], " Delta_0 lead:", str(lead(Delta[0])))
            # hand arithmetic, graph 4 colouring (0,0): (1/2)(1/24)^2 [E11 + 2 Dil(2) E01 + Dil(2)^2 E00] at s^-10
            e11 = lead(E(0, 0, 1, 1)).rational(); e01 = lead(E(0, 0, 0, 1)).rational(); e00 = lead(E(0, 0, 0, 0)).rational(); d2 = lead(Dil(0, 2)).rational()
            hand = Fr(1, 2) * Fr(1, 24) ** 2 * (e11 + 2 * d2 * e01 + d2 ** 2 * e00)
            log("  graph 4 (0,0) s^-10 by hand from leading coefficients:", hand)
            # theta graph colour 0: (1/8) Delta^-1 [2 (E10 + E01) E00 + Dil(2) E00^2]
            e10 = lead(E(0, 0, 1, 0)).rational()
            hand3 = Fr(1, 8) * (2 * (e10 + e01) * e00 + d2 * e00 ** 2) / lead(Delta[0]).rational()
            log("  theta graph colour 0 s^-10 by hand from leading coefficients:", hand3)
        # graphs
        def W(g, d): return WK[(g, tuple(sorted(d)))]
        vals = {}
        # G1 vertex g = 2
        g1 = {}
        for i in range(5):
            t = Dil(i, 4).scale(W(2, (4,))) + (Dil(i, 2) * Dil(i, 3)).scale(W(2, (2, 3))) + (Dil(i, 2) * Dil(i, 2) * Dil(i, 2)).scale(W(2, (2, 2, 2)) / 6)
            g1[i] = Delta[i] * t
        vals[1] = g1
        # G2 g=1 with loop
        g2 = {}
        for i in range(5):
            t = (E(i, i, 2, 0) + E(i, i, 1, 1) + E(i, i, 0, 2)).scale(W(1, (0, 2)))
            t = t + (Dil(i, 2) * (E(i, i, 1, 0) + E(i, i, 0, 1))).scale(W(1, (0, 1, 2)))
            t = t + (Dil(i, 3) * E(i, i, 0, 0)).scale(W(1, (0, 0, 3)))
            t = t + (Dil(i, 2) * Dil(i, 2) * E(i, i, 0, 0)).scale(W(1, (0, 0, 2, 2)) / 2)
            g2[i] = t.scale(Fr(1, 2))
        vals[2] = g2
        # G3 theta
        g3 = {}
        for i in range(5):
            t = ((E(i, i, 1, 0) + E(i, i, 0, 1)) * E(i, i, 0, 0)).scale(2 * W(0, (0, 0, 0, 1)))
            t = t + (Dil(i, 2) * E(i, i, 0, 0) * E(i, i, 0, 0)).scale(W(0, (0, 0, 0, 0, 2)))
            g3[i] = (Dinv[i] * t).scale(Fr(1, 8))
        vals[3] = g3
        # G4 (1,1)
        g4 = {}
        for i in range(5):
            for j in range(5):
                t = E(i, j, 1, 1) + Dil(i, 2) * E(i, j, 0, 1) + Dil(j, 2) * E(i, j, 1, 0) + Dil(i, 2) * Dil(j, 2) * E(i, j, 0, 0)
                g4[(i, j)] = t.scale(W(1, (1,)) ** 2 / 2)
        vals[4] = g4
        if tag == "a=0,b=0":
            for key in ((0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 1), (1, 0)):
                log("  graph 4 colouring %s: %s" % (key, g4[key].fmt(-5, 3)))
            orb = g4[(0, 0)]
            for k in range(1, 5): orb = orb + g4[(0, k)]
            log("  graph 4: 5 * [sum over (0,k), k=0..4] s^-10 coefficient =", orb.coef(-10) * 5)
        # G5 (0,1): loop at the g=0 vertex (colour i), edge to the g=1 vertex (colour j)
        g5 = {}
        for i in range(5):
            for j in range(5):
                t = E(i, i, 0, 0) * (E(i, j, 0, 1) + Dil(j, 2) * E(i, j, 0, 0))
                g5[(i, j)] = (Dinv[i] * t).scale(W(1, (1,)) / 2)
        vals[5] = g5
        # G6 triple edge
        g6 = {}
        for i in range(5):
            for j in range(5):
                e = E(i, j, 0, 0)
                g6[(i, j)] = (Dinv[i] * Dinv[j] * e * e * e).scale(Fr(1, 12))
        vals[6] = g6
        # G7
        g7 = {}
        for i in range(5):
            for j in range(5):
                g7[(i, j)] = (Dinv[i] * Dinv[j] * E(i, i, 0, 0) * E(j, j, 0, 0) * E(i, j, 0, 0)).scale(Fr(1, 8))
        vals[7] = g7
        total = S(0, [], BIG)
        for gi in range(1, 8):
            acc = S(0, [], BIG)
            vs = []
            for key, v in vals[gi].items():
                acc = acc + v; vs.append(v.v)
            log("[%s] graph %d: colour-summed Cf series: %s   per-colouring valuations: min %d max %d" % (tag, gi, acc.fmt(10, 5), min(vs), max(vs)))
            total = total + acc
        log("[%s] TOTAL F_2^cf = %s" % (tag, total.fmt(20, 8)))
        Fu = s_to_u(total)
        T = polar_in_t(Fu, tcu)
        log("[%s] in t_c' = g/varpi0: (%s)/t^2 + (%s)/t + (%s) + (%s) t + ..." % (tag, T.coef(-2), T.coef(-1), T.coef(0), T.coef(1)))
        return total, vals, T

    tot0, vals0, T0 = run_frame(Fr(0), Fr(0), "a=0,b=0")
    # hand check of the gap from the printed series
    c10 = tot0.coef(-10).rational(); c5 = tot0.coef(-5).rational(); c0 = tot0.coef(0).rational()
    log("gap check by hand: [s^-5]/[s^-10] =", c5 / c10, " (= -2 gamma = -7/5 iff c' = 0)")
    if "--all" in sys.argv:
        tot1, _, T1 = run_frame(Fr(1), Fr(0), "a=1,b=0")
        tot2, _, T2 = run_frame(Fr(0), Fr(1), "a=0,b=1")
        tot3, _, T3 = run_frame(Fr(2, 3), Fr(-1, 7), "a=2/3,b=-1/7")
    open(sys.argv[1] if len(sys.argv) > 1 else "ref_check_out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
