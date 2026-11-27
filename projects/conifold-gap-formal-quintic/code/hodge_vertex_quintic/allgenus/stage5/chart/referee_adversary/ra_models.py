# -*- coding: utf-8 -*-
"""ra_models.py -- stage5/chart/referee_adversary: the integrands tested by ra_run.py (as specs for ra_sp.py).

Colour point sigma* = 1 (q = 0, L = 1) unless stated.  The quintic summand Gamma(1 + 5s) prod_(j=0..4) Gamma(zeta^j/z - s) (-q)^s,
sigma = z s, is:  Sing(alpha0 = -1) for Gamma((1 - sigma)/z);  group 'y' (alpha 5, M = 1) for Gamma(1 + 5 sigma/z);
group 'x' (alpha -1, four members a_j = zeta^j - sigma, M = 0) for j = 1..4.
Power sums of the group 'x':  route 'm' (exact, any field): the e_j = zeta^j - 1 are the roots of
m(s) = s^4 + 5 s^3 + 10 s^2 + 10 s + 5 = ((1+s)^5 - 1)/s, and Sum_r t^(r-1) Sum_j (e_j - eps)^(-r) = -m'(eps + t)/m(eps + t);
route 'zeta' (F_p only): explicit zeta in F_p, Sum_j (zeta^j - 1 - eps)^(-r) computed directly.  The two routes share nothing.
"""
from fractions import Fraction as Fr
from math import comb
from ra_core import FP, RSym, RNum
from ra_sp import Group, Sing

M_POLY = [Fr(5), Fr(10), Fr(10), Fr(5), Fr(1)]           # m(s), coefficients of s^0..s^4


def _div_series(F, num, den, D):
    """num/den as a power series to degree D (den[0] invertible), lists over F"""
    i0 = F.inv(den[0]); out = []
    for n in range(D + 1):
        s = num[n] if n < len(num) else F.zero
        for i in range(1, min(n, len(den) - 1) + 1):
            s -= den[i] * out[n - i]
        out.append(F.red(s * i0))
    return out


def _shift_poly(F, p, v):
    """coefficients of p(v + t) in t"""
    out = [F.zero] * len(p)
    for k, a in enumerate(p):
        for i in range(k + 1):
            out[i] = F.red(out[i] + a * comb(k, i) * (v ** (k - i) if not isinstance(F, FP) else pow(v, k - i, F.p)))
    return out


class XPsumM:
    """route 'm' for Sum_j (e_j - eps)^(-r):  RSym via G(s) = -m'(s)/m(s) around s = 0, RNum via the shifted quotient"""

    def __init__(self):
        self.cache = {}

    def G(self, F, v, D):
        key = (F.name, str(v), D)
        if key not in self.cache:
            vv = F.c(v)
            m = _shift_poly(F, [F.c(a) for a in M_POLY], vv)
            dm = [F.red(-(k + 1) * m[k + 1]) for k in range(len(m) - 1)]   # -m'(v + t) in t
            self.cache[key] = _div_series(F, dm, m, D)
        return self.cache[key]

    def pe(self, R, r):
        F = R.F
        if isinstance(R, RSym):
            G = self.G(F, 0, r + R.K + 2)
            return [F.red(G[r - 1 + e] * comb(r - 1 + e, r - 1)) for e in range(R.K)]
        G = self.G(F, R.v if isinstance(F, FP) else R.v, r + 2)
        return G[r - 1]

    def ps(self, F, r):
        return self.G(F, 0, r + 2)[r - 1]


class XPsumZeta:
    """route 'zeta' (F_p): explicit zeta^j"""

    def __init__(self, F, colour=0, members=None):
        self.F = F; self.z = F.zeta5(); self.colour = colour
        self.members = members if members is not None else [j for j in range(5) if j != colour]
        self.cache = {}

    def astar(self, j):
        F = self.F
        return F.red(pow(self.z, j, F.p) - pow(self.z, self.colour, F.p))

    def pe(self, R, r):
        key = (R.name, r)
        if key not in self.cache:
            tot = R.zero()
            for j in self.members:
                if isinstance(R, RSym):                # the R-element astar_j - eps
                    b = R.inv([self.astar(j)] + R.sub(R.zero(), R.epsp(1))[1:])
                else:
                    b = R.inv(R.F.red(self.astar(j) - R.v))
                p = R.const(1)
                for _ in range(r):
                    p = R.mul(p, b)
                tot = R.add(tot, p)
            self.cache[key] = tot
        return self.cache[key]

    def ps(self, F, r):
        return F.red(sum(pow(F.inv(self.astar(j)), r, F.p) for j in self.members))


def single(alpha, astar):
    """power sums for one factor a(sigma) with a(sigma*) = astar (rational), a(sigma* + eps) = astar + alpha eps"""
    cache = {}

    def pe(R, r):
        key = (R.name, r)
        if key not in cache:
            if isinstance(R, RSym):
                b = R.inv(R.add(R.const(astar), R.epsp(1, alpha)))
            else:
                b = R.inv(R.add(R.const(astar), R.scal(alpha, R.v)))
            p = R.const(1)
            for _ in range(r):
                p = R.mul(p, b)
            cache[key] = p
        return cache[key]

    def ps(F, r):
        return F.c(Fr(1) / Fr(astar) ** r)
    return pe, ps


def quintic(route='m', F=None, My=1, M0=0, Mx=0, d=5, alpha0=-1, sing_bk='std', y_bk='std', x_bk='std',
            sing_expo=Fr(-1, 2), y_expo=Fr(-1, 2), count=1, gstar=None, fextra=None, extra=None):
    """the (possibly modified) quintic summand at colour 0"""
    xp = XPsumM() if route == 'm' else XPsumZeta(F)
    groups = []
    if d != 0:
        pe, ps = single(d, d)
        groups.append(Group('y', d, pe, ps, M=My, bk=y_bk, expo=y_expo))
    groups.append(Group('x', -1, xp.pe, xp.ps, M=Mx, bk=x_bk))
    groups += (extra or [])
    return {'sing': Sing(alpha0=alpha0, M0=M0, bk=sing_bk, expo=sing_expo, count=count), 'groups': groups,
            'gstar': gstar or {}, 'fextra': fextra or {}}


def colour_i(F, i):
    """the quintic summand expanded at the colour-i point sigma* = zeta^i (F_p, explicit zeta); singular factor j = i"""
    z = F.zeta5(); zi = pow(z, i, F.p)
    xp = XPsumZeta(F, colour=i)
    cache = {}

    def pe(R, r):                                      # a = 5 sigma, a(sigma* + eps) = 5 zeta^i + 5 eps
        key = (R.name, r)
        if key not in cache:
            if isinstance(R, RSym):
                b = R.inv([F.red(5 * zi)] + R.epsp(1, 5)[1:])
            else:
                b = R.inv(F.red(5 * zi + 5 * R.v))
            p = R.const(1)
            for _ in range(r):
                p = R.mul(p, b)
            cache[key] = p
        return cache[key]

    def ps(F_, r):
        return pow(F.inv(5 * zi), r, F.p)
    # a_0 = zeta^i - sigma = -(sigma - zeta^i): alpha0 = -1 (the shift of the colour point does not enter the tau-data)
    return {'sing': Sing(alpha0=-1), 'groups': [Group('y', 5, pe, ps, M=1), Group('x', -1, xp.pe, xp.ps, M=0)],
            'gstar': {}, 'fextra': {}}


def rational_weights(lams, d=5, My=1):
    """Gamma(My + d sigma/z) prod_j Gamma((lam_j - sigma)/z) with lam_0 = 1 = sigma* singular, lam_j (j >= 1) rational"""
    groups = []
    if d != 0:
        pe, ps = single(d, d)
        groups.append(Group('y', d, pe, ps, M=My))
    for k, lam in enumerate(lams):
        pe, ps = single(-1, Fr(lam) - 1)
        groups.append(Group('x%d' % k, -1, pe, ps, M=0))
    return {'sing': Sing(alpha0=-1), 'groups': groups, 'gstar': {}, 'fextra': {}}


class XPsumOrb:
    """five members a_j = zeta^j - sigma at the orbifold point sigma* = 0: Sum_j 1/(zeta^j - s) = 5 s^4/(1 - s^5)"""

    def G(self, F, v, D):
        vv = F.c(v)
        num = _shift_poly(F, [F.zero] * 4 + [F.c(5)], vv)
        den = _shift_poly(F, [F.one] + [F.zero] * 4 + [F.c(-1)], vv)
        return _div_series(F, num, den, D)

    def pe(self, R, r):
        F = R.F
        if isinstance(R, RSym):
            G = self.G(F, 0, r + R.K + 2)
            return [F.red(G[r - 1 + e] * comb(r - 1 + e, r - 1)) for e in range(R.K)]
        return self.G(F, R.v, r + 2)[r - 1]

    def ps(self, F, r):
        return self.G(F, 0, r + 2)[r - 1]


def orbifold():
    """the quintic summand at the orbifold point L = 0 (q = oo): the SINGULAR factor is now Gamma(1 + 5 sigma/z)
    (alpha0 = 5, M0 = 1); all five Gamma((zeta^j - sigma)/z), j = 0..4, are regular there"""
    xo = XPsumOrb()
    return {'sing': Sing(alpha0=5, M0=1), 'groups': [Group('x5', -1, xo.pe, xo.ps, M=0)], 'gstar': {}, 'fextra': {}}


def core_lemma(alpha0, M0, gstar, fextra, count=1, bk='std'):
    """only the singular factor, an arbitrary extra phase g (no linear term) and an arbitrary amplitude F"""
    return {'sing': Sing(alpha0=alpha0, M0=M0, count=count, bk=bk), 'groups': [], 'gstar': gstar, 'fextra': fextra}
