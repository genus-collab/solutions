# -*- coding: utf-8 -*-
"""rf_tr.py -- referee's own topological recursion engine (independent of tr_rec.py), genus-0 curves with y = p, B = dp1dp2/(p1-p2)^2.

Route: at every branch point a (no symmetry used) the local coordinate zeta is defined by x(a + t) - x(a) = l2 zeta^2, zeta = t + O(t^2)
(t(zeta) by series inversion), so the deck transformation is zeta -> -zeta.  Kernel (arXiv:1509.06954 (2-28)):
  K = (1/2) Int_{-zeta}^{zeta} B(p0, .) / ((y(zeta) - y(-zeta)) dx) = Sum_{M>=2} e_{a,M}(p0) kap_M(zeta)/dzeta,
  kap_M = odd(t^(M-1)) / (4 l2 zeta odd(t)),  e_{j,M}(p) = dp/(p - c_j)^M  (M = pole order >= 2).
omega_{g,n} is stored as {((j1,M1),..,(jn,Mn)): coefficient of prod e_{j,M}}.  Pull-backs to zeta at a: e_{j,M} -> t'/(a - c_j + t)^M,
by the deck map: -(that)(-zeta).  B(q, p_i) = Sum_M (M-1) t^(M-2) t' dzeta e_{a,M}(p_i);  B(q, sigma q) = -t'(z) t'(-z)/(t(z) - t(-z))^2.
F_g = (1/(2-2g)) Sum_a Res Phi omega_{g,1}, dPhi = y dx:  Res_{c_a} Phi e_{a,M} = (1/(M-1)) [t^(M-2)] (y x')(a + t).
Checks built in: every residue asserts sufficient precision; pole orders > 6g-4+2n must vanish (two extra orders computed).
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import Ls, ps_power


def Mmax(g, n): return 6 * g - 4 + 2 * n


class Curve(object):
    """F: field; c: branch points; xt(a, N): Taylor coefficients of x(a+t), t^0..t^N; yxt(a, N): of y x' (y = p)"""
    def __init__(self, F, c, xt, yxt, name):
        self.F = F; self.c = c; self.xt = xt; self.yxt = yxt; self.name = name


class Local(object):
    def __init__(self, C, ia, N):
        F = C.F; self.F = F; a = C.c[ia]; self.a = a; self.ia = ia; self.C = C
        l = C.xt(a, N + 4); assert F.iszero(l[1]) and not F.iszero(l[2]), "not a simple branch point"
        self.l2 = l[2]; h = [x / l[2] for x in l[2:N + 4]]               # h(t) = (x(a+t) - x(a))/(l2 t^2)
        phi = ps_power(h, Fr(-1, 2), N + 1, F)                            # zeta = t h^(1/2)  =>  t = zeta phi(t), phi = h^(-1/2)
        T = [F.zero]; P = [F.one] + [F.zero] * (N + 1)
        for n in range(1, N + 2):                                         # [zeta^n] t = (1/n) [t^(n-1)] phi^n
            P = [sum((P[i] * phi[k - i] for i in range(k + 1)), F.zero) for k in range(N + 2)]
            T.append(P[n - 1] / F.of(n))
        self.u = Ls(0, T[1:N + 2], N)                                     # t = zeta u(zeta)
        self.t = Ls(1, T[1:N + 2], N + 1)
        self.tp = self.t.deriv(F)                                         # t'(zeta), valuation 0, hi N
        odd = Ls(1, [x if (1 + k) % 2 == 1 else F.zero for k, x in enumerate(self.t.c)], self.t.hi)
        den = Ls(2, odd.c, odd.hi + 1).scale(F.of(4) * self.l2)           # 4 l2 zeta odd(t)
        self.dinv = den.inv(F)
        self.uinv = self.u.inv(F)
        self._pw = {}; self._E = {}; self._kap = {}; self._q = {}
        two_odd = Ls(1, [x * 2 for x in odd.c], odd.hi)
        self.Bqs = self.tp.mul(self.tp.flip(), F).scale(-F.one).mul(two_odd.mul(two_odd, F).inv(F), F)

    def upow(self, k):
        """u^k (k may be negative)"""
        key = ('u', k); r = self._pw.get(key)
        if r is None:
            if k == 0: r = Ls(0, [self.F.one] + [self.F.zero] * self.u.hi, self.u.hi)
            elif k > 0: r = self.upow(k - 1).mul(self.u, self.F)
            else: r = self.upow(k + 1).mul(self.uinv, self.F)
            self._pw[key] = r
        return r

    def wpow(self, j, k):
        """(a - c_j + t)^(-k), j != ia"""
        key = ('w', j, k); r = self._pw.get(key)
        if r is None:
            F = self.F
            if k == 1:
                d = self.a - self.C.c[j]
                r = Ls(0, [d] + self.t.c, self.t.hi).inv(F)
            else: r = self.wpow(j, k - 1).mul(self.wpow(j, 1), F)
            self._pw[key] = r
        return r

    def E(self, e):
        """pull-back of an element at zeta (as a series in zeta, dzeta stripped): ('P', j, M) or ('B', M)"""
        r = self._E.get(e)
        if r is None:
            F = self.F
            if e[0] == 'B':
                M = e[1]; tm = self.upow(M - 2); tm = Ls(M - 2, tm.c, tm.hi + M - 2)
                r = tm.mul(self.tp, F).scale(F.of(M - 1))
            else:
                j, M = e[1], e[2]
                if j == self.ia:
                    um = self.upow(-M); r = Ls(-M, um.c, um.hi - M).mul(self.tp, F)
                else: r = self.wpow(j, M).mul(self.tp, F)
            self._E[e] = r
        return r

    def Em(self, e):
        """pull-back by the deck transformation: -E(e)(-zeta)"""
        key = ('m',) + e; r = self._E.get(key)
        if r is None: r = self.E(e).flip().scale(-self.F.one); self._E[key] = r
        return r

    def kap(self, M):
        r = self._kap.get(M)
        if r is None:
            F = self.F; tm = self.upow(M - 1); tm = Ls(M - 1, tm.c, tm.hi + M - 1)
            odd = Ls(tm.v, [x if (tm.v + k) % 2 == 1 else F.zero for k, x in enumerate(tm.c)], tm.hi)
            r = odd.mul(self.dinv, F); self._kap[M] = r
        return r

    def res(self, M0, e1, e2):
        """[zeta^-1] kap_M0 E(e1) Em(e2); e2 = None means the B(q, sigma q) term (e1 ignored)"""
        F = self.F
        if e2 is None: return _res(self.kap(M0).mul(self.Bqs, F), None, F)
        key = (M0, e2); q = self._q.get(key)
        if q is None: q = self.kap(M0).mul(self.Em(e2), F); self._q[key] = q
        return _res(self.E(e1), q, F)


def _res(A, B, F):
    if B is None:
        return A.coef(-1, F) if A.v <= -1 else F.zero
    lo = A.v; hi = -1 - B.v
    if hi < lo: return F.zero
    assert A.hi >= hi and B.hi >= -1 - lo, ("precision too low for residue", A.v, A.hi, B.v, B.hi)
    s = F.zero
    for d in range(lo, hi + 1):
        x = A.c[d - A.v]
        if not F.iszero(x): s = s + x * B.c[-1 - d - B.v]
    return s


class TR(object):
    def __init__(self, C, N):
        self.C = C; self.F = C.F; self.N = N
        self.loc = [Local(C, ia, N) for ia in range(len(C.c))]
        self.om = {}

    def omega(self, g, n):
        key = (g, n)
        if key not in self.om: self.om[key] = self._rec(g, n - 1)
        return self.om[key]

    def _factor(self, ia, g1, k, cap):
        if g1 == 0 and k == 1:
            return [(('B', M), ((ia, M),), self.F.one) for M in range(2, cap + 1)]
        return [(('P',) + key[0], key[1:], c) for key, c in self.omega(g1, k + 1).items()]

    def _rec(self, g, n):
        """omega_{g,n+1}(p0, p_1..p_n)"""
        F = self.F; mtop = Mmax(g, n + 1); cap = mtop + 2; out = {}
        def add(key, v):
            if not F.iszero(v): out[key] = out[key] + v if key in out else v
        idx = list(range(n))
        for ia, L in enumerate(self.loc):
            if g >= 1:
                if g == 1 and n == 0:
                    for M0 in range(2, cap + 1): add(((ia, M0),), L.res(M0, None, None))
                else:
                    for key, c in self.omega(g - 1, n + 2).items():
                        e1 = ('P',) + key[0]; e2 = ('P',) + key[1]
                        for M0 in range(2, cap + 1): add(((ia, M0),) + key[2:], c * L.res(M0, e1, e2))
            for g1 in range(g + 1):
                g2 = g - g1
                for r in range(n + 1):
                    for I in itertools.combinations(idx, r):
                        J = [x for x in idx if x not in I]
                        if (g1 == 0 and r == 0) or (g2 == 0 and len(J) == 0): continue
                        for e1, k1, c1 in self._factor(ia, g1, r, cap):
                            for e2, k2, c2 in self._factor(ia, g2, len(J), cap):
                                rest = [None] * n
                                for pos, x in zip(I, k1): rest[pos] = x
                                for pos, x in zip(J, k2): rest[pos] = x
                                rest = tuple(rest); cc = c1 * c2
                                for M0 in range(2, cap + 1): add(((ia, M0),) + rest, cc * L.res(M0, e1, e2))
        bad = [k for k, v in out.items() if any(e[1] > mtop for e in k)]
        assert not bad, ("pole order above 6g-4+2n", g, n + 1, bad[:3])
        return out

    def phi_res(self, ia, M):
        """Res_{p = c_ia} Phi(p) dp/(p - c_ia)^M = (1/(M-1)) [t^(M-2)] (y x')"""
        key = ('phi', ia, M)
        if key not in self.om:
            yx = self.C.yxt(self.C.c[ia], M); self.om[key] = yx[M - 2] / self.F.of(M - 1)
        return self.om[key]

    def dilaton(self, W):
        out = {}
        for key, c in W.items():
            v = c * self.phi_res(*key[0]); rest = key[1:]
            out[rest] = out[rest] + v if rest in out else v
        return {k: v for k, v in out.items() if not self.F.iszero(v)}

    def Fg(self, g):
        return self.dilaton(self.omega(g, 1)).get((), self.F.zero) / self.F.of(2 - 2 * g)

    def xmoment(self, W):
        """Sum_a Res x(p) omega(p) for a one-point differential:  Res x e_{a,M} = [t^(M-1)] x(a + t)"""
        tot = self.F.zero
        for key, c in W.items():
            ia, M = key[0]; tot = tot + c * self.C.xt(self.C.c[ia], M)[M - 1]
        return tot

    def symmetric(self, g, n):
        W = self.omega(g, n); bad = 0
        for k, c in W.items():
            for perm in itertools.permutations(range(n)):
                o = W.get(tuple(k[i] for i in perm))
                if o is None or not self.F.iszero(o - c): bad += 1
        return bad
