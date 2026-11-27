# -*- coding: utf-8 -*-
"""rf_tr.py -- referee's own direct topological recursion, EXACT over Q, for y = z and x a Laurent polynomial in z whose
critical points (zeros of x') are all simple and RATIONAL (so no number field is needed).  Written from scratch.

Convention (standard EO, arXiv:math-ph/0702045): omega_{0,1} = y dx, omega_{0,2} = B = dz dz'/(z - z')^2,
  omega_{g,n+1}(z0, J) = Sum_a Res_{z->a} K_a(z0, z) [omega_{g-1,n+2}(z, s z, J) + Sum' omega(z, I) omega(s z, J minus I)],
  K_a(z0, z) = (1/2) (Int_{s z}^{z} B(z0, .)) / ((y(z) - y(s z)) dx(z)),  s = local deck involution at a.
Representation (2g-2+n > 0): omega = Sum C[key] prod_i dz_i/(z_i - a_{j_i})^(k_i + 1), key = ((j_1,k_1),..), k_i >= 1.
F_g (g >= 2) by the dilaton formula F_g = (1/(2-2g)) Sum_a Res_a Phi omega_{g,1}, dPhi = y dx.
"""
import itertools
from fractions import Fraction as Fr
from rf_ser import (gbinom, pmul, pinv, lmul, coeff, acc_add, acc_to_laurent, taylor_laurentpoly, lp_deriv)


def lp_eval(L, z):
    z = Fr(z)
    return sum((c * z ** e for e, c in L.items()), Fr(0))


class TR(object):
    def __init__(self, X, bps, kmax):
        """X: Laurent polynomial {e: c}; bps: ALL zeros of x' (rational, simple); kmax: largest k needed"""
        self.X = {e: Fr(c) for e, c in X.items() if c}
        self.a = [Fr(b) for b in bps]
        dX = lp_deriv(self.X)
        lo = min(dX)
        q = {e - lo: c for e, c in dX.items()}                   # polynomial, q(0) != 0
        assert len(set(self.a)) == len(self.a) and max(q) == len(self.a), "branch point list incomplete"
        for b in self.a:
            assert b != 0 and lp_eval(dX, b) == 0 and lp_eval(lp_deriv(dX), b) != 0
        self.kmax = kmax
        self.N = 2 * kmax + 10
        self.top = kmax + 3                       # absolute order kept for single-factor expansions
        self.loc = [self._local(b) for b in self.a]
        self.om = {}
        self._cp, self._cs = {}, {}

    def _local(self, b):
        N = self.N
        xs = taylor_laurentpoly(self.X, b, N + 3)
        assert xs[1] == 0 and xs[2] != 0
        u = [xs[k + 2] / xs[2] for k in range(N + 1)]            # (x(b+t) - x(b))/(xs2 t^2)
        sq = [Fr(1)] + [Fr(0)] * N                                # sqrt(u)
        for k in range(1, N + 1):
            sq[k] = (u[k] - sum((sq[i] * sq[k - i] for i in range(1, k)), Fr(0))) / 2
        zeta = [Fr(0)] + sq[:N]                                   # zeta(t) = t sqrt(u), zeta^2 xs2 = x - x(b)
        # reversion T(w): zeta(T(w)) = w, Lagrange: [w^n] T = (1/n) [t^(n-1)] (t/zeta)^n
        h = pinv(sq, N)
        T = [Fr(0)] * (N + 1)
        hp = [Fr(1)] + [Fr(0)] * N
        for n in range(1, N + 1):
            hp = pmul(hp, h, N)
            T[n] = hp[n - 1] / n
        mz = [-c for c in zeta]
        s = [Fr(0)] * (N + 1)                                     # s(t) = T(-zeta(t)), Horner
        for c in reversed(T):
            s = pmul(s, mz, N)
            s[0] += c
        assert s[0] == 0 and s[1] == -1
        # check x(b + s(t)) = x(b + t) mod t^(N+1)
        F = [Fr(0)] * (N + 1)
        for c in reversed(xs[: N + 1]):
            F = pmul(F, s, N)
            F[0] += c
        assert F == xs[: N + 1], "involution check failed"
        ds = [k * s[k] for k in range(1, N + 1)] + [Fr(0)]
        su = s[1:] + [Fr(0)]                                      # s = t su
        spow = [[Fr(1)] + [Fr(0)] * N]
        for n in range(1, N + 1):
            spow.append(pmul(spow[-1], s, N))
        dxu = [(k + 2) * xs[k + 2] for k in range(N + 1)]        # x'(b+t) = t dxu(t)
        inv2 = [c / 2 for c in pinv(dxu, N)]                      # 1/(2x') = t^-1 inv2
        tms = [Fr(1) - su[0]] + [-c for c in su[1:]]              # (t - s)/t
        w = pinv(tms, N)
        Bss = (-2, pmul(pmul(w, w, N), ds, N))                    # s'/(t-s)^2
        spds = [pmul(spow[m], ds, self.top) for m in range(self.top + 1)]
        return dict(b=b, s=s, ds=ds, su=su, spow=spow, spds=spds, inv2=inv2, Bss=Bss, xs=xs, sun={})

    def _su_neg(self, L, k):
        if k not in L["sun"]:
            L["sun"][k] = pinv(L["su"], self.N) if k == 1 else pmul(self._su_neg(L, k - 1), self._su_neg(L, 1), self.N)
        return L["sun"][k]

    def e_t(self, j0, i, k):
        """dz/(z - a_i)^(k+1) at z = a_j0 + t, per dt, Laurent to exponent self.top"""
        key = (j0, i, k)
        if key not in self._cp:
            if i == j0:
                self._cp[key] = (-(k + 1), [Fr(1)])
            else:
                d = self.a[j0] - self.a[i]
                self._cp[key] = (0, [gbinom(-(k + 1), n) * d ** (-(k + 1) - n) for n in range(self.top + 1)])
        return self._cp[key]

    def e_s(self, j0, i, k):
        """dz/(z - a_i)^(k+1) at z = a_j0 + s(t), per dt (includes s'), Laurent to exponent self.top"""
        key = (j0, i, k)
        if key not in self._cs:
            L = self.loc[j0]
            n = self.top + k + 2
            if i == j0:
                ser = pmul(L["ds"], self._su_neg(L, k + 1), self.N)
                self._cs[key] = (-(k + 1), ser[: n + 1])
            else:
                d = self.a[j0] - self.a[i]
                acc = [Fr(0)] * (self.top + 1)
                for m in range(self.top + 1):
                    cf = gbinom(-(k + 1), m) * d ** (-(k + 1) - m)
                    sp = L["spds"][m]
                    for e in range(self.top + 1):
                        acc[e] += cf * sp[e]
                self._cs[key] = (0, acc)
        return self._cs[key]

    def expand(self, j0, g, n, sig):
        """omega_{g,n}(first var at a_j0 + t or a_j0 + s(t); rest) as {rest_key: Laurent}"""
        L = self.loc[j0]
        out = {}
        if (g, n) == (0, 2):
            for m in range(self.top + 1):
                if sig:
                    out[((j0, m + 1),)] = (0, [(m + 1) * c for c in L["spds"][m]])
                else:
                    out[((j0, m + 1),)] = (m, [Fr(m + 1)])
            return out
        acc = {}
        for key, c in self.om[(g, n)].items():
            ser = self.e_s(j0, *key[0]) if sig else self.e_t(j0, *key[0])
            acc_add(acc.setdefault(key[1:], {}), ser, c)
        return {r: acc_to_laurent(d, self.top) for r, d in acc.items()}

    def step(self, g, n1):
        n = n1 - 1
        S = tuple(range(n))
        res = {}
        for j0, L in enumerate(self.loc):
            R = {}
            if g >= 1:
                if (g - 1, n + 2) == (0, 2):
                    acc_add(R.setdefault((), {}), (L["Bss"][0], L["Bss"][1][:3]))
                else:
                    for key, c in self.om[(g - 1, n + 2)].items():
                        acc_add(R.setdefault(key[2:], {}), lmul(self.e_t(j0, *key[0]), self.e_s(j0, *key[1]), 0), c)
            for g1 in range(g + 1):
                for sz in range(n + 1):
                    for I in itertools.combinations(S, sz):
                        J = tuple(x for x in S if x not in I)
                        if (g1, len(I)) == (0, 0) or (g - g1, len(J)) == (0, 0):
                            continue
                        E1 = self.expand(j0, g1, len(I) + 1, False)
                        E2 = self.expand(j0, g - g1, len(J) + 1, True)
                        for r1, s1 in E1.items():
                            for r2, s2 in E2.items():
                                rest = [None] * n
                                for p, v in enumerate(I):
                                    rest[v] = r1[p]
                                for p, v in enumerate(J):
                                    rest[v] = r2[p]
                                acc_add(R.setdefault(tuple(rest), {}), lmul(s1, s2, 0))
            for rest, d in R.items():
                Rl = acc_to_laurent(d, 0)
                if not Rl[1]:
                    continue
                mmax = 1 - Rl[0]
                A = [Fr(1)] + [Fr(0)] * self.N                     # A_m = (t^m - s^m)/(t - s)
                for m in range(1, mmax + 1):
                    if m > 1:
                        A = [Fr(0)] + A[:-1]
                        A = [x + y for x, y in zip(A, L["spow"][m - 1])]
                    kap = (-1, pmul(A, L["inv2"], self.N))
                    c = sum((a * coeff(kap, -1 - (Rl[0] + jj)) for jj, a in enumerate(Rl[1]) if a), Fr(0))
                    if c:
                        k2 = ((j0, m),) + rest
                        res[k2] = res.get(k2, Fr(0)) + c
        res = {k: c for k, c in res.items() if c}
        self.om[(g, n1)] = res
        return res

    def asym(self, g, n):
        """number of keys violating the symmetry of omega_{g,n}"""
        om = self.om[(g, n)]
        bad = 0
        for key, c in om.items():
            for p in itertools.permutations(range(n)):
                if om.get(tuple(key[i] for i in p), Fr(0)) != c:
                    bad += 1
        return bad

    def value(self, g, n, pts):
        tot = Fr(0)
        for key, c in self.om[(g, n)].items():
            t = c
            for (j, k), p in zip(key, pts):
                t *= Fr(1) / (Fr(p) - self.a[j]) ** (k + 1)
            tot += t
        return tot

    def pair(self, g, n, legs):
        """Sum Res .. Res prod Phi_i omega_{g,n}; legs = functions (a, K) -> [t^k] Phi(a + t), k = 0..K"""
        K = max(k for key in self.om[(g, n)] for (_, k) in key) if self.om[(g, n)] else 1
        tabs = [[leg(b, K) for b in self.a] for leg in legs]
        tot = Fr(0)
        for key, c in self.om[(g, n)].items():
            t = c
            for i, (j, k) in enumerate(key):
                t *= tabs[i][j][k]
            tot += t
        return tot
