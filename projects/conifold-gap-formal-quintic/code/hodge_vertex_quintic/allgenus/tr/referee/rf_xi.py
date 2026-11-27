# -*- coding: utf-8 -*-
"""rf_xi.py -- referee's xi-function basis at the special point x = kappa (p^4/4 + 1/p) (branch points c_j = w^-j), independent code.

Derivation (REFEREE.md s.3): with x - x(c_j) = -w_j^2/2 (arXiv:1509.06954 (2-30)-(2-35)) and x''(c_j) = 5 kappa c_j^2,
w_j = sqrt(Lam) c_j (p - c_j) + .., Lam = -5 kappa;  xi^j = Int B(w_j, .)/dw_j = xit^j/sqrt(Lam), xit^j = -1/(c_j (p - c_j));
Delta_j^(1/2) = sqrt(Lam) c_j.  Hence (TFT scaling e_j = sqrt(Lam) et_j, metric /Lam, Omega_g x Lam^(g-1), and R_TR(z) = M_0(z/kappa)):
     omega_{g,n} = Lam^(g-1) Sum kappa^(-D) <prod tau_{d_l}(et_{i_l})>^frozen prod d((d/dx)^(d_l) xit^{i_l}),  D = 3g-3+n-Sum d.
Rational functions vanishing at infinity with poles at the c_j only: {(j, k): a} = Sum a/(p - c_j)^k;  d/dx = (1/x') d/dp with
1/x' = p^2/(kappa (p^5 - 1)) = Sum_j r_j/(p - c_j),  r_j = 1/(5 kappa c_j^2).  Differentials: {(j, M): a} = Sum a dp/(p - c_j)^M.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)

NS = (1, -1, -1, 1, 1)


class XiSpecial(object):
    def __init__(self, F, kappa, dmax):
        self.F = F; self.kappa = Fr(kappa); self.dmax = dmax
        self.c = [F.pw(-j) for j in range(5)]; kap = F.of(self.kappa)
        self.r = [F.one / (F.of(5) * kap * cj * cj) for cj in self.c]
        p0 = F.of(Fr(3, 7)); lhs = p0 * p0 / (kap * (p0 * p0 * p0 * p0 * p0 - F.one))
        rhs = sum((self.r[j] / (p0 - self.c[j]) for j in range(5)), F.zero)
        assert F.iszero(lhs - rhs), "partial fractions of 1/x'"
        self.xit = {}; self.dxi = {}
        for i in range(5):
            f = {(i, 1): -(F.one / self.c[i])}
            for d in range(dmax + 1):
                self.xit[(i, d)] = f
                self.dxi[(i, d)] = {(j, k + 1): a * F.of(-k) for (j, k), a in f.items() if not F.iszero(a)}
                f = self.ddx(f)

    def ddx(self, f):
        F = self.F; c = self.c; out = {}
        def add(key, v):
            out[key] = out[key] + v if key in out else v
        for (j, k), a in f.items():                               # derivative: -k a/(p - c_j)^(k+1)
            m = k + 1; b = a * F.of(-k)
            for i in range(5):                                    # times r_i/(p - c_i)
                if i == j: add((j, m + 1), b * self.r[i]); continue
                dd = F.one / (c[i] - c[j])                        # 1/((p-a)(p-b)^m) = A/(p-a) - Sum_l 1/((a-b)^(m-l+1) (p-b)^l)
                pw = [F.one]
                for _ in range(m + 1): pw.append(pw[-1] * dd)
                add((i, 1), b * self.r[i] * pw[m])
                for l in range(1, m + 1): add((j, l), -(b * self.r[i] * pw[m - l + 1]))
        return {k: v for k, v in out.items() if not F.iszero(v)}

    def convert_vec(self, V):
        F = self.F; V = dict(V); out = {}
        for d in range(self.dmax, -1, -1):
            for i in range(5):
                a = V.get((i, 2 * d + 2))
                if a is None or F.iszero(a): continue
                X = self.dxi[(i, d)]; w = a / X[(i, 2 * d + 2)]; out[(i, d)] = w
                for key, x in X.items(): V[key] = V[key] - w * x if key in V else -(w * x)
        return out, {k: v for k, v in V.items() if not F.iszero(v)}

    def convert(self, W, n):
        """e-basis tensor -> coefficients of prod d xit^{i_l}_{d_l}; returns (dict, number of non-zero remainders)"""
        cur = dict(W); nbad = 0
        for slot in range(n):
            groups = {}
            for k, v in cur.items(): groups.setdefault(k[:slot] + k[slot + 1:], {})[k[slot]] = v
            new = {}
            for rest, V in groups.items():
                co, rem = self.convert_vec(V)
                nbad += 1 if rem else 0
                for key, v in co.items(): new[rest[:slot] + (key,) + rest[slot:]] = v
            cur = new
        return {k: v for k, v in cur.items() if not self.F.iszero(v)}, nbad

    def frozen(self, W, g, n, Lam=None, kap=None):
        """frozen colour correlators {((i,d),..): value} = Lam^(1-g) kappa^D x coefficient"""
        F = self.F; Lam = Fr(-5) * self.kappa if Lam is None else Fr(Lam); kap = self.kappa if kap is None else Fr(kap)
        T, nbad = self.convert(W, n); out = {}
        for key, v in T.items():
            D = 3 * g - 3 + n - sum(d for (_, d) in key)
            out[key] = v * F.of(Lam ** (1 - g) * kap ** D)
        return out, nbad

    def flat(self, col, avec, dvec):
        """<prod tau_{d_l}(H^{a_l})> from colour correlators, H^a = N_a Sum_i zeta^(i(a+1)) et_i (zeta = F.w)"""
        F = self.F; tot = F.zero; n = len(avec)
        import itertools
        for idx in itertools.product(range(5), repeat=n):
            v = col.get(tuple((i, d) for i, d in zip(idx, dvec)))
            if v is None: continue
            w = F.one
            for i, a in zip(idx, avec): w = w * F.of(NS[a]) * F.pw(i * (a + 1))
            tot = tot + w * v
        return tot
