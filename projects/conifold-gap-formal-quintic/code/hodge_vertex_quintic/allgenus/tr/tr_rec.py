# -*- coding: utf-8 -*-
"""tr_rec.py -- Eynard-Orantin topological recursion on a genus-0 curve with y = p (tr_core curves), exact residue calculus.

Convention: arXiv:1509.06954 (DNOPS) eqs. (2-26)-(2-29), quoted in TR_MODEL.md s.1:
  omega_{0,1} = y dx, omega_{0,2} = B = dp1 dp2/(p1-p2)^2,
  omega_{g,n+1}(p0, J) = Sum_b Res_{p -> c_b} [Int_p^{sigma_b(p)} B(., p0)] / [2 (omega01(sigma_b p) - omega01(p))]
                          x [ omega_{g-1,n+2}(p, sigma_b p, J) + Sum' omega_{g1,1+|I|}(p, I) omega_{g2,1+|J-I|}(sigma_b p, J-I) ],
  Sum' over g1 + g2 = g, I subset J, excluding (g1, I) = (0, {}) and (g2, J-I) = (0, {}) (omega_{0,2}-factors are included).
With p0 - c_b = T, p - c_b = t:  kernel = Sum_{m>=1} e_{b,m}(p0) K_m(t)/dt,  K_m = (t^m - sigma^m)/(2 (t - sigma) x'(c_b + t)).
Basis: omega_{g,n} (2g-2+n > 0) = Sum C[((j1,k1),..,(jn,kn))] prod e_{j,k}(p),  e_{j,k}(p) = dp/(p - c_j)^(k+1), k >= 1.
B(p, q) near p = c_b:  Sum_{m>=0} (m+1) t^m dt e_{b,m+1}(q)  (element ('R', m) of the first slot).
Symmetric curves (ModelCurve, nbp = N): residues at c_0 only, then the Z/N rotation (zeta = F.z(1) of order N)  C[((r,m), (j+r,k)..)] = zeta^(-r(m + Sum k)) C[((0,m), (j,k)..)]
(pull-back by p -> zeta^-r p; TR is equivariant since y dx and B are invariant).  full=True computes all branch points directly.
Checks: pole orders <= 6g-4+2n (two margin orders computed and asserted zero); symmetry in the arguments; dilaton; string.
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_core import smul, sinv, lmul, lcoef, binom


class TR(object):
    def __init__(self, curve, full=False, log=None):
        self.C = curve; self.F = curve.F; self.N = curve.N
        self.sym = curve.symmetric and not full
        self.bps = [0] if self.sym else list(range(curve.nbp))
        self.om = {}; self._S = {}; self._K = {}; self._L1 = {}; self._L2 = {}; self._Q = {}; self._rec = {}
        self.log = log
        self.T = self.N - 3                                                # truncation degree of Laurent products

    # ------------------------------------------------ local series
    def S(self, b, n):
        """sigma^n sigma' at c_b (Laurent series, n any integer)"""
        key = (b, n); r = self._S.get(key)
        if r is not None: return r
        F = self.F; L = self.C.local(b); N = self.N
        su = L.sigma[1:N + 1] + [F.zero]                                  # sigma = t su(t)
        base = su if n >= 0 else sinv(su, N, F)
        u = [F.one] + [F.zero] * N
        for _ in range(abs(n)): u = smul(u, base, N, F)
        r = (n, smul(u, L.dsig, N, F)); self._S[key] = r
        return r

    def K(self, b, m):
        key = (b, m); r = self._K.get(key)
        if r is not None: return r
        F = self.F; L = self.C.local(b); N = self.N; sig = L.sigma
        sm = [F.one] + [F.zero] * N
        for _ in range(m): sm = smul(sm, sig, N, F)
        num = [-c for c in sm]; num[m] = num[m] + F.one                       # t^m - sigma^m
        tms = [-c for c in sig]; tms[1] = tms[1] + F.one                      # t - sigma = t (2 + ..)
        den = smul(tms[1:] + [F.zero], L.dx[1:] + [F.zero], N, F)              # (t - sigma) x' = t^2 den
        den = [c * 2 for c in den]
        r = (-2, smul(num, sinv(den, N, F), N, F)); self._K[key] = r
        return r

    def L1(self, b, e):
        key = (b, e); r = self._L1.get(key)
        if r is not None: return r
        F = self.F; N = self.N
        if e[0] == 'R': r = (e[1], [F.one])
        else:
            j, k = e[1], e[2]
            if j == b: r = (-k - 1, [F.one])
            else:
                di = F.one / (self.C.c[b] - self.C.c[j]); base = _pow(di, k + 1, F); c = []
                for l in range(N + 1):
                    c.append(base * F.of(binom(-k - 1, l))); base = base * di
                r = (0, c)
        self._L1[key] = r
        return r

    def L2(self, b, e):
        """element pulled back by sigma_b:  f(sigma(t)) sigma'(t)"""
        key = (b, e); r = self._L2.get(key)
        if r is not None: return r
        F = self.F; v, c = self.L1(b, e); acc = None
        for k, x in enumerate(c):
            if x.iszero() or v + k > self.N: continue
            sv, sc = self.S(b, v + k)
            term = (sv, [x * y for y in sc])
            acc = term if acc is None else _ladd(acc, term, F)
        r = _ltrunc(acc, self.N); self._L2[key] = r
        return r

    def rec(self, b, m, e1, e2):
        key = (b, m, e1, e2); r = self._rec.get(key)
        if r is not None: return r
        F = self.F; qk = (b, m, e2); q = self._Q.get(qk)
        if q is None:
            q = lmul(self.K(b, m), self.L2(b, e2), self.T, F); self._Q[qk] = q
        v1, c1 = self.L1(b, e1); s = F.zero
        for k, x in enumerate(c1):
            if not x.iszero():
                y = lcoef(q, -1 - v1 - k, F)
                if not y.iszero(): s = s + x * y
        self._rec[key] = s
        return s

    # ------------------------------------------------ the recursion
    @staticmethod
    def kmax(g, n): return 6 * g - 5 + 2 * n

    def omega(self, g, n):
        key = (g, n); W = self.om.get(key)
        if W is None:
            W = self._compute(g, n); self.om[key] = W
        return W

    def _factor(self, b, g1, nI, slot, mcut):
        """terms (element of the p- or sigma(p)-slot, elements for the I-variables, coefficient) of omega_{g1, 1+nI}"""
        F = self.F
        if g1 == 0 and nI == 1:
            return [(('R', m), ((b, m + 1),), F.of(m + 1)) for m in range(mcut + 1)]
        W = self.omega(g1, 1 + nI)
        return [(('P',) + k[0], k[1:], c) for k, c in W.items()]

    def _compute(self, g, n):
        F = self.F; S = n - 1; kmx = self.kmax(g, n); ms = list(range(1, kmx + 3)); mcut = kmx + 2
        assert 2 * kmx + 6 <= self.N, "series order N = %d too small for (g,n) = (%d,%d)" % (self.N, g, n)
        W = {}
        for b in self.bps:
            out = {}
            def add(m, eS, val):
                if val.iszero(): return
                k = (m, eS); out[k] = out[k] + val if k in out else val
            if g >= 1:
                if (g - 1, n + 1) == (0, 2):                                   # B(p, sigma p) = sigma' dt^2/(t - sigma)^2
                    L = self.C.local(b); N = self.N; sig = L.sigma
                    u = [-c for c in sig[1:]] + [F.zero]; u[0] = u[0] + F.one        # (t - sigma)/t
                    G = (-2, smul(L.dsig, sinv(smul(u, u, N, F), N, F), N, F))
                    for m in ms: add(m, (), lcoef(lmul(self.K(b, m), G, 0, F), -1, F))
                else:
                    for k, c in self.omega(g - 1, n + 1).items():
                        e1 = ('P',) + k[0]; e2 = ('P',) + k[1]; eS = k[2:]
                        for m in ms: add(m, eS, self.rec(b, m, e1, e2) * c)
            idx = list(range(S))
            for g1 in range(g + 1):
                g2 = g - g1
                for r in range(S + 1):
                    for I in itertools.combinations(idx, r):
                        J = tuple(x for x in idx if x not in I)
                        if (g1 == 0 and len(I) == 0) or (g2 == 0 and len(J) == 0): continue
                        f1 = self._factor(b, g1, len(I), 1, mcut); f2 = self._factor(b, g2, len(J), 2, mcut)
                        for e1, eI, c1 in f1:
                            for e2, eJ, c2 in f2:
                                eS = [None] * S
                                for pos, x in zip(I, eI): eS[pos] = x
                                for pos, x in zip(J, eJ): eS[pos] = x
                                eS = tuple(eS); cc = c1 * c2
                                for m in ms: add(m, eS, self.rec(b, m, e1, e2) * cc)
            bad = [k for k, v in out.items() if (k[0] > kmx or any(e[1] > kmx for e in k[1])) and not v.iszero()]
            assert not bad, ("pole order exceeds 6g-4+2n", (g, n), bad[:4])
            for (m, eS), c in out.items():
                if c.iszero() or m > kmx: continue
                if self.sym:
                    tot = m + sum(e[1] for e in eS)
                    nb = self.C.nbp
                    for r in range(nb):
                        key = ((r, m),) + tuple(((j + r) % nb, k) for (j, k) in eS)
                        W[key] = c * F.z(-r * tot)
                else:
                    W[((b, m),) + eS] = c
        if self.log: self.log("omega_{%d,%d}: %d coefficients" % (g, n, len(W)))
        return W

    # ------------------------------------------------ checks and operators
    def symmetric_failures(self, g, n):
        W = self.omega(g, n); bad = 0
        for k, c in W.items():
            for perm in itertools.permutations(range(n)):
                o = W.get(tuple(k[i] for i in perm))
                if o is None or not (o - c).iszero(): bad += 1
        return bad

    def dilaton(self, W, slot=0):
        """Sum_b Res_{c_b} Phi(p) W(.., p, ..) in the given slot, dPhi = y dx:  Res_{c_j} Phi e_{j,k} = (1/k) [t^(k-1)] (y x')(c_j + t)"""
        F = self.F; out = {}; tab = {}
        for k, c in W.items():
            j, kk = k[slot]
            if (j, kk) not in tab: tab[(j, kk)] = self.C.ydx_taylor(j, kk)[kk - 1] * F.of(Fr(1, kk))
            key = k[:slot] + k[slot + 1:]; v = c * tab[(j, kk)]
            out[key] = out[key] + v if key in out else v
        return {k: v for k, v in out.items() if not v.iszero()}

    def string(self, W, slot=0):
        """Sum_b Res_{c_b} y(p) W(.., p, ..):  Res_{c_j} p e_{j,k} = delta_{k,1}"""
        out = {}
        for k, c in W.items():
            if k[slot][1] != 1: continue
            key = k[:slot] + k[slot + 1:]
            out[key] = out[key] + c if key in out else c
        return {k: v for k, v in out.items() if not v.iszero()}

    def Fg(self, g):
        """F_g = (1/(2-2g)) Sum_b Res_{c_b} Phi omega_{g,1}"""
        d = self.dilaton(self.omega(g, 1))
        return d.get((), self.F.zero) * self.F.of(Fr(1, 2 - 2 * g))


def _pow(a, e, F):
    r = F.one
    for _ in range(e): r = r * a
    return r


def _ladd(A, B, F):
    va, ca = A; vb, cb = B; v = min(va, vb); n = max(va + len(ca), vb + len(cb)) - v
    out = [F.zero] * n
    for k, x in enumerate(ca): out[va - v + k] = out[va - v + k] + x
    for k, x in enumerate(cb): out[vb - v + k] = out[vb - v + k] + x
    return (v, out)


def _ltrunc(A, T):
    v, c = A
    return (v, c[:max(T - v + 1, 0)])
