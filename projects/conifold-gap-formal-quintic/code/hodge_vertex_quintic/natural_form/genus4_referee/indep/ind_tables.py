# -*- coding: utf-8 -*-
"""ind_tables.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): edge, leaf and vertex tables of the hatted formal quintic.

Input: only the from-scratch hatted rows of natural_form/referee (ref_fp, ref_rows, ref_hat: Qh[a][i][n] as u-series, xhat = C0,
t_c) and ind_wk.  No fq_* / fq3_* / fq4_* module is imported.
Definitions (point frame, colours = normalised canonical basis, eta^{ab} = (1/5) delta_{a+b = 3 mod 5}):
  N_ij(z, w) = Sum_ab Q_a^(i)(-z) eta^ab Q_b^(j)(-w),   E_R(z, w) = (1 - N)/(z + w),   E_R(a, b)[i][j] = [z^a w^b];
     the division by z + w is done by the recursion e(m, n) = d(m+1, n) - e(m+1, n-1) and its exactness (d(0, n) = e(0, n-1),
     i.e. symplecticity) is asserted;
  T_R(e)[i] = (-1)^e Q_{0, e-1}^(i)   (e >= 2; T(z) = z (1 - R^-1(z) 1));
  M_0[n][j][i] = (1/5) Sum_{a+b=3} Q_{b,0}^(j) Q_{a,n}^(i) at u = 0 (Teleman form), M_0^-1 by series inversion (and asserted to be
     M_0(-z)^T), E_M0 = (1 - M_0^-1(z) M_0^-1(w)^T)/(z + w) by the same division; asserted: E_M0 = [u^0] E_R;
  M = M_0(tau zhat): E_M(a, b) = tau^(a+b+1) E_M0(a, b) (asserted against the division applied to M^-1 = M_0^-1(tau z));
     variants: 'nat' tau = 1/Sum 4 u^n/(5n+4); 'M0' tau = 1; 'id' E_M = 0; 'nat_t<k>' tau_k shifted by 1 (control);
  E_out = E_R - E_M;
  vertex (g, half-edge exponents hp, hs = len(hp)), colour c:
     V[c] = Sum_{partitions lam of r = 3g-3+hs-|hp|} WK(g; hp, lam_i + 1)/|Aut lam| (zeta^(4c) xhat)^(2g-2+hs+len(lam)) prod_i T_R(lam_i+1)[c].
Fourier transform in the colour (Z/5 covariance is NOT assumed: all components are computed and the sparsity is asserted):
  vector  Vh[k] = (1/5) Sum_c zeta^(-kc) V[c]: exactly one k = kappa is non-zero;
  matrix  Eh[k][l] = (1/25) Sum_ij zeta^(-ki-lj) E[i][j]: non-zero only for k + l = sigma(a, b) mod 5.
Series: lists of K residues (index = power of u), None = 0.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "referee")))
from ref_fp import P89, P61, of, inv
from ref_hat import HatRows
from ind_wk import wk

PRIMES = {"89": P89, "61": P61}


def partitions(r, cap=None):
    """partitions of r as non-increasing tuples"""
    if cap is None: cap = r
    if r == 0: yield (); return
    for x in range(min(cap, r), 0, -1):
        for rest in partitions(r - x, x): yield (x,) + rest


class Tables(object):
    def __init__(self, prime, g, K, variant="nat", ab=(0, 0)):
        p = self.p = PRIMES[prime]; self.g = g; self.K = K; self.variant = variant; NZ = self.NZ = 3 * g - 3
        H = self.H = HatRows(p, NZ, K, alpha=Fr(ab[0]), beta=Fr(ab[1]))
        assert all(H.checks.values()), H.checks
        self.zeta = z = H.zeta; Q = H.Q; i5 = inv(5, p); self.zp = [pow(z, k, p) for k in range(5)]
        one = [1] + [0] * (K - 1)
        # ---------------- E_R
        def dR(M, N):
            sg = (p - 1) if (M + N) % 2 == 0 else 1                      # -(-1)^(M+N)
            out = [[None] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    acc = None
                    for a in range(5): acc = self.add(acc, self.mul(Q[a][i][M], Q[(3 - a) % 5][j][N]))
                    acc = self.scale(acc, sg * i5 % p)
                    if M == 0 and N == 0 and i == j: acc = self.add(acc, one)
                    out[i][j] = acc
            return out
        self.ER = self.divide(dR, NZ)
        # ---------------- leaves
        self.TR = {e: [self.clean(self.scale(Q[0][c][e - 1], 1 if e % 2 == 0 else p - 1)) for c in range(5)] for e in range(2, NZ + 2)}
        # ---------------- M_0, its inverse, E_M0
        M0 = [[[sum(Q[(3 - a) % 5][j][0][0] * Q[a][i][n][0] for a in range(5)) * i5 % p for i in range(5)] for j in range(5)] for n in range(NZ + 1)]
        assert M0[0] == [[1 if i == j else 0 for j in range(5)] for i in range(5)], "M_0(0) != 1"
        assert all(M0[n] == H.Rt.get((n, 0), [[0] * 5 for _ in range(5)]) for n in range(NZ + 1))
        Mi = [M0[0]]
        for n in range(1, NZ + 1):
            Mi.append([[(-sum(M0[m][i][k] * Mi[n - m][k][j] for m in range(1, n + 1) for k in range(5))) % p for j in range(5)] for i in range(5)])
        assert all(Mi[n][i][j] == M0[n][j][i] * (1 if n % 2 == 0 else p - 1) % p for n in range(NZ + 1) for i in range(5) for j in range(5)), "M_0 not symplectic"
        self.M0 = M0; self.M0inv = Mi
        def dM0(M, N):
            return [[[((1 if (M == 0 and N == 0 and i == j) else 0) - sum(Mi[M][i][k] * Mi[N][j][k] for k in range(5))) % p] + [0] * (K - 1) for j in range(5)] for i in range(5)]
        EM0 = self.divide(lambda M, N: [[self.clean(x) for x in row] for row in dM0(M, N)], NZ)
        self.E0 = {key: [[(x[0] if x else 0) for x in row] for row in tab] for key, tab in EM0.items()}        # scalars
        self.frozen_ok = all((self.ER[key][i][j] or [0])[0] == self.E0[key][i][j] for key in self.ER for i in range(5) for j in range(5))
        assert self.frozen_ok, "[u^0] E_R != E_M0"
        # ---------------- tau, E_M, E_out
        if variant.startswith("nat"):
            ser = [Fr(4, 5 * k + 4) for k in range(K)]; tq = [Fr(1)]
            for k in range(1, K): tq.append(-sum(ser[j] * tq[k - j] for j in range(1, k + 1)))
            tau = [of(x, p) for x in tq]
            if variant.startswith("nat_t"): kk = int(variant[5:]); tau[kk] = (tau[kk] + 1) % p
        else: tau = list(one)
        self.tau = tau; self.tp = [list(one)]
        for n in range(1, NZ + 1): self.tp.append(self.mul(self.tp[-1], tau))
        if variant == "id":
            self.EM = {key: [[None] * 5 for _ in range(5)] for key in self.ER}
        else:
            self.EM = {(a, b): [[self.scale(self.tp[a + b + 1], self.E0[(a, b)][i][j]) for j in range(5)] for i in range(5)] for (a, b) in self.ER}
            def dM(M, N):                                                 # M^-1(z) = M_0^-1(tau z)
                return [[self.scale(self.tp[M + N], dM0(M, N)[i][j][0]) for j in range(5)] for i in range(5)]
            assert self.divide(dM, NZ) == self.EM, "E_M != tau^(a+b+1) E_M0"
        self.EO = {key: [[self.sub(self.ER[key][i][j], self.EM[key][i][j]) for j in range(5)] for i in range(5)] for key in self.ER}
        assert all(self.ER[(a, b)][i][j] == self.ER[(b, a)][j][i] for (a, b) in self.ER for i in range(5) for j in range(5)), "E_R(a,b)^T != E_R(b,a)"
        # ---------------- Fourier transforms of the edges
        self.FO = {}; self.FI = {}; self.sig = {}; self.valO = {}
        for key in self.ER:
            so, fo = self.fourier2(self.EO[key]); si, fi = self.fourier2([[([x] + [0] * (K - 1)) if x else None for x in row] for row in self.E0[key]])
            assert so is None or si is None or so == si, "sigma differs between E_out and E_M"
            s = so if so is not None else si; self.sig[key] = s if s is not None else 0
            assert s is None or s == (-(key[0] + key[1] + 1)) % 5, "unexpected Z/5 charge of an edge"
            self.FI[key] = [0] * 5 if variant == "id" else [(x[0] if x else 0) for x in fi]
            v = min([self.val(x) for x in fo]); self.valO[key] = v
            self.FO[key] = [(x[v:] if x is not None else None) for x in fo]          # divided by u^v, length K - v
        # ---------------- xhat powers, vertex tables (filled on demand)
        self.xp = [list(one)]
        for n in range(1, 6 * g): self.xp.append(self.mul(self.xp[-1], H.C0[:K]))
        self._vt = {}

    # ------------------------------------------------------------ series
    def clean(self, a): return a if (a is not None and any(a)) else None

    def mul(self, a, b, L=None):
        if a is None or b is None: return None
        L = min(len(a), len(b)) if L is None else L; out = [0] * L
        for i in range(min(L, len(a))):
            x = a[i]
            if x:
                for j in range(min(L - i, len(b))):
                    y = b[j]
                    if y: out[i + j] += x * y
        p = self.p; out = [c % p for c in out]
        return out if any(out) else None

    def scale(self, a, c):
        c %= self.p
        if a is None or not c: return None
        p = self.p; return [x * c % p for x in a]

    def add(self, a, b):
        if a is None: return b
        if b is None: return a
        p = self.p; out = [(x + y) % p for x, y in zip(a, b)]
        return out if any(out) else None

    def sub(self, a, b):
        return self.add(a, self.scale(b, self.p - 1))

    def val(self, a):
        if a is None: return self.K
        return next((k for k, x in enumerate(a) if x), self.K)

    # ------------------------------------------------------------ (1 - N)/(z + w)
    def divide(self, d, NZ):
        """d(M, N) -> 5x5 of series, M + N <= NZ; returns {(m, n): 5x5} with m + n + 1 <= NZ and asserts exact divisibility"""
        D = {(M, N): d(M, N) for M in range(NZ + 1) for N in range(NZ + 1 - M)}
        assert all(x is None for row in D[(0, 0)] for x in row), "N(0, 0) != 1"
        E = {}
        for t in range(NZ):                                              # t = m + n
            for n in range(t + 1):
                m = t - n; cur = D[(m + 1, n)]
                if n > 0: cur = [[self.sub(cur[i][j], E[(m + 1, n - 1)][i][j]) for j in range(5)] for i in range(5)]
                E[(m, n)] = cur
            assert D[(0, t + 1)] == E[(0, t)], "(1 - N) is not divisible by z + w at degree %d" % (t + 1)
        return E

    # ------------------------------------------------------------ Fourier
    def fourier2(self, E):
        """(sigma, [Eh[k][sigma - k] for k in 0..4]) for a 5x5 matrix of series; asserts the selection rule k + l = sigma"""
        p = self.p; zp = self.zp; i25 = inv(25, p); comp = {}
        for k in range(5):
            for l in range(5):
                acc = None
                for i in range(5):
                    for j in range(5):
                        acc = self.add(acc, self.scale(E[i][j], zp[(-k * i - l * j) % 5]))
                acc = self.scale(acc, i25)
                if acc is not None: comp[(k, l)] = acc
        sig = set((k + l) % 5 for (k, l) in comp)
        assert len(sig) <= 1, "Z/5 selection rule violated on an edge"
        if not sig: return None, [None] * 5
        s = sig.pop()
        return s, [comp.get((k, (s - k) % 5)) for k in range(5)]

    # ------------------------------------------------------------ vertices
    def vertex_colours(self, g, hp):
        """[V[c] for c in 0..4] for the sorted exponent tuple hp"""
        p = self.p; hs = len(hp); r = 3 * g - 3 + hs - sum(hp); out = [None] * 5
        if r < 0: return out
        for lam in partitions(r):
            w = wk(g, tuple(hp) + tuple(x + 1 for x in lam))
            if w == 0: continue
            sym = 1
            for x in set(lam):
                m = lam.count(x)
                for q in range(2, m + 1): sym *= q
            c0 = of(w / sym, p); n = 2 * g - 2 + hs + len(lam)
            for c in range(5):
                t = self.scale(self.xp[n], c0 * self.zp[(4 * c * n) % 5])
                for x in lam: t = self.mul(t, self.TR[x + 1][c])
                out[c] = self.add(out[c], t)
        return out

    def vertex_table(self, g, hs):
        """{sorted hp (non-increasing): (kappa, Vh[kappa], [V[c]])}, only non-zero entries"""
        key = (g, hs); tab = self._vt.get(key)
        if tab is not None: return tab
        p = self.p; zp = self.zp; i5 = inv(5, p); tab = {}
        def tuples(n, cap, tot):
            if n == 0: yield (); return
            for x in range(min(cap, tot), -1, -1):
                for rest in tuples(n - 1, x, tot - x): yield (x,) + rest
        for hp in tuples(hs, 3 * g - 3 + hs, 3 * g - 3 + hs):
            V = self.vertex_colours(g, hp)
            if all(x is None for x in V): continue
            Vh = []
            for k in range(5):
                acc = None
                for c in range(5): acc = self.add(acc, self.scale(V[c], zp[(-k * c) % 5]))
                Vh.append(self.scale(acc, i5))
            nz = [k for k in range(5) if Vh[k] is not None]
            assert len(nz) == 1, "Z/5 selection rule violated on a vertex"
            assert Vh[nz[0]] == V[0]
            tab[hp] = (nz[0], Vh[nz[0]], V)
        self._vt[key] = tab
        return tab
