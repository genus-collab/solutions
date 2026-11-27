# -*- coding: utf-8 -*-
"""rf_s15_outer.py -- referee (own code): the OUTER data of the flat P-action graph sum, mapped into the frozen idempotent frame.

Only project input: natural_form/theory/ur_gauge.gauge (read-only import; the stored Theorem-U gauge p_{k,b}(u) over Q, and the
Conifold object for t_c), exactly the source named in REDUCTION_DRAFT s.1.4 / UNIFORM_REDUCTION s.3.  Everything else is own code.
  P (lowered, P_k = Sum_b p_{k,b} E_{b+k,b});  canonical: P_can = P c(z)^-1, c(z) = P_11(u = 0, z);
  vector form A = J P J (J the involution a <-> 3 - a, eta = 5 J), A^-1 = J P^-1 J (model flat -> quintic flat, inverse);
  d = (A_0)_00 = (P_0)_33;  leaf point q = A^-1(z)(-z 1) = -(z/d) 1 + Sum_k z^(k+1) T_k,  T_k = -A^-1_k e_0;
  flat edge E_A = (J/5 - A^-1(psi) (J/5) A^-1(psi')^T)/(psi + psi');
  frozen scaling: <prod tau_k(H^a)>_h(v) = v^(2-2h) prod v^(w(k,a)) <..>_h(1),  w(k, a) = (c_a + 4k - 4)/5  (asserted integral on every
  non-zero edge entry and leaf);  vertex dilaton factor d^(2h-2+n) (outer leaves carry d T).
Variants: 'can', 'stored' (P instead of P_can), 'nodil' (no d at all), 'noT' (no outer leaves), 'nosign' (quotient without (-1)^t).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
from rf_s15_base import CEXP, smul, sinv, spow_q, shift, fr, inv


def wgt(k, a): return Fr(CEXP[a % 5] + 4 * k - 4, 5)


def load_gauge(K, Nu, a=0, b=0):
    """stored gauge {k: 5x5 list of Fraction lists (u-series, length Nu)} and t_c (Fraction list)"""
    for q in (os.path.join(PROJ, "natural_form", "theory"), PROJ):
        if q not in sys.path: sys.path.insert(0, q)
    import ur_gauge
    NC = Nu + 3 * K + 8
    if (Fr(a), Fr(b)) == (0, 0): res = ur_gauge.gauge(K, NC)
    else:
        import fq_conifold
        Orig = fq_conifold.Conifold
        class Cf(Orig):
            def __init__(self, R, Nu=10, **kw): Orig.__init__(self, R, Nu=Nu, alpha=Fr(a), beta=Fr(b))
        ur_gauge.Conifold = Cf
        try: res = ur_gauge.gauge(K, NC)
        finally: ur_gauge.Conifold = Orig
    P = {}
    for k in range(K + 1):
        M = [[[Fr(0)] * Nu for _ in range(5)] for _ in range(5)]
        for bb in range(5):
            S = res['p'][k][bb]; assert S.prec >= Nu, ("gauge precision", k, bb, S.prec, Nu)
            M[(bb + k) % 5][bb] = [Fr(S.coef(n)) for n in range(Nu)]
        P[k] = M
    tc = [Fr(res['C'].tc.coef(n)) for n in range(Nu)]
    return P, tc


class Outer(object):
    def __init__(self, Fz, P, N, variant='can'):
        p = Fz.p; self.p = p; self.N = N; K = max(P); self.K = K; self.var = variant
        if variant != 'stored':                    # canonical normalisation P_can = P c(z)^-1
            c = {k: P[k][1][1][0] for k in range(0, K + 1, 5)}
            assert c[0] == 1
            ci = {0: Fr(1)}
            for k in range(5, K + 1, 5): ci[k] = -sum(c[j] * ci[k - j] for j in range(5, k + 1, 5))
            Pc = {}
            for k in range(K + 1):
                Pc[k] = [[[sum(ci[j] * P[k - j][x][y][n] for j in range(0, k + 1, 5)) for n in range(N)] for y in range(5)] for x in range(5)]
            self.cz = c; P = Pc
        Pm = {k: [[[fr(P[k][x][y][n], p) for n in range(N)] for y in range(5)] for x in range(5)] for k in P}
        # z-inverse of P (P_0 diagonal)
        P0i = [[[0] * N for _ in range(5)] for _ in range(5)]
        for i in range(5):
            assert all(not any(Pm[0][i][j]) for j in range(5) if j != i)
            P0i[i][i] = sinv(Pm[0][i][i], p, N)
        Pi = {0: P0i}
        for n in range(1, K + 1):
            acc = [[[0] * N for _ in range(5)] for _ in range(5)]
            for k in range(1, n + 1):
                T = self.mm(Pm[k], Pi[n - k])
                acc = [[[(a + b) % p for a, b in zip(acc[i][j], T[i][j])] for j in range(5)] for i in range(5)]
            Pi[n] = self.mm(P0i, [[[(-a) % p for a in acc[i][j]] for j in range(5)] for i in range(5)])
        J = lambda x: (3 - x) % 5
        self.Ai = {k: [[Pi[k][J(x)][J(y)] for y in range(5)] for x in range(5)] for k in Pi}
        self.d = sinv(self.Ai[0][0][0], p, N)
        # v, its powers
        itau = [fr(Fr(4, 5 * n + 4), p) for n in range(N + 2)]
        self.V = spow_q(itau, Fr(5, 4), p, N)            # v = u V,  V = tau^(-5/4)
        self.Vinv = sinv(self.V, p, N)
        self.tau = sinv(itau, p, N)

    def mm(self, A, B):
        p = self.p; N = self.N; out = [[[0] * N for _ in range(5)] for _ in range(5)]
        for i in range(5):
            for k in range(5):
                if not any(A[i][k]): continue
                for j in range(5):
                    if any(B[k][j]):
                        s = smul(A[i][k], B[k][j], p, N); out[i][j] = [(a + b) % p for a, b in zip(out[i][j], s)]
        return out

    def vpow(self, e):
        """v^e for an integer e >= 0 as a plain power series"""
        assert e >= 0
        s = [1] + [0] * (self.N - 1)
        for _ in range(e): s = smul(s, self.V, self.p, self.N)
        return shift(s, e, self.N)

    def flat_edges(self):
        """{((m, x), (n, y)): series} of v^(2 + w + w') E_A (vector form, flat model labels)"""
        p = self.p; N = self.N; K = self.K; Ai = self.Ai; out = {}
        f = {}
        for r in range(K + 1):
            for s in range(K + 1 - r):
                if r == 0 and s == 0: continue
                M = [[[0] * N for _ in range(5)] for _ in range(5)]
                for x in range(5):
                    for y in range(5):
                        acc = [0] * N
                        for xx in range(5):
                            a = Ai[r][x][xx]; yy = (3 - xx) % 5; b = Ai[s][y][yy]
                            if any(a) and any(b): acc = [(c1 + c2) % p for c1, c2 in zip(acc, smul(a, b, p, N))]
                        M[x][y] = [(-c * inv(5, p)) % p for c in acc]
                f[(r, s)] = M
        for Mt in range(1, K + 1):          # divisibility by psi + psi' (symplecticity of P)
            for x in range(5):
                for y in range(5):
                    tot = [0] * N
                    for t in range(Mt + 1): tot = [(a + (-1) ** t * b) % p for a, b in zip(tot, f[(Mt - t, t)][x][y])]
                    assert not any(tot), ("E_A not divisible: P not symplectic", Mt, x, y)
        for m in range(K):
            for n in range(K - m):
                for x in range(5):
                    for y in range(5):
                        acc = [0] * N
                        for t in range(n + 1):
                            sg = 1 if self.var == 'nosign' else (-1) ** t
                            acc = [(a + sg * b) % p for a, b in zip(acc, f[(m + 1 + t, n - t)][x][y])]
                        if not any(acc): continue
                        e = 2 + wgt(m, x) + wgt(n, y); assert e.denominator == 1 and e >= 0, ("edge weight", m, x, n, y, e)
                        out[((m, x), (n, y))] = smul(acc, self.vpow(int(e)), p, N)
        return out

    def flat_leaves(self):
        """{(k, c): series} of d v^w T (k >= 2 = power of z, c flat label)"""
        p = self.p; N = self.N; out = {}
        if self.var == 'noT': return out
        for k in range(1, self.K + 1):
            for c in range(5):
                s = [(-a) % p for a in self.Ai[k][c][0]]
                if not any(s): continue
                e = wgt(k + 1, c); assert e.denominator == 1 and e >= 1, ("leaf weight", k + 1, c, e)
                if self.var != 'nodil': s = smul(s, self.d, p, N)
                out[(k + 1, c)] = smul(s, self.vpow(int(e)), p, N)
        return out

    def to_idempotent(self, Fz, Kmax):
        """outer edges Eo[((p,i),(q,j))] and outer leaves XO[(i,k)] in the frozen idempotent frame (powers <= Kmax)"""
        p = self.p; N = self.N; Bi = Fz.Ri; Pin = Fz.Psi_inv
        L = [[[sum(Bi[r][i][l] * Pin[l][x] for l in range(5)) % p for x in range(5)] for i in range(5)] for r in range(Fz.K + 1)]
        fe = self.flat_edges(); fl = self.flat_leaves()
        XO = {}
        for (k, c), s in fl.items():
            for r in range(0, Kmax - k + 1):
                for i in range(5):
                    a = L[r][i][c]
                    if a: XO[(i, k + r)] = [(x + a * y) % p for x, y in zip(XO.get((i, k + r), [0] * N), s)]
        Eo = {}
        for ((m, x), (n, y)), s in fe.items():
            for r in range(0, Kmax - m + 1):
                for t in range(0, Kmax - n + 1):
                    for i in range(5):
                        a = L[r][i][x]
                        if not a: continue
                        for j in range(5):
                            b = L[t][j][y]
                            if not b: continue
                            key = ((m + r, i), (n + t, j)); ab = a * b % p
                            Eo[key] = [(u + ab * w) % p for u, w in zip(Eo.get(key, [0] * N), s)]
        Eo = {k: v for k, v in Eo.items() if any(v)}
        return Eo, XO

    def global_factor(self, g):
        """(offset, series) of v^(2-2g) d^(2g-2) (d dropped for 'nodil')"""
        p = self.p; N = self.N; s = [1] + [0] * (N - 1)
        for _ in range(2 * g - 2):
            s = smul(s, self.Vinv, p, N)
            if self.var != 'nodil': s = smul(s, self.d, p, N)
        return 2 - 2 * g, s
