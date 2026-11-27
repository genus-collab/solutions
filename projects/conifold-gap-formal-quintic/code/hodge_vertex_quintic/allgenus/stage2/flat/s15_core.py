# -*- coding: utf-8 -*-
"""s15_core.py -- the FLAT P-action graph sum of P.Omega^mod_{v(u)} (REDUCTION_DRAFT.md s.1.5, Theorem R's form) for ANY stable
leg-free graph, modulo a prime p (2^89-1 or 2^61-1), as truncated Laurent series in u.

Ingredients (exact over Q, then reduced mod p):
  P      = canonical Theorem-U gauge (theory/referee/rf_flat_core: ur_gauge read-only + P c(z)^-1), lowered form, z-order <= K;
  Pvi    = vector form of P^-1;  d = (P_0)_33 (dilaton factor);  leaf point q_T = Pv^-1(-z phi_0) = -(z/d) phi_0 + T,
           dT[(k+1, c)] = -d (Pvi_k)_{c0}, k >= 1;
  E      = flat edge (eta^-1 - Pv^-1(psi) eta^-1 Pv^-1(psi')^T)/(psi + psi'), entries E[((m,x),(n,y))] (rf_g2flat convention);
  vertex = d^(2h-2+n) Sum_{leaf multisets L} prod_{l in L} dT_l / mult(L)!  <labels, L>_h(v(u)),
           <prod tau_k(phi_a)>_h(v) = C v^e,  e = 2-2h+Sum w(k,a),  C = flat model correlator at v = 1 (MVertex6 fed with the
           referee's LG M_0, rf_mla_core.lg_M0; flat <-> idempotent vectors c[a][i] = NSG[a] zeta^(i(a+1)), as rf_mla_core.Corr).
Series: (offset o, list of N coefficients) = u^o (c_0 + .. + c_{N-1} u^{N-1}) + O(u^{o+N}); offsets are never raised, so every
result is exact through u^{o_min+N-1}.  Every vertex term has o >= 2-2h-n (w >= -1, leaves w > 0), so a genus-g graph is exact
through u^{N-(2g-2)-1}.
Variants (failing controls): 'nodil' = drop the dilaton factor d^(2h-2+n) and use the unscaled leaf T (= dT/d);
'noT' = leaves removed; 'nosign' = edge quotient without the alternating sign (-1)^t; 'stored' = the stored (non-canonical)
Theorem-U gauge P_stored = P_can c(z) instead of the canonical one (the (N)-failure mode of S15_PROOF.md s.5).
Reused read-only: rf_flat_core, rf_mla_core (theory referee), fq_fp, fq3_mvertex.MVertex6/MConst, fq3_graphs.outer_graphs.
"""
import os, sys, itertools, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
REF = os.path.join(PROJ, "allgenus", "theory", "referee")
for _p in (PROJ, REF):
    if _p not in sys.path: sys.path.insert(0, _p)
import rf_flat_core as C
from rf_flat_core import yc_series, spow, smul, sadd, sscale, zinv, eta_conj, mT, mmul, mscale, madd, mzero, tau_v, load_gauge, canonical
from rf_mla_core import lg_M0, w, NSG
from fq_fp import FpRing
from fq3_mvertex import MVertex6, MConst
from fq3_graphs import outer_graphs


class LS(object):
    """truncated Laurent series mod p: u^o (c_0 + .. + c_{N-1} u^{N-1})"""
    __slots__ = ("o", "c")
    def __init__(self, o, c): self.o = o; self.c = c


class Flat(object):
    def __init__(self, p, N, K=8, a=0, b=0, variant=None, log=print):
        self.p = p; self.N = N; self.K = K; self.var = variant; R = FpRing(p); self.R = R
        C.NU[0] = N + 1           # one extra order: W = v/u loses one (the last entry of v[1:] + [0] is not a coefficient)
        tau, v, itau = tau_v()
        Ps, tc, var0, res = load_gauge(K, Fr(a), Fr(b))
        P, self.cz = canonical(Ps, K)
        if variant == 'stored': P = Ps          # the stored Theorem-U normalisation P_stored = P_can c(z), c = 1 + m5' z^5 + ..
        Y = yc_series(v, K + 2); te = [Fr(0)] * (N + 1)
        for k in range(K + 1):
            for j in range(5): te = sadd(te, smul(P[k][0][j], Y[k + 1][j]))
        self.teff = sscale(te, Fr(-1)); self.tc = tc; self.tau = tau
        Pi = zinv(P, K); Pvi = {k: eta_conj(Pi[k]) for k in Pi}
        d = P[0][3][3]; self.dQ = d
        nd = variant == 'nodil'
        self.dT = {}
        for k in range(1, K + 1):
            for c in range(5):
                s = sscale(Pvi[k][c][0], Fr(-1))
                if any(s): self.dT[(k + 1, c)] = self.ser(s if nd else smul(d, s))
        if variant == 'noT': self.dT = {}
        for (k, c) in self.dT: assert w(k, c).denominator == 1, ("leaf weight not integral", k, c)
        EI = mzero()
        for i in range(5): EI[i][(3 - i) % 5] = C.sc(Fr(1, 5))
        Num = {}
        for m in range(K + 1):
            for n in range(K + 1 - m): Num[(m, n)] = mscale(mmul(mmul(Pvi[m], EI), mT(Pvi[n])), Fr(-1))
        Num[(0, 0)] = madd(Num[(0, 0)], EI)
        self.E = {}
        for m in range(K):
            for n in range(K - m):
                acc = mzero()
                for t in range(n + 1): acc = madd(acc, Num[(m + 1 + t, n - t)], (-1) ** t if variant != 'nosign' else 1)
                for x in range(5):
                    for y in range(5):
                        if any(acc[x][y]): self.E[((m, x), (n, y))] = self.ser(acc[x][y])
        W = v[1:] + [Fr(0)]
        self.Wp = {e: self.ser(spow(W, Fr(e))) for e in range(-8, 12)}
        self.dp = {e: self.ser([Fr(1)] + [Fr(0)] * (N - 1) if nd else spow(d, Fr(e))) for e in range(-3, 9)}
        M, Mi = lg_M0(R, 6)
        self.mv = MVertex6(R, MConst(R, M, Mi, check=True), kmax=6)
        self.cv = [[(NSG[aa] * R.z(i * (aa + 1))) % p for i in range(5)] for aa in range(5)]
        self.vcache = {}; self.leafkeys = sorted(self.dT)
        log("  Flat setup p = %d, N = %d, K = %d, frame (%s, %s), variant %s: d = %s.., leaf keys %s, edge entries %d"
            % (p.bit_length(), N, K, a, b, variant, [str(x) for x in d[:3]], self.leafkeys, len(self.E)))

    # ---------------- series
    def ser(self, s):
        p = self.p; return LS(0, [(x.numerator * pow(x.denominator, -1, p)) % p for x in s[:self.N]] + [0] * (self.N - len(s)))

    def mul(self, a, b):
        N = self.N; p = self.p; A = a.c; B = b.c; c = [0] * N
        for i in range(N):
            x = A[i]
            if x:
                for j in range(N - i):
                    y = B[j]
                    if y: c[i + j] += x * y
        return LS(a.o + b.o, [x % p for x in c])

    def scal(self, a, s): p = self.p; return LS(a.o, [x * s % p for x in a.c])

    def add(self, a, b):
        if a is None: return b
        o = min(a.o, b.o); N = self.N; c = [0] * N
        for s in (a, b):
            sh = s.o - o
            for i in range(N - sh): c[i + sh] += s.c[i]
        return LS(o, [x % self.p for x in c])

    def coeffs(self, a, lo, hi):
        """[u^lo .. u^hi] of a, asserting validity (hi <= o + N - 1)"""
        if a is None: return [0] * (hi - lo + 1)
        assert hi <= a.o + self.N - 1, ("precision", a.o, self.N, hi)
        return [a.c[k - a.o] if 0 <= k - a.o < self.N else 0 for k in range(lo, hi + 1)]

    # ---------------- vertex tensors in flat labels
    def leaf_multisets(self, budget):
        out = []; keys = self.leafkeys
        def rec(i, cur, b):
            out.append(tuple(cur))
            for j in range(i, len(keys)):
                cost = keys[j][0] - 1
                if cost <= b: cur.append(keys[j]); rec(j, cur, b - cost); cur.pop()
        rec(0, [], budget)
        return out

    def _contract_last(self, T, m, vec):
        s = 5 ** (m - 1); p = self.p; out = [0] * s
        for k in range(5):
            x = vec[k]
            if not x: continue
            blk = T[k * s:(k + 1) * s]
            out = [(u + x * y) % p for u, y in zip(out, blk)]
        return out

    def _to_flat(self, T, n):
        p = self.p; cv = self.cv
        for l in range(n):
            s = 5 ** l; new = [0] * len(T)
            for hi in range(5 ** (n - l - 1)):
                for lo in range(s):
                    base = hi * 5 * s + lo; vals = [T[base + i * s] for i in range(5)]
                    if not any(vals): continue
                    for a in range(5): new[base + a * s] = sum(cv[a][i] * vals[i] for i in range(5)) % p
            T = new
        return T

    def vertex(self, h, ks):
        """{flat colour tuple: LS} for labels with psi-powers ks (half-edge order), leaves summed"""
        key = (h, tuple(ks))
        if key in self.vcache: return self.vcache[key]
        n = len(ks); budget = 3 * h - 3 + n - sum(ks); out = {}; p = self.p
        if budget >= 0:
            for L in self.leaf_multisets(budget):
                if 2 * h - 2 + n + len(L) <= 0: continue
                T = self.mv.tensor(h, tuple(ks) + tuple(k for k, c in L))
                if not any(T): continue
                m = n + len(L)
                for (k, c) in reversed(L): T = self._contract_last(T, m, self.cv[c]); m -= 1
                if not any(T): continue
                T = self._to_flat(T, n)
                lser = None; mult = 1
                for kc in set(L): mult *= math.factorial(L.count(kc))
                for kc in L: lser = self.dT[kc] if lser is None else self.mul(lser, self.dT[kc])
                wl = sum(w(k, c) for k, c in L)
                for idx, x in enumerate(T):
                    if not x: continue
                    cols = []; r = idx
                    for _ in range(n): cols.append(r % 5); r //= 5
                    e = 2 - 2 * h + sum(w(k, a) for k, a in zip(ks, cols)) + wl
                    assert e.denominator == 1, ("non-integral exponent with non-zero correlator", h, ks, L, cols)
                    e = int(e)
                    term = self.mul(self.Wp[e], self.dp[2 * h - 2 + n])
                    term = LS(term.o + e, term.c)
                    if lser is not None: term = self.mul(term, lser)
                    term = self.scal(term, x * pow(mult, -1, p) % p)
                    tc = tuple(cols); out[tc] = self.add(out.get(tc), term)
        self.vcache[key] = out
        return out
