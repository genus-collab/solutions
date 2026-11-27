# -*- coding: utf-8 -*-
"""fq3_mvertex.py -- the colour-mixing vertex CohFT OmegaHat^M = M.omega^ of the formal quintic to M-ORDER 6 (genus-3 programme).

Same conventions as fq_mvertex.MVertex (genus-2 code, M-orders <= 3), generalised to psi-deficit D <= 6, i.e. M(zhat) mod zhat^7:
  OmegaHat_{g,n}(i_1..i_n; d) = Sum_{legged stable graphs Gamma' (fq3_graphs.legged_graphs(g, n, D))} (1/|Aut|) Sum_{decorations}
     prod_v [ WK(g_v; exponents at v)/m_v! zeta^{4 j_v (2g_v-2+n_v+m_v)} prod_{kappa-leaves} That[e][j_v] ]
     prod_edges Ehat[(a, b)][j][j']  prod_legs (M^{-1}_{k_l})_{j_{v(l)} i_l},
  Ehat[(a,b)] = Sum_{t=0}^{b} (-1)^t Gh[(a+1+t, b-t)],  Gh[(Mo,No)][j][j'] = -Sum_k (M^{-1}_Mo)_{jk} (M^{-1}_No)_{j'k}  ((0,0): 0),
  That[e][j] = -Sum_k (M^{-1}_{e-1})_{jk} zeta^{-4k}  (e >= 2),  Delta^hat_j^{1/2} = zeta^{4j}, unit 1^hat_k = zeta^{-4k}.
Total M-weight of a term = D automatically (dimension constraints at every vertex); inner graphs have <= D edges.
Implementation: per inner graph, the vertex options are GROUPED by (ks_v, hp_v) (M-orders of the legs at v, psi-exponents of the
half-edges at v) with the sum over (m_v, leaves) done once per group as a 5-vector over the colour j_v; for every combination of
groups the colour sum over the V vertices is accumulated into a core tensor per leg-order tuple ks, and the legs are contracted
with the rows of M^{-1}_{k_l} once per ks (vertex by vertex).  Entries are F_p integers; M is any object with .Minv = list of
>= 7 matrices (5 x 5, F_p) -- e.g. fq3_m.MConst -- only M^{-1} enters.  tensor(g, d) is cached per (g, d).
Checks (self-test): equality with fq_mvertex.MVertex (random M of order 3 extended to order 6) on all (g, d) with D <= 3, g <= 2;
diagonal M = diag(exp(Sum_{k odd} a_k zeta^{-ki} zhat^k)) gives zeta^{4i(2g-2+n)} Int exp(Sum s_k ch_k) prod psi^d with
s_1 = 12 a_1 zeta^{-i}, s_3 = -720 a_3 zeta^{-3i}, s_5 = 30240 a_5 zeta^{-5i} (fq_phi conventions) and vanishes off the diagonal,
through genus 3 (D <= 6); string and dilaton equations.
"""
import itertools, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_hodge import wk
from fq_graphs import compositions, tuples_sum_le
from fq3_graphs import legged_graphs

KMAX = 6


class MVertex6(object):
    def __init__(self, R, M, kmax=KMAX):
        self.R = R; p = R.p; self.p = p; self.kmax = kmax
        self.Minv = [[[x % p for x in row] for row in M.Minv[k]] for k in range(kmax + 1)]
        self.zh = [R.z(4 * j) for j in range(5)]
        self.unit = [R.z(-4 * k) for k in range(5)]
        self.cache = {}
        Gh = {}
        for Mo in range(kmax + 1):
            for No in range(kmax + 1 - Mo):
                if Mo == 0 and No == 0: Gh[(0, 0)] = [[0] * 5 for _ in range(5)]; continue
                Gh[(Mo, No)] = [[(-sum(self.Minv[Mo][j][k] * self.Minv[No][jj][k] for k in range(5))) % p for jj in range(5)] for j in range(5)]
        self.Eh = {}
        for a in range(kmax):
            for b in range(kmax - a):
                tab = [[0] * 5 for _ in range(5)]
                for t in range(b + 1):
                    Gt = Gh[(a + 1 + t, b - t)]; sg = 1 if t % 2 == 0 else p - 1
                    for j in range(5):
                        for jj in range(5): tab[j][jj] = (tab[j][jj] + sg * Gt[j][jj]) % p
                self.Eh[(a, b)] = tab
        self.Th = {e: [(-sum(self.Minv[e - 1][j][k] * self.unit[k] for k in range(5))) % p for j in range(5)] for e in range(2, kmax + 2)}
        self._zpow = {}

    def _zp(self, j, e):
        key = (j, e); r = self._zpow.get(key)
        if r is None: r = pow(self.zh[j], e, self.p); self._zpow[key] = r
        return r

    def vertex_groups(self, g, hs, legd, D):
        """{(ks, hp): 5-vector over j} = Sum_{m, leaves} WK(g; legd+ks, hp, leaves)/m! zeta^{4j(2g-2+hs+nl+m)} prod That[e][j],
        with Sum(legd+ks) + Sum hp + Sum leaves = 3g-3+hs+nl+m and M-weight Sum ks + Sum(e-1) <= D"""
        R = self.R; p = self.p; nl = len(legd); out = {}
        dim0 = 3 * g - 3 + hs + nl
        for m in range(0, D + 1):
            dim = dim0 + m; rem = dim - sum(legd)
            if rem < 2 * m or rem < 0: continue
            leaf_opts = [()] if m == 0 else [c for s in range(2 * m, min(rem, 2 * m + D) + 1) for c in compositions(s, m, 2)]
            fm = R.inv(math.factorial(m))
            for leaves in leaf_opts:
                wl = sum(e - 1 for e in leaves)
                if wl > D: continue
                left = rem - sum(leaves)
                for ks in tuples_sum_le(nl, min(left, D - wl)):
                    left2 = left - sum(ks)
                    for hp in tuples_sum_le(hs, left2):
                        if sum(hp) != left2: continue
                        w = wk(g, tuple(d + k for d, k in zip(legd, ks)) + tuple(hp) + tuple(leaves))
                        if w == 0: continue
                        c = R.of(w) * fm % p
                        key = (tuple(ks), tuple(hp)); vec = out.get(key)
                        if vec is None: vec = [0] * 5; out[key] = vec
                        for j in range(5):
                            val = c * self._zp(j, 2 * g - 2 + hs + nl + m) % p
                            for e in leaves: val = val * self.Th[e][j] % p
                            vec[j] = (vec[j] + val) % p
        return {k: v for k, v in out.items() if any(v)}

    def tensor(self, g, d):
        """OmegaHat_{g,n}(i_1..i_n; d): flat list of length 5^n, index Sum_l i_l 5^l"""
        d = tuple(d); key = (g, d)
        r = self.cache.get(key)
        if r is not None: return r
        n = len(d); D = 3 * g - 3 + n - sum(d)
        assert D <= self.kmax, "psi-deficit %d > %d" % (D, self.kmax)
        out = [0] * 5 ** n
        if D >= 0 and 2 * g - 2 + n > 0:
            for Gm in legged_graphs(g, n, D):
                self._add_graph(Gm, d, D, out)
        self.cache[key] = out
        return out

    def _add_graph(self, Gm, d, D, out):
        p = self.p; n = len(d)
        V = Gm['V']; gv = Gm['gv']; ed = Gm['ed']; hsAt = Gm['hsAt']; legsAt = Gm['legsAt']; E = len(ed)
        groups = [self.vertex_groups(gv[v], len(hsAt[v]), tuple(d[l] for l in legsAt[v]), D) for v in range(V)]
        if any(not gr for gr in groups): return
        inv_aut = self.R.inv(Gm['aut'])
        cols = list(itertools.product(range(5), repeat=V))
        cores = {}                                                 # ks_all (tuple over legs 0..n-1) -> {col: value}
        keys = [list(gr.keys()) for gr in groups]
        for combo in itertools.product(*keys):
            hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = combo[v][1][r]
            etabs = [self.Eh[(hp[2 * e], hp[2 * e + 1])] for e in range(E)]
            vecs = [groups[v][combo[v]] for v in range(V)]
            ks_all = [0] * n
            for v in range(V):
                for r, l in enumerate(legsAt[v]): ks_all[l] = combo[v][0][r]
            core = cores.setdefault(tuple(ks_all), {})
            for col in cols:
                val = inv_aut
                for v in range(V):
                    val = val * vecs[v][col[v]] % p
                    if not val: break
                if not val: continue
                for e in range(E):
                    val = val * etabs[e][col[ed[e][0]]][col[ed[e][1]]] % p
                    if not val: break
                if not val: continue
                core[col] = (core.get(col, 0) + val) % p
        for ks_all, core in cores.items():
            self._contract_legs(core, ks_all, Gm, n, out)

    def _contract_legs(self, core, ks_all, Gm, n, out):
        """out[i_1..i_n] += Sum_col core[col] prod_l Minv[k_l][col[v(l)]][i_l], vertex by vertex"""
        p = self.p; V = Gm['V']; legsAt = Gm['legsAt']
        # current tensor: dict  (col of unprocessed vertices v..V-1, leg indices so far in leg order) -> value
        cur = {col: val for col, val in core.items() if val}
        done_legs = []
        for v in range(V):
            L = legsAt[v]; mats = [self.Minv[ks_all[l]] for l in L]
            new = {}
            for key, val in cur.items():
                c = key[0]; rest = key[1:]
                if not L:
                    new[rest] = (new.get(rest, 0) + val) % p; continue
                rows = [mats[r][c] for r in range(len(L))]
                for idx in itertools.product(range(5), repeat=len(L)):
                    x = val
                    for r in range(len(L)):
                        x = x * rows[r][idx[r]] % p
                        if not x: break
                    if not x: continue
                    k2 = rest + idx
                    new[k2] = (new.get(k2, 0) + x) % p
            cur = new; done_legs += list(L)
        if n == 0:
            out[0] = (out[0] + cur.get((), 0)) % p; return
        pos_of = {l: 5 ** l for l in range(n)}
        for key, val in cur.items():
            if not val: continue
            pos = sum(key[r] * pos_of[done_legs[r]] for r in range(n))
            out[pos] = (out[pos] + val) % p


# ------------------------------------------------------------------ helpers for tests
class MConst(object):
    """constant M(zhat) = Sum_k M[k] zhat^k mod zhat^{kmax+1} with inverse; duck-typed for MVertex6 / rows (.M, .Minv, .R, .p, .is_identity)"""
    def __init__(self, R, M, Minv, check=True):
        from fq_sser import mat_mul, mat_T, mat_id, mat_zero, mat_add, mat_eq
        p = R.p; self.R = R; self.p = p; self.M = M; self.Minv = Minv; self.kmax = len(M) - 1
        self.is_identity = all(not any(x % p for r in M[k] for x in r) for k in range(1, self.kmax + 1))
        if check:
            for k in range(self.kmax + 1):
                S1 = mat_zero(5); S2 = mat_zero(5)
                for a in range(k + 1):
                    S1 = mat_add(S1, mat_mul(mat_T(M[a]), M[k - a], p), p, (1 if (k - a) % 2 == 0 else p - 1))
                    S2 = mat_add(S2, mat_mul(M[a], Minv[k - a], p), p)
                target = mat_id(5, p) if k == 0 else mat_zero(5)
                assert mat_eq(S1, target, p), "MConst: not symplectic at zhat^%d" % k
                assert mat_eq(S2, target, p), "MConst: M M^-1 != 1 at zhat^%d" % k

    @staticmethod
    def from_log(R, mu, kmax=KMAX, check=True):
        """M = exp(mu(zhat)), mu = [0, mu_1, .., mu_kmax] (mu_k^T = (-1)^{k+1} mu_k)"""
        from fq_symbol import mat_series_exp
        p = R.p
        Z = [[0] * 5 for _ in range(5)]
        mu = [(mu[k] if (k < len(mu) and mu[k] is not None) else Z) for k in range(kmax + 1)]
        neg = [[[(-x) % p for x in r] for r in m] for m in mu]
        return MConst(R, mat_series_exp(mu, R, kmax), mat_series_exp(neg, R, kmax), check)

    @staticmethod
    def identity(R, kmax=KMAX):
        I = [[1 if i == j else 0 for j in range(5)] for i in range(5)]; Z = [[0] * 5 for _ in range(5)]
        return MConst(R, [I] + [Z] * kmax, [I] + [Z] * kmax, False)

    @staticmethod
    def diagonal_odd(R, a, kmax=KMAX, check=True):
        """M = diag(exp(Sum_{k odd} a[k][i] zhat^k)); a = {k: list of 5}"""
        mu = [[[0] * 5 for _ in range(5)] for _ in range(kmax + 1)]
        for k, lst in a.items():
            for i in range(5): mu[k][i][i] = R.of(lst[i])
        return MConst.from_log(R, mu, kmax, check)


def hodge_integral(R, g, d, s):
    """Int_{Mbar_{g,n}} exp(Sum_{k odd} s[k] ch_k) prod psi^d  (F_p), via fq_hodge.hi"""
    from fq_hodge import hi, partitions_odd, dmax
    from fractions import Fraction as Fr
    p = R.p; n = len(d); D = 3 * g - 3 + n - sum(d)
    if D < 0 or 2 * g - 2 + n <= 0: return 0
    if D == 0: return R.of(wk(g, d)) if n else 0
    if D > dmax(g): return 0
    tot = 0
    for K in partitions_odd(D, max(2 * g - 1, 1)):
        x = hi(g, d, K)
        if x == 0: continue
        c = R.of(x); cnt = {}
        for k in K: cnt[k] = cnt.get(k, 0) + 1
        for k, m in cnt.items(): c = c * pow(s.get(k, 0), m, p) % p * R.inv(math.factorial(m)) % p
        tot = (tot + c) % p
    return tot


if __name__ == "__main__":
    import random, time
    from fractions import Fraction as Fr
    from fq_fp import FpRing, P2
    from fq_sser import SympM
    from fq_mvertex import MVertex
    R = FpRing(P2); p = R.p; rng = random.Random(7)
    # (1) agreement with the genus-2 MVertex on D <= 3 (random order-3 M extended to order 6)
    Ms = SympM.random(R, rng); M6 = MConst.from_log(R, [None, Ms.m1, Ms.m2, Ms.m3])
    assert all(M6.M[k] == Ms.M[k] and M6.Minv[k] == Ms.Minv[k] for k in range(4))
    mv2 = MVertex(R, Ms); mv6 = MVertex6(R, M6); t = time.time(); cnt = 0
    for g in range(3):
        for n in range(0, 5):
            if 2 * g - 2 + n <= 0: continue
            for d in set(tuple(sorted(t_)) for t_ in itertools.product(range(5), repeat=n)):
                D = 3 * g - 3 + n - sum(d)
                if D < 0 or D > 3: continue
                assert mv2.tensor(g, d) == mv6.tensor(g, d), (g, d); cnt += 1
    print("(1) MVertex6 == MVertex on %d tensors with D <= 3, n <= 4 (%.1f s)" % (cnt, time.time() - t))
    # (2) diagonal M: Hodge-type integrals through genus 3 (D <= 6), off-diagonal entries vanish
    a1 = [rng.randrange(p) for _ in range(5)]; a3 = [rng.randrange(p) for _ in range(5)]; a5 = [rng.randrange(p) for _ in range(5)]
    Md = MConst.diagonal_odd(R, {1: a1, 3: a3, 5: a5}); mvd = MVertex6(R, Md); t = time.time(); cnt = 0
    tests = [(3, ()), (3, (2,)), (3, (1,)), (3, (3,)), (3, (1, 2)), (3, (2, 2)), (3, (1, 1, 1)), (3, (2, 1, 1)), (3, (3, 2)),
             (2, (0, 0)), (2, (1, 0)), (2, (0,)), (2, (0, 0, 0)), (2, (1, 1)), (2, ()), (1, (0, 0, 0)), (1, (0, 0)), (1, (0,)), (1, (0, 0, 0, 0)),
             (0, (0, 0, 0, 0, 0)), (0, (0, 0, 0, 0, 0, 0)), (0, (1, 0, 0, 0, 0, 0)), (0, (0, 0, 0, 0))]
    for g, d in tests:
        n = len(d); D = 3 * g - 3 + n - sum(d)
        T = mvd.tensor(g, d); t1 = time.time() - t
        for idx in itertools.product(range(5), repeat=n):
            val = T[sum(i * 5 ** l for l, i in enumerate(idx))]
            if n and len(set(idx)) > 1: assert val == 0, (g, d, idx); continue
            i = idx[0] if n else None
            if n == 0:
                # colour sum of the leafless vertex: Sum_i zeta^{4i(2g-2)} I_i
                exp_ = sum(R.z(4 * i * (2 * g - 2)) * hodge_integral(R, g, d, {1: 12 * a1[i] % p, 3: (-720) * a3[i] % p, 5: 30240 * a5[i] % p}) for i in range(5)) % p
            else:
                exp_ = R.z(4 * i * (2 * g - 2 + n)) * hodge_integral(R, g, d, {1: 12 * a1[i] % p, 3: (-720) * a3[i] % p, 5: 30240 * a5[i] % p}) % p
            assert val == exp_, (g, d, idx, R.fmt(val), R.fmt(exp_))
        cnt += 1
        print("   (2) genus %d d=%s D=%d: Hodge-type check ok (%d inner graphs, %.1f s)" % (g, d, D, len(legged_graphs(g, n, D)), time.time() - t)); t = time.time()
    print("(2) diagonal M: %d tensors agree with Int exp(s1 ch1 + s3 ch3 + s5 ch5) prod psi^d" % cnt)
    # (3) string and dilaton equations at a random full M (order 6)
    mu = [[[0] * 5 for _ in range(5)] for _ in range(7)]
    for k in range(1, 7):
        for i in range(5):
            for j in range(i, 5):
                if k % 2 == 1: mu[k][i][j] = mu[k][j][i] = rng.randrange(p)
                elif i != j: x = rng.randrange(p); mu[k][i][j] = x; mu[k][j][i] = (-x) % p
    Mr = MConst.from_log(R, mu); mvr = MVertex6(R, Mr); unit = [R.z(-4 * k) for k in range(5)]
    def contract_unit(T, n, pos):
        out = [0] * 5 ** (n - 1)
        for idx in itertools.product(range(5), repeat=n):
            rest = idx[:pos] + idx[pos + 1:]
            out[sum(i * 5 ** l for l, i in enumerate(rest))] = (out[sum(i * 5 ** l for l, i in enumerate(rest))] + T[sum(i * 5 ** l for l, i in enumerate(idx))] * unit[idx[pos]]) % p
        return out
    t = time.time()
    for g, d in [(1, (1,)), (1, (2, 0)), (2, (1, 1)), (2, (3,)), (2, (2, 1)), (3, (2, 1)), (3, (3,)), (3, (2, 2, 1)), (0, (1, 1, 0, 0)), (1, (0, 1, 1))]:
        n = len(d) + 1
        # string: Omega_{g,n+1}(.., 1; d, 0) = Sum_j Omega_{g,n}(..; d with d_j - 1)
        T = mvr.tensor(g, d + (0,)); lhs = contract_unit(T, n, n - 1)
        rhs = [0] * 5 ** (n - 1)
        for j in range(len(d)):
            if d[j] == 0: continue
            d2 = list(d); d2[j] -= 1; T2 = mvr.tensor(g, tuple(d2))
            rhs = [(x + y) % p for x, y in zip(rhs, T2)]
        assert lhs == rhs, ("string", g, d)
        # dilaton: Omega_{g,n+1}(.., 1; d, 1) = (2g-2+n-1) Omega_{g,n}(..; d)
        T = mvr.tensor(g, d + (1,)); lhs = contract_unit(T, n, n - 1)
        T0 = mvr.tensor(g, d); rhs = [(2 * g - 2 + len(d)) * x % p for x in T0]
        assert lhs == rhs, ("dilaton", g, d)
    print("(3) string and dilaton equations at a random order-6 M: ok (%.1f s)" % (time.time() - t))
    print("fq3_mvertex ok")
