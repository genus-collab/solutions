# -*- coding: utf-8 -*-
"""rf_mla_core.py -- referee re-derivation of the ML_Lambda test (REDUCTION_DRAFT.md s.2.4): own code for the LG rows (closed recursion
of LG_MODEL.md s.3), the Teleman matrix M_0 (colour matrix C_n[j][i] = (1/5) Sum_a zeta^(j(4-a) + i(a+1-n)) g^a_n), the frames of T_v at
GENERAL v (homogeneity checked), the flat contraction and the list enumeration.  Read-only project modules: fq_fp.FpRing (field with a
fifth root of unity), fq3_mvertex.MVertex6 and MConst (the vertex graph sum, fed with the referee's own M_0).  No worker module.
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
from rf_flat_core import yc_monomials
CEXP = (0, -1, 3, 2, 1); NSG = (1, -1, -1, 1, 1)


def w(k, a): return Fr(4 * k + CEXP[a] - 4, 5)          # v-weight of tau_k(phi_a)


def lg_rows(N):
    """g^a_n over Q, n <= N (LG_MODEL.md s.3 closed recursion)"""
    g = [[Fr(1)] for _ in range(5)]
    for n in range(1, N + 1):
        d = [Fr(0)] * 5
        for a in range(4): d[a + 1] = d[a] + Fr(4 * n - 3 - CEXP[a], 5) * g[a][n - 1]
        assert d[4] + Fr(4 * n - 3 - CEXP[4], 5) * g[4][n - 1] == 0, "cycle does not close"
        t = -sum((4 * n + 1 - CEXP[a]) * d[a] for a in range(5)) / (20 * n)
        for a in range(5): g[a].append(t + d[a])
    return g


def lg_M0(R, N):
    """(M, Minv) mod p: M[n][j][i] = C_n[j][i], Minv[n][i][j] = (-1)^n C_n[j][i]"""
    p = R.p; g = lg_rows(N); inv5 = R.inv(5); M = []; Mi = []
    for n in range(N + 1):
        C = [[inv5 * sum(R.z(j * (4 - a) + i * (a + 1 - n)) * R.of(g[a][n]) for a in range(5)) % p for i in range(5)] for j in range(5)]
        M.append(C); Mi.append([[(C[j][i] * (1 if n % 2 == 0 else p - 1)) % p for j in range(5)] for i in range(5)])
    return M, Mi


# ---------------- frames of T_v at general v ----------------
def pairing(vec, Y=None, sign=True):
    """non-negative part of <Y^C(-z), n(z)>:  {J >= 0: (coefficient, v-exponent)} for n = {(k, c): (x, e)} (x v^e z^k phi_c)"""
    Y = Y if Y is not None else yc_monomials(60)
    out = {}
    for (k, c), (x, e) in vec.items():
        for (b, j, n, cf) in Y:
            if b != c or n > k: continue
            J = k - n; val = x * cf * ((-1) ** n if sign else 1); E = e + j
            if J in out:
                assert out[J][1] == E or out[J][0] == 0 or val == 0, "inhomogeneous pairing"
                out[J] = (out[J][0] + val, E)
            else: out[J] = (val, E)
    return {J: t for J, t in out.items() if t[0] != 0}


def frame(k, a, Y=None, sign=True):
    """frame(k, a) = z^k phi_a + Sum_{j<k} f_j v^(e_j) z^j phi_1 in T_v: {(k', c): (x, e)}"""
    assert a != 1
    vec = {(k, a): (Fr(1), Fr(0))}
    for _ in range(4 * k + 10):
        P = pairing(vec, Y, sign)
        if not P: break
        J = max(P); x, E = P[J]
        # <Y^C(-z), z^J phi_1> = -z^J + lower:  add x v^E z^J phi_1
        old = vec.get((J, 1), (Fr(0), E))
        assert old[0] == 0 or old[1] == E
        vec[(J, 1)] = (old[0] + x, E)
    assert not pairing(vec, Y, sign), "elimination did not terminate"
    for (kk, c), (x, e) in vec.items():          # homogeneity: weight of every term = w(k, a)
        assert e + w(kk, c) == w(k, a) and e >= 0 and e.denominator == 1, ("not homogeneous", k, a, kk, c, e)
    return {kc: t for kc, t in vec.items() if t[0] != 0}


def frames_upto(wmax, Y=None, sign=True):
    out = {}
    for a in (0, 2, 3, 4):
        k = 0
        while w(k, a) <= wmax:
            out[(k, a)] = frame(k, a, Y, sign); k += 1
    return out


# ---------------- correlators ----------------
class Corr(object):
    """flat ancestor correlators <prod tau_k(phi_a)>_h of the frozen model (v = 1) via MVertex6 with the referee's own M_0"""
    def __init__(self, p, kmax=6):
        sys.path.insert(0, PROJ)
        from fq_fp import FpRing
        from fq3_mvertex import MVertex6, MConst
        self.R = FpRing(p); self.p = p; R = self.R
        M, Mi = lg_M0(R, kmax)
        self.M = MConst(R, M, Mi, check=True)
        self.mv = MVertex6(R, self.M, kmax=kmax); self.kmax = kmax
        self.c = [[(NSG[a] * R.z(i * (a + 1))) % p for i in range(5)] for a in range(5)]
        self.cache = {}

    def deficit(self, h, ins): return 3 * h - 3 + len(ins) - sum(k for k, a in ins)

    def flat(self, h, ins):
        ins = tuple(sorted(ins)); key = (h, ins)
        if key in self.cache: return self.cache[key]
        n = len(ins); D = self.deficit(h, ins); p = self.p
        if 2 * h - 2 + n <= 0 or D < 0: val = 0
        elif D > self.kmax: raise ValueError("deficit")
        else:
            T = self.mv.tensor(h, tuple(k for k, a in ins)); val = 0
            for idx in range(5 ** n):
                t = T[idx]
                if not t: continue
                r = idx; cf = t
                for l in range(n):
                    cf = cf * self.c[ins[l][1]][r % 5] % p; r //= 5
                val = (val + cf) % p
        self.cache[key] = val
        return val

    def multi(self, h, vecs, stats=False):
        """<n_1 .. n_r>_h (multilinear, v = 1) and the number of non-zero flat terms (with stats: also the number of flat terms
        that are dimensionally allowed, deficit >= 0, and the number of those that are not killed by the string equation alone)"""
        p = self.p; tot = 0; nzt = 0; nallow = 0; nnontriv = 0
        for combo in itertools.product(*[list(v.items()) for v in vecs]):
            cf = Fr(1)
            for (_, (x, e)) in combo: cf *= x
            ins = [kc for (kc, _) in combo]
            if stats and 2 * h - 2 + len(ins) > 0 and self.deficit(h, ins) >= 0:
                nallow += 1
                # string-trivial: a tau_0(1) insertion while every other insertion has psi-degree 0 (and the list is not (0, 3))
                if not ((0, 0) in ins and all(k == 0 for k, a in ins) and (h, len(ins)) != (0, 3)): nnontriv += 1
            val = self.flat(h, ins)
            if val:
                term = self.R.of(cf) * val % p
                if term: nzt += 1; tot = (tot + term) % p
        if stats: return tot, nzt, nallow, nnontriv
        return tot, nzt

    def fmt(self, x):
        r = self.R.ratrec(x % self.p)
        return str(r) if r is not None else "[%d mod p]" % (x % self.p)


def singular_lists(h, NT, frames, emax=1):
    """[(frame keys, m, e)] with n >= 1 frames, m tau_0(H), n + m <= NT, integer e = 2 - 2h - m + Sum w <= emax"""
    keys = sorted(frames); out = []
    for ntot in range(1, NT + 1):
        if 2 * h - 2 + ntot <= 0: continue
        for m in range(0, ntot):
            n = ntot - m
            for combo in itertools.combinations_with_replacement(keys, n):
                e = 2 - 2 * h - m + sum(w(*c) for c in combo)
                if e.denominator != 1 or e > emax: continue
                out.append((combo, m, int(e)))
    return out
