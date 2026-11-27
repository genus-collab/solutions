# -*- coding: utf-8 -*-
"""xy_tr.py -- direct Eynard-Orantin topological recursion (own code) for x = p^4 + 4/p, y = p, B = dp dp'/(p - p')^2.

Convention (EO, arXiv:math-ph/0702045): omega_{0,1} = y dx, omega_{0,2} = B,
   omega_{g,n+1}(p0, S) = Sum_i Res_{p -> a_i} K_i(p0, p) [omega_{g-1,n+2}(p, sig p, S) + Sum' omega(p, I) omega(sig p, J)],
   K_i(p0, p) = (Int_{sig p}^{p} B(p0, .)) / (2 (y(p) - y(sig p)) dx(p)).
Representation (2g-2+n > 0): omega_{g,n} = Sum C[(i_1,k_1),..,(i_n,k_n)] prod_j e_{i_j,k_j}(p_j),  e_{i,k}(p) = dp/(p - a_i)^(k+1),
a_i = zeta^i, k >= 1 (poles at the branch points only, no residues).  Residue at a_0 only; the other branch points by the
symmetry p -> zeta p (x -> zeta^4 x, y -> zeta y):  C[(i_j + r, k_j)] = zeta^(r Sum_j k_j) C[(i_j, k_j)].
Field F = xy_field.QZ (exact, Q(zeta_5)) or xy_field.make_fp(p).  Nothing here writes a file.
"""
import itertools
from fractions import Fraction as Fr
from math import comb
from xy_local import LocalQ, lconv, lmul, ladd_into, acc_to_l


class TR(object):
    def __init__(self, F, gmax, nmax):
        self.F = F; self.zero = F.zero()
        kmax = max(6 * g - 5 + 2 * n for g in range(gmax + 1) for n in range(1, nmax + 1) if 2 * g - 2 + n > 0)
        self.kmax = kmax
        self.Ptop = kmax + 3                       # absolute order to which first-variable expansions are kept
        self.N = 2 * kmax + 12
        self.L = LocalQ(self.N, kmax + 3)
        self.om = {}
        self._ep = {}; self._es = {}; self._pp = {}
        self.zpow = [F.zeta(j) for j in range(5)]

    # ---- expansions of e_{i,k} at p = 1 + t and at sigma(p), divided by dt, truncated at t^Ptop
    def e_p(self, i, k):
        if (i, k) not in self._ep:
            F = self.F; top = self.Ptop
            if i == 0: s = (-(k + 1), [F.one()])
            else:
                w = F.inv1mz(i); wp = F.one()
                for _ in range(k + 1): wp = wp * w
                out = []
                for n in range(top + 1):
                    out.append(wp.scale(Fr((-1) ** n * comb(k + n, n)))); wp = wp * w
                s = (0, out)
            self._ep[(i, k)] = s
        return self._ep[(i, k)]

    def e_s(self, i, k):
        if (i, k) not in self._es:
            F = self.F; L = self.L; top = self.Ptop + k + 2
            if i == 0:
                from xy_local import ps_mul
                s = lconv(F, (-(k + 1), ps_mul(L.dsig, L.su_neg(k + 1), L.N)), self.Ptop)
            else:
                w = F.inv1mz(i); wp = F.one()
                for _ in range(k + 1): wp = wp * w
                acc = [self.zero] * (self.Ptop + 1)
                for n in range(self.Ptop + 1):
                    cf = wp.scale(Fr((-1) ** n * comb(k + n, n))); wp = wp * w
                    sp = L.sp[n]
                    for j in range(n, self.Ptop + 1):
                        if sp[j]: acc[j] = acc[j] + cf.scale(sp[j])
                s = (0, acc)
            self._es[(i, k)] = s
        return self._es[(i, k)]

    def expand(self, g, n, at_sigma):
        """omega_{g,n}(p or sigma p, rest) as {rest_key: Laurent}; (g, n) = (0, 2) means B(., s)"""
        F = self.F; L = self.L; out = {}
        if (g, n) == (0, 2):
            for m in range(self.Ptop + 1):
                if at_sigma: out[((0, m + 1),)] = lconv(F, (0, [Fr(m + 1) * a for a in L.sp[m]]), self.Ptop)
                else: out[((0, m + 1),)] = (m, [F.frac(m + 1)])
            return out
        acc = {}
        for key, cf in self.om[(g, n)].items():
            ser = self.e_s(*key[0]) if at_sigma else self.e_p(*key[0])
            d = acc.setdefault(key[1:], {})
            ladd_into(d, ser, cf)
        for r, d in acc.items(): out[r] = acc_to_l(d, self.zero, self.Ptop)
        return out

    def pair(self, k1, k2):
        if (k1, k2) not in self._pp:
            self._pp[(k1, k2)] = lmul(self.e_p(*k1), self.e_s(*k2), 0, self.zero)
        return self._pp[(k1, k2)]

    def step(self, g, n1):
        """compute omega_{g, n1} (n1 = n + 1 >= 1)"""
        F = self.F; L = self.L; n = n1 - 1; S = tuple(range(n)); R = {}
        if g >= 1:
            if (g - 1, n + 2) == (0, 2):
                ladd_into(R.setdefault((), {}), lconv(F, L.Bss, 0))
            else:
                for key, cf in self.om[(g - 1, n + 2)].items():
                    ladd_into(R.setdefault(key[2:], {}), self.pair(key[0], key[1]), cf)
        for g1 in range(g + 1):
            g2 = g - g1
            for sz in range(n + 1):
                for I in itertools.combinations(S, sz):
                    J = tuple(s for s in S if s not in I)
                    if (g1, len(I)) == (0, 0) or (g2, len(J)) == (0, 0): continue
                    E1 = self.expand(g1, len(I) + 1, False); E2 = self.expand(g2, len(J) + 1, True)
                    for r1, s1 in E1.items():
                        for r2, s2 in E2.items():
                            rest = [None] * n
                            for pos, var in enumerate(I): rest[var] = r1[pos]
                            for pos, var in enumerate(J): rest[var] = r2[pos]
                            ladd_into(R.setdefault(tuple(rest), {}), lmul(s1, s2, 0, self.zero))
        out = {}; mmax = 6 * g - 4 + 2 * n1 + 1
        for rest, d in R.items():
            Rl = acc_to_l(d, self.zero, 0)
            for m in range(1, mmax + 1):
                kv, kc = L.kap[m]
                v = self.zero
                for j, a in enumerate(Rl[1]):
                    e = Rl[0] + j; jj = -1 - e - kv
                    if 0 <= jj < len(kc) and kc[jj]: v = v + a.scale(kc[jj])
                if not v.is_zero(): out[((0, m),) + rest] = v
        assert all(k[0][1] <= 6 * g - 5 + 2 * n1 for k in out), "pole order above 6g-4+2n"
        full = {}
        for key, cf in out.items():
            ksum = sum(k for (_, k) in key)
            for r in range(5):
                nk = tuple(((i + r) % 5, k) for (i, k) in key)
                full[nk] = cf * self.zpow[(r * ksum) % 5]
        self.om[(g, n1)] = full
        return full

    def run(self, order):
        for (g, n1) in order: self.step(g, n1)

    def value(self, g, n1, pts):
        """omega_{g,n}/(dp_1..dp_n) at rational points (Fractions, p^5 != 1)"""
        F = self.F; tot = self.zero
        cache = {}
        for key, cf in self.om[(g, n1)].items():
            t = cf
            for (i, k), p in zip(key, pts):
                if (i, k, p) not in cache:
                    b = F.inv_a_minus_zeta(p, i); x = F.one()
                    for _ in range(k + 1): x = x * b
                    cache[(i, k, p)] = x
                t = t * cache[(i, k, p)]
            tot = tot + t
        return tot

    def free_energy(self, g):
        """F_g = 1/(2-2g) Sum_i Res_{a_i} Phi omega_{g,1}, Phi = (4/5) p^5 - 4 log p  (g >= 2)"""
        tot = self.zero
        for key, cf in self.om[(g, 1)].items():
            (i, k), = key
            if i == 0:
                ph = Fr(4, 5) * comb(5, k) - Fr(4 * (-1) ** (k + 1), k)
                tot = tot + cf.scale(ph)
        return tot.scale(Fr(5, 2 - 2 * g))

    def symmetric_failures(self, g, n1):
        om = self.om[(g, n1)]; bad = 0
        for key, cf in om.items():
            for perm in itertools.permutations(range(n1)):
                k2 = tuple(key[j] for j in perm)
                if k2 != key and not (om.get(k2, self.zero) - cf).is_zero(): bad += 1
        return bad
