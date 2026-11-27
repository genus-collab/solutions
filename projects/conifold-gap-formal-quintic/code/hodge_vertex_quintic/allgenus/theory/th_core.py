# -*- coding: utf-8 -*-
"""th_core.py -- flat-frame ancestor correlators of the frozen vertex theory (= the LG model x^4/4 + v/x at v = 1, UNIFORM_REDUCTION
Prop. S) and the frame of the subspace T_v of REDUCTION_DRAFT.md s.2 (exact; read-only use of the project modules).

Model at v = 1, flat frame phi_a = H^a (a in Z/5), eta(phi_a, phi_b) = 5 delta_{a+b=3}, H * phi_a = bmod_a phi_{a+1},
bmod = (1, -1/v, 1, -1, -1).  In the normalised canonical frame e_i of the project (colour i):
    phi_a = Sum_i c[a][i] e_i,   c[a][i] = N_a zeta^(i(a+1)),  N = (1, -1, -1, 1, 1)      (zhat^0 part of the conifold symbol),
and <tau_k1(phi_a1) .. tau_kn(phi_an)>_g = Sum_i prod_l c[a_l][i_l] * MVertex6(M_0).tensor(g, (k_1..k_n))[i_1..i_n].
Vanishing solution (residue at the pole), evaluated at -z, v = 1:
    Y^C_a(-z) = SGN_a I_{MEXP_a}(-z),  I_m(-z) = Sum_{i>=0} z^-(m+5i+1) / (4^i i! (m+4i+1)!)    (terms with m+4i+1 < 0 absent).
T_v = { n in H[z] : (Y^C(-z), n(z)) has no non-negative power of z } = common kernel of l_k = Omega(z^-k Y^C, .), k >= 1.
A vector of H[z] is a dict {(k, a): Fraction} = Sum coefficient z^k phi_a.
Nothing here writes a file.
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_par import setup
from fq3_mvertex import MVertex6

MEXP = (0, -1, 3, 2, 1)
SGN = (1, -1, -1, 1, 1)
NSIGN = (1, -1, -1, 1, 1)
CEXP = (0, -1, 3, 2, 1)          # hatted weights c_a


def fact(n):
    out = 1
    for k in range(2, n + 1): out *= k
    return out


def yc_minus(a, nmax):
    """{n: coefficient of z^-n} of Y^C_a(-z) at v = 1, n <= nmax"""
    m = MEXP[a]; out = {}
    i = 0
    while True:
        j = m + 4 * i + 1
        if j >= 0:
            n = m + 5 * i + 1
            if n > nmax: break
            out[n] = Fr(SGN[a], 4 ** i * fact(i) * fact(j))
        i += 1
        if i > nmax + 2: break
    return out


def pair_nonneg(vec):
    """{j >= 0: coefficient of z^j} of (Y^C(-z), n(z)) for n = vec"""
    out = {}
    for (k, a), cf in vec.items():
        for n, y in yc_minus(a, k).items():
            j = k - n
            if j >= 0: out[j] = out.get(j, Fr(0)) + cf * y
    return {j: x for j, x in out.items() if x != 0}


def proj_T(vec):
    """n -> n - Sum_j l_j(n) (z^(j-1) H + ..): the projection of n to T_1 along H C[z] (adds multiples of z^j H only)"""
    vec = dict(vec)
    for _ in range(200):
        P = pair_nonneg(vec)
        if not P: return {k: x for k, x in vec.items() if x != 0}
        j = max(P)
        # (Y^C(-z), z^j H) = -z^j + lower: add P[j] z^j H
        vec[(j, 1)] = vec.get((j, 1), Fr(0)) + P[j]
    raise RuntimeError("proj_T did not terminate")


def frame(k, a):
    """the frame vector of T_1 with leading term z^k phi_a (a != 1)"""
    assert a != 1
    return proj_T({(k, a): Fr(1)})


def weight5(k, a):
    """5 x (u-exponent of an insertion tau_k(phi_a)) = 4k + c_a - 4"""
    return 4 * k + CEXP[a] - 4


def vec_str(vec):
    names = ("1", "H", "H^2", "H^3", "H^4")
    return " + ".join("(%s) z^%d %s" % (cf, k, names[a]) for (k, a), cf in sorted(vec.items()))


class Model(object):
    def __init__(self, p=P2):
        st = setup((p, "0", "0", 34, 6))
        self.R = st['R']; self.p = self.R.p; R = self.R
        self.mv = MVertex6(R, st['M0'])
        self.c = [[(NSIGN[a] * R.z(i * (a + 1))) % self.p for i in range(5)] for a in range(5)]
        self.cache = {}

    def corr(self, g, ins):
        """<prod tau_k(phi_a)>_g for ins = [(k, a), ..] (F_p); 0 if unstable or of the wrong dimension"""
        ins = tuple(sorted(ins)); key = (g, ins)
        r = self.cache.get(key)
        if r is not None: return r
        n = len(ins); p = self.p
        d = tuple(k for k, a in ins); D = 3 * g - 3 + n - sum(d)
        if 2 * g - 2 + n <= 0 or D < 0:
            self.cache[key] = 0; return 0
        if D > self.mv.kmax: raise ValueError("psi-deficit %d out of range for %s" % (D, key))
        T = self.mv.tensor(g, d)
        cur = T
        for l in range(n - 1, -1, -1):          # contract the last index first (index = Sum i_l 5^l)
            a = ins[l][1]; size = 5 ** l
            cur = [sum(self.c[a][i] * cur[i * size + r] for i in range(5)) % p for r in range(size)]
        val = cur[0] % p
        self.cache[key] = val
        return val

    def mcorr(self, g, vecs):
        """multilinear extension: <n_1, .., n_r>_g for vectors n_i (dicts)"""
        p = self.p; R = self.R; tot = 0
        for combo in itertools.product(*[list(v.items()) for v in vecs]):
            cf = Fr(1)
            for (_, x) in combo: cf *= x
            val = self.corr(g, [ka for (ka, _) in combo])
            tot = (tot + R.of(cf) * val) % p
        return tot

    def fmt(self, x):
        r = self.R.ratrec(x % self.p)
        return str(r) if r is not None else "[%d mod p]" % (x % self.p)
