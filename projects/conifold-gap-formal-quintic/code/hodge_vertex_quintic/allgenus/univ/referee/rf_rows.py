# -*- coding: utf-8 -*-
"""rf_rows.py -- REFEREE (univ): symbol rows of the LG model from the closed recursion (LG_MODEL.md s.3, own implementation),
hatted rows Qhat = Ghat(u, zhat) Rhat(tau zhat) from a flat-frame gauge p[k][b] (rf_gauge), left dressings in hatted variables.

Rhat_a^(i)[m] = N_a zeta^(i(a+1-m)) g^a_m,  N = (1,-1,-1,1,1),  c = (0,-1,3,2,1);
Ghat_{b+k, b}[zhat^k] = tau^((1-c_b)/4) u^mu(b,k) p[k][b](u),  mu(b,k) = (4k + c_b - c_{b+k})/5;
Qhat_a^(i)[n] = Sum_{k <= n} Ghat_{a, a-k}[zhat^k] tau^(n-k) Rhat_{a-k}^(i)[n-k].
Rows are returned mod p as Q[a][i][n] = list of K residues (u^0..u^(K-1)).  The only project import is fq_fp (the primes and the
project's choice of zeta, so that rows can be compared residue by residue with the project's rows).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, HERE)
from rf_series import LS, tau_series

CEXP = (0, -1, 3, 2, 1); NS = (1, -1, -1, 1, 1)


def mu(b, k):
    v = 4 * k + CEXP[b % 5] - CEXP[(b + k) % 5]; assert v % 5 == 0; return v // 5


def prime_zeta(name):
    sys.path.insert(0, PROJ)
    from fq_fp import FpRing, P1, P2
    p = {"89": P2, "61": P1}[name]; R = FpRing(p)
    z = R.z(1) % p; assert pow(z, 5, p) == 1 and z != 1
    return p, z


def symbol_g(nmax):
    g = [[Fr(1)] for _ in range(5)]
    for n in range(1, nmax + 1):
        d = [Fr(0)]
        for a in range(4): d.append(d[a] + Fr(4 * n - 3 - CEXP[a], 5) * g[a][n - 1])
        assert d[4] + Fr(4 * n - 3 - CEXP[4], 5) * g[4][n - 1] == 0, "cycle does not close"
        t = -sum((4 * n + 1 - CEXP[a]) * d[a] for a in range(5)) / (20 * n)
        for a in range(5): g[a].append(t + d[a])
    return g


def red(x, p): return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def lsred(s, K, p):
    assert s.v() >= 0 or s.iszero(), "pole"
    assert s.prec >= K, "precision %d < %d" % (s.prec, K)
    return [red(s.coef(n), p) for n in range(K)]


def smul(a, b, p, K):
    out = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i): out[i + j] += a[i] * b[j]
    return [x % p for x in out]


def hatted_gauge(pg, tau, K, KZ, p):
    """{(a, b, k): list of K residues} for k <= KZ"""
    tp = {b: tau.powr(Fr(1 - CEXP[b], 4)) for b in range(5)}
    out = {}
    for k in range(KZ + 1):
        for b in range(5):
            m = mu(b, k)
            if m >= K: continue
            s = tp[b].mul(pg[k][b]).trunc(K - m + 1).shift(m)
            r = lsred(s, K, p)
            if any(r): out[((b + k) % 5, b, k)] = r
    return out


def rows(Gh, tau, K, KZ, p, zeta):
    """Q[a][i][n] from hatted left factor Gh {(a, b, k): series} and the symbol at tau zhat"""
    g = symbol_g(KZ)
    t1 = lsred(tau, K, p); tpw = [[1] + [0] * (K - 1)]
    for m in range(1, KZ + 1): tpw.append(smul(tpw[-1], t1, p, K))
    Q = [[[[0] * K for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (a, b, k), h in Gh.items():
        for n in range(k, KZ + 1):
            ht = smul(h, tpw[n - k], p, K)
            for i in range(5):
                c = NS[b] * pow(zeta, (i * (b + 1 - (n - k))) % 5, p) * red(g[b][n - k], p) % p
                if c:
                    q = Q[a][i][n]
                    for j in range(K): q[j] = (q[j] + ht[j] * c) % p
    return Q


def dress(Q, D, K, KZ, p):
    """rows -> D rows, D = {(a, b, k): series (list of K)} INCLUDING the identity"""
    out = [[[[0] * K for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (a, b, k), d in D.items():
        for i in range(5):
            for n in range(k, KZ + 1):
                t = smul(d, Q[b][i][n - k], p, K); q = out[a][i][n]
                for j in range(K): q[j] = (q[j] + t[j]) % p
    return out


def project_rows(prime, g, K, ab=(0, 0)):
    """the project's hatted rows (fq_rows -> fq_conifold -> fq_sser -> fq3_msum.HatRowsU), read-only import"""
    sys.path.insert(0, PROJ)
    from fq_fp import FpRing, P1, P2
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq3_msum import HatRowsU
    p = {"89": P2, "61": P1}[prime]; KZ = 3 * g - 3; R = FpRing(p)
    rows_ = Rows(R, Nz=KZ); C = Conifold(R, Nu=K + 4, alpha=Fr(ab[0]), beta=Fr(ab[1])); S = SRows(rows_, C, prec_abs=5 * K + 4)
    H = HatRowsU(S, K); assert H.val_ok and H.min_prec >= K
    return [[[[H.Q[a][i][n].coef(k) % p for k in range(K)] for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
