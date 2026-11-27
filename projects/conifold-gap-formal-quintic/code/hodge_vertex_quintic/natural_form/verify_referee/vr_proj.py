# -*- coding: utf-8 -*-
"""vr_proj.py -- VERIFIER (adversarial check of natural_form/referee): the factorisation R~ = G M_0(tau zhat) recomputed from the
PROJECT row machinery (fq_rows.Rows -> fq_conifold -> fq_sser.SRows -> fq3_msum.HatRowsU; ring-N WKB, s-series) pushed BEYOND the
orders of the original test (Nz > 6), with the verifier's own matrix-series arithmetic.  Nothing is imported from referee/.
The referee's rows come from a different derivation (hatted Picard-Fuchs recursion in (zhat, u)); agreement of the two at Nz > 6
is an independent confirmation of the referee's new statements (deg G_5 = 6, deg G_6 = 8, rule W, chain relation).

Checks printed:
  (1) degrees of G_n; number of vanishing coefficients above floor((5n+2)/4)
  (2) every Fourier component of G against the referee's rational table G_fourier_NZ16_K12.json (read-only), both directions
  (3) rule W as an IFF: every Z/5-allowed component with weight >= 0 is nonzero, every one with weight < 0 vanishes
  (4) chain relation at weight 0, counted from the verifier's own G
  (5) M_0 on the left: candidates for c_1 at the two top zhat-orders, per Fourier component
  (6) A = M_0 Z M_0^-1 is not a polynomial within the truncation (lemma behind the uniqueness proof)
Run: py -I -B vr_proj.py --prime 89 --Nz 10 --K 8 --prec 46 --Nu 14
"""
import os, sys, time, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(os.path.dirname(HERE)); sys.path.insert(0, PROJ)

W = (0, 4, 3, 2, 1); PSI = (1, -1, -1, 1, 1)


def mmul(A, B, p): return [[sum(A[i][k] * B[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
def madd(A, B, p, c=1): return [[(A[i][j] + c * B[i][j]) % p for j in range(5)] for i in range(5)]
def mzero(): return [[0] * 5 for _ in range(5)]
def mid(): return [[int(i == j) for j in range(5)] for i in range(5)]
def nz(A): return any(x for r in A for x in r)


def bmul(X, Y, p, NZ, K):
    out = {}
    for (n1, k1), A in X.items():
        for (n2, k2), B in Y.items():
            if n1 + n2 > NZ or k1 + k2 >= K: continue
            key = (n1 + n2, k1 + k2); C = mmul(A, B, p)
            out[key] = madd(out[key], C, p) if key in out else C
    return {k: v for k, v in out.items() if nz(v)}


def tau_series(K, p):
    F = [Fr(4, 5 * k + 4) for k in range(K)]; t = [Fr(1)]
    for k in range(1, K): t.append(-sum(F[j] * t[k - j] for j in range(1, k + 1)))
    return t, [x.numerator % p * pow(x.denominator, p - 2, p) % p for x in t]


def fourier(X, z, p):
    i25 = pow(25, p - 2, p); out = {}
    for m in range(5):
        for mp in range(5):
            v = sum(pow(z, (m * i + mp * j) % 5, p) * X[i][j] for i in range(5) for j in range(5)) % p
            if v: out[(m, mp)] = v * i25 % p
    return out


def build(p, Nz, K, prec, Nu, alpha=0, beta=0):
    from fq_fp import FpRing
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq3_msum import HatRowsU
    R = FpRing(p); rows = Rows(R, Nz=Nz); C = Conifold(R, Nu=Nu, alpha=Fr(alpha), beta=Fr(beta)); S = SRows(rows, C, prec_abs=prec)
    H = HatRowsU(S, K)
    assert H.val_ok and H.min_prec >= K, (H.val_ok, H.min_prec)
    Rt = {}
    for n in range(Nz + 1):
        for k in range(K):
            Mx = mzero()
            for j in range(5):
                for i in range(5):
                    acc = 0
                    for (a, b), e in H.eta_inv.items():
                        s = H.Q[b][j][0].mul(H.Q[a][i][n]); assert s.prec > k
                        acc += e * s.coef(k)
                    Mx[j][i] = acc % p
            if nz(Mx): Rt[(n, k)] = Mx
    return R, H, Rt
