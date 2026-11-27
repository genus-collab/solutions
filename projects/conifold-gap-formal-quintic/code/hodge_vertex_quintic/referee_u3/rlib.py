# -*- coding: utf-8 -*-
"""referee helpers: Fourier (flat-basis) coordinates of Z/5-covariant sp series, adjoint actions, loading of the stored data"""
import sys, os, json
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); sys.path.insert(0, PROJ)
HERE = os.path.dirname(os.path.abspath(__file__))
from fq_symbol import mat_series_mul, mat_series_exp, mat_series_log
from fq3_m import sp_mats, sp_params, in_sp, NSP, KMAX
from fq3_circsys import CircSys
from fq3_taut import circulant_coords


def load(tag="61"):
    S = CircSys(tag); d = json.load(open(os.path.join(HERE, "m0_%s.json" % tag)))
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    return S, d, cs


def zeros(): return [[[0] * 5 for _ in range(5)] for _ in range(KMAX + 1)]


def ser_from_circ(S, v15):
    """15 circulant coordinates -> sp series (no N_1^taut added)"""
    p = S.p
    return sp_mats([sum(S.cols[j][k] * v15[j] for j in range(15) if v15[j]) % p for k in range(NSP)], p)


def circ_of(S, N):
    t, names = circulant_coords(S.R, N)
    return t


def mul(X, Y, p): return mat_series_mul(X, Y, p, KMAX)


def Ad(M, Mi, X, p):
    """M X M^-1"""
    return mul(mul(M, X, p), Mi, p)


def add(X, Y, p, c=1): return [[[(X[k][i][j] + c * Y[k][i][j]) % p for j in range(5)] for i in range(5)] for k in range(KMAX + 1)]
def scale(X, c, p): return [[[x * c % p for x in r] for r in Xk] for Xk in X]
def is_zero(X, p): return not any(x % p for Xk in X for r in Xk for x in r)


def fourier(S, N):
    """hat c_k(m): N_k phi_m = hat c_k(m) phi_{m-k}, (phi_m)_j = zeta^{m j}; returns ({(k, m): value}, covariant?)"""
    p = S.p; z = S.R.z; out = {}; ok = True
    for k in range(KMAX + 1):
        for m in range(5):
            vec = [sum(N[k][i][j] * z(m * j) for j in range(5)) % p for i in range(5)]
            c = vec[0]; out[(k, m)] = c
            if any((vec[i] - c * z((m - k) * i)) % p for i in range(5)): ok = False
    return out, ok


def from_fourier(S, F):
    """inverse of fourier: (N_k)_{ij} = (1/5) Sum_m hat c_k(m) zeta^{(m-k) i} zeta^{-m j}"""
    p = S.p; z = S.R.z; inv5 = pow(5, p - 2, p); N = zeros()
    for (k, m), c in F.items():
        if not c % p: continue
        for i in range(5):
            for j in range(5): N[k][i][j] = (N[k][i][j] + c * z((m - k) * i - m * j) * inv5) % p
    return N


def show(S, name, N, rec=None, kmax=6):
    F, ok = fourier(S, N); p = S.p
    print("%s  (covariant: %s; in sp: %s)" % (name, ok, in_sp(N, p)))
    for k in range(1, kmax + 1):
        row = []
        for m in range(5):
            c = F[(k, m)]
            row.append("0" if c == 0 else (rec.fmt(c) if rec else "*"))
        print("   k=%d: m=0..4 (H^a -> H^(a-k), a = m-1): %s" % (k, row))
    return F
