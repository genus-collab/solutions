# -*- coding: utf-8 -*-
"""un_ann.py -- construction of finite (non-infinitesimal) left dressings D = exp(X) of the natural R' with prescribed properties.

All objects are in HATTED variables, entries are u-series mod u^K (un_core.Arith with L = 1):
   X = {(a, b, k): series},  a = b + k mod 5 (Z/5 grading);   flat-regular  <=>  X[(a, b, k)] is divisible by u^mu(b, k);
   infinitesimally symplectic  <=>  X_{k,b} = -(-1)^k X_{k,b*}, b* = 3 - b - k  (entry (k, col) := X[(col + k, col, k)]).
Hatted vanishing covector of the quintic:  Yhat^C = Sum_{n>=0} zhat^-n yhat_n(u) e_{(1-n) mod 5},  yhat_n = y_n / u^(n - floor(n/5))
(un_lin.yc_series), yhat_0 = -1; all yhat_n are units.  For a hatted matrix D:
   (D Yhat^C) at zhat^m, component e_{m+1}:   delta_m(D) = Sum_{n >= max(0,-m)} D[(m + 1, (1 - n) mod 5, m + n)] yhat_n .
Hypothesis (A) for D:  delta_m(D) = 0 for all m >= 1.  The flat-frame transported covector D Y^C_Q has
   z^0 coefficient  delta_0 e_1   and   z^-1 coefficient  u delta_{-1} e_0,   so   T_eff := -u delta_{-1}(D)   (= t_c for D = 1).

random_X(...)   random flat-regular symplectic X with every hatted entry O(u) (so exp(X) is a polynomial mod u^K);
solve_ann(...)  adjusts the column-1 entries (and their symplectic partners) of X until delta_m(exp X) = 0 mod u^K for m = 1..KZ
                (u-adic Newton iteration; the corrections are flat-regular automatically, which is asserted).
"""
import os, sys, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from un_core import mu, mat_exp_left
from un_lin import orbit_of


def yhat_series(ys, R, K):
    """yhat_n = y_n / u^nu_n as lists of K residues, nu_n = n - floor(n/5)"""
    out = []
    for n, y in enumerate(ys):
        nu = n - n // 5
        assert all(c == 0 for c in y[:nu]) and len(y) >= nu + K, "y_%d: valuation or precision" % n
        out.append([R.of(c) for c in y[nu:nu + K]])
    return out


def set_entry(ar, X, k, col, ser):
    """set the entry (k, col) and its symplectic partner"""
    bs = (3 - col - k) % 5
    if bs == col:
        assert k % 2 == 1 or ser is None, "self-dual entry at even k must vanish"
        X[((col + k) % 5, col, k)] = ser
    else:
        X[((col + k) % 5, col, k)] = ser
        X[((bs + k) % 5, bs, k)] = ser if k % 2 == 1 else ar.neg(ser)


def random_X(ar, R, KZ, rng, flat_regular=True, with_k0=False, kmin_val=1, exclude_col1=True, deg=None):
    """random symplectic graded X; coefficient of u^j of the entry (k, b): a random small rational, for
    j >= max(kmin_val, mu(b, k)) if flat_regular, j >= kmin_val otherwise (hatted-regular only).  Column-1 orbits are left out
    if exclude_col1 (they are then fixed by solve_ann)."""
    K = ar.KU; X = {}
    for k in range(0 if with_k0 else 1, KZ + 1):
        seen = set()
        for b in range(5):
            if b in seen: continue
            bs = (3 - b - k) % 5; seen.add(b); seen.add(bs)
            if bs == b and k % 2 == 0: continue
            if exclude_col1 and (b == 1 or bs == 1): continue
            j0 = max(kmin_val, mu(b, k)) if flat_regular else kmin_val
            if j0 >= K: continue
            coefs = [0] * K
            for j in range(j0, K if deg is None else min(K, j0 + deg + 1)):
                coefs[j] = R.of(Fr(rng.randint(-9, 9), rng.randint(1, 9)))
            set_entry(ar, X, k, b, ar.from_u(coefs))
    return {key: v for key, v in X.items() if v is not None}


def delta(ar, D, yh, m, KZ):
    """delta_m(D) for D = {(a, b, k): series} INCLUDING the identity (use with_identity)"""
    acc = None
    for n in range(max(0, -m), len(yh)):
        k = m + n
        if k > KZ: break
        s = D.get(((m + 1) % 5, (1 - n) % 5, k))
        if s is None: continue
        acc = ar.add(acc, ar.mul(s, ar.from_u(yh[n])))
    return acc


def with_identity(ar, Xm1):
    D = dict(Xm1)
    for a in range(5): D[(a, a, 0)] = ar.add(D.get((a, a, 0)), ar.one())
    return D


def expX(ar, X, KZ):
    return with_identity(ar, mat_exp_left(ar, X, KZ, ar.KU + KZ + 2))


def solve_ann(ar, X, yh, KZ, log=None):
    """returns X' = X + (column-1 corrections) with delta_m(exp X') = 0 mod u^K for 1 <= m <= KZ; asserts flat-regularity of X'"""
    K = ar.KU; X = dict(X)
    for it in range(4 * K + 4):
        D = expX(ar, X, KZ); dl = {m: delta(ar, D, yh, m, KZ) for m in range(1, KZ + 1)}
        if all(v is None for v in dl.values()): break
        dX = {}
        for m in range(KZ, 0, -1):
            acc = dl[m]
            for n in range(1, len(yh)):
                k = m + n
                if k > KZ: break
                s = dX.get(((m + 1) % 5, (1 - n) % 5, k))
                if s is not None: acc = ar.add(acc, ar.mul(s, ar.from_u(yh[n])))
            if acc is None: continue
            if (2 - m) % 5 == 1 and m % 2 == 0:
                raise AssertionError("delta_%d != 0 cannot be corrected (self-dual column-1 entry at even order)" % m)
            set_entry(ar, dX, m, 1, acc)
        for key, v in dX.items(): X[key] = ar.add(X.get(key), v)
        X = {key: v for key, v in X.items() if v is not None}
    else:
        raise AssertionError("solve_ann did not converge")
    if log: log("      solve_ann: converged after %d iteration(s)" % it)
    bad = [key for key, v in X.items() if ar.val(v) < mu(key[1], key[2])]
    return X, bad


def flat_irregular(ar, X):
    """entries violating divisibility by u^mu(b, k)"""
    return [key for key, v in X.items() if v is not None and ar.val(v) < mu(key[1], key[2])]


def lie_symplectic_defect(ar, X):
    bad = []
    for (a, b, k), v in X.items():
        bs = (3 - b - k) % 5; w = X.get(((bs + k) % 5, bs, k))
        tgt = v if k % 2 == 1 else ar.neg(v)
        if bs == b:
            if k % 2 == 0 and v is not None: bad.append((a, b, k))
        elif w != tgt: bad.append((a, b, k))
    return bad


def teff(ar, D, yh, KZ, R):
    """T_eff/u = -delta_{-1}(D) as a list of K residues, and delta_0(D) (the z^0 coefficient of the transported covector)"""
    d1 = delta(ar, D, yh, -1, KZ); d0 = delta(ar, D, yh, 0, KZ)
    return [(-x) % ar.p for x in ar.part(d1, 0)], ar.part(d0, 0)
