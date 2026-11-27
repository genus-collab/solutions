# -*- coding: utf-8 -*-
"""rf_e3core.py -- REFEREE (univ): own construction of left dressings D of R' (hatted variables, entries mod p as lists of K
residues), different from the worker's: D = CAYLEY(X) = (1 - X/2)^-1 (1 + X/2) (symplectic for infinitesimally symplectic X), the
column-1 orbits fixed by a Jacobi iteration; every property is then CHECKED on D itself, independently of the construction:
  symplectic      D(-zhat)^T eta D(zhat) = eta  (eta_{a,3-a} = 1; conjugation by diag(s^(1-c_a)) preserves eta since c_a + c_{3-a} = 2)
  flat-regular    D[(a, b, k)] divisible by u^mu(b, k)
  annihilating    Delta_m = Sum_n D[(m+1, 1-n, m+n)] yhat_n = 0 mod u^K, m = 1..KZ, yhat_n = y_n/u^(n - n//5) (flat:
                  u^-mu(1,m) Delta_m is the z^m coefficient of D Y^C_Q, using mu(1-n, m+n) = mu(1, m) + n - n//5, checked here)
  T_eff = -u Delta_{-1},  delta_0 = Delta_0 (z^0 coefficient / e_1; non-constant => no z-independent target connection).
"""
import os, sys, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_rows import mu, smul, red


def madd(A, B, p):
    out = dict(A)
    for key, v in B.items():
        w = out.get(key)
        out[key] = v if w is None else [(x + y) % p for x, y in zip(w, v)]
    return {k: v for k, v in out.items() if any(v)}


def mscale(A, c, p): return {k: [x * c % p for x in v] for k, v in A.items() if c % p}


def mmul(A, B, KZ, p, K):
    out = {}
    for (a, b, k1), s in A.items():
        for (b2, c, k2), t in B.items():
            if b2 != b or k1 + k2 > KZ: continue
            key = (a, c, k1 + k2); pr = smul(s, t, p, K); w = out.get(key)
            out[key] = pr if w is None else [(x + y) % p for x, y in zip(w, pr)]
    return {k: v for k, v in out.items() if any(v)}


def ident(K): return {(a, a, 0): [1] + [0] * (K - 1) for a in range(5)}


def cayley(X, KZ, p, K):
    h = mscale(X, pow(2, p - 2, p), p); inv = ident(K); pw = ident(K)
    for _ in range(K + KZ + 2):
        pw = mmul(pw, h, KZ, p, K)
        if not pw: break
        inv = madd(inv, pw, p)
    return mmul(inv, madd(ident(K), h, p), KZ, p, K)


def set_orbit(X, k, b, ser, p):
    bs = (3 - b - k) % 5
    if bs == b:
        assert k % 2 == 1 or not any(ser); X[((b + k) % 5, b, k)] = ser; return
    X[((b + k) % 5, b, k)] = ser
    X[((bs + k) % 5, bs, k)] = ser if k % 2 == 1 else [(-x) % p for x in ser]


def random_X(rng, KZ, K, p, with_k0=False, minval=1, flat_regular=True, skip_col1=True):
    X = {}
    for k in range(0 if with_k0 else 1, KZ + 1):
        done = set()
        for b in range(5):
            bs = (3 - b - k) % 5
            if b in done: continue
            done.update((b, bs))
            if bs == b and k % 2 == 0: continue
            if skip_col1 and 1 in (b, bs): continue
            j0 = max(minval, mu(b, k)) if flat_regular else minval
            ser = [0] * K
            for j in range(j0, K): ser[j] = red(Fr(rng.randint(-9, 9), rng.randint(1, 9)), p)
            if any(ser): set_orbit(X, k, b, ser, p)
    return X


def Delta(D, yh, m, KZ, p, K):
    acc = [0] * K
    for n in range(max(0, -m), len(yh)):
        k = m + n
        if k > KZ: break
        s = D.get(((m + 1) % 5, (1 - n) % 5, k))
        if s: acc = [(x + y) % p for x, y in zip(acc, smul(s, yh[n], p, K))]
    return acc


def solve_ann(X, yh, KZ, p, K, linear=False, iters=None):
    """adjust the column-1 orbits of X until Delta_m(D) = 0 mod u^K for m = 1..KZ except the self-dual even m (left to check);
    D = cayley(X) (or D = 1 + X if linear)"""
    skip = [m for m in range(1, KZ + 1) if (2 - m) % 5 == 1 and m % 2 == 0]
    X = dict(X)
    for it in range(iters or (3 * (K + KZ) + 6)):
        D = madd(ident(K), X, p) if linear else cayley(X, KZ, p, K)
        dl = {m: Delta(D, yh, m, KZ, p, K) for m in range(1, KZ + 1) if m not in skip}
        if not any(any(v) for v in dl.values()): return X, D, it, skip
        for m, v in dl.items():                                   # Delta_m = -X[(m+1,1,m)] + ..  (yhat_0 = -1)
            if any(v):
                cur = X.get(((m + 1) % 5, 1, m), [0] * K)
                set_orbit(X, m, 1, [(x + y) % p for x, y in zip(cur, v)], p)
    raise AssertionError("annihilation iteration did not converge")


def checks(D, yh, KZ, p, K):
    """(symplectic defects, flat-irregular entries, {m: Delta_m != 0}, delta_0, T_eff/u)"""
    sd = []
    for n in range(KZ + 1):
        for c in range(5):
            cp = (3 - c - n) % 5; acc = [0] * K
            for j in range(n + 1):
                s = D.get(((c + j) % 5, c, j)); t = D.get(((cp + n - j) % 5, cp, n - j))
                if s and t:
                    pr = smul(s, t, p, K); sg = 1 if j % 2 == 0 else p - 1
                    acc = [(x + sg * y) % p for x, y in zip(acc, pr)]
            if n == 0: acc[0] = (acc[0] - 1) % p
            if any(acc): sd.append((n, c))
    irr = [key for key, v in D.items() if any(v[:min(mu(key[1], key[2]), K)])]
    ann = {m: v for m in range(1, KZ + 1) for v in [Delta(D, yh, m, KZ, p, K)] if any(v)}
    d0 = Delta(D, yh, 0, KZ, p, K); dm1 = Delta(D, yh, -1, KZ, p, K)
    return sd, irr, ann, d0, [(-x) % p for x in dm1]
