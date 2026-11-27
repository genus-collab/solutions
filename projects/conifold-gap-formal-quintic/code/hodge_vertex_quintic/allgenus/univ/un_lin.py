# -*- coding: utf-8 -*-
"""un_lin.py -- linear algebra mod p and the bookkeeping of the first-order analysis (E3).

Variables.  A flat-regular, graded, infinitesimally symplectic left dressing is p(u, z) = Sum_{orbits (k, b)} p_{k,b}(u) H(k, b),
H(k, b) = z^k [E_{b+k,b} - (-1)^k E_{b*+k,b*}];  in hatted variables H(k, b) = u^mu(b,k) h(k, b) (mu(b*, k) = mu(b, k)).  Its visible part
at genus g (window u^0..u^W, W = 2g - 3) are the Taylor coefficients p_{k,b,j} with mu(b, k) + j <= W: the variable list VARS(g).
  entry(k, col) of p: the coefficient of z^k E_{col+k,col}:  + p_{k,b} if col = b,  -(-1)^k p_{k,b} if col = b* != b.

Response equations (un_e3_resp):  [u^n] Sum p_{k,b}(u) u^mu rho_{k,b}(Gamma; u) = 0,  n <= W, Gamma an outer graph with edges.
Annihilation equations (hypothesis (A)):  (p Y^C_Q)_+ = 0 with Y^C_Q = Sum_{n>=0} z^-n y_n(u) e_{(1-n) mod 5}, y_0 = -1,
  y_n' = b_{(1-n) mod 5} y_{n-1}, y_n(0) = 0 (n >= 1):   C_m := Sum_{n>=0} entry(m + n, (1 - n) mod 5)(u) y_n(u) = 0  (m >= 1),
  visible through u^(W - mu(1, m)) (the windows match: mu((1-n) mod 5, m + n) = mu(1, m) + ord_u y_n).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from un_core import mu


def rref(rows, p, ncols):
    """reduced row echelon form mod p; returns (pivot columns, reduced rows as dicts col -> value)"""
    M = [[x % p for x in r] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return piv, M[:r]


def rank(rows, p, ncols): return len(rref(rows, p, ncols)[0]) if rows else 0


def nullspace(rows, p, ncols):
    """basis of {x : rows x = 0} as lists of length ncols"""
    piv, M = rref(rows, p, ncols) if rows else ([], [])
    free = [c for c in range(ncols) if c not in piv]; out = []
    for f in free:
        v = [0] * ncols; v[f] = 1
        for i, c in enumerate(piv): v[c] = (-M[i][f]) % p
        out.append(v)
    return out


def in_span(vec, basis, p, ncols):
    return rank(basis + [vec], p, ncols) == rank(basis, p, ncols)


def orbit_of(k, col):
    """(representative b, sign) of the entry (k, col), or None if the entry is forced to vanish (self-dual, k even)"""
    bs = (3 - col - k) % 5
    if bs == col: return (col, 1) if k % 2 == 1 else None
    b = min(col, bs)
    # representative = the first of the pair met when scanning b = 0..4 (as in un_e3_resp.orbits)
    return (b, 1) if col == b else (b, -1 if k % 2 == 0 else 1)


def variables(g, kmax=None, with_k0=True):
    """[(k, b, j)] visible at genus g"""
    W = 2 * g - 3; kmax = 3 * g - 3 if kmax is None else kmax; out = []
    for k in range(0 if with_k0 else 1, kmax + 1):
        seen = set()
        for b in range(5):
            if b in seen: continue
            bs = (3 - b - k) % 5; seen.add(b); seen.add(bs)
            if bs == b and k % 2 == 0: continue
            for j in range(0, W - mu(b, k) + 1): out.append((k, b, j))
    return out


def response_rows(resp_g, outer_g, g, VARS, p, graphs=None):
    """rows of the response map: one row per (outer graph with edges, n <= W) that is not identically zero; returns (rows, labels)"""
    W = 2 * g - 3; rows = []; labels = []
    idxs = sorted(set(int(i) for d in resp_g.values() for i in d))
    for i in idxs:
        if i == 0: continue
        if graphs is not None and i not in graphs: continue
        for n in range(W + 1):
            row = []
            for (k, b, j) in VARS:
                c = resp_g.get("%d,%d" % (k, b), {}).get(str(i)); d = n - j - mu(b, k)
                row.append(c[d] % p if (c is not None and 0 <= d < len(c)) else 0)
            if any(row): rows.append(row); labels.append((i, n))
    return rows, labels


def vertex_rows(resp_g, g, VARS, p):
    """the response of the vertex graph: rows n = 0..W"""
    W = 2 * g - 3; rows = []
    for n in range(W + 1):
        row = []
        for (k, b, j) in VARS:
            c = resp_g.get("%d,%d" % (k, b), {}).get("0"); d = n - j - mu(b, k)
            row.append(c[d] % p if (c is not None and 0 <= d < len(c)) else 0)
        rows.append(row)
    return rows


def yc_series(bq, nmax, NU):
    """y_n(u) (n = 0..nmax) of the pure-1/z solution of th Y_a = b_a Y_{a+1}: lists of Fractions [u^0..u^(NU-1)]"""
    ys = [[Fr(-1)] + [Fr(0)] * (NU - 1)]
    for n in range(1, nmax + 1):
        a = (1 - n) % 5; b = bq[a]; prev = ys[-1]; out = [Fr(0)] * NU
        # integrand = b * y_{n-1}; b may have a simple pole (a = 1)
        vb = b.v if b.c else 0
        for i in range(-1, NU - 1):                       # coefficient of u^i of the integrand
            s = Fr(0)
            for l in range(NU):
                if prev[l] == 0: continue
                e = i - l
                if e < vb or e >= b.prec: continue
                s += prev[l] * b.coef(e)
            if i == -1: assert s == 0, "logarithm in y_%d" % n
            else: out[i + 1] = s / (i + 1)
        ys.append(out)
    return ys


def annihilation_rows(g, VARS, ys, R, kmax=None):
    """rows of the conditions C_m [u^i] = 0, m = 1..kmax, i = 0..W - mu(1, m)"""
    W = 2 * g - 3; kmax = 3 * g - 3 if kmax is None else kmax; pos = {v: n for n, v in enumerate(VARS)}; rows = []; labels = []
    for m in range(1, kmax + 1):
        for i in range(0, W - mu(1, m) + 1):
            row = [0] * len(VARS)
            for n in range(0, len(ys)):
                k = m + n; col = (1 - n) % 5
                if k > kmax: break
                ob = orbit_of(k, col)
                if ob is None: continue
                b, sg = ob
                for j in range(0, i + 1):
                    key = (k, b, j)
                    if key not in pos: continue
                    c = ys[n][i - j] if i - j < len(ys[n]) else None
                    assert c is not None
                    if c: row[pos[key]] = (row[pos[key]] + sg * R.of(c)) % R.p
            if any(row): rows.append(row); labels.append((m, i))
    return rows, labels


def teff_row(g, VARS, ys, R, n_order):
    """coefficient row of [u^n_order] of  delta T := -[z^-1 e_0]((p Y^C_Q)) = -Sum_{k>=0} entry(k, -k mod 5)(u) y_{k+1}(u)"""
    pos = {v: n for n, v in enumerate(VARS)}; row = [0] * len(VARS)
    for k in range(0, len(ys) - 1):
        col = (-k) % 5; ob = orbit_of(k, col)
        if ob is None: continue
        b, sg = ob
        for j in range(0, n_order + 1):
            key = (k, b, j)
            if key not in pos: continue
            c = ys[k + 1][n_order - j]
            if c: row[pos[key]] = (row[pos[key]] - sg * R.of(c)) % R.p
    return row
