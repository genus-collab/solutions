# -*- coding: utf-8 -*-
"""ref_unique_core.py -- referee code (key gaps/G2/referee): Theorem U's uniqueness clause tested on RANDOM one-pole systems
(not the quintic, not the model), mod p = 2^61 - 1, by explicit null-space bases (reduced row echelon form), independent of
GAPS_G2's scripts.  Weights are truncated Laurent series mod p, stored as dicts {exponent: coefficient}, exponents >= -1.
Stdlib only.
"""
import random

PR = (1 << 61) - 1


def inv(x):
    return pow(x % PR, PR - 2, PR)


def s_mul(f, g, N):
    """product of dict-series, truncated at exponent < N"""
    out = {}
    for i, a in f.items():
        for j, b in g.items():
            if i + j < N:
                out[i + j] = (out.get(i + j, 0) + a * b) % PR
    return {k: v for k, v in out.items() if v}


def s_inv_unit(f, N):
    """inverse of a unit power series (f[0] != 0), truncated at exponent < N"""
    a0 = f.get(0, 0); assert a0 % PR
    ia = inv(a0); out = {0: ia}
    for n in range(1, N):
        s = sum(f.get(j, 0) * out.get(n - j, 0) for j in range(1, n + 1)) % PR
        if s:
            out[n] = (-s * ia) % PR
    return out


def s_shift(f, k):
    return {e + k: v for e, v in f.items()}


def s_deriv(f):
    return {e - 1: (e * v) % PR for e, v in f.items() if e != 0 and (e * v) % PR}


def s_pow_unit(f, num, den, N):
    """f^(num/den) for f = 1 + O(u), via the binomial recursion n g_n = Sum_k (k(alpha+1) - n) f_k g_(n-k) (f_0 = 1)"""
    assert f.get(0, 0) == 1
    al = num * inv(den) % PR
    g = {0: 1}
    for n in range(1, N):
        s = 0
        for k in range(1, n + 1):
            fk = f.get(k, 0)
            if fk:
                s += (k * (al + 1) - n) % PR * fk % PR * g.get(n - k, 0)
        s %= PR
        if s:
            g[n] = s * inv(n) % PR
    return g


def rnd_unit(rng, N):
    return {e: (rng.randrange(1, PR) if e == 0 else rng.randrange(PR)) for e in range(N)}


def systems(rng, N, pole=True, equal=True, same=False):
    """weights (b, bt) with (H) [pole: b_1, bt_1 have a simple pole] and (EP) [equal: prod b = prod bt]; same: bt = b"""
    b = [rnd_unit(rng, N) for _ in range(5)]
    if pole:
        b[1] = s_shift(b[1], -1)
    if same:
        return b, [dict(x) for x in b]
    bt = [rnd_unit(rng, N) for _ in range(4)]
    if pole:
        bt[1] = s_shift(bt[1], -1)
    num = {0: 1}
    for x in b:
        num = s_mul(num, s_shift(x, 1) if x.get(-1) else x, N + 1)
    den = {0: 1}
    for x in bt:
        den = s_mul(den, s_shift(x, 1) if x.get(-1) else x, N + 1)
    bt4 = s_mul(num, s_inv_unit(den, N + 1), N)                  # prod b / (bt_0..bt_3): a unit
    if not equal:
        bt4 = s_mul(bt4, {0: 1, 1: 1}, N)                       # bt_4 -> (1 + u) bt_4: products differ
    bt.append(bt4)
    return b, bt


def theta_frakd(b, bt, N):
    """the diagonal frak d Theta_U (entries frak d d^U_a) as dict-series, computed from the weights"""
    d = [{0: 1}]
    for a in range(4):
        ra = s_mul(s_shift(bt[a], 1) if bt[a].get(-1) else bt[a],
                   s_inv_unit(s_shift(b[a], 1) if b[a].get(-1) else b[a], N), N)
        d.append(s_mul(d[a], ra, N))
    det = {0: 1}
    for x in d:
        det = s_mul(det, x, N)
    c0 = det[0]; detn = {e: v * inv(c0) % PR for e, v in det.items()}
    fd = s_pow_unit(detn, -1, 5, N)
    return [s_mul(fd, x, N) for x in d]


def equations(b, bt, K, M, delta):
    """degree-delta part of z P' = A P - P At: unknowns (k, bb, m) = [u^m] of the entry (k + bb - delta, bb) of P_k"""
    rows = []
    for k in range(K + 1):
        for bb in range(5):
            a = (k - 1 + bb - delta) % 5
            for j in range(-1, M):
                r = {}
                if k >= 1 and j + 1 <= M and j + 1 > 0:
                    r[(k - 1, bb, j + 1)] = (j + 1) % PR           # (P_(k-1)')_(a,bb) at u^j
                for e, c in b[a].items():                         # - b_a (P_k)_(a+1, bb)
                    m = j - e
                    if 0 <= m <= M:
                        r[(k, bb, m)] = (r.get((k, bb, m), 0) - c) % PR
                for e, c in bt[(bb - 1) % 5].items():             # + (P_k)_(a, bb-1) bt_(bb-1)
                    m = j - e
                    if 0 <= m <= M:
                        key = (k, (bb - 1) % 5, m); r[key] = (r.get(key, 0) + c) % PR
                r = {c: v for c, v in r.items() if v}
                if r:
                    rows.append(r)
    return rows


def nullspace(rows, cols):
    """explicit basis of the null space (reduced row echelon form mod p); cols = list of all column keys"""
    order = {c: i for i, c in enumerate(cols)}
    piv = {}                                                      # pivot column -> reduced row (dict), leading coefficient 1
    for r0 in rows:
        r = dict(r0)
        for pc, pr_ in piv.items():                               # reduce by existing pivots
            v = r.get(pc)
            if v:
                for c, w in pr_.items():
                    x = (r.get(c, 0) - v * w) % PR
                    if x:
                        r[c] = x
                    else:
                        r.pop(c, None)
        if not r:
            continue
        lead = min(r, key=lambda c: order[c]); il = inv(r[lead])
        r = {c: v * il % PR for c, v in r.items()}
        for pc in list(piv):                                      # keep the pivot rows fully reduced
            v = piv[pc].get(lead)
            if v:
                row = piv[pc]
                for c, w in r.items():
                    x = (row.get(c, 0) - v * w) % PR
                    if x:
                        row[c] = x
                    else:
                        row.pop(c, None)
        piv[lead] = r
    free = [c for c in cols if c not in piv]
    basis = []
    for f in free:
        vec = {f: 1}
        for pc, row in piv.items():
            v = row.get(f)
            if v:
                vec[pc] = (-v) % PR
        basis.append(vec)
    return basis


def rank_of(vectors, coords):
    """rank mod p of the projections of the given sparse vectors onto the coordinate list coords"""
    idx = {c: i for i, c in enumerate(coords)}
    rows = []
    for v in vectors:
        r = {idx[c]: x for c, x in v.items() if c in idx and x}
        if r:
            rows.append(r)
    piv = {}
    for r in rows:
        r = dict(r)
        while r:
            lead = min(r)
            if lead not in piv:
                il = inv(r[lead]); piv[lead] = {c: x * il % PR for c, x in r.items()}
                break
            pv = piv[lead]; f = r[lead]
            for c, x in pv.items():
                y = (r.get(c, 0) - f * x) % PR
                if y:
                    r[c] = y
                else:
                    r.pop(c, None)
    return len(piv)
