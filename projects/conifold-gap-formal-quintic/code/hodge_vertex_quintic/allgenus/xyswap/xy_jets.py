# -*- coding: utf-8 -*-
"""xy_jets.py -- the closed formula (7.9) of arXiv:2212.00320 for n >= 1 points, evaluated exactly over Q at rational points
z_i = b_i (x = z^4 + v/z, v a rational number) by Taylor jets in delta_i = z_i - b_i.  Returns omega^{ABDKS}_{g,n}/(dz_1..dz_n);
multiply by (-1)^n for the EO convention.  Graphs: connected simple graphs on n numbered vertices, weight hbar^(2 g(Gamma)).
"""
import itertools
from fractions import Fraction as Fr
from math import comb
from xy_closed import exponent_terms


def vertex_table(g, v):
    """-(1/w) exp(E) as {(ew, eh): {ez: c}} with eh <= 2g (v numeric)"""
    terms = [(tw, th, tz, tc * Fr(v) ** tv) for (tw, th, tz, tv, tc) in exponent_terms(g)]
    cur = {(0, 0): {0: Fr(1)}}; tot = {}; n = 0
    while cur:
        for (ew, eh), pol in cur.items():
            d = tot.setdefault((ew - 1, eh), {})
            for ez, c in pol.items(): d[ez] = d.get(ez, 0) - c
        n += 1; nxt = {}
        for (ew, eh), pol in cur.items():
            for (tw, th, tz, tc) in terms:
                if eh + th > 2 * g: continue
                d = nxt.setdefault((ew + tw, eh + th), {})
                for ez, c in pol.items(): d[ez + tz] = d.get(ez + tz, 0) + c * tc / n
        cur = {k: p for k, p in nxt.items() if any(p.values())}
    return tot


def binom_gen(m, j):
    """binomial(m, j) for any integer m, j >= 0"""
    r = Fr(1)
    for i in range(j): r = r * (m - i) / (i + 1)
    return r


def taylor_pow(b, m, R):
    """(b + d)^m, m integer, as list up to d^R"""
    return [binom_gen(m, j) * Fr(b) ** (m - j) for j in range(R + 1)]


def taylor_laurent(pol, b, R):
    out = [Fr(0)] * (R + 1)
    for ez, c in pol.items():
        for j, x in enumerate(taylor_pow(b, ez, R)): out[j] += c * x
    return out


def taylor_invxp(b, v, R):
    """1/x'(z) = z^2/(4 z^5 - v) at z = b + d"""
    num = taylor_pow(b, 2, R)
    den = [4 * x for x in taylor_pow(b, 5, R)]; den[0] -= v
    inv = [Fr(0)] * (R + 1); inv[0] = 1 / den[0]
    for k in range(1, R + 1):
        inv[k] = -sum(den[i] * inv[k - i] for i in range(1, k + 1)) / den[0]
    return [sum(num[i] * inv[k - i] for i in range(k + 1)) for k in range(R + 1)]


def jmul(a, b, R):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            if all(x <= r for x, r in zip(e, R)): out[e] = out.get(e, 0) + ca * cb
    return out


def uni(i, n, ser):
    return {tuple(j if q == i else 0 for q in range(n)): c for j, c in enumerate(ser) if c}


def edge_jet(i, j, n, m, b, R):
    """(z_i - z_j)^(-m) at (b_i + d_i) - (b_j + d_j)"""
    d0 = Fr(b[i]) - Fr(b[j]); out = {}
    for l in range(R[i] + R[j] + 1):
        c = binom_gen(-m, l) * d0 ** (-m - l)
        for s in range(l + 1):                       # (d_i - d_j)^l
            if s > R[i] or l - s > R[j]: continue
            e = [0] * n; e[i] = s; e[j] = l - s
            out[tuple(e)] = out.get(tuple(e), 0) + c * comb(l, s) * (-1) ** (l - s)
    return out


def apply_D(jet, i, r, u, R):
    """apply D_{x_i}^r: f -> d/dz_i (f / x'(z_i)) r times"""
    n = len(R); cur = jet; Rc = list(R)
    for _ in range(r):
        prod = jmul(cur, uni(i, n, u[: Rc[i] + 1]), Rc)
        Rc[i] -= 1; nxt = {}
        for e, c in prod.items():
            if e[i] >= 1 and e[i] - 1 <= Rc[i]:
                e2 = list(e); e2[i] -= 1; nxt[tuple(e2)] = nxt.get(tuple(e2), 0) + c * e[i]
        cur = nxt
    return cur, Rc


def graphs(n):
    pairs = list(itertools.combinations(range(n), 2)); out = []
    for k in range(len(pairs) + 1):
        for E in itertools.combinations(pairs, k):
            adj = {i: set() for i in range(n)}
            for (a, c) in E: adj[a].add(c); adj[c].add(a)
            seen = {0}; st = [0]
            while st:
                x = st.pop()
                for y in adj[x]:
                    if y not in seen: seen.add(y); st.append(y)
            if len(seen) == n: out.append(E)
    return out


def closed_value(g, b, v):
    """omega^{ABDKS}_{g,n}/(dz_1..dz_n) at z = b (tuple of rationals)"""
    n = len(b); VT = vertex_table(g, v); total = Fr(0)
    for E in graphs(n):
        gG = len(E) - n + 1
        budget = 2 * g - 2 * gG
        if budget < 0: continue
        # edge choices: k_e >= 0, a_e in 0..2k_e : w_i^(a+1) w_j^(2k-a+1) hbar^(2k) (1/4)^k binom(2k,a) / (z_i - z_j)^(2k+2)
        ranges = [[(k, a) for k in range(budget // 2 + 1) for a in range(2 * k + 1)] for _ in E]
        for ch in itertools.product(*ranges):
            hb = sum(2 * k for (k, a) in ch)
            if hb > budget: continue
            wpow = [0] * n; cf = Fr(1)
            for (e, (k, a)) in zip(E, ch):
                wpow[e[0]] += a + 1; wpow[e[1]] += 2 * k - a + 1; cf *= Fr(comb(2 * k, a), 4 ** k)
            vb = budget - hb
            vkeys = [[(ew, eh) for (ew, eh) in VT if eh <= vb and ew + wpow[i] >= 0] for i in range(n)]
            for vc in itertools.product(*vkeys):
                if sum(eh for (_, eh) in vc) != vb: continue
                R = [ew + wpow[i] for i, (ew, eh) in enumerate(vc)]
                jet = {tuple([0] * n): cf}
                for i, key in enumerate(vc): jet = jmul(jet, uni(i, n, taylor_laurent(VT[key], b[i], R[i])), R)
                for (e, (k, a)) in zip(E, ch): jet = jmul(jet, edge_jet(e[0], e[1], n, 2 * k + 2, b, R), R)
                Rc = list(R)
                for i in range(n):
                    jet, Rc = apply_D(jet, i, R[i], taylor_invxp(b[i], v, R[i]), Rc)
                total += jet.get(tuple([0] * n), 0)
    return total
