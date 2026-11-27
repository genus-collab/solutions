# -*- coding: utf-8 -*-
"""tr_v0.py -- behaviour of omega_{g,n} of x = W_u(p) = p^4/4 + u/p, y = p as u -> 0 at fixed p (TR_MODEL.md s.5).

Homogeneity (PROVED, TR_MODEL.md s.5): with s^5 = u,  omega^(u)_{g,n}(p) = s^(5(2-2g-n)) f(p/s) prod dp_l/s  where omega^(1) = f(p) prod dp.
Expanding omega^(1) at p_l = infinity, omega^(1) = Sum_m c_m prod dp_l/p_l^(m_l+1) (|p_l| > 1; m_l >= 1 since omega has no residues
and only poles at the branch points), gives for |p_l| > |u|^(1/5) the convergent expansion
      omega^(u)_{g,n}(p) = Sum_{K >= 0} u^(2-2g-n+K) Omega_K(p),   Omega_K = Sum_{|m| = 5K} c_m prod dp_l/p_l^(m_l+1),
(|m| = 0 mod 5 by the Z/5 symmetry, asserted).  omega is singular at u = 0 iff Omega_K != 0 for some K < 2g-2+n.
e_{j,k}(p) = dp/(p - c_j)^(k+1) = Sum_{m >= k} binom(m, k) c_j^(m-k) dp/p^(m+1).
"""
import itertools, os, sys
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)


def at_infinity(F, c, W, n, mtot):
    """{m-tuple: c_m} for |m| <= mtot from the e-basis tensor W of an n-point differential (branch points c);
    slot by slot: e_{j,k} -> Sum_{m >= k} binom(m, k) c_j^(m-k) [m], keys pruned by the bound |m| <= mtot"""
    cp = {}
    def cpow(j, e):
        key = (j, e)
        if key not in cp:
            r = F.one
            for _ in range(e): r = r * c[j]
            cp[key] = r
        return cp[key]
    cur = dict(W)
    for slot in range(n):
        new = {}
        for key, v in cur.items():
            j, k = key[slot]
            fixed = sum((e if isinstance(e, int) else e[1]) for l, e in enumerate(key) if l != slot)
            for m in range(k, mtot - fixed + 1):
                w = v * cpow(j, m - k) * comb(m, k)
                nk = key[:slot] + (m,) + key[slot + 1:]
                new[nk] = new[nk] + w if nk in new else w
        cur = {kk: x for kk, x in new.items() if not x.iszero()}
    return cur


def v0_profile(F, c, W, g, n):
    """returns (dict K -> {m: c_m}, Kmin, bad) for |m| <= 5 (2g-2+n); bad = number of non-zero c_m with |m| != 0 mod 5"""
    mt = 5 * (2 * g - 2 + n)
    E = at_infinity(F, c, W, n, mt)
    byK = {}; bad = 0
    for m, x in E.items():
        if sum(m) % 5: bad += 1; continue
        byK.setdefault(sum(m) // 5, {})[m] = x
    Kmin = min(byK) if byK else None
    return byK, Kmin, bad
