# -*- coding: utf-8 -*-
"""chk3.py -- referee NOTE: consistency of s.6c (O(u^2) on the hyperplane m_1 = 0) with s.6d on the locus c_1(0) = 0.
For a low part on {m_1 = 0, genus-2 rows} the O(u^2) conditions on the 7 high circulant coordinates (c4, c5, c6) are affine-linear:
rows J3 (2) and the five quadrics.  Report rank, consistency and fibre dimension at generic points and at points with c_1(0) = 0."""
import sys, json, os, random
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq3_circsys import CircSys
from fq3_quad import IDX_2EDGE


def rank(M, p):
    M = [list(r) for r in M]; rk = 0
    for c in range(len(M[0])):
        piv = None
        for r in range(rk, len(M)):
            if M[r][c] % p: piv = r; break
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], p - 2, p); M[rk] = [x * inv % p for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] % p:
                f = M[r][c]; M[r] = [(x - f * y) % p for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def main(tag):
    S = CircSys(tag); p = S.p; rng = random.Random(7)
    print("=" * 30, "prime", tag)
    # q_HH = 0 : the quadrics are affine in the high part for a fixed low part
    hh = [0] * 8 + [rng.randrange(p) for _ in range(7)]
    print("quadrics on a random pure-high vector (q_HH):", [S.quad(i, hh) for i in IDX_2EDGE])
    for label, c10 in (("generic c_1(0)", None), ("c_1(0) = 0", 0)):
        for trial in range(3):
            c1 = [rng.randrange(p) if c10 is None else 0, rng.randrange(p), rng.randrange(p)]
            low = S.low(c1, rng.randrange(p), rng.randrange(p)); v = low + [0] * 7
            rows = []; rhs = []
            for r in range(2):
                rows.append([S.J3[r][j] for j in range(8, 15)]); rhs.append((-sum(S.J3[r][j] * v[j] for j in range(8))) % p)
            for i in IDX_2EDGE:
                g = S.quad_grad(i, v); rows.append([g[j] for j in range(8, 15)]); rhs.append((-S.quad(i, v)) % p)
            rk = rank(rows, p); rka = rank([r + [b] for r, b in zip(rows, rhs)], p)
            zero_rows = [nm for nm, r, b in zip(['J3a', 'J3b'] + ['q%d' % i for i in IDX_2EDGE], rows, rhs) if not any(r) and not b]
            # m_2 : q_17 = m_1 m_2, so on m_1 = 0 grad q_17 = m_2(Delta) grad m_1 ; m_2(Delta) = 0 iff the full gradient vanishes
            g17 = S.quad_grad(17, v)
            print("  %-15s trial %d: 7 equations on 7 high coordinates: rank %d, augmented rank %d (%s), fibre dimension %s; identically satisfied rows %s; grad q_17 = 0 at (low, 0): %s"
                  % (label, trial, rk, rka, "consistent" if rk == rka else "INCONSISTENT", (7 - rk) if rk == rka else "-", zero_rows, not any(g17)))


if __name__ == "__main__":
    for tag in ("61", "89"): main(tag)
