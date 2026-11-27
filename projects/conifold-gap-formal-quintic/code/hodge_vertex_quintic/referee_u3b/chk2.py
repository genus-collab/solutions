# -*- coding: utf-8 -*-
"""chk2.py -- referee NOTE: the N_2 system (2 genus-3 rows + 5 polar forms) at the solution points: rank, left kernel, and the
consistency condition lambda . rhs = 0 read off the stored P1 passes (g3_u3_<p>_<cls>.json, key A1), when present."""
import sys, json, os, random
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq3_circsys import CircSys
from fq3_quad import IDX_2EDGE

ROWS = ['J3(2)', 'J3(5)'] + ['gradq_%d' % i for i in IDX_2EDGE]


def nullspace_left(M, p):
    """basis of {lam : lam M = 0} by elimination on [M | I]"""
    n = len(M); m = len(M[0]); W = [list(M[r]) + [1 if r == c else 0 for c in range(n)] for r in range(n)]; rk = 0
    for c in range(m):
        piv = None
        for r in range(rk, n):
            if W[r][c] % p: piv = r; break
        if piv is None: continue
        W[rk], W[piv] = W[piv], W[rk]; inv = pow(W[rk][c], p - 2, p); W[rk] = [x * inv % p for x in W[rk]]
        for r in range(n):
            if r != rk and W[r][c] % p:
                f = W[r][c]; W[r] = [(x - f * y) % p for x, y in zip(W[r], W[rk])]
        rk += 1
    return rk, [W[r][m:] for r in range(rk, n)]


def main(tag):
    S = CircSys(tag); p = S.p; cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    print("=" * 30, "prime", tag)
    A = cs['Pi_A']; B = cs['Pi_B']; star = cs['circ_star']; sol = cs['circ_b_solutions'][0]; extra = cs['circ_b_extra']
    AB2 = [(u + 2 * w) % p for u, w in zip(A, B)]
    mat = lambda v: [S.J3[0], S.J3[1]] + [S.quad_grad(i, v) for i in IDX_2EDGE]
    # the pure-high point of the 'circ' pass: c5 = (1, 2, -1), c6 = (3, 1), c4 from the two genus-3 rows
    hi = [0] * 8 + [0, 0, 1, 2, p - 1, 3, 1]
    a, b, c, d = S.J3[0][8], S.J3[0][9], S.J3[1][8], S.J3[1][9]; det = (a * d - b * c) % p
    r0 = (-sum(S.J3[0][j] * hi[j] for j in range(10, 15))) % p; r1 = (-sum(S.J3[1][j] * hi[j] for j in range(10, 15))) % p
    idet = pow(det, p - 2, p); hi[8] = (r0 * d - b * r1) * idet % p; hi[9] = (a * r1 - c * r0) * idet % p
    print("pure-high point: J3 rows involve c4:", det != 0, " residuals:", [k for k, x in S.residuals(hi).items() if x])
    pts = [("pure-high (circ)", hi, 'circ'), ("A+2B (circpi)", AB2, 'circpi'), ("genuine c1_0=0 (circb)", sol, 'circb'), ("Delta* (circstar)", star, 'circstar'),
           ("B of Pi", B, None), ("random point of Pi", [(3 * u + 5 * w) % p for u, w in zip(A, B)], None)] + [("extra[%d]" % n, e, None) for n, e in enumerate(extra)]
    for name, v, cls in pts:
        rk, ker = nullspace_left(mat(v), p)
        line = "%-24s rank %d" % (name, rk)
        for lam in ker:
            line += "; left-kernel support " + str([ROWS[r] for r in range(7) if lam[r]])
        f = os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls)) if cls else None
        if f and os.path.exists(f):
            d = json.load(open(f))
            if 'A1' in d:
                A1 = d['A1']
                rhs = [(-A1['2'][2]) % p, (-A1['5'][2]) % p] + [(-A1[str(i)][3]) % p for i in IDX_2EDGE]
                line += "; P1 stored: rhs non-zero rows " + str([ROWS[r] for r in range(7) if rhs[r]])
                line += "; lambda.rhs = " + str([sum(l * x for l, x in zip(lam, rhs)) % p for lam in ker])
                # the three-edge [u^3] of the P1 pass (cubic conditions by the full hatted theory)
                line += "; three-edge [u^0..u^3] all zero: " + str(all(x == 0 for i in (4, 9, 10, 11, 12, 18, 19, 20, 32) for x in A1[str(i)][:4]))
                line += "; one-edge [u^0,u^1] zero: " + str(all(x == 0 for i in (2, 5) for x in A1[str(i)][:2]))
                line += "; two-edge [u^0..u^2] zero: " + str(all(x == 0 for i in IDX_2EDGE for x in A1[str(i)][:3]))
            else:
                line += "; P1 not stored yet (keys: %s)" % sorted(d.keys())
            if 'n2_rank' in d: line += "; script n2_rank %s" % d['n2_rank']
            if 'failed' in d: line += "; script failed list %s" % d['failed']
        print(line)
    # is the left kernel the same combination of rows along the scaling line s * Delta ?  (rows gradq scale linearly, J3 rows do not)
    for name, v in (("Delta*", star), ("genuine", sol)):
        k1 = nullspace_left(mat(v), p)[1][0]; v2 = [2 * x % p for x in v]; k2 = nullspace_left(mat(v2), p)[1][0]
        print("   %s: left kernel has zero J3 components: %s ; proportional at Delta and 2 Delta: %s" % (
            name, k1[0] == 0 and k1[1] == 0, all((k1[i] * k2[6] - k2[i] * k1[6]) % p == 0 for i in range(7))))


if __name__ == "__main__":
    for tag in ("61", "89"): main(tag)
