# -*- coding: utf-8 -*-
"""ref3b_pi.py -- the neighbourhood of Pi: fibre of the high system over the low part of points of Pi, kernel of the Jacobian, and whether
the tangent directions beyond Pi + fibre integrate (second-order test).   Run: py -I ref3b_pi.py --prime 61|89"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *


def nullspace(M, p):
    """basis of the kernel of M (rows) over F_p"""
    M = [[x % p for x in r] for r in M]; n = len(M[0]); piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    out = []
    for f in [c for c in range(n) if c not in piv]:
        x = [0] * n; x[f] = 1
        for k, c in enumerate(piv): x[c] = (-M[k][f]) % p
        out.append(x)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); rng = random.Random(3); I = IRing(p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    print("=" * 100); print("ref3b_pi, p = 2^%s-1" % a.prime)
    to12 = lambda v: [v[0], v[1], v[2], v[3], v[7]] + list(v[8:])
    A12, B12 = to12(A), to12(B)
    conds = [k for k in CONDS if P[k]]
    def jac(pt): return [[mevalR(mdiff(P[k], i, p), pt, I) for i in range(NV)] for k in conds]
    for label, (s, t) in (("random point of Pi", (rng.randrange(p), rng.randrange(p))), ("random point of Pi", (rng.randrange(p), rng.randrange(p))), ("A", (1, 0)), ("B", (0, 1))):
        pt = [(s * A12[j] + t * B12[j]) % p for j in range(NV)]
        rows, rhs, names = high_system(P, pt[:5], I, p); sol, rk, free = solve_field(I, rows, rhs)
        K = nullspace(jac(pt), p)
        # split the kernel: span of Pi (A12, B12), pure-high directions (first 5 entries zero), the rest
        hi = nullspace([r[:] for r in jac(pt)] + [[1 if j == i else 0 for j in range(NV)] for i in range(5)], p)
        def rank(M): return rank_field(I, M) if M else 0
        rPi = rank([A12, B12] + hi); rAll = rank(K + [A12, B12] + hi)
        print("   %s: high system rank %d (free high directions %d); Jacobian kernel dimension %d; pure-high kernel directions %d; dim(Pi + pure-high) = %d; kernel contains it: %s" % (label, rk, len(free), len(K), len(hi), rPi, rAll == len(K)))
        # second-order test for a kernel vector k: is there a curve pt + e k + e^2 k2 + ... in the solution set?  need for all conditions f:
        #   (1/2) D^2 f[k, k] + Df[k2] = 0 solvable for k2, for k in the kernel (a quadratic condition on the kernel coordinates)
        if len(K) > rPi:
            J = jac(pt); nk = len(K)
            def second(kv):
                # exact second-order coefficient of f(pt + e kv) : coefficient of e^2 by three evaluations e = 1, -1, 2 would need degree <= 3; do it exactly
                out = []
                for k_ in conds:
                    f0 = mevalR(P[k_], pt, I); v1 = mevalR(P[k_], [(x + y) % p for x, y in zip(pt, kv)], I); vm = mevalR(P[k_], [(x - y) % p for x, y in zip(pt, kv)], I)
                    out.append((v1 + vm - 2 * f0) * pow(2, p - 2, p) % p)                      # e^2 coefficient (odd terms cancel; degree <= 3)
                return out
            # obstruction space: left kernel of J (functionals l with l J = 0); l . second(kv) must vanish
            Lk = nullspace([list(c) for c in zip(*J)], p)
            # quadratic forms on the kernel: q_l(c) = l . second(sum c_i K_i); test on random kernel vectors outside Pi + pure-high
            tests = []
            for _ in range(6):
                c = [rng.randrange(p) for _ in range(nk)]; kv = [sum(c[i] * K[i][j] for i in range(nk)) % p for j in range(NV)]
                sec = second(kv); tests.append(any(sum(l[i] * sec[i] for i in range(len(sec))) % p for l in Lk))
            print("      second-order obstruction for random kernel vectors (True = obstructed, the direction does not integrate): %s  [left kernel dimension %d]" % (tests, len(Lk)))
            # which kernel vectors are unobstructed?  the quadratic forms q_l restricted to the kernel: compute their common zero cone dimension roughly
            # by polarisation: Gram matrices of each q_l on the kernel
            def qmat(l):
                M = [[0] * nk for _ in range(nk)]; sec1 = [second(K[i]) for i in range(nk)]
                for i in range(nk):
                    M[i][i] = sum(l[c_] * sec1[i][c_] for c_ in range(len(conds))) % p
                    for j in range(i + 1, nk):
                        s2 = second([(x + y) % p for x, y in zip(K[i], K[j])])
                        M[i][j] = M[j][i] = (sum(l[c_] * (s2[c_] - sec1[i][c_] - sec1[j][c_]) for c_ in range(len(conds))) * pow(2, p - 2, p)) % p
                return M
            Qs = [qmat(l) for l in Lk]; span = rank([[x for r in Q for x in r] for Q in Qs])
            # radical common to all forms: kernel vectors c with Q c = 0 for all Q
            rad = nullspace([r for Q in Qs for r in Q], p)
            print("      obstruction quadrics on the %d-dimensional kernel: %d forms spanning %d; ranks %s; common radical dimension %d" % (nk, len(Qs), span, sorted(set(rank(Q) for Q in Qs)), len(rad)))
    print("REF3B done")


if __name__ == "__main__":
    main()
