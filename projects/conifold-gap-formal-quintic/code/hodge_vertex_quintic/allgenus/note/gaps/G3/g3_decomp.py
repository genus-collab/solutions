# -*- coding: utf-8 -*-
"""g3_decomp.py -- exact checks over Q[[u]] (mod u^N) of the u-adic statements of gap cluster G3 (lem:tra-frames step 4,
lem:tra-conv (a), (b)), with v(u) = u tau^(-5/4) and the control v(u) = 1 + u (v(0) != 0).

T1 Theta = 1 + N is invertible mod u^N (Neumann series), both compositions = identity.
T2 pr(q) := q - sigma(q)(z) H, sigma = Theta^-1 (Pi_i(q))_i, is annihilated by every Pi_i (random q, z-degree K).
T3 pr(q) = Sum_(k <= K, a != 1) [z^k] q^a frame(k,a)(v(u)) coefficient by coefficient (frame coordinates).
T4 ord_u of the H-coefficient f_j of frame(k,a) at v(u) equals (4(k-j) + nu_a + 1)/5 exactly (k <= KF).
T5 the infinite element Sum_(k >= 0) z^k H^a: the coefficients [u^m z^j H] of pr(truncations) stabilise at the predicted K;
   CONTROL v(u) = 1 + u: the exact partial sums of [u^0 z^0 H] never stabilise (each new frame changes them).
T6 a toy homogeneous bilinear form mu with the dimension bound (h = 2, n = 2, slot degrees <= 5) and random coefficients
   (singular exponents included): mu(X^(K), X^(K)) for frame truncations X^(K) of an element of K with all frame
   coordinates 1 converges u-adically, with ord(mu(X^(K+5)) - mu(X^(K))) >= 4K/5 - 3; CONTROL v(u) = 1 + u: it does not.
Usage: py -I -B g3_decomp.py N K   (log to stdout; no file written)
"""
import os
import sys
import random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g3_series as S
from g3_frames import pi0, frame_theta, NU

AS = (0, 2, 3, 4)


def Pi(i, q, VP, N):
    """Pi_i(q) mod u^N for q = {(k, b): series}; Pi_i(z^k H^b) = Pi_0(z^(k-i) H^b), Pi_0(z^m H^b) = pi0(m,b) v^(m - m//5)."""
    out = S.zero(N)
    for (k, b), c in q.items():
        m = k - i
        if m < 0:
            continue
        p = pi0(m, b)
        if p:
            out = S.add(out, S.scal(p, S.mul(c, VP.trunc(m - m // 5, N))))
    return out


def theta(s, VP, N):
    L = len(s)
    out = []
    for i in range(L):
        t = s[i][:]
        for q in range(1, (L - 1 - i) // 5 + 1):
            t = S.add(t, S.scal(pi0(5 * q, 1), S.mul(VP.trunc(4 * q, N), s[i + 5 * q])))
        out.append(t)
    return out


def theta_inv(t, VP, N, maxit=10 ** 4):
    L = len(t)
    sig, term, r = [S.zero(N) for _ in range(L)], [x[:] for x in t], 0
    while any(any(x) for x in term):
        sig = [S.add(a, b) for a, b in zip(sig, term)]
        nt = theta(term, VP, N)
        term = [S.scal(-1, S.add(a, S.scal(-1, b))) for a, b in zip(nt, term)]   # -(Theta - 1) term
        r += 1
        if r > maxit:
            raise RuntimeError("Neumann series does not terminate")
    return sig, r


def rand_series(N, rng, lo=-3, hi=3):
    return [Fr(rng.randint(lo, hi)) for _ in range(N)]


def frame_at(k, a, VP, N):
    """frame(k,a) at v(u): {(j, b): series}."""
    out = {(k, a): S.one(N)}
    for j, c in frame_theta(k, a).items():
        e = (4 * (k - j) + NU[a] + 1) // 5
        out[(j, 1)] = S.scal(c, VP.trunc(e, N))
    return out


def pr(q, K, VP, N):
    t = [Pi(i, q, VP, N) for i in range(K + 1)]
    sig, r = theta_inv(t, VP, N)
    out = {kb: c[:] for kb, c in q.items()}
    for j, c in enumerate(sig):
        out[(j, 1)] = S.add(out.get((j, 1), S.zero(N)), S.scal(-1, c))
    return out, r


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    rng = random.Random(20261115)
    VP = S.VPowers(N + 2)
    print("g3_decomp.py  N = %d (mod u^N), K = %d (z-degree); v(u) = u tau^(-5/4); exact over Q" % (N, K))
    ok = True
    # T1
    s = [rand_series(N, rng) for _ in range(K + 1)]
    a1, _ = theta_inv(theta(s, VP, N), VP, N)
    a2 = theta(theta_inv(s, VP, N)[0], VP, N)
    t1 = (a1 == s) and (a2 == s)
    print("[T1] Theta^-1 Theta = Theta Theta^-1 = 1 on a random vector of length %d: %s" % (K + 1, t1))
    ok &= t1
    # T2, T3
    q = {(k, b): rand_series(N, rng) for k in range(K + 1) for b in range(5)}
    p, r = pr(q, K, VP, N)
    t2 = all(not any(Pi(i, p, VP, N)) for i in range(K + 4))
    print("[T2] Pi_i(pr q) = 0 mod u^%d for 0 <= i <= %d (Neumann terms used: %d): %s" % (N, K + 3, r, t2))
    rec = {}
    for k in range(K + 1):
        for a in AS:
            for kb, c in frame_at(k, a, VP, N).items():
                rec[kb] = S.add(rec.get(kb, S.zero(N)), S.mul(q[(k, a)], c))
    keys = set(rec) | set(p)
    t3 = all(rec.get(kb, S.zero(N)) == p.get(kb, S.zero(N)) for kb in keys)
    print("[T3] pr(q) = Sum [z^k]q^a frame(k,a)(v(u)) on all %d coefficients: %s" % (len(keys), t3))
    ok &= t2 and t3
    # T4
    KF, nchk, nbad = 3 * K, 0, 0
    VPb = S.VPowers(4 * KF // 5 + 4)
    for k in range(KF + 1):
        for a in AS:
            for j, c in frame_theta(k, a).items():
                e = (4 * (k - j) + NU[a] + 1) // 5
                ser = S.scal(c, VPb.trunc(e, VPb.N))
                nchk += 1
                if S.ordu(ser) != e:
                    nbad += 1
    print("[T4] ord_u f_j(v(u)) = (4(k-j)+nu_a+1)/5 for all %d non-zero frame coefficients, k <= %d: failures %d" % (nchk, KF, nbad))
    ok &= nbad == 0
    # T5
    ok &= t5(N, VP)
    # T6
    ok &= t6(N, rng)
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


def t5(N, VP):
    a = 0
    good = True
    for (m, j) in ((0, 0), (N - 1, 0), (N - 1, 3)):
        hist = []
        for KK in range(j, j + 5 * N // 4 + 12):
            qq = {(k, a): S.one(N) for k in range(KK + 1)}
            p, _ = pr(qq, KK, VP, N)
            hist.append(p.get((j, 1), S.zero(N))[m])
        last = max([i for i in range(1, len(hist)) if hist[i] != hist[i - 1]] or [0]) + j
        pred = j + (5 * m - 1) // 4             # largest k with (4(k-j)+1)/5 <= m: later frames vanish mod u^(m+1)
        good &= last <= max(pred, j)
        print("[T5] Sum_k z^k 1: [u^%d z^%d H] of pr(truncation at K) last changes at K = %d (predicted bound %d); value %s"
              % (m, j, last, pred, hist[-1]))
    # control v(u) = 1 + u: frames at v(0) = 1; coefficient [u^0 z^0 H] of Sum_(k <= K) frame(k,0)
    VPc = S.VPowers(N, kind="shift")
    ps, changes = S.zero(N)[0], 0
    for KK in range(0, 61):
        fr = frame_at(KK, a, VPc, N)
        new = ps + fr.get((0, 1), S.zero(N))[0]
        if new != ps:
            changes += 1
        ps = new
    print("[T5-C] CONTROL v = 1 + u: the exact partial sums of [u^0 z^0 H] of Sum_(k<=K) frame(k,0) change %d times for K <= 60 "
          "(at every K = 1 mod 5): no u-adic limit (expected)" % changes)
    good &= changes >= 12
    return good


def t6(N, rng):
    h, n, D = 2, 2, 5
    C = {}
    for k1 in range(D + 1):
        for k2 in range(D + 1 - k1):
            for a1 in range(5):
                for a2 in range(5):
                    num = 4 * (k1 + k2) + NU[a1] + NU[a2] - 8
                    if num % 5 == 0:
                        C[(k1, a1, k2, a2)] = (Fr(rng.randint(1, 9)) * rng.choice((1, -1)), 2 - 2 * h + num // 5)
    good = True
    for kind in ("model", "shift"):
        M = 2 * N
        VP = S.VPowers(M + 8, kind=kind)

        def flat(KK):
            comp = {}
            for k in range(KK + 1):
                for a in AS:
                    for (j, b), c in frame_at(k, a, VP, M + 8).items():
                        if j <= D:
                            comp[(j, b)] = S.add(comp.get((j, b), S.zero(M + 8)), c)
            comp[(0, 1)] = S.add(comp.get((0, 1), S.zero(M + 8)), S.one(M + 8))     # + 1 * H
            return comp

        def mu(comp):
            tot = (-8, S.zero(M + 8))
            for (k1, a1, k2, a2), (c, e) in C.items():
                x, y = comp.get((k1, a1)), comp.get((k2, a2))
                if x is None or y is None:
                    continue
                term = S.lmul((0, S.scal(c, S.mul(x, y))), VP.get(e))
                tot = S.ladd(tot, term, M)
            return tot
        vals = [mu(flat(KK)) for KK in range(5, 41, 5)]
        ords = []
        for A, B in zip(vals, vals[1:]):
            d = S.ladd(B, (A[0], [-x for x in A[1]]), M)
            o = S.ordu(d[1])
            ords.append(None if o is None else o + d[0])
        if kind == "model":
            ok = all(o is None or o >= 4 * K0 / 5 - 3 for o, K0 in zip(ords, range(5, 41, 5)))
            print("[T6] toy form, v = u tau^(-5/4): ord_u(mu(X^(K+5)) - mu(X^(K))) for K = 5..35: %s (bound 4K/5 - 3; None = 0 "
                  "mod u^%d): %s" % (ords, M, ok))
            good &= ok
        else:
            bad = sum(1 for o in ords if o is not None and o <= 0)
            print("[T6-C] CONTROL v = 1 + u: the same differences have ord_u %s: no convergence (expected)" % ords)
            good &= bad == len(ords)          # every difference is non-zero already at u^0
    return good


if __name__ == "__main__":
    main()
