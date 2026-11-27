"""Stage-2 critic: the leg dictionary LD, D(ii) once more, and an end-to-end check on the frames of K_0.
Conventions and the pairing: see cr2_legcore.py (imported read-only from this folder).

CLAIM LD: q_{Phi_k^b} := Sum_{a,d} <Phi_k^b, d xif_{a,d}> H^a z^d = -c_b tau_k(H^b), c = (-1, 1, -1, -3, -2), for all b, k
(proof: STAGE2_SYNTHESIS.md s.2, step 7.2).  So the TR leg of an insertion s in S_w is -Lambda(s), Lambda as in DICT_PROOF s.1.
Checks:
  L1 LD for b <= 4, k <= KMAX, a <= 4, d <= KMAX + 3, mod 2^61-1 and 2^89-1, at v0 = 1 and v0 = 3^5;
  L2 Res_{c_j} f(x) d xif_{a,d} = 0 at EACH branch point (f = log coefficient of each leg): the branch of log p is irrelevant;
  L3 D(ii), exact over Q: log coefficient of Phi_k^b = c_b Sum_i PiC(tau_(k-i) H^b) (-x)^i/i!, PiC from Lemma V;
  L4 frames frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H in K_0 (a != 1, k <= 12), f_j from the PiC kernel: the leg
     -Lambda(frame) is log-free, every term v^m p^j has j >= -m (so it lies in Q[p,q]), f_j has a v-power >= 0, and the
     ancestor vector of the leg, computed by residues, IS the frame.
Controls (must fail): C1 step +Int dx; C2 Phi_0^2 := x; C3 xif_1 without its v-power at v0 = 3^5; C4 log p as a constant.
Usage: py -I -B cr2_legdict.py [KMAX]      (stdlib only; writes nothing; about 75 s for KMAX = 24)
"""
import os
import sys
from fractions import Fraction as Fr
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cr2_legcore as C  # noqa: E402

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
DMAX = KMAX + 3
PRIMES = (2 ** 61 - 1, 2 ** 89 - 1)


def run_ld(P, svals, ctrl=None):
    T = C.legs(KMAX, ctrl)
    bad, l2tot, n = 0, 0, 0
    for s in svals:
        v0, data = C.branch_data(s, P, DMAX, ctrl)
        for (b, k), F in T.items():
            q, l2 = C.ancestor(F, data, v0, P, DMAX, ctrl)
            l2tot += l2
            for (a, d), val in q.items():
                n += 1
                bad += val != ((-C.CB[b]) % P if (a, d) == (b, k) else 0)
    print(" P = 2^%d-1, v0 = s^5 for s in %s, ctrl %s: LD mismatches %d of %d pairings; non-zero L2 residues %d"
          % (P.bit_length(), svals, ctrl, bad, n, l2tot))
    return bad, l2tot


def frame_coeffs(a, k):
    """f[j] = (coefficient, v-power) with Pi_i(z^k H^a + Sum_j f_j z^j H) = 0 for all i (PiC(tau_0 H) = 1)."""
    f = {}
    for i in range(k - 1, -1, -1):
        val = {}
        for (bb, kk, (co, e)) in [(a, k, (Fr(1), 0))] + [(1, j, f[j]) for j in range(i + 1, k)]:
            pc = C.pic(bb, kk - i)
            if pc:
                val[pc[1] + e] = val.get(pc[1] + e, 0) + pc[0] * co
        val = {e: x for e, x in val.items() if x}
        C.check(len(val) <= 1, "frame(%d,%d) not homogeneous" % (k, a))
        f[i] = (-val[min(val)], min(val)) if val else (Fr(0), 0)
    return f


def frames(P, s, kf):
    v0, data = C.branch_data(s, P, DMAX)
    T, n_ok = C.legs(KMAX), 0
    for a in (0, 2, 3, 4):
        for k in range(kf + 1):
            f = frame_coeffs(a, k)
            leg = {}
            for (bb, kk, (co, e)) in [(a, k, (Fr(1), 0))] + [(1, j, f[j]) for j in range(k)]:
                for (j, l, m), x in T[(bb, kk)].items():
                    leg[(j, l, m + e)] = leg.get((j, l, m + e), 0) - x * co / C.CB[bb]
            leg = {key: x for key, x in leg.items() if x}
            ok = C.check(all(l == 0 for (_, l, _) in leg), "frame(%d,%d): leg not log-free" % (k, a))
            ok &= C.check(all(j >= -m for (j, _, m) in leg), "frame(%d,%d): leg outside Q[p,q]" % (k, a))
            ok &= C.check(all(e >= 0 for (co, e) in f.values() if co), "frame(%d,%d): negative v-power" % (k, a))
            q, _ = C.ancestor(leg, data, v0, P, DMAX)
            want = {(a, k): 1}
            for j, (co, e) in f.items():
                want[(1, j)] = C.fr(co, P) * pow(v0, e, P) % P
            ok &= C.check(all(q[key] == want.get(key, 0) % P for key in q), "frame(%d,%d): ancestor vector != frame" % (k, a))
            n_ok += ok
            if (a, k) in ((0, 1), (4, 2)) and P == PRIMES[0]:
                print("  frame(%d,%d): f = %s; leg %s" % (k, a, {j: (str(co), e) for j, (co, e) in f.items() if co},
                                                          {key: str(x) for key, x in sorted(leg.items())}))
    print(" L4 P = 2^%d-1, v0 = %d^5: %d of %d frames (a in 0,2,3,4; k <= %d) pass" % (P.bit_length(), s, n_ok, 4 * (kf + 1), kf))


def main():
    print("cr2_legdict.py KMAX = %d, DMAX = %d" % (KMAX, DMAX))
    T, nd = C.legs(KMAX), 0
    for (b, k), F in T.items():
        f = C.as_poly_in_x(C.split(F)[1])
        want = {}
        for i in range(k + 1):
            pc = C.pic(b, k - i)
            if pc:
                want[(i, pc[1])] = want.get((i, pc[1]), 0) + C.CB[b] * pc[0] * Fr((-1) ** i, factorial(i))
        nd += C.check(f is not None and {x: y for x, y in f.items() if y} == want, "D(ii) at b=%d k=%d" % (b, k))
    print(" L3 D(ii) exact over Q: %d of %d legs" % (nd, len(T)))
    for P in PRIMES:
        bad, l2 = run_ld(P, (1, 3))
        C.check(bad == 0 and l2 == 0, "LD or L2 at P = 2^%d-1" % P.bit_length())
        frames(P, 3, min(KMAX, 12))
    for ctrl, s in (("C1", (1,)), ("C2", (1,)), ("C3", (3,)), ("C4", (1,))):
        bad, _ = run_ld(PRIMES[0], s, ctrl)
        C.check(bad > 0, "control %s did not fail" % ctrl)
        print("  control %s fails as it must: %s" % (ctrl, bad > 0))
    print("RESULT:", "ALL CHECKS PASSED" if not C.FAILS else "%d FAILURES" % len(C.FAILS))


if __name__ == "__main__":
    main()
