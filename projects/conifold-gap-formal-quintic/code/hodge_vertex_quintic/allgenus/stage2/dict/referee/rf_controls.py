"""Referee failing controls for D(iii): each control must FAIL the test it targets.  Exact over Q; rf_core only.
Usage: py -I -B rf_controls.py
 K1 Phi_0^2 := x (a function of x): the leg map S_w -> N_w loses rank, first at w = 4 (Lemma 3 needs q = -p^4/4 mod Q[v,x]).
 K2 Phi_0^2 := v^2/p^6 (outside Q[p,q]): legs leave M; some log-free combination leaves Q[p,q] and L.
 K3 Phi_0^0 := p + v/p^4 (outside Q[p,q], weight 1): the same.
 K4 integration that drops the logarithm (Int dp/p := 0): the log map vanishes, 'log-free' = S_w, dim != d(w).
 K5 L-test anchors: 1/p (leg of tau_0(H^2)/v, weight -1), v^2/p^5 (w=5), v/p^2 (w=3) rejected; v/p, v^2/p^2, p^3 accepted.
 K6 a genuine K_0 leg plus one monomial v^m p^j: rejected iff j < -m (detection power of the L test).
 K7 coefficients in v^-1 Q[v] (m >= -1 allowed): log-free combinations with legs outside L exist (sharpness of Q[v]).
"""
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import (Ar, PHI0, W0, XP, mul, lin, legs_iter, integ, I_step, wt, in_Qpq, in_L, in_M, nf_vx,  # noqa: E402
                     x_coords, triples, rank, left_kernel)

ar = Ar(None)
FAILS = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def dd(w):
    return w // 5 + 1 - (1 if w % 5 == 0 else 0)


def logfree_legs(T, w, mmin=0):
    """Kernel of the log map (raw log-part dicts as rows: representation-free) on triples with m >= mmin."""
    S = [(b, k, (w - W0[b] - 4 * k) // 5) for b in range(5) for k in range(max(0, w) // 4 + 3)
         if (w - W0[b] - 4 * k) % 5 == 0 and (w - W0[b] - 4 * k) // 5 >= mmin and (b, k) in T]
    kern = left_kernel(ar, [T[(b, k)][1] for (b, k, m) in S])
    return S, kern, [lin(ar, *[(c, T[S[i][:2]][0]) for i, c in y.items()]) for y in kern]


def safe_legs(K, b0, phi_b0):
    """Genuine legs for b != b0; for b0 the modified start, iterated until a (log p)^2 would appear."""
    T = {key: val for key, val in legs_iter(ar, K).items() if key[0] != b0}
    A0, A1 = ar.lift(phi_b0[0]), ar.lift(phi_b0[1])
    for k in range(K + 1):
        T[(b0, k)] = (A0, A1)
        try:
            A0, A1 = I_step(ar, A0, A1)
        except ValueError:
            break
    return T


def k1():
    phi = {**dict(enumerate(PHI0)), 2: ({4: Fr(1, 4), -1: 1}, {})}
    T = legs_iter(ar, 12, phi)
    bad = []
    for w in range(41):
        S = triples(w)
        rows = [dict([(("p", j), c) for j, c in nf_vx(ar, T[(b, k)][0], wt(b, k))[0].items()] +
                     [(("L", n), c) for n, c in x_coords(ar, T[(b, k)][1], wt(b, k))[0].items()]) for (b, k, m) in S]
        if rank(ar, rows) < len(S):
            bad.append(w)
    check(bool(bad) and bad[0] == 4, "K1 Phi_0^2 := x: leg map S_w -> N_w not injective at w = %s (expected first 4)" % bad[:6])


def k2_k3(label, b0, phi_b0):
    T = safe_legs(12, b0, phi_b0)
    outM = sum(1 for (b, k), (A0, A1) in T.items() if not in_M(ar, A0, A1, wt(b, k)))
    nq = nl = ntot = 0
    for w in range(41):
        S, kern, legs = logfree_legs(T, w)
        for A0 in legs:
            ntot += 1
            nq += not in_Qpq(A0, w)
            nl += not in_L(ar, A0, w)
    check(outM > 0 and nq > 0 and nl > 0, "%s: %d legs outside M; of %d log-free combinations (w <= 40) %d outside Q[p,q], "
          "%d outside L (expected > 0)" % (label, outM, ntot, nq, nl))


def k4():
    def I_nolog(A0, A1):
        xp = ar.lift(XP)
        P0, P1 = mul(ar, A0, xp), mul(ar, A1, xp)
        P0 = {a: c for a, c in P0.items() if a != -1}            # Int dp/p := 0
        P1 = {a: c for a, c in P1.items() if a != -1}
        B0, B1 = integ(ar, P0, P1)
        return lin(ar, (-1, B0)), lin(ar, (-1, B1))
    T = {}
    for b in range(5):
        A0, A1 = ar.lift(PHI0[b][0]), {}                        # no log p anywhere
        for k in range(13):
            T[(b, k)] = (A0, A1)
            A0, A1 = I_nolog(A0, A1)
    bad = [w for w in range(41) if len(logfree_legs(T, w)[1]) != dd(w)]
    check(bool(bad) and bad[0] == 0, "K4 log-dropping integration: dim 'K_0' != d(w) at w = %s (expected from w = 0)" % bad[:6])


def k5():
    tests = [({-1: 1}, -1, False), ({-5: 1}, 5, False), ({-2: 1}, 3, False), ({-1: 1}, 4, True), ({-2: 1}, 8, True),
             ({3: 1}, 3, True), ({-6: 1}, 4, False), ({-3: 1}, 12, True)]
    res = [(in_L(ar, ar.lift(A), w), exp) for A, w, exp in tests]
    check(all(r == e for r, e in res), "K5 L-test anchors (1/p, v^2/p^5, v/p^2, v^2/p^6 rejected; v/p, v^2/p^2, p^3, "
          "v^3/p^3 accepted): %s" % [r for r, e in res])


def k6():
    T = legs_iter(ar, 12)
    nrej = nacc = nbad = ngood = 0
    for w in (13, 21, 29, 36):
        S, kern, legs = logfree_legs(T, w)
        for A0 in legs:
            for m in range(0, w // 5 + w // 4 + 3):
                j = w - 5 * m
                B = lin(ar, (1, A0), (1, {j: Fr(7, 3)}))
                if j >= -m:
                    ngood += 1
                    nacc += in_L(ar, B, w) and in_Qpq(B, w)
                else:
                    nbad += 1
                    nrej += (not in_L(ar, B, w)) and (not in_Qpq(B, w))
    check(nrej == nbad > 0 and nacc == ngood > 0, "K6 K_0 legs + (7/3) v^m p^j, w in (13, 21, 29, 36): %d/%d with j < -m "
          "rejected by in_L and in_Qpq, %d/%d with j >= -m accepted by both" % (nrej, nbad, nacc, ngood))


def k7():
    T = legs_iter(ar, 12)
    nbad = ntot = 0
    first = None
    for w in range(-1, 41):
        S, kern, legs = logfree_legs(T, w, mmin=-1)
        for A0 in legs:
            ntot += 1
            if not in_L(ar, A0, w):
                nbad += 1
                first = first if first is not None else (w, sorted(A0.items())[:3])
    check(nbad > 0, "K7 m >= -1 allowed: %d of %d log-free combinations (w <= 40) have legs outside L; first %s"
          % (nbad, ntot, first))


if __name__ == "__main__":
    print("rf_controls.py (exact over Q): every control below must detect its defect")
    k1()
    k2_k3("K2 Phi_0^2 := v^2/p^6", 2, ({-6: 1}, {}))
    k2_k3("K3 Phi_0^0 := p + v/p^4", 0, ({1: 1, -4: 1}, {}))
    k4()
    k5()
    k6()
    k7()
    print("RESULT:", "ALL CONTROLS BEHAVE AS EXPECTED" if not FAILS else "%d UNEXPECTED" % len(FAILS))
