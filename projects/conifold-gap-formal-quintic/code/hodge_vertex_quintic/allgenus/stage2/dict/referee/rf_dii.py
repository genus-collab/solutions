"""Referee: D(ii) in the exact form used by ../DICT_PROOF.md, re-derived independently (rf_core only).
Usage: py -I -B rf_dii.py FIELD KD N1 N2 WD        FIELD in {Q, 61, 89}
 D1 log coefficient of every iterated leg = Sum_i R_(k-i)^b (-x)^i/i!,  R_n^b := Res_{p=0} W^n dPhi_0^b / n!
    (generating-function lemma, REFEREE.md s.3), k <= KD.
 D2 R_n^b = c_b PiC(tau_n(H^b)), c = (-1, 1, -1, -3, -2), with PiC = (-1)^(n+1) [z^-n] Y^C_b,
    Y^C = (I_0, -I_-1, -I_3, I_2, I_1), I^C_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!  (Stage-1 closed form);
    (a) W^n by repeated multiplication, n <= N1; (b) binomial expansion, n <= N2; anchors = the seven values listed
    in model/MODEL_SINGULAR.md s.1.
 D3 K_0 := {Pi_i = 0 for all i >= 0}, Pi_i(v^m tau_k(H^b)) = v^m PiC(tau_(k-i)(H^b)), computed from the closed form ONLY,
    equals ker(log map) computed from the legs, w <= WD (same subspace of S_w; dimension d(w)).
 Controls (must fail): c_3 = +3; the sign (-1)^n instead of (-1)^(n+1); Pi_i shifted by one in D3; D3 without c_b.
"""
import os
import sys
import time
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import Ar, P61, P89, DPHI0, W0, mul, legs_iter, x_coords, triples, rank, left_kernel  # noqa: E402

FIELD = sys.argv[1] if len(sys.argv) > 1 else "Q"
KD = int(sys.argv[2]) if len(sys.argv) > 2 else 60
N1 = int(sys.argv[3]) if len(sys.argv) > 3 else 200
N2 = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
WD = int(sys.argv[5]) if len(sys.argv) > 5 else 160
ar = Ar({"Q": None, "61": P61, "89": P89}[FIELD])
CB = (-1, 1, -1, -3, -2)
YC = ((0, 1), (-1, -1), (3, -1), (2, 1), (1, 1))      # Y^C_b = sgn * I_m
FAILS = []
T0 = time.time()


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg, flush=True)


def pic(b, n, sign_exp=1):
    """PiC(tau_n(H^b)) as (coefficient, v-exponent); (0, None) if zero."""
    m, sgn = YC[b]
    if (n - m - 1) % 5 or n - m - 1 < 0:
        return ar.q(0), None
    i = (n - m - 1) // 5
    ell = m + 4 * i + 1
    if ell < 0:
        return ar.q(0), None
    return ar.q(sgn * (-1) ** (i + ell + n + sign_exp), 4 ** i * factorial(i) * factorial(ell)), ell


def eq(a, b):
    return not ar.nz(a - b)


def res_rep(N):
    """R_n^b, n <= N, with W^n by repeated multiplication."""
    out = {}
    Wd = ar.lift({4: ar.q(1, 4), -1: 1})
    for b in range(5):
        d, Wn = ar.lift(DPHI0[b]), {0: ar.q(1)}
        for n in range(N + 1):
            out[(b, n)] = mul(ar, Wn, d).get(-1, 0) * ar.q(1, factorial(n))
            Wn = mul(ar, Wn, Wd)
    return out


def res_binom(b, n):
    (s, c), = DPHI0[b].items()
    t5 = 4 * n + s + 1
    if t5 % 5 or not 0 <= t5 // 5 <= n:
        return ar.q(0)
    t = t5 // 5
    return ar.q(c * comb(n, t), 4 ** (n - t) * factorial(n))


def d2(cb=CB, sign_exp=1, label="D2"):
    bad = 0
    for (b, n), r in RREP.items():
        pc, e = pic(b, n, sign_exp)
        bad += (not eq(r, cb[b] * pc)) or (e is not None and 5 * e != W0[b] + 4 * n)
    for b in range(5):
        for n in range(N2 + 1):
            pc, e = pic(b, n, sign_exp)
            bad += not eq(res_binom(b, n), cb[b] * pc)
    return bad


def d1(T):
    bad = 0
    for (b, k), (A0, A1) in T.items():
        xc, ok = x_coords(ar, A1, W0[b] + 4 * k)
        want = ar.clean({i: RREP[(b, k - i)] * ar.q((-1) ** i, factorial(i)) for i in range(k + 1)})
        bad += not ok or ar.clean(xc) != want
    return bad


def d3(T, shift=0, scale=True):
    bad = 0
    for w in range(WD + 1):
        S = triples(w)
        prow = []
        for (b, k, m) in S:
            r = {}
            for i in range(k + 2):
                pc, e = pic(b, k - i + shift) if k - i + shift >= 0 else (0, None)
                if e is not None and ar.nz(pc):
                    r[i] = pc
            prow.append(r)
        kp = left_kernel(ar, prow)
        kraw = left_kernel(ar, [ar.clean(x_coords(ar, T[(b, k)][1], W0[b] + 4 * k)[0]) for (b, k, m) in S])
        # raw-leg coordinates y_t (Sum y_t v^m Phi_k^b log-free) -> insertion coordinates a_t = c_b y_t (Lambda = c_b^-1 Phi)
        kl = [ar.clean({i: c * (CB[S[i][0]] if scale else 1) for i, c in y.items()}) for y in kraw]
        d = w // 5 + 1 - (1 if w % 5 == 0 else 0)
        same = len(kp) == len(kl) == rank(ar, kp + kl) == d
        bad += not same
    return bad


if __name__ == "__main__":
    print("rf_dii.py FIELD=%s KD=%d N1=%d N2=%d WD=%d" % (FIELD, KD, N1, N2, WD), flush=True)
    anchors = [((1, 0), 1, 0), ((0, 1), -1, 1), ((4, 2), ar.q(-1, 2), 2), ((3, 3), ar.q(-1, 6), 3),
               ((2, 4), ar.q(1, 24), 4), ((1, 5), ar.q(1, 96), 4), ((0, 6), ar.q(-1, 480), 5)]
    for (b, n), c, e in anchors:
        pc, pe = pic(b, n)
        check(eq(pc, ar.q(c)) and pe == e, "anchor PiC(tau_%d(H^%d))" % (n, b))
    print("anchors (model note s.1): 7 values of PiC checked", flush=True)
    RREP = res_rep(max(N1, KD))
    nb = d2()
    check(nb == 0, "D2 mismatches %d" % nb)
    print("D2 R_n^b = c_b PiC(tau_n(H^b)): repeated products n <= %d, binomial n <= %d, mismatches %d [%.0fs]"
          % (max(N1, KD), N2, nb, time.time() - T0), flush=True)
    T = legs_iter(ar, max(KD, WD // 4 + 1))
    nb = d1({key: val for key, val in T.items() if key[1] <= KD})
    check(nb == 0, "D1 mismatches %d" % nb)
    print("D1 log coefficient of Phi_k^b = Sum_i R_(k-i)(-x)^i/i!, k <= %d: mismatches %d [%.0fs]" % (KD, nb, time.time() - T0),
          flush=True)
    nb = d3(T)
    check(nb == 0, "D3 failing weights %d" % nb)
    print("D3 common kernel of the Pi_i (closed form) = ker(log map) (legs), dim d(w), w <= %d: failing weights %d [%.0fs]"
          % (WD, nb, time.time() - T0), flush=True)
    print("== Controls (must fail)")
    c1 = d2(cb=(-1, 1, -1, 3, -2))
    c2 = d2(sign_exp=0)
    c3 = d3(T, shift=1)
    c4 = d3(T, scale=False)
    print(" c_3 = +3: %d mismatches; sign (-1)^n: %d mismatches; Pi_i shifted by one: %d failing weights; raw legs"
          " without the c_b rescaling: %d failing weights (of %d)" % (c1, c2, c3, c4, WD + 1))
    check(c1 > 0 and c2 > 0 and c3 > 0 and c4 > 0, "a control did not fail")
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS), "[%.0fs]" % (time.time() - T0))
