# -*- coding: utf-8 -*-
"""rf_misc.py -- REFEREE: small exact identities used by GAMMA.md / GAMMA_PROOFS.md (own code, exact over Q).
(M1) Stirling functional equation Gamma_St(y+1) = y Gamma_St(y) in Q[[v]], v = 1/y:
     (1/v + 1/2) log(1 + v) - 1 + Sum_k b_k v^(2k-1) [(1 + v)^(1-2k) - 1] = 0   (input of (S4)-PF, (R3)-PF and the single-Gamma lemma).
(M2) cyclotomic Stirling terms: Sum_j (zeta^j x - 1)^(-m) = (-1)^m 5 Sum_(5|r) C(m-1+r, r) x^r (GAMMA_PROOFS s.3 for odd m) against
     the referee's Newton recurrence for the roots of (1 - x^5) a^5 + 5a^4 + 10a^3 + 10a^2 + 5a + 1.
(M3) the saddle-prefactor identity e^(S_1) (-S_0'')^(-1/2) = const * sigma at sigma = L (gm_fs docstring), via the x-form:
     (-W''(s))^(-1/2) A(s; 0) = 5^(-1/2) L  with W''(s) = -5 s^2/(1 - u): checked at 5 rational s.
(M4) cross-check of the referee's saddle code (rf_sp.fs_series, Newton sums, no roots of unity) against a second route with explicit
     5th roots of unity in F_p61 and Horner composition (rf_sp.fs_generic_modp); control: a non-root in place of zeta.
Run: py -I -B rf_misc.py [N]   log: logs/rf_misc_<N>.log
"""
import os, sys
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_ring import FQ, FP, P61, add, log1, stirling_b
from rf_res import A_x
from rf_sp import fs_series, fs_generic_modp


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    Q = FQ(); os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "rf_misc_%d.log" % N), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(nm)
    lg = log1(Q, [Fr(1), Fr(1)] + [Fr(0)] * (N + 1), N + 2)
    ser = [lg[k + 1] for k in range(N + 1)]                               # (1/v) log(1 + v)
    ser = add(Q, ser, [x / 2 for x in lg[:N + 1]]); ser[0] -= 1
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        for r in range(1, N + 1 - m):                                     # (1+v)^(-m) - 1
            ser[m + r] += stirling_b(k) * (-1) ** r * comb(m + r - 1, r)
    check("(M1) Gamma_St(y+1) = y Gamma_St(y) as formal series in 1/y, through v^%d" % N, not any(ser[:N + 1]))
    ser2 = list(ser); ser2[3] += Fr(1, 7)
    check("(M1) control: perturbed series is detected", any(ser2))
    _, p = A_x(N + 1)
    ok = True
    for m in range(1, 30):
        cf = [Fr(0)] * (N + 1)
        for r in range(0, N + 1, 5): cf[r] = Fr((-1) ** m * 5 * comb(m - 1 + r, r))
        ok &= (p[m][:N + 1] == cf)
    check("(M2) Sum_j (zeta^j x - 1)^(-m) = (-1)^m 5 Sum_(5|r) C(m-1+r, r) x^r, m < 30, through x^%d (Newton recurrence)" % N, ok)
    okp = True
    for s in (Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(-1, 4), Fr(3, 2)):
        u = s ** 5; W2 = -5 * s * s / (1 - u); A0sq = 1 / (1 - u)        # A(s;0)^2
        okp &= (A0sq / (-W2) == Fr(1, 5) / s ** 2)                        # ((-W'')^(-1/2) A(s;0))^2 = L^2/5
    check("(M3) ((-W''(s))^(-1/2) A(s;0))^2 = L^2/5 at s = 1/2, 1/3, 2/3, -1/4, 3/2 (so the saddle prefactor is 5^(-1/2) L)", okp)
    Fp = FP(P61); g = 2                                                   # (M4) own cross-check of rf_sp.fs_series
    while pow(g, (P61 - 1) // 5, P61) == 1: g += 1
    zeta = pow(g, (P61 - 1) // 5, P61); okc = True
    for s in (Fr(1, 2), Fr(5, 7)):
        okc &= fs_series(Fp, s, 14) == fs_generic_modp(Fp, s, 14, zeta)
    bad = fs_generic_modp(Fp, Fr(1, 2), 8, 12345) != fs_series(Fp, Fr(1, 2), 8)
    check("(M4) saddle series: Newton-sum route = explicit 5th roots of unity in F_p61 + composition, L = 2, 7/5, through z^14;"
          " a non-root 'zeta' differs", okc and bad)
    log("RF_MISC %d: %s" % (N, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL))
    OUT.close()


if __name__ == "__main__":
    main()
