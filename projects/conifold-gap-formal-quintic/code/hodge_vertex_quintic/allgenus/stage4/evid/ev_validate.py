# -*- coding: utf-8 -*-
"""ev_validate.py -- Stage 4 evidence: validation of the colour-0 machinery (ev_core.py, used by ev_reduced.py) against the referee's
from-scratch modules natural_form/referee/ref_rows.py, ref_hat.py and the direct actual gauge of gauge/gi_actual_norm.py (all imported
read-only), mod p, in the window of the direct test (EVID.md s.3).
  V1 periods: C0 = varpi_0, C1 = D t_c equal ref_hat.HatRows(..).C0, .C1 in the frame given
  V2 QRR constants h_n (exact over Q via Newton power sums) equal ref_rows.hcoef (colour 0, n <= NZ + 4)
  V3 WKB polynomials ph_n (incremental recursion, QRR free coefficients) equal ref_rows.wkb (colour 0), n <= NZ + 4
  V4 hatted rows Qh_a (colour 0) equal HatRows.Q[a][0], all a, zhat^n (n <= NZ), u^k (k < K)
  V5 Prop. S: the model rows F_b (ev_core.model_F, from the model alone) equal the quintic's own symbol Psi0 M_0 = Qhat(u = 0) for all
     b, all colours i, n <= NZ:  [zhat^n] R^_b^(i) = zeta^(i(b+1-n)) F_{b,n}
  V6 the colour-reduced pairing  Sum_{n<=N} Qh_{a,n}(u) (-tau)^(N-n) F_{3-c,N-n}  [N = a - c mod 5, else 0]  equals gi_actual_norm.actual's
     Ghat^act = Qhat M_0(tau zhat)^-1 Psi0^-1 (five colours, matrix inverse) in all 25 entries, all N <= NZ, all u^k, k < K
  V7 the same computation over Q (exact) reduced mod p equals the mod-p computation (ph_n, rows, F_b): the two code paths agree
Run (cwd outside the project): py -I -B ev_validate.py [alpha beta NZ K primes]  (defaults 0 0 35 29 89,61)
   log: logs/ev_validate_<frame>_<NZ>_<K>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(PROJ, "allgenus", "gauge"))
import ev_core as E
import gi_actual_norm as GA                                   # read-only: actual (direct Ghat^act); it puts ref_* on sys.path
from ref_rows import wkb, hcoef                               # read-only
from ref_hat import HatRows, M0_of                            # read-only
from fq_fp import FpRing, P1, P2
EPS = (1, -1, -1, 1, 1)


def main():
    a = sys.argv[1:]
    al = Fr(a[0]) if len(a) > 0 else Fr(0); be = Fr(a[1]) if len(a) > 1 else Fr(0)
    NZ = int(a[2]) if len(a) > 2 else 35; K = int(a[3]) if len(a) > 3 else 29
    primes = [{"89": (P2, "2^89-1"), "61": (P1, "2^61-1")}[x] for x in (a[4] if len(a) > 4 else "89,61").split(",")]
    fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    name = "ev_validate_%s_%d_%d.log" % (fr, NZ, K)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)

    t0 = time.time()
    log("ev_validate: frame (%s, %s), window zhat <= %d, u < %d; log file logs/%s" % (al, be, NZ, K, name))
    QF = E.QF(); hq = E.qrr_h(QF, NZ + 4); Wq = E.Wkb(QF); rq = E.qrr_free(QF, hq)
    for n in range(1, NZ + 5): Wq.step(rq)
    C0q, C1q, _, _ = E.periods(QF, K, al, be); Rq = E.Rows(QF, C0q, C1q, K, amax=4)
    for n in range(NZ + 1): Rq.set(n, Wq.ph[n])
    Fq = E.model_F(QF, NZ)
    log("   exact colour-0 data over Q in %.1f s" % (time.time() - t0))
    for (p, pn) in primes:
        t1 = time.time(); F = E.PF(p); zeta = FpRing(p).zeta
        H = HatRows(p, NZ, K, al, be, zeta=zeta)
        C0, C1, _, _ = E.periods(F, K, al, be)
        check("[%s] V1 periods C0, C1 (u < %d) = ref_hat" % (pn, K), C0 == [x % p for x in H.C0] and C1 == [x % p for x in H.C1])
        h = E.qrr_h(F, NZ + 4)
        check("[%s] V2 h_n = ref_rows.hcoef, n <= %d" % (pn, NZ + 4), h == [x % p for x in hcoef(p, zeta, 0, NZ + 4)])
        W = E.Wkb(F); rule = E.qrr_free(F, h)
        for n in range(1, NZ + 5): W.step(rule)
        Phi, info = wkb(p, zeta, 0, NZ + 4)
        check("[%s] V3 ph_n = ref_rows.wkb, n <= %d (degrees %s..)" % (pn, NZ + 4, info['deg'][:12]),
              all(E.padd(F, W.ph[n], []) == E.padd(F, Phi[n], []) for n in range(NZ + 5)))
        R = E.Rows(F, C0, C1, K, amax=4)
        for n in range(NZ + 1): R.set(n, W.ph[n])
        check("[%s] V4 rows Qh_a (colour 0) = HatRows.Q[a][0], a = 0..4, n <= %d, u < %d" % (pn, NZ, K),
              all(R.Q[a_][n] == [x % p for x in H.Q[a_][0][n]] for a_ in range(5) for n in range(NZ + 1)))
        Fm = E.model_F(F, NZ); M0, _ = M0_of(H.Rt, NZ, p)
        bad5 = []
        for n in range(NZ + 1):
            for b in range(5):
                for i in range(5):
                    rh = sum(EPS[b] * pow(zeta, (j * (b + 1)) % 5, p) * M0[n][j][i] for j in range(5)) % p   # (Psi0 M0)[b][i]
                    if rh != pow(zeta, (i * (b + 1 - n)) % 5, p) * Fm[b][n] % p: bad5.append((n, b, i))
        check("[%s] V5 Prop. S: model rows F_b = quintic symbol Psi0 M_0, all b, i, n <= %d" % (pn, NZ), not bad5, bad5[:4])
        _, G, psi0_ok, chk = GA.actual(p, NZ, K, al, be, lambda *x: None)
        tau = E.tau_series(F, K); tp = [[1] + [0] * (K - 1)]
        for m in range(1, NZ + 1): tp.append(E.smul(F, tp[-1], tau, K))
        bad6 = []; nz6 = 0
        for a_ in range(5):
            for c in range(5):
                for N in range(NZ + 1):
                    red = E.pair(F, R.Q[a_], Fm[(3 - c) % 5], tp, N, K) if (N - a_ + c) % 5 == 0 else [0] * K
                    if red != [x % p for x in G[(N, a_, c)]]: bad6.append((N, a_, c))
                    nz6 += sum(1 for x in red if x)
        check("[%s] V6 colour-reduced pairing = direct Ghat^act (gi_actual_norm.actual): 25 entries x %d orders x %d u-terms (%d non-zero)"
              % (pn, NZ + 1, K, nz6), not bad6 and psi0_ok, bad6[:4])
        ok7 = all(E.padd(F, [F.of(x) for x in Wq.ph[n]], []) == E.padd(F, W.ph[n], []) for n in range(NZ + 5))
        ok7 = ok7 and all([F.of(x) for x in Rq.Q[a_][n]] == R.Q[a_][n] for a_ in range(5) for n in range(NZ + 1))
        ok7 = ok7 and all([F.of(x) for x in Fq[b]] == Fm[b] for b in range(5))
        check("[%s] V7 exact-Q computation reduced mod p = mod-p computation (ph_n, rows, F_b)" % pn, ok7)
        log("   [%s] done in %.0f s" % (pn, time.time() - t1))
    log("EV_VALIDATE %s NZ=%d K=%d: %s   (%.0f s)" % (fr, NZ, K, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED,
                                                    time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
