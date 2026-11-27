# -*- coding: utf-8 -*-
"""ev_direct.py -- Stage 4 evidence, part (1): the direct test E4 (gauge/gi_actual_norm.py) extended to a larger window (EVID.md s.1).

Actual side (no gauge input): the referee's from-scratch hatted rows natural_form/referee/ref_hat.HatRows (QRR normalisation at q = 0,
all five colours, zhat <= NZ, u < K, own checks h_auto/psi_orth/sympl/covar/circ), the symbol M_0 from the same rows, and
        Ghat^act = Qhat M_0(tau zhat)^-1 Psi0^-1          (gauge/gi_actual_norm.actual, imported read-only).
Canonical side (no R-matrix input): natural_form/theory/ur_gauge.gauge over Q via gauge/gi_core.gauge_frame (read-only), normalised by
P_11(0, z) = 1, hatted by gi_actual_norm.hatted_can.  Checks per prime:
  (0) Ghat^act is Z/5-graded and [zhat^k] Ghat^act_{b+k,b} = O(u^mu(b,k))  [= (H) + Theorem U in the window];
  (1) Ghat^act = Ghat^can in every slot (k <= NZ, N < K, N >= mu(b,k));
  (C) control: the stored normalisation (p^s_{5q,0}(0) = 0) differs, first at zhat^5;
  (2) normalisation slots c^act_{5q} = [u^4q zhat^5q] Ghat^act_11 = [..] Ghat^act_22 for all q with 5q <= NZ, 4q < K, and
      log c^act = Sum a_j z^(5j): a_j for odd j are the INDEPENDENT tests (a_1, a_3, a_5 = z^25, a_7 = z^35), even j forced by
      symplecticity (c(z) c(-z) = 1).
Run (working directory outside the project):  py -I -B ev_direct.py [alpha beta NZ K primes NU]
      defaults 0 0 35 29 89,61 NZ+K+10;  log: logs/ev_direct_<frame>_<NZ>_<K>.log  (paths printed relative to this folder)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(PROJ, "allgenus", "gauge"))
from gi_core import gauge_frame, p11_at_zero, inv_poly, normalise       # read-only
import gi_actual_norm as GA                                               # read-only: actual, hatted_can, mu
from fq_fp import FpRing, P1, P2


def series_log(c, Q, p=None):
    """log(1 + Sum_q c[q] x^q) mod x^(Q+1) over F_p (p given) or exactly over Q (p None, c Fractions): list over q"""
    L = [0] * (Q + 1)
    for n in range(1, Q + 1):            # n L_n = n c_n - Sum_{k<n} k L_k c_(n-k)   (from c' = c L')
        s = n * c[n] - sum(k * L[k] * c[n - k] for k in range(1, n))
        L[n] = s % p * pow(n, p - 2, p) % p if p else Fr(s, 1) / n
    return L


def main():
    a = sys.argv[1:]
    al = Fr(a[0]) if len(a) > 0 else Fr(0); be = Fr(a[1]) if len(a) > 1 else Fr(0)
    NZ = int(a[2]) if len(a) > 2 else 35; K = int(a[3]) if len(a) > 3 else 29
    primes = [{"89": (P2, "2^89-1"), "61": (P1, "2^61-1")}[x] for x in (a[4] if len(a) > 4 else "89,61").split(",")]
    NU = int(a[5]) if len(a) > 5 else NZ + K + 10
    fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    name = "ev_direct_%s_%d_%d%s.log" % (fr, NZ, K, "" if len(primes) == 2 else "_" + primes[0][1][2:4])
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)

    t0 = time.time()
    log("ev_direct: frame (%s, %s); referee rows zhat <= %d, u < %d; canonical gauge over Q, k <= %d, working precision O(u^%d)"
        % (al, be, NZ, K, NZ, NU))
    log("   log file logs/%s" % name)
    res = gauge_frame(NZ, NU, al, be); ps = res['p']; tau = res['tau']
    c11 = p11_at_zero(ps); pc = normalise(ps, inv_poly(c11, NZ))
    can = GA.hatted_can(pc, tau, NZ, K); sto = GA.hatted_can(ps, tau, NZ, K)
    log("   canonical gauge in %.0f s; stored normalisation P^s_11(0, z) = %s + ..." % (time.time() - t0,
        " + ".join("(%s) z^%d" % (c11[k], k) for k in sorted(c11) if k <= 10)))
    Qmax = min(NZ // 5, (K - 1) // 4)
    for (p, pname) in primes:
        t1 = time.time()
        R, G, psi0_ok, chk = GA.actual(p, NZ, K, al, be, log)
        log("-- p = %s: actual Ghat from the referee rows in %.0f s" % (pname, time.time() - t1))
        check("[%s] referee rows pass their own checks (h_auto = (H) for 5 not | n <= %d, psi_orth, sympl, covar, circ); Psi0 ok"
              % (pname, NZ), all(chk.values()) and psi0_ok)
        grad = all(not any(G[(n, a_, b)]) for (n, a_, b) in G if (a_ - b - n) % 5)
        vals = all(not any(G[(n, (b + n) % 5, b)][:min(GA.mu(b, n), K)]) for n in range(NZ + 1) for b in range(5))
        check("[%s] (0) Ghat^act Z/5-graded, [zhat^k] Ghat^act_{b+k,b} = O(u^mu(b,k)) for all k <= %d" % (pname, NZ), grad and vals)
        tested = nonzero = 0; bad = []; dst = []
        for n in range(NZ + 1):
            for b in range(5):
                a_ = (b + n) % 5; g = G[(n, a_, b)]
                for N in range(GA.mu(b, n), K):
                    tested += 1; c = R.of(can[(n, b)][N]); s = R.of(sto[(n, b)][N])
                    if g[N]: nonzero += 1
                    if g[N] != c: bad.append((n, a_, b, N))
                    if g[N] != s: dst.append((n, a_, b, N))
        check("[%s] (1) Ghat^act = Ghat^can in all %d slots (k <= %d, N < %d; %d non-zero)" % (pname, tested, NZ, K, nonzero),
              not bad, bad[:5])
        check("[%s] (C) control: stored normalisation differs from the actual gauge in %d slots, first at (k, a, b, N) = %s"
              % (pname, len(dst), min(dst) if dst else None), bool(dst) and min(x[0] for x in dst) == 5)
        cs = [1] + [G[(5 * q, 1, 1)][4 * q] for q in range(1, Qmax + 1)]
        c22 = [1] + [G[(5 * q, 2, 2)][4 * q] for q in range(1, Qmax + 1)]
        L = series_log(cs, Qmax, p); Ls = series_log([Fr(1)] + [Fr(c11[5 * q]) for q in range(1, Qmax + 1)], Qmax)
        log("   [%s] normalisation slots c^act_{5q} = [u^4q zhat^5q] Ghat^act_11 (= Ghat^act_22), log c^act = Sum a_q z^(5q):" % pname)
        for q in range(1, Qmax + 1):
            sr = str(Ls[q]) if len(str(Ls[q])) < 40 else "%s.. (%d/%d digits)" % (str(Ls[q])[:12], len(str(Ls[q].numerator)),
                                                                                 len(str(Ls[q].denominator)))
            log("      z^%-3d: c^act = %s, Ghat_22 slot = %s, a_%d = %s   [stored gauge, exact: %s]   %s" % (
                5 * q, cs[q], c22[q], q, L[q], sr,
                "INDEPENDENT (odd)" if q % 2 else "forced by symplecticity"))
        check("[%s] (2) c^act = 1 through z^%d: a_1 = a_3 = ... = a_%d = 0 (independent), even orders forced and also 0"
              % (pname, 5 * Qmax, Qmax if Qmax % 2 else Qmax - 1), all(x == 0 for x in cs[1:] + c22[1:]))
    log("EV_DIRECT %s NZ=%d K=%d: %s   (%.0f s)" % (fr, NZ, K, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED,
                                                    time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
