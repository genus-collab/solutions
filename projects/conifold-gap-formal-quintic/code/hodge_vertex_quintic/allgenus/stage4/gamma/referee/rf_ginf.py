# -*- coding: utf-8 -*-
"""rf_ginf.py -- REFEREE: the free constants x_m = R_5m(L = 0) (q = infinity values of the QRR-normalised colour-0 WKB) against
Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))), the Stirling factor of prod_j Gamma(zeta^j/z + 1/5), at high order
mod two primes (2^127-1, 2^107-1).  Ginf is computed twice: Bernoulli-polynomial closed form, and the direct re-expansion of
log Gamma(y + a) - [(y + a - 1/2) log y - y + log(2 pi)/2] in 1/y (rf_run_gamma.ginf_direct), both exactly over Q.
Run: py -I -B rf_ginf.py N      log: logs/rf_ginf_<N>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_ring import FP, P127
from rf_wkb import wkb, qrr_h
from rf_run_gamma import ginf_bern, ginf_direct, P107


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "rf_ginf_%d.log" % N), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("rf_ginf: free constants vs Ginf through z^%d (m <= %d)" % (N, N // 5))
    gb = ginf_bern(N); log("     closed form over Q: %.1f s" % (time.time() - t0))
    gd = ginf_direct(N); log("     direct re-expansion over Q: %.1f s" % (time.time() - t0))
    check("(G0) closed form = direct Stirling re-expansion of prod_j Gamma(zeta^j/z + 1/5), exact over Q, z^0..z^%d" % N, gb == gd)
    for p in (P127, P107):
        F = FP(p); R, info = wkb(F, N, qrr_h(F, N), checks=False)
        bad = [m for m, x in info['free'].items() if x != F.of(gb[5 * m])]
        check("[%s] (G1) x_m = R_5m(0) = [z^5m] Ginf for m = 1..%d" % (F.name, N // 5), not bad, bad[:4])
        bad2 = [m for m in info['free'] if F.of(gb[5 * m] + Fr(1, 10 ** 6)) == info['free'][m]]
        check("[%s] (C) control: Ginf + 10^-6 at every z^5m matches no x_m" % F.name, not bad2)
    log("RF_GINF %d: %s   (%.0f s)" % (N, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
