# -*- coding: utf-8 -*-
"""gm_compare.py -- identity (Gamma): the full-Stirling Mellin-Barnes saddle series c(z; L) (gm_fs) equals the QRR-normalised
colour-0 WKB series R(z; L) = Sum_n z^n L^(4n) ph_n(L^-5) (gm_qrr), and the free constants are the q = infinity Gamma constants.

Checks (each printed PASS/FAIL):
  (Q)  exact over Q, n <= NQ: c_n(L) = R_n(L) for L in LQ; every R_n(L) a polynomial with exponents = 4n mod 5 ((H): h_auto);
  (P)  mod p (2^61-1 and 2^89-1), n <= NP: D(z) := c(z; L)/R(z; L) = 1 + O(z^(NP+1)) for L in LP (D is L-independent in any case
       if c is a WKB series; D = 1 is the claim);
  (G)  free constants x_m = R_5m(L = 0) = [z^5m] g(z), g = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))), m <= NP/5;
  (C)  controls: perturbing ONE free constant (x_1 or x_5) of the QRR series makes (P) fail, at the predicted order.
Run: py -I -B gm_compare.py NQ NP [tag]      log: logs/gm_compare_<NQ>_<NP>.log (paths relative to this folder)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from gm_series import Rq, Rp, P61, P89, bernoulli_list, bernoulli_poly, sexp0, smul, sinv
from gm_qrr import wkb0, R_of_L
from gm_fs import fs_series

LQ = (Fr(2), Fr(1, 2), Fr(7, 3))
LP = (Fr(2), Fr(3, 5))


def gseries(R, N):
    B = bernoulli_list(N + 2); lg = [Fr(0)] * (N + 1)
    for l in range(1, N // 5 + 1):
        lg[5 * l] = Fr((-1) ** (l + 1)) * bernoulli_poly(5 * l + 1, Fr(1, 5), B) / (l * (5 * l + 1))
    return sexp0(R, [R.of(x) for x in lg])


def first_nonzero(R, s):
    for k, x in enumerate(s):
        if k > 0 and R.norm(x) != 0: return k
    return None


def main():
    NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 16; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    name = "gm_compare_%d_%d.log" % (NQ, NP); OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii")
    FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + (": " + str(det) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("gm_compare: NQ = %d (exact over Q), NP = %d (mod p)" % (NQ, NP))
    R = Rq(); ph, info = wkb0(R, NQ)
    check("(Q0) QRR rows over Q: (H) h_auto for all n <= %d with 5 not | n, resonance coefficients zero" % NQ,
          all(info['h_auto']) and all(info['res_ok']))
    check("(Q1) deg ph_n <= floor(4n/5) for n <= %d, so R_n(L) in Q[L] with exponents 4n - 5k >= 0" % NQ,
          all(len(P) - 1 <= 4 * n // 5 for n, P in enumerate(ph)), info['deg'])
    for L in LQ:
        c = fs_series(R, L, NQ); Rl = R_of_L(R, ph, L)
        bad = [n for n in range(NQ + 1) if c[n] != Rl[n]]
        check("(Q2) L = %s: c_n(L) = R_n(L) exactly over Q for n <= %d" % (L, NQ), not bad, bad[:5])
    log("   e.g. c_5(2) = %s" % fs_series(R, Fr(2), 5)[5])
    log("   (Q) done in %.1f s" % (time.time() - t0))
    for p in (P61, P89):
        Rm = Rp(p); t1 = time.time()
        php, infop = wkb0(Rm, NP)
        check("[%s] (P0) (H) h_auto for n <= %d, resonance coefficients zero" % (Rm.name, NP), all(infop['h_auto']) and all(infop['res_ok']))
        g = gseries(Rm, NP)
        badg = [m for m, x in infop['free'].items() if x != g[5 * m]]
        check("[%s] (G) free constants x_m = R_5m(0) = [z^5m] g(z) for m = 1..%d" % (Rm.name, NP // 5), not badg, badg[:5])
        for L in LP:
            c = fs_series(Rm, L, NP); Rl = R_of_L(Rm, php, L)
            D = smul(Rm, c, sinv(Rm, Rl, NP + 1), NP + 1)
            fn = first_nonzero(Rm, D)
            check("[%s] (P1) L = %s: D(z) = c(z;L)/R(z;L) = 1 + O(z^%d)" % (Rm.name, L, NP + 1), fn is None and D[0] == 1, fn)
        for (mpert, val) in ((1, 1), (NP // 5, 3)):
            if mpert < 1: continue
            fr = dict(infop['free']); fr[mpert] = Rm.norm(fr[mpert] + val)
            phc, _ = wkb0(Rm, NP, free=fr)
            c = fs_series(Rm, LP[0], NP); Rl = R_of_L(Rm, phc, LP[0])
            D = smul(Rm, c, sinv(Rm, Rl, NP + 1), NP + 1)
            fn = first_nonzero(Rm, D)
            check("[%s] (C) control: x_%d shifted by %d -> D - 1 first non-zero at z^%d (expected z^%d)" % (Rm.name, mpert, val, fn if fn else -1, 5 * mpert),
                  fn == 5 * mpert)
        log("   [%s] done in %.1f s" % (Rm.name, time.time() - t1))
    log("GM_COMPARE %d %d: %s   (%.0f s)" % (NQ, NP, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
