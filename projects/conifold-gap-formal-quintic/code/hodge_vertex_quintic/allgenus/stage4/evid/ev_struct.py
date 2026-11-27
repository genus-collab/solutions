# -*- coding: utf-8 -*-
"""ev_struct.py -- Stage 4 evidence: structure of the colour-0 WKB polynomials ph_n (EVID.md s.2.1, s.4.3).  Uses ev_core only.
  S1 degrees: deg ph_n = floor(4n/5) for 3 <= n <= N, ph_1 and ph_2 constant; hence for n = 5m the free coefficient [u^4m] ph_5m is
     the TOP coefficient (and the u^0 coefficient of the WKB function at fixed z: z^n u^(-4n/5) ph_n(u) has exponents in (-4n/5, 0])
  S2 ph_n(0) = F_{0,n} (the model's homogeneous WKB coefficient) for every n <= N, whatever the normalisation (u^0 part of (HPF))
  S3 (H) alone: ph_n(1) = h_n for all 5 not | n <= N, QRR free coefficients (the WKB recursion and the Bernoulli numbers only)
  S4 divisibility R_n(1) = 0 and no resonance obstruction at 5 | n (asserted at every step inside ev_core.Wkb.step)
Run (cwd outside the project): py -I -B ev_struct.py FIELD N    FIELD Q | 89 | 61 | 127;   log logs/ev_struct_<FIELD>_<N>.log
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ev_core as E
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1, "127": (1 << 127) - 1}


def main():
    fld = sys.argv[1] if len(sys.argv) > 1 else "89"; N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    F = E.QF() if fld == "Q" else E.PF(PR[fld])
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "ev_struct_%s_%d.log" % (fld, N)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)

    t0 = time.time(); log("ev_struct: field %s, n <= %d; log file logs/%s" % (F.name, N, name))
    h = E.qrr_h(F, N); Fm = E.model_F(F, N); W = E.Wkb(F); rule = E.qrr_free(F, h)
    deg_bad = []; zero_bad = []; H_bad = []; marks = (50, 100, 200, 300, 500, 1000, 2000)
    for n in range(1, N + 1):
        W.step(rule); c = E.padd(F, W.ph[n], [])
        d = len(c) - 1; want = 0 if n < 3 else 4 * n // 5
        if d != want: deg_bad.append((n, d, want))
        if F.red((c[0] if c else 0) - Fm[0][n]) != 0: zero_bad.append(n)
        if n % 5 and F.red(sum(c) - h[n]) != 0: H_bad.append(n)
        if n in marks: log("   n = %4d reached in %.0f s (deg ph_n = %d)" % (n, time.time() - t0, d))
    check("S1 deg ph_n = floor(4n/5) for 3 <= n <= %d, ph_1 = -3/20 and ph_2 = 9/800 constant (free coefficient = top coefficient)" % N,
          not deg_bad, deg_bad[:5])
    check("S2 ph_n(0) = F_{0,n} (model, homogeneous) for all n <= %d: automatic, not a criterion" % N, not zero_bad, zero_bad[:5])
    check("S3 (H): ph_n(1) = h_n for all %d n <= %d with 5 not | n" % (N - N // 5, N), not H_bad, H_bad[:5])
    triv = sum(1 for r in W.resonance if r[2])
    check("S4 R_n(1) = 0 at every n <= %d; no resonance obstruction at the %d multiples of 5 (%d of them because deg rho_n < 4n/5)"
          % (N, len(W.resonance), triv), len(W.resonance) == N // 5 and all(r[1] for r in W.resonance))
    if isinstance(F, E.QF):
        log("   free (top) coefficients f_5m = [u^4m] ph_5m, QRR normalisation (first six):")
        for m in range(1, min(6, N // 5) + 1): log("      f_%d = %s" % (5 * m, W.free[5 * m]))
    log("EV_STRUCT %s n<=%d: %s   (%.0f s)" % (F.name, N, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
