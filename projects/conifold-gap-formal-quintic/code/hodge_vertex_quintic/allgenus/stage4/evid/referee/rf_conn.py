# -*- coding: utf-8 -*-
"""rf_conn.py -- REFEREE: Lemma-N-free test of the CONNECTION FORM of EVID.md s.2.4, exact over Q (or mod p), colour 0 only.
  ph^can_N(u) := C0(u) [zhat^N] (Ghat^can R^(tau zhat))_(0, colour 0) = C0 Sum_k Ghat^can_(0,-k)[zhat^k](u) tau^(N-k) F_(-k mod 5, N-k)
with Ghat^can from rf_gauge (own Theorem-U solver from the periods, canonical normalisation P_11(0, z) = 1, regularity ASSERTED, not
assumed) and F_b the model rows (rf_core.model).  No WKB recursion, no K_q, no Lemma N enters ph^can.  Checks, for every N <= NZ:
  C1 ph^can_N is a polynomial: its u-coefficients above floor(4N/5) vanish (u < K), and the degree is exactly floor(4N/5) (N >= 3)
  C2 ph^can_N(1) = h_N (QRR constants, rf_core.qrr), i.e. Phi^can(zhat, 1) = h(zhat): the connection form of D = 1, at EVERY N
  C3 ph^can_N = the QRR WKB polynomial ph_N of rf_core.Wkb (equivalently D = 1 through z^N);  C4 control: the stored gauge fails C2 at N = 5
Run (cwd outside the project): py -I -B rf_conn.py FIELD NZ K [alpha beta NU]      log logs/rf_conn_<FIELD>_<frame>_<NZ>_<K>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
import rf_gauge as G
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def phcan(F, gh, Fm, tp, C0, NZ, K):
    out = []
    for N in range(NZ + 1):
        s = [0] * K
        for k in range(N + 1):
            b = (-k) % 5; c = Fm[b][N - k]
            if c: s = [F.red(x + c * y) for x, y in zip(s, C.mul(F, gh[(k, b)], tp[N - k], K))]
        out.append(C.mul(F, C0, s, K))
    return out


def main():
    a = sys.argv[1:]; fld = a[0]; NZ = int(a[1]); K = int(a[2]); al = Fr(a[3]) if len(a) > 3 else Fr(0); be = Fr(a[4]) if len(a) > 4 else Fr(0)
    NU = int(a[5]) if len(a) > 5 else NZ + K + 10; F = C.QF() if fld == "Q" else C.PF(PR[fld])
    fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "rf_conn_%s_%s_%d_%d.log" % (fld, fr, NZ, K)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)
    t0 = time.time(); log("rf_conn (referee, independent): field %s, frame (%s, %s), N <= %d, u < %d, gauge precision NU = %d; log logs/%s"
                          % (F.name, al, be, NZ, K, NU, name))
    C0 = C.periods(F, K, al, be)[0]; tau = C.inv(F, C.itau(F, K), K); Fm = C.model(F, NZ); h = C.qrr(F, NZ)
    tp = [[F.of(1)] + [0] * (K - 1)]
    for m in range(1, NZ + 1): tp.append(C.mul(F, tp[-1], tau, K))
    res = G.gauge(F, NZ, NU, al, be, "can"); log("   canonical gauge k <= %d in %.0f s (regular: asserted at every order)" % (NZ, time.time() - t0))
    pc = phcan(F, G.hatted(F, res, K), Fm, tp, C0, NZ, K)
    W = C.Wkb(F, F.of(1)); rule = C.qrr_rule(F, h); ph = [[F.of(1)]]
    for n in range(1, NZ + 1): W.step(rule); ph.append(W.ph(n))
    rows = []
    for N in range(NZ + 1):
        d = 4 * N // 5 if N >= 3 else 0; p = pc[N]
        rows.append((N, not any(p[d + 1:]) and p[d] != 0, F.red(sum(p[:d + 1]) - h[N]) == 0, p[:d + 1] == ph[N] + [0] * (d + 1 - len(ph[N]))))
    ok = [r for r in rows if 4 * r[0] // 5 + 1 < K]; Nm = max(r[0] for r in ok)
    check("C1 ph^can_N is a polynomial of degree exactly floor(4N/5) (N >= 3; constant for N <= 2), N <= %d" % Nm, all(r[1] for r in ok),
          [r[0] for r in ok if not r[1]][:5])
    check("C2 connection form ph^can_N(1) = h_N at EVERY N <= %d (5 | N: (N); 5 not | N: (H))" % Nm, all(r[2] for r in ok),
          [r[0] for r in ok if not r[2]][:5])
    check("C3 ph^can_N = ph_N (QRR WKB polynomial from the PF operator), N <= %d" % Nm, all(r[3] for r in ok), [r[0] for r in ok if not r[3]][:5])
    if isinstance(F, C.QF):
        log("   f^can_5m = [u^4m] ph^can_5m: " + "; ".join("f_%d = %s" % (5 * m, pc[5 * m][4 * m]) for m in range(1, min(3, Nm // 5) + 1)))
    ps = phcan(F, G.hatted(F, G.gauge(F, min(NZ, 10), min(NU, 40), al, be, "stored"), min(K, 12)), Fm, tp, C0, min(NZ, 10), min(K, 12))
    dev = F.red(sum(ps[5][:5]) - h[5])
    check("C4 control: the stored gauge gives ph^s_5(1) - h_5 = m5 = -43397/807975000 (fails C2 by exactly m5)", dev == F.of(Fr(-43397, 807975000)), dev)
    log("RF_CONN %s %s N<=%d K=%d: %s   (%.0f s)" % (F.name, fr, NZ, K, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
