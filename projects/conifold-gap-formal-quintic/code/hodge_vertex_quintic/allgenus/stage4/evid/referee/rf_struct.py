# -*- coding: utf-8 -*-
"""rf_struct.py -- REFEREE: independent re-derivation of the structural claims of EVID.md s.2.1 and s.8.2 (uses rf_core only).
  T1 (H) alone: Phi(zhat, u = 1) = h_n at every 5 not | n <= N (WKB from the PF operator, QRR free coefficients)
  T2 deg ph_n = floor(4n/5) for 3 <= n <= N, ph_1 and ph_2 constant
  T3 ph_n(0) = F_(0,n) (model rows of rf_core.model) for n <= N
  T4 every resonance (5 | n) trivially unobstructed: the right-hand side never reaches u^(4n/5)
  T5 (exact field only) p_m = Sum_(j=1..4) (1 - zeta^j)^-m = 0 for every odd multiple m of 5, m <= N, and p_m != 0 for every other
     m <= N with m odd (Galois trace in Q(zeta5)); hence [z^m] log h_(w_0) = B_(m+1)/(m(m+1) 5^m) at those m
Run (cwd outside the project): py -I -B rf_struct.py FIELD N     FIELD Q | 89 | 61;   log logs/rf_struct_<FIELD>_<N>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def main():
    fld = sys.argv[1]; N = int(sys.argv[2]); F = C.QF() if fld == "Q" else C.PF(PR[fld])
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "rf_struct_%s_%d.log" % (fld, N)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)
    t0 = time.time(); log("rf_struct (referee, independent): field %s, n <= %d; log logs/%s" % (F.name, N, name))
    h = C.qrr(F, N); Fm = C.model(F, N); W = C.Wkb(F, F.of(1)); rule = C.qrr_rule(F, h)
    log("   h_n (n <= %d) and model rows in %.0f s" % (N, time.time() - t0))
    Hbad, dbad, zbad = [], [], []
    for n in range(1, N + 1):
        W.step(rule); ph = W.ph(n)
        if n % 5 and F.red(sum(ph) - h[n]) != 0: Hbad.append(n)
        if len(ph) - 1 != (4 * n // 5 if n >= 3 else 0): dbad.append((n, len(ph) - 1))
        if F.red(ph[0] - Fm[0][n]) != 0: zbad.append(n)
        if n in (100, 200, 300, 500, 1000): log("   n = %d in %.0f s" % (n, time.time() - t0))
    check("T1 (H): Phi(zhat, 1) = h at all %d orders n <= %d with 5 not | n" % (N - N // 5, N), not Hbad, Hbad[:5])
    check("T2 deg ph_n = floor(4n/5) (3 <= n <= %d), ph_1 = %s, ph_2 = %s constant" % (N, W.ph(1), W.ph(2)), not dbad, dbad[:5])
    check("T3 ph_n(0) = F_(0,n) of the model, n <= %d" % N, not zbad, zbad[:5])
    check("T4 resonances: %d of %d multiples of 5 had a right-hand side reaching u^(4n/5) (all unobstructed)" % (len(W.res), N // 5),
          all(ok for (_, ok) in W.res))
    E = min(N, 300); badT = [(e, j) for e in range(-1, -E - 1, -1) for j in range(1, 6) if W.T(e)[j] and max(W.T(e)[j]) > -1]
    check("T6 orbifold lemma: for -%d <= e <= -1 every T_j(s^e), j = 1..5, has only s-exponents <= -1 (P maps x Q[x][z] to itself, "
          "x = 1/s); with T_1 this forces deg ph_n <= 4n/5 and trivial resonances for ALL n" % E, not badT, badT[:3])
    low = []
    for m in range(1, 6):                                       # Pi(x^m) = (1/3125) prod_k (5M + kz) x^m, M x^a = x^(a+1) - (a z/5)(x^a - x^(a+5))
        t = {(0, m): F.of(1)}
        for k in range(1, 6):
            nx = {}
            for (j, a), c in t.items():
                for key, v in (((j, a + 1), 5 * c), ((j + 1, a), (k - a) * c), ((j + 1, a + 5), a * c)): nx[key] = nx.get(key, 0) + v
            t = {key: F.red(v) for key, v in nx.items() if F.red(v)}
        low.append(min(a for (j, a) in t))
    check("T7 prod_(k=1..5)(5M + kz) x^m has x-valuation >= 6 for m = 1..5 (the orbifold indicial roots 1..5): valuations %s" % low,
          all(v >= 6 for v in low))
    if isinstance(F, C.QF):
        log("   free (top) coefficients f_5m = [u^4m] ph_5m (QRR): " + ", ".join("f_%d = %s" % (5 * m, W.psi[5 * m].get(-1)) for m in (1, 2, 3)))
        p = C.powersums_Q(N)
        z5 = [m for m in range(5, N + 1, 10) if p[m] != 0]; nz = [m for m in range(1, N + 1, 2) if m % 5 and p[m] == 0]
        check("T5 p_m = 0 for all odd multiples m of 5 (m <= %d), and p_m != 0 for the other odd m" % N, not z5 and not nz, (z5[:3], nz[:3]))
    log("RF_STRUCT %s n<=%d: %s   (%.0f s)" % (F.name, N, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
