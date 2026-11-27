# -*- coding: utf-8 -*-
"""rf_t4_frozen.py -- referee check of TR_MODEL.md s.2-3: TR on x = kappa W (referee engine rf_tr + referee xi basis rf_xi) against the
project's frozen theory MVertex6(M_0) (fq3_mvertex, fq_fp imported READ-ONLY; M_0 rebuilt from natural_form/lg/lg_symbol_rows_Q.json
by LG_MODEL.md s.2 (D): C_n[j][i] = (1/5) Sum_a zeta^(j(4-a) + i(a+1-n)) g^a_n, M_0^-1 = M_0(-z)^T).
kappa = 1, 4, -1 with the derived normalisation Lam^(1-g) kappa^D, Lam = -5 kappa; controls that must FAIL:
  (C1) kappa = -1 read with factor 5^(1-g) and no kappa^D;  (C2) kappa = 1 with Lam = +5;  (C3) MVertex6 built from M_0(-z).
Then flat values (kappa = 1) by rational reconstruction at both primes.  Writes rf_t4_frozen.log.  Run: py -I -B rf_t4_frozen.py
"""
import os, sys, time, json, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
for d in (HERE, PROJ):
    if d not in sys.path: sys.path.insert(0, d)
from rf_ser import FpField, P61, P89
from rf_tr import TR, Mmax
from rf_curves import special
from rf_xi import XiSpecial
from fq_fp import FpRing
from fq3_mvertex import MVertex6, MConst

LIST = ((0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1))


def m0(R, nmax, sign=1):
    g = json.load(open(os.path.join(PROJ, "natural_form", "lg", "lg_symbol_rows_Q.json")))['g']
    p = R.p; M = []; Minv = []
    for n in range(nmax + 1):
        C = [[0] * 5 for _ in range(5)]
        for j in range(5):
            for i in range(5):
                s = sum(R.z(j * (4 - a) + i * (a + 1 - n)) * R.of(Fr(g[a][n])) for a in range(5))
                C[j][i] = s * R.of(Fr(1, 5)) * (sign ** n) % p
        M.append(C); Minv.append([[C[j][i] * (1 if n % 2 == 0 else p - 1) % p for j in range(5)] for i in range(5)])
    return MConst(R, M, Minv, check=True)


def compare(F, col, mv, g, n):
    ncmp = nnz = 0; bad = []
    for d in itertools.product(range(3 * g - 3 + n + 1), repeat=n):
        D = 3 * g - 3 + n - sum(d)
        if D < 0 or D > 6: continue
        T = mv.tensor(g, d)
        for idx in itertools.product(range(5), repeat=n):
            val = T[sum(i * 5 ** l for l, i in enumerate(idx))]
            mine = col.get(tuple(zip(idx, d)), F.zero).v
            ncmp += 1; nnz += (val % F.p != 0)
            if (mine - val) % F.p: bad.append((idx, d))
    return ncmp, nnz, bad


def main():
    log = open(os.path.join(HERE, "rf_t4_frozen.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t4_frozen.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    flatvals = {}
    for p in (P61, P89):
        R = FpRing(p); F = FpField(p); tag = "2^%d-1" % p.bit_length()
        assert F.wv == R.zeta, "zeta conventions differ"
        mv = MVertex6(R, m0(R, 6)); mvneg = MVertex6(R, m0(R, 6, -1))
        for kappa in (1, 4, -1):
            t = time.time(); tr = TR(special(F, kappa), 2 * Mmax(2, 2) + 6); XB = XiSpecial(F, kappa, 7); tot = totnz = 0
            for g, n in LIST:
                W = tr.omega(g, n); col, nbad = XB.frozen(W, g, n)
                ncmp, nnz, bad = compare(F, col, mv, g, n); tot += ncmp; totnz += nnz
                check(not bad and nbad == 0, "%s kappa = %s: (%d,%d) TR x Lam^(1-g) kappa^D == MVertex6(M_0): %d colour correlators (D <= 6), %d non-zero, %d mismatches, %d remainders" % (
                    tag, kappa, g, n, ncmp, nnz, len(bad), nbad))
                if kappa == 1:
                    flatvals.setdefault((g, n), []).append((F, col, XB))
                    if (g, n) in ((0, 3), (1, 1), (2, 1)):
                        b2 = compare(F, col, mvneg, g, n)[2]          # (0,3) has only D = 0 (TFT level): R cannot enter there
                        check(len(b2) > 0 or (g, n) == (0, 3), "%s control C3 (MVertex6 from M_0(-z)) on (%d,%d): %d mismatches (expected > 0 unless D = 0 only)" % (tag, g, n, len(b2)))
                        colw, _ = XB.frozen(W, g, n, Lam=5)
                        b3 = compare(F, colw, mv, g, n)[2]
                        check(len(b3) > 0 or g == 1, "%s control C2 (Lam = +5) on (%d,%d): %d mismatches (expected > 0 unless g = 1)" % (tag, g, n, len(b3)))
                if kappa == -1 and (g, n) in ((1, 1), (2, 1), (1, 2), (0, 4)):
                    colw, _ = XB.frozen(W, g, n, Lam=5, kap=1)
                    b1 = compare(F, colw, mv, g, n)[2]
                    check(len(b1) > 0, "%s control C1 (x = -W read with 5^(1-g), no kappa^D) on (%d,%d): %d mismatches (expected > 0)" % (tag, g, n, len(b1)))
            for g in (2, 3):
                fz = mv.tensor(g, ())[0]; mine = tr.Fg(g) * F.of((Fr(-5) * kappa) ** (1 - g) * Fr(kappa) ** (3 * g - 3))
                check((mine.v - fz) % p == 0, "%s kappa = %s: F_%d(TR) Lam^(1-g) kappa^(3g-3) = %s == MVertex6 vertex value %s" % (tag, kappa, g, F.rat(mine), F.rat(fz)))
            say("  %s kappa = %s: %d correlators compared, %d non-zero (%.1f s)" % (tag, kappa, tot, totnz, time.time() - t))
    # flat values (H^a = N_a Sum zeta^(i(a+1)) et_i), both primes must reconstruct the same rational
    say("flat correlators of the frozen theory from the referee's TR (kappa = 1):")
    wanted = [((1, 1), (1,), (0,)), ((1, 1), (0,), (1,)), ((1, 2), (1, 1), (0, 0)), ((1, 3), (1, 1, 1), (0, 0, 0)),
              ((2, 1), (1,), (0,)), ((2, 1), (0,), (1,)), ((2, 1), (4,), (2,)), ((2, 1), (3,), (3,)), ((2, 1), (2,), (4,)),
              ((2, 2), (1, 1), (0, 0)), ((3, 1), (0,), (1,)), ((3, 1), (1,), (0,)), ((3, 1), (4,), (2,)), ((3, 1), (3,), (3,)),
              ((3, 1), (2,), (4,)), ((3, 1), (1,), (5,)), ((3, 1), (0,), (6,)), ((3, 1), (4,), (7,))]
    for gn, avec, dvec in wanted:
        vals = [Fv.rat(XB.flat(col, avec, dvec)) for (Fv, col, XB) in flatvals[gn]]
        name = " ".join("t%d(H^%d)" % (d, a) for a, d in zip(avec, dvec))
        check(vals[0] is not None and vals[0] == vals[1], "<%s>_%d = %s (2^61-1: %s, 2^89-1: %s)" % (name, gn[0], vals[1], vals[0], vals[1]))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
