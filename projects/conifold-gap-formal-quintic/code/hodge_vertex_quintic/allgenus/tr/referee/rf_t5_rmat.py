# -*- coding: utf-8 -*-
"""rf_t5_rmat.py -- referee check of TR_MODEL.md s.2.3 (R_TR(z) = M_0(z/kappa), dilaton leaf), EXACT over Q(zeta5), own code.

arXiv:1509.06954 (2-31) with x - x(c_j) = -w_j^2/2 and the Gaussian moments Int w^(2n) e^(-w^2/(2z)) = sqrt(2 pi z) (2n-1)!! z^n gives
   X_n[j][i] := [z^n] R^-1(z)^j_i = -(2n-1)!! [w_j^(2n-1)] xi^i     (n >= 0, (-1)!! = 1; the n = 0 term is delta_ij)
and (2-32):   Sum_k X_n[i][k] Delta_k^(-1/2) = (2n+1)!! [w_i^(2n+1)] y.
With the referee's local coordinate zeta_j (x - x(c_j) = (5/2) kappa c_j^2 zeta_j^2): w_j = sqrt(Lam) c_j zeta_j, Lam = -5 kappa,
xi^i = xit^i/sqrt(Lam), Delta_k^(-1/2) = 1/(sqrt(Lam) c_k).  Compared with C_n of M_0 (LG_MODEL.md s.2 (D), rows from
lg_symbol_rows_Q.json): reading A: X_n[j][i] = (-1/kappa)^n C_n[i][j] (i.e. R_TR = M_0(z/kappa)); reading B (transposed).
Writes rf_t5_rmat.log.  Run: py -I -B rf_t5_rmat.py
"""
import os, sys, time, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import Q5Field, Ls
from rf_tr import Local
from rf_curves import special

NM = 10


def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r


def main():
    log = open(os.path.join(HERE, "rf_t5_rmat.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t5_rmat.py  %s  (exact Q(zeta5), z^0..z^%d)" % (time.strftime("%Y-%m-%d %H:%M"), NM))
    F = Q5Field(); g = json.load(open(os.path.join(PROJ, "natural_form", "lg", "lg_symbol_rows_Q.json")))['g']
    C = []
    for n in range(NM + 1):
        C.append([[sum((F.pw(j * (4 - a) + i * (a + 1 - n)) * Fr(g[a][n]) for a in range(5)), F.zero) * Fr(1, 5) for i in range(5)] for j in range(5)])
    for kappa in (1, 4, -1):
        cur = special(F, kappa); Lam = Fr(-5 * kappa); c = cur.c
        locs = [Local(cur, j, 2 * NM + 6) for j in range(5)]
        X = [[[None] * 5 for _ in range(5)] for _ in range(NM + 1)]
        for j in range(5):
            L = locs[j]
            for i in range(5):
                if i == j:
                    s = L.uinv; xi = Ls(-1, s.c, s.hi - 1).scale(-(F.one / c[i]))          # -1/(c_j t)
                else:
                    xi = L.wpow(i, 1).scale(-(F.one / c[i]))                                # -1/(c_i (c_j - c_i + t))
                for n in range(NM + 1):
                    X[n][j][i] = xi.coef(2 * n - 1, F) * Fr(-dfact(2 * n - 1)) / (Lam ** n) / _pw(c[j], 2 * n - 1, F)
        okA = [all(X[n][j][i] == C[n][i][j] * Fr(-1, kappa) ** n for i in range(5) for j in range(5)) for n in range(NM + 1)]
        okB = [all(X[n][j][i] == C[n][j][i] * Fr(-1, kappa) ** n for i in range(5) for j in range(5)) for n in range(NM + 1)]
        say("kappa = %s (x = %s W, Lam = %s):" % (kappa, kappa, Lam))
        check(all(okA), "(2-31) Laplace transform of B, row = branch point: X_n = (-1/kappa)^n C_n^T for n = 0..%d (all 25 entries); orders passing: %s" % (NM, [n for n in range(NM + 1) if okA[n]]))
        say("     transposed reading X_n = (-1/kappa)^n C_n holds only for n = %s" % [n for n in range(NM + 1) if okB[n]])
        ok32 = True; okM = True; T0 = []
        for n in range(NM + 1):
            for i in range(5):
                lhs = sum((X[n][i][k] / c[k] for k in range(5)), F.zero)
                rhs = locs[i].t.coef(2 * n + 1, F) * Fr(dfact(2 * n + 1)) / (Lam ** n) / _pw(c[i], 2 * n + 1, F)
                ok32 = ok32 and lhs == rhs
                okM = okM and lhs == sum((C[n][k][i] / c[k] for k in range(5)), F.zero) * Fr(-1, kappa) ** n
            if kappa == 1 and 1 <= n <= 4: T0.append(-sum((X[n][0][k] / c[k] for k in range(5)), F.zero))
        check(ok32, "(2-32) (Laplace transform of dy) = R^-1(z) applied to the unit, n = 0..%d, all colours" % NM)
        check(okM, "(2-32) = M_0(z/kappa)^-1 1^ with 1^_k = 1/c_k (frozen unit), n = 0..%d" % NM)
        if kappa == 1:
            check(all(x.a[1] == x.a[2] == x.a[3] == 0 for x in T0), "dilaton leaf T(z) = z(1^ - M_0^-1 1^), colour 0, z^2..z^5: %s" % [str(x.a[0]) for x in T0])
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


def _pw(a, e, F):
    r = F.one
    if e < 0: a = F.one / a; e = -e
    for _ in range(e): r = r * a
    return r


if __name__ == "__main__":
    main()
