# -*- coding: utf-8 -*-
"""critic_gauge_vs_table.py -- COMPLETENESS CRITIC, one light cross-check between two reports that were never compared directly.

(theory)   the regular gauge P(u, z) of UNIFORM_REDUCTION.md s.3, computed over Q from the periods alone
           (here: the verifier's proof-route implementation verify_theory/v2_gauge.py, imported read-only), was compared with the
           actual G only for zhat <= 6, u <= 5;
(referee)  the table referee/G_fourier_NZ16_K12.json holds all 397 non-zero Fourier components of the actual G for zhat <= 16, u <= 11
           (from-scratch rows; confirmed on the project rows by verify_referee).
This script predicts every Fourier component of [u^n zhat^k] G, n <= 11, from the gauge with the CANONICAL normalisation
P_11(0, z) = 1 (imposed at z^5 and z^10; z^15 lies at u^12, outside the table) and compares with the table in both directions.
Dictionary (verify_theory/v3_fourier.py):  Xhat(m, m') = (1/5) eps_a eps_b tau^((c_a - c_b)/4) u^mu(b,k) p_{k,b}(u) / p_{0,a}(u),
a = b + k = m - 1, b = -m' - 1 (mod 5), eps = (1,-1,-1,1,1), c = (0,-1,3,2,1), mu(b,k) = (4k + c_b - c_a)/5.
Frame (0, 0) only.  Exact rationals; single process.  Run: py -I -B critic_gauge_vs_table.py [K] [N]
"""
import os, sys, json, time
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "verify_theory")))
from v2_gauge import gauge
CEXP = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)
OUT = open(os.path.join(HERE, "critic_gauge_vs_table.log"), "w", encoding="utf-8")


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 36
    UMAX = 11
    t0 = time.time(); res = gauge(K, N); P = res['p']; tau = res['tau']
    log("gauge (proof route, periods only), K = %d, N = %d: %.0f s; precisions of p_k: %s"
        % (K, N, time.time() - t0, [min(x.prec for x in P[k]) for k in range(K + 1)]))
    # canonical normalisation: P -> P c(z^5) with P_11(0, z) = 1
    cc = []
    for m in range(1, K // 5 + 1):
        c = P[5 * m][1].coef(0) / P[0][1].coef(0); cc.append(-c)
        for k in range(K, 5 * m - 1, -1):
            P[k] = [P[k][b] - P[k - 5 * m][b].scale(c) for b in range(5)]
    log("canonical normalisation P_11(0,z) = 1: scalar factor c(z^5) = 1 + (%s) z^5 + (%s) z^10" % tuple(str(x) for x in (cc + [None])[:2]))
    log("   expected from the other reports: c_5 = 43397/807975000 (fit to the R-matrix); c_10 = 153074543599907/261129440250000000 (symplecticity)")
    ok5 = cc[0] == Fr(43397, 807975000); ok10 = len(cc) > 1 and cc[1] == Fr(153074543599907, 261129440250000000)
    log(("PASS " if ok5 else "FAIL ") + "c_5 agrees"); log(("PASS " if ok10 else "FAIL ") + "c_10 agrees (canonical = symplectic normalisation at z^10)")
    r1 = all(P[5 * m][2].coef(0) == 0 for m in range(1, K // 5 + 1))
    log(("PASS " if r1 else "FAIL ") + "P_22(0,z) = P_11(0,z) = 1 at z^5, z^10 after the normalisation ((R1))")
    # prediction of the Fourier components
    pred = {}; untested = 0
    for k in range(0, K + 1):
        for b in range(5):
            a = (b + k) % 5; mu = (4 * k + CEXP[b] - CEXP[a]) // 5
            if mu > UMAX: continue
            f = tau.powq(Fr(CEXP[a] - CEXP[b], 4)) * P[k][b] * P[0][a].inv()
            for n in range(mu, UMAX + 1):
                if n - mu >= f.prec: untested += 1; continue
                val = f.coef(n - mu)
                if val: pred[(n, k, (a + 1) % 5, (-b - 1) % 5)] = Fr(EPS[a] * EPS[b], 5) * val
    tab = json.load(open(os.path.join(HERE, "..", "referee", "G_fourier_NZ16_K12.json"), encoding="utf-8"))["G_fourier"]
    table = {}
    for key, val in tab.items():
        un, zk, mm = key.split(); m, mp = mm.strip("()").split(",")
        table[(int(un[1:]), int(zk[1:]), int(m), int(mp))] = Fr(val)
    table = {k: v for k, v in table.items() if k[0] <= UMAX and k[1] <= K}
    both = set(pred) & set(table); bad = [k for k in both if pred[k] != table[k]]
    only_pred = sorted(set(pred) - set(table)); only_tab = sorted(set(table) - set(pred))
    log("table entries (u <= %d, zhat <= %d): %d;  predicted non-zero components: %d;  coefficients skipped for lack of precision: %d"
        % (UMAX, K, len(table), len(pred), untested))
    log(("PASS " if not bad else "FAIL ") + "all %d common components agree exactly (rationals)%s" % (len(both), "" if not bad else ": first bad " + str(bad[:3])))
    log(("PASS " if not only_pred else "FAIL ") + "no predicted component is missing from the table%s" % ("" if not only_pred else ": " + str(only_pred[:6])))
    log(("PASS " if not only_tab else "NOTE ") + "table components not reproduced: %d %s" % (len(only_tab), only_tab[:8]))
    kmax = max(k[1] for k in both); nmax = max(k[0] for k in both)
    log("   highest orders compared: zhat^%d, u^%d;  components with zhat >= 7: %d" % (kmax, nmax, sum(1 for k in both if k[1] >= 7)))
    log("CRITIC_GAUGE_VS_TABLE " + ("ALL PASS" if ok5 and ok10 and r1 and not bad and not only_pred and not only_tab else "SEE ABOVE")); OUT.close()


if __name__ == "__main__":
    main()
