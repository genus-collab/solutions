# -*- coding: utf-8 -*-
"""ur_gauge_run.py -- driver for check (C) of UNIFORM_REDUCTION.md.

1. Computes the regular gauge P(u, z) (ur_gauge.gauge) over Q for k <= K from the periods alone and checks that every p_{k,b} is a
   POWER series in u (no pole): the regularity theorem, order by order.
2. Reads off the first u-order n(k, b) at which the entry (b + k, b) of [zhat^k] Ghat appears, hence the predicted zhat-degree of
   [u^n] G for n beyond the range of the production data.
3. Compares with the project's R-matrix: Ghat^proj(u, zhat) := Qhat(u, zhat) R^(tau zhat)^-1 (hatted rows of fq3_msum.HatRowsU, symbol
   of fq_symbol.Symbol; orders zhat^k, k <= 6, and u^n, n <= 5; both primes) must equal the gauge of step 1 up to the one constant
   c_5 of P -> P (1 + c_5 z^5) (it enters at u^4 zhat^5 and u^4 zhat^6 u^(mu) only).
Run: py -I ur_gauge_run.py [K] [NU]      (defaults 12, 34; single process)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import P1, P2
from fq_series import Ser
from ur_gauge import gauge, hatted, first_orders, mu, CEXP
OUT = open(os.path.join(HERE, "ur_gauge_run.log"), "w", encoding="utf-8")
FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def project_ghat(st, tau, KZ=6, KU=6):
    """{(a, b, k): [u^0..u^(KU-1)] mod p} of Qhat(u, zhat) R^(tau zhat)^-1,  (R^(t)^-1)_{i b} = Sum_b' R^_{b'}^{(i)}(-t) eta^{b' b}"""
    R = st['R']; p = R.p; H = st['H']; sym = st['sym']
    tauR = Ser(R, 0, [R.of(x) for x in tau.c[:KU]])
    tp = [Ser.const(R, 1, KU)]
    for m in range(1, KZ + 1): tp.append(tp[-1].mul(tauR))
    out = {}
    for a in range(5):
        for b in range(5):
            for k in range(KZ + 1):
                acc = Ser.zero(R, KU)
                for (b1, b2), e in H.eta_inv.items():
                    if b2 != b: continue
                    for i in range(5):
                        for m in range(k + 1):
                            c = sym.Rh[b1][i][m] * e % p * (1 if m % 2 == 0 else p - 1) % p
                            if c: acc = acc.add(H.Q[a][i][k - m].mul(tp[m]).scale(c))
                out[(a, b, k)] = [acc.coef(n) if n < acc.prec else None for n in range(KU)]
    return out


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    NU = int(sys.argv[2]) if len(sys.argv) > 2 else 34
    t0 = time.time(); log("regular gauge over Q: K = %d, working precision O(u^%d)" % (K, NU))
    res = gauge(K, NU, log); log("   (%.0f s)" % (time.time() - t0))
    tau = res['tau']
    log("tau(u) = " + " + ".join("(%s) u^%d" % (tau.coef(n), n) for n in range(5)))
    v = res['v']; log("v(u) = u tau^(-5/4) = " + " + ".join("(%s) u^%d" % (v.coef(n), n) for n in range(1, 5)))
    hol = all(x is None or x >= 0 for inf in res['info'] for x in inf['vals'])
    check("REGULARITY: every p_{k,b}(u), k <= %d, is a power series in u (no pole at the conifold point)" % K, hol)
    check("resonances (free constants) occur only for k = 0 mod 5, at u^0", all((not inf['resonances']) or (inf['k'] % 5 == 0 and inf['resonances'] == [0]) for inf in res['info']),
          [(inf['k'], inf['resonances']) for inf in res['info'] if inf['resonances']])
    fo = first_orders(res)
    log("first u-order n(k, b) of the entry (a, b) = (b + k, b) of [zhat^k] Ghat   [a priori bound mu(b, k) in brackets]:")
    for k in range(K + 1):
        log("   k = %2d: %s" % (k, "  ".join("b=%d: %s [%d]" % (b, fo[k][b], mu(b, k)) for b in range(5))))
    nk = {k: min(x for x in fo[k].values() if x is not None) for k in range(K + 1)}
    log("first u-order at which zhat^k appears: %s" % nk)
    nmax = nk[K] - 1
    deg = {n: max(k for k in range(K + 1) if nk[k] <= n) for n in range(nmax + 1)}
    log("PREDICTED zhat-degree of [u^n] G, n = 0..%d: %s" % (nmax, [deg[n] for n in range(nmax + 1)]))
    check("degrees 0, 1, 3, 4, 5 for n = 0..4 (the production finding) and >= 6 for n = 5", [deg[n] for n in range(5)] == [0, 1, 3, 4, 5] and deg[5] >= 6)
    # the H (x) H entries (a, b) = (2, 1): k = 1 mod 5
    log("the H (x) H entries p_{k,1}, k = 1 mod 5 (a = 2, b = 1): valuations %s" % {k: res['p'][k][1].v for k in range(1, K + 1, 5)})
    p1 = res['p'][1]
    log("p_{1,b}(u) lowest terms: %s" % ["b=%d: %s" % (b, p1[b].fmt(p1[b].v + 2) if p1[b].c else "0") for b in range(5)])
    check("p_{1,1} = O(u^2): no z H (x) H term at orders u^0, u^1 of the flat gauge (hence none in G_0, G_1)", p1[1].v >= 2)
    hat = hatted(res)
    # ---- comparison with the project's R-matrix
    from fq3_par import setup
    KZ = 6; KU = 6
    for prime, tag in ((P2, "2^89-1"), (P1, "2^61-1")):
        st = setup((prime, "0", "0", 34, KU)); R = st['R']; p = R.p
        G = project_ghat(st, tau, KZ, KU)
        log("-" * 100); log("p = %s" % tag)
        offgrade = [key for key, val in G.items() if (key[0] - key[1] - key[2]) % 5 and any(val)]
        check("Z/5 grading of Ghat^proj: entries (a, b) of zhat^k vanish unless a - b = k mod 5", not offgrade, offgrade[:5])
        rec = {(k, b): [R.of(hat[(k, b)].coef(n)) for n in range(KU)] for k in range(KZ + 1) for b in range(5)}
        bad = [(k, b, n) for k in range(5) for b in range(5) for n in range(KU) if G[((b + k) % 5, b, k)][n] != rec[(k, b)][n]]
        check("Ghat^proj = gauge of step 1 for zhat^0..zhat^4, all entries, u^0..u^5 (no free constant)", not bad, bad[:6])
        bad = [(k, b, n) for k in (5, 6) for b in range(5) for n in range(4) if G[((b + k) % 5, b, k)][n] != rec[(k, b)][n]]
        check("Ghat^proj = gauge of step 1 for zhat^5, zhat^6 through u^3", not bad, bad[:6])
        c5 = (G[(0, 0, 5)][4] - rec[(5, 0)][4]) % p            # rec_{0,0}[u^0] = 1
        try: c5s = st['R'].fmt(c5)
        except Exception: c5s = str(c5)
        log("   c_5 (the constant of P -> P (1 + c_5 z^5), fitted on the entry (0, 0) at u^4 zhat^5): %s" % c5s)
        p51 = res['p'][5][1].coef(0)
        check("CANONICAL NORMALISATION: c_5 = -p_{5,1}(0) = -p_{5,2}(0) (= %s), i.e. the gauge of the actual R-matrix has P_11(0, z) = P_22(0, z) = 1 + O(z^10)" % (-p51),
              (c5 + R.of(p51)) % p == 0 and res['p'][5][2].coef(0) == p51)
        bad = []
        for k in (5, 6):
            for b in range(5):
                for n in range(KU):
                    corr = rec[(k - 5, b)][n - 4] * c5 % p if n >= 4 else 0
                    if (rec[(k, b)][n] + corr) % p != G[((b + k) % 5, b, k)][n]: bad.append((k, b, n))
        check("Ghat^proj = gauge of step 1 times (1 + c_5 z^5) for zhat^5, zhat^6, all entries, u^0..u^5 (one fitted constant, %d identities)" % (2 * 5 * 2 - 1), not bad, bad[:6])
    log("UR_GAUGE_RUN " + ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); OUT.close()


if __name__ == "__main__":
    main()
