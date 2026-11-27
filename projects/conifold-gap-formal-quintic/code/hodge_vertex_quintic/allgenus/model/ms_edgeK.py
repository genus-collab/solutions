# -*- coding: utf-8 -*-
"""ms_edgeK.py -- PREDICTION test (not used in any fit): the edges of R' lie in K_0 (x) K_0 and the leaves in K.

Lemma K (MODEL_SINGULAR.md s.4, proved there from (G) + symplecticity): with D(x) := <phi(1/x), B(x, psi')> the contraction of the
vanishing covector into the first slot of the edge B of a symplectic gauge P with P Y^C = Y^C_Q,  (x + psi') D(x) has no positive power
of x, hence D has no non-negative power: EVERY shifted covector Pi_k (k >= 0; Pi_k(tau_m) = PiC(tau_(m-k)), the descendant time
t^H_k) is contracted to zero, not only Pi_0 (the only one checked in gauge/gi_hat_edge.py).  Leaves: Pi_k(T) = 0 for k >= 1.
Hatted form: phi^(k)_{i,m} = tau^-k phi_{i,m-k} (common factors dropped), phi_{i,m} the hatted vanishing covector of gauge/gi_hat_core.py.
Rows: prod (fq3_par.setup, zhat <= 6, u < 6) or ref (natural_form/referee/ref_hat.HatRows, zhat <= NZ, u < K); production EdgeLeaf.
Imports read-only: project modules, natural_form/referee/ref_hat.py, gauge/gi_core.py and gauge/gi_hat_core.py (helpers only).
Run: py -I -B ms_edgeK.py ref|prod 89|61 [alpha beta] [NZ K] [--c5 x]   (--c5: control R' -> R' exp(x z^5), FAILs expected)
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
GAUGE = os.path.normpath(os.path.join(HERE, "..", "gauge")); sys.path.insert(0, GAUGE)
from gi_core import PROJ
from gi_hat_core import EPS, y_coeffs, edge_degree_bound
sys.path.insert(0, os.path.join(PROJ, "natural_form", "referee"))
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_sser import EdgeLeaf


class Duck(object):
    pass


def setup(src, p, al, be, NZ, K, C5=None):
    R = FpRing(p)
    if src == "ref":
        from ref_hat import HatRows, M0_of
        H = HatRows(p, NZ, K, al, be, zeta=R.zeta); assert all(H.checks.values())
        Qh = [[[Ser(R, 0, list(H.Q[a][i][n])) for n in range(NZ + 1)] for i in range(5)] for a in range(5)]
        M0, M0inv = M0_of(H.Rt, NZ, p)
    else:
        from fq3_par import setup as st_
        st = st_((p, str(al), str(be), 34, K)); H = st['H']; assert H.Nz >= NZ
        Qh = [[[H.Q[a][i][n].truncate(K) for n in range(NZ + 1)] for i in range(5)] for a in range(5)]
        M0inv = st['M0'].Minv
    tauQ = [Fr(1)]; inv_ = [Fr(4, 5 * k + 4) for k in range(K)]
    for k in range(1, K): tauQ.append(-sum(inv_[j] * tauQ[k - j] for j in range(1, k + 1)))
    tau = Ser(R, 0, [R.of(x) for x in tauQ]); taui = tau.inv()
    tp = [Ser.const(R, 1, K)]; tn = [Ser.const(R, 1, K)]
    for k in range(1, NZ + 3): tp.append(tp[-1].mul(tau)); tn.append(tn[-1].mul(taui))
    Qp = [[[None] * (NZ + 1) for i in range(5)] for a in range(5)]
    for a in range(5):
        for i in range(5):
            for n in range(NZ + 1):
                acc = Ser.zero(R, K)
                for k in range(n + 1):
                    for j in range(5):
                        c = M0inv[k][j][i]
                        if c: acc = acc.add(Qh[a][j][n - k].mul(tp[k]).scale(c))
                Qp[a][i][n] = acc
    if C5 is not None:                                   # control: R' -> R' exp(c5 z^5) (flat z^5 = u^4 zhat^5), symplectic, breaks (G)
        c5 = R.of(C5); Qn = [[[None] * (NZ + 1) for i in range(5)] for a in range(5)]
        for a in range(5):
            for i in range(5):
                for n in range(NZ + 1):
                    acc = Ser.zero(R, K); cj = 1; j = 0
                    while 5 * j <= n and 4 * j < K:
                        acc = acc.add(Qp[a][i][n - 5 * j].shift(4 * j).truncate(K).scale(cj)); j += 1; cj = cj * c5 % p * R.inv(j) % p
                    Qn[a][i][n] = acc
        Qp = Qn
    i5 = R.inv(5); S = Duck(); S.R = R; S.Nz = NZ
    S.eta_inv = {(0, 3): i5, (1, 2): i5, (2, 1): i5, (3, 0): i5, (4, 4): i5}
    EL = EdgeLeaf(S, Qp, None)
    z = R.zeta
    psi0inv = [[EPS[b] * pow(z, (-i * (b + 1)) % 5, p) * i5 % p for b in range(5)] for i in range(5)]
    y = y_coeffs(NZ + 2)
    phi = [[tn[m].scale(sum(psi0inv[i][b] * R.of(y[b].get(m, 0)) for b in range(5)) * (p - 1 if m % 2 else 1) % p)
            for m in range(NZ + 2)] for i in range(5)]
    return R, EL, phi


def main():
    src = sys.argv[1]; p = P2 if sys.argv[2] == "89" else P1
    al = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); be = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    NZ = int(sys.argv[5]) if len(sys.argv) > 5 else (16 if src == "ref" else 6)
    K = int(sys.argv[6]) if len(sys.argv) > 6 else (12 if src == "ref" else 6)
    C5 = None
    if "--c5" in sys.argv: C5 = Fr(sys.argv[sys.argv.index("--c5") + 1]); sys.argv = sys.argv[:sys.argv.index("--c5")]
    t0 = time.time(); R, EL, phi = setup(src, p, al, be, NZ, K, C5)
    if C5 is not None: print("CONTROL RUN: R' replaced by R' exp(c5 z^5), c5 = %s (symplectic, violates (G)); FAILs are EXPECTED" % C5)
    print("ms_edgeK: rows %s, p = 2^%s-1, frame (%s, %s), zhat <= %d, u < %d (setup %.0f s)" % (src, sys.argv[2], al, be, NZ, K, time.time() - t0))
    E = {(m, n): EL.E(m, n) for m in range(NZ) for n in range(NZ - m)}
    okN = [N for N in range(K) if edge_degree_bound(N) <= NZ]
    rnd = random.Random(7); allok = True
    ctrl = [[Ser.const(R, rnd.randrange(1, p), K) for m in range(NZ + 2)] for i in range(5)]
    nev = 0
    for kk in range(0, NZ):
        res = {}
        for name, cov in (("control (random covector)", ctrl), ("vanishing", phi)):
            tested = nonzero = nontriv = 0
            for n in range(NZ):
                for j in range(5):
                    acc = Ser.zero(R, K); cnt = 0
                    for m in range(kk, NZ - n):
                        for i in range(5):
                            if cov[i][m - kk].c and E[(m, n)][i][j].c: acc = acc.add(E[(m, n)][i][j].mul(cov[i][m - kk])); cnt += 1
                    for N in okN:
                        tested += 1; nonzero += 1 if acc.coef(N) else 0; nontriv += 1 if cnt >= 2 else 0
            res[name] = nonzero
            if tested == 0: continue
            if name != "vanishing":
                tag = "PASS" if nonzero > 0 else "VACUOUS (the control vanishes too at this truncation; k not counted)"
            else:
                if res["control (random covector)"] == 0: continue
                ok = nonzero == 0; allok = allok and ok; tag = "PASS" if ok else "FAIL"; nev += nontriv
            print("%s (E_%d) %s covector shifted by %d: %d coefficients (u^N, N <= %d), %d non-zero, %d from sums of >= 2 entries"
                  % (tag, kk, name, kk, tested, max(okN), nonzero, nontriv))
    print("   non-trivial vanishing sums over all counted shifts: %d" % nev)
    # leaves: Sum_{i, d >= 2} phi_{i, d-k} T(d)[i] = 0 for k >= 1
    q, r = divmod(NZ + 1, 5); Ncmp = min(min(K, 4 * q + r) + 1, K)
    for kk in range(1, NZ - 3):
        acc = Ser.zero(R, K)
        for d in range(max(2, kk), NZ + 2):
            T = EL.T(d)
            for i in range(5): acc = acc.add(T[i].mul(phi[i][d - kk]))
        acc2 = Ser.zero(R, K)
        for d in range(max(2, kk), NZ + 2):
            T = EL.T(d)
            for i in range(5): acc2 = acc2.add(T[i].mul(ctrl[i][d - kk]))
        ok = acc.truncate(Ncmp).iszero(); allok = allok and ok
        print("%s (T_%d) leaf: Sum phi_{i,d-%d} T(d)[i] = 0 mod u^%d   [control, random covector: %s]"
              % ("PASS" if ok else "FAIL", kk, kk, Ncmp, "non-zero" if not acc2.truncate(Ncmp).iszero() else "ZERO (vacuous)"))
    print("MS_EDGEK %s %s (%s,%s): %s  (%.0f s)" % (src, sys.argv[2], al, be, "ALL PASS" if allok else "FAILURES", time.time() - t0))


if __name__ == "__main__":
    main()
