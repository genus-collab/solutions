# -*- coding: utf-8 -*-
"""gi_hat_edge.py -- the vanishing covector against the ACTUAL edge and leaf tables of R' = Psi(u) G(u, zhat) (mod p).

Rows: (ref)  natural_form/referee/ref_hat.HatRows (from-scratch hatted rows, zhat <= NZ, u < K; any conifold frame), or
      (prod) fq3_par.setup -> fq3_msum.HatRowsU + the production symbol M_0 (zhat <= 6, u < 6).
R' rows:  Q'[a][i][n] = Sum_{k,j} Qhat[a][j][n-k] tau^k (M_0^-1_k)[j][i]      (M = M_0(tau zhat), right placement).
Edge and leaf tables: the PRODUCTION class fq_sser.EdgeLeaf (E(m, n)[i][j] = [psi^m psi'^n] E'_ij, T(d)[i] = [z^d] T'_i), unchanged.
Checks (all exact mod p):
 (L)  Lemma E (edge contraction), regression with RANDOM covectors phi (support m <= 3):
        Sum_{i,m} phi_{i,m} E(m, n)[i][j] = [psi'^n] Sum_ab eta^ab Q'^(j)_b(-psi') h_a(psi')/psi',   h = [R'(x) phi(-x)]_(>0)
 (E)  the vanishing covector phi_{i,m} = (-1)^m [x^-m] (Psi0^-1 y(tau x))_i is contracted to ZERO by the edge tables
        (all n, j; all u^N with edge_degree_bound(N) <= NZ; the bound is also checked against the tables themselves)
 (T)  leaf:  u tau^(-1/2) Sum_{i,d>=2} phi_{i,d} T(d)[i] = -t_c + u tau^(-3/2)/xhat;   Q'^(i)_{0,0} = zeta^i/xhat
 (C)  control: with the covector of the wrong sign convention (phi_{i,m} without (-1)^m) the contraction does NOT vanish.
Run: py -I -B gi_hat_edge.py ref|prod 89|61 [alpha beta] [NZ K]      log: gi_hat_edge_<src>_<prime>_<frame>.log
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import PROJ
from gi_hat_core import CEXP, EPS, y_coeffs, edge_degree_bound
sys.path.insert(0, os.path.join(PROJ, "natural_form", "referee"))
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_sser import EdgeLeaf
from fq_conifold import Conifold


class Duck(object):
    pass


def main():
    src = sys.argv[1]; prime = sys.argv[2]; p = P2 if prime == "89" else P1
    al = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); be = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    NZ = int(sys.argv[5]) if len(sys.argv) > 5 else (16 if src == "ref" else 6)
    K = int(sys.argv[6]) if len(sys.argv) > 6 else (12 if src == "ref" else 6)
    tag = "%s_%s_ab_%s_%s" % (src, prime, str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "gi_hat_edge_%s.log" % tag), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(name)

    t0 = time.time(); R = FpRing(p); z = R.zeta
    log("gi_hat_edge: rows = %s, p = 2^%s-1, frame (alpha, beta) = (%s, %s), zhat <= %d, u < %d" % (src, prime, al, be, NZ, K))
    # ---------------------------------------------------------------- rows and symbol
    if src == "ref":
        from ref_hat import HatRows, M0_of
        H = HatRows(p, NZ, K, al, be, zeta=z)
        check("referee rows: own checks (h_auto, psi_orth, sympl, covar, circ)", all(H.checks.values()), H.checks)
        Qh = [[[Ser(R, 0, list(H.Q[a][i][n])) for n in range(NZ + 1)] for i in range(5)] for a in range(5)]
        M0, M0inv = M0_of(H.Rt, NZ, p)
    else:
        from fq3_par import setup
        st = setup((p, str(al), str(be), 34, K)); H = st['H']
        assert H.Nz >= NZ
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
    i5 = R.inv(5); S = Duck(); S.R = R; S.Nz = NZ
    S.eta_inv = {(0, 3): i5, (1, 2): i5, (2, 1): i5, (3, 0): i5, (4, 4): i5}
    EL = EdgeLeaf(S, Qp, None)
    psi0 = [[Qh[b][i][0].coef(0) for i in range(5)] for b in range(5)]
    check("Psi0[b][i] = Qhat_{b,0}^(i)(u = 0) = eps_b zeta^(i(b+1))", all(psi0[b][i] == EPS[b] * pow(z, (i * (b + 1)) % 5, p) % p for b in range(5) for i in range(5)))
    psi0inv = [[EPS[b] * pow(z, (-i * (b + 1)) % 5, p) * i5 % p for b in range(5)] for i in range(5)]
    y = y_coeffs(NZ + 2)
    # phi[i][m]: Ser in u (without the common factor u tau^(-1/2))
    def covector(signed=True):
        phi = [[None] * (NZ + 2) for i in range(5)]
        for i in range(5):
            for m in range(NZ + 2):
                c = sum(psi0inv[i][b] * R.of(y[b].get(m, 0)) for b in range(5)) % p
                if signed and m % 2: c = (-c) % p
                phi[i][m] = tn[m].scale(c)
        return phi
    E = {(m, n): EL.E(m, n) for m in range(NZ) for n in range(NZ - m)}
    log("   rows, R' and %d edge tables E(m, n), m + n + 1 <= %d: %.0f s" % (len(E), NZ, time.time() - t0))
    # empirical degree of the edge tables: largest m + n + 1 with a non-zero [u^N] entry
    top = {}
    for (m, n), tab in E.items():
        for i in range(5):
            for j in range(5):
                s = tab[i][j]
                for N in range(min(K, s.prec)):
                    if s.coef(N): top[N] = max(top.get(N, 0), m + n + 1)
    log("   largest zhat-order m + n + 1 present in the edge tables at u^N, N = 0..%d: %s" % (K - 1, [top.get(N, 0) for N in range(K)]))
    log("   bound edge_degree_bound(N) (two factors of G, Theorem D):                %s" % [edge_degree_bound(N) for N in range(K)])
    okN = [N for N in range(K) if edge_degree_bound(N) <= NZ]
    check("the edge tables respect the degree bound wherever it is below the truncation (N <= %d)" % max(okN), all(top.get(N, 0) <= edge_degree_bound(N) for N in okN))
    # ---------------------------------------------------------------- (L) Lemma E with random covectors
    rnd = random.Random(20261113); bad = 0; tested = 0
    for trial in range(3):
        ph = [[rnd.randrange(p) for m in range(4)] for i in range(5)]
        h = {}                                               # h[a][r], r >= 1
        for a in range(5):
            for r in range(1, NZ + 1):
                acc = Ser.zero(R, K)
                for m in range(4):
                    if r + m > NZ: continue
                    for i in range(5):
                        acc = acc.add(Qp[a][i][r + m].scale(ph[i][m] * (-1) ** m))
                h[(a, r)] = acc
        for n in range(NZ - 4):
            for j in range(5):
                lhs = Ser.zero(R, K)
                for m in range(4):
                    for i in range(5): lhs = lhs.add(E[(m, n)][i][j].scale(ph[i][m]))
                rhs = Ser.zero(R, K)
                for (a, b), e in S.eta_inv.items():
                    for r in range(1, n + 2):
                        k = n + 1 - r
                        rhs = rhs.add(Qp[b][j][k].mul(h[(a, r)]).scale(e * (-1) ** k))
                tested += 1
                if not lhs.sub(rhs).truncate(K).iszero(): bad += 1
    check("(L) Lemma E: contraction of a covector into an edge = R'^-1(psi') [R'(x) phi(-x)]_(>0)(psi')/psi'  (3 random covectors, %d series identities mod u^%d)" % (tested, K), bad == 0, bad)
    # ---------------------------------------------------------------- (E) the vanishing covector
    for signed, name in ((True, "(E) vanishing covector"), (False, "(C) control (no (-1)^m)")):
        phi = covector(signed); tested = 0; nonzero = []; nontriv = 0
        for n in range(NZ):
            for j in range(5):
                acc = Ser.zero(R, K); terms = 0
                for m in range(NZ - n):
                    for i in range(5):
                        if phi[i][m].c and E[(m, n)][i][j].c: acc = acc.add(E[(m, n)][i][j].mul(phi[i][m])); terms += 1
                for N in okN:
                    tested += 1
                    if acc.coef(N): nonzero.append((n, j, N))
                    cnt = sum(1 for m in range(NZ - n) for i in range(5) if phi[i][m].c and E[(m, n)][i][j].c and E[(m, n)][i][j].truncate(N + 1).c)
                    if cnt >= 2: nontriv += 1
        if signed:
            check("%s: Sum_{i,m} phi_{i,m} E(m, n)[i][j] = 0 for all n < %d, j, and u^N, N <= %d (%d coefficients, %d of them sums of >= 2 non-zero table entries)"
                  % (name, NZ, max(okN), tested, nontriv), not nonzero, nonzero[:5])
        else:
            check("%s: the contraction does NOT vanish (%d non-zero of %d)" % (name, len(nonzero), tested), len(nonzero) > 0)
    # ---------------------------------------------------------------- (T) leaf
    C = Conifold(R, Nu=K + 4, alpha=al, beta=be)
    tc = Ser(R, C.tc.v, [R.of(c) for c in C.tc.c]).truncate(K + 1)
    xh = Ser(R, C.varpi0.v, [R.of(c) for c in C.varpi0.c]).truncate(K)
    check("Q'^(i)_{0,0} = zeta^i/xhat (the unit of the hatted TFT), all colours, mod u^%d" % K,
          all(Qp[0][i][0].sub(xh.inv().scale(pow(z, i, p))).truncate(K).iszero() for i in range(5)))
    phi = covector(True); acc = Ser.zero(R, K)
    for d in range(2, NZ + 2):
        T = EL.T(d)
        for i in range(5): acc = acc.add(T[i].mul(phi[i][d]))
    half = Ser(R, 0, [R.of(x) for x in tauQ])               # tau^(-1/2) via the Q-series
    from fq_fp import QRing
    QQ = QRing(); tq = Ser(QQ, 0, tauQ); lt = tq.log()
    tm12 = lt.scale(Fr(-1, 2)).exp(); tm32 = lt.scale(Fr(-3, 2)).exp()
    tm12 = Ser(R, 0, [R.of(x) for x in tm12.c]); tm32 = Ser(R, 0, [R.of(x) for x in tm32.c])
    lhs = acc.mul(tm12).shift(1)
    rhs = tm32.mul(xh.inv()).shift(1).sub(tc)
    # completeness: row 0 of zhat^k enters at u^(4q + r), k = 5q + r; the first missing order is k = NZ + 1
    q, r = divmod(NZ + 1, 5); Nmax = min(K, 4 * q + r)       # lhs known mod u^(Nmax + 1) (the factor u), and mod u^K
    Ncmp = min(Nmax + 1, K)
    check("(T) LEAF: u tau^(-1/2) Sum_{i, d >= 2} phi_{i,d} T(d)[i] = -t_c + u tau^(-3/2)/xhat  mod u^%d" % Ncmp, lhs.sub(rhs).truncate(Ncmp).iszero(),
          "" if lhs.sub(rhs).truncate(Ncmp).iszero() else lhs.sub(rhs).fmt(4))
    w1 = R.of(C.varpi0.coef(1)); g2 = R.of(C.tc.coef(2))
    check("    leading term: Pi(T') = -(gamma_t + w_1 - 2/3) u^2 = (p_{1,4}(0)/2) u^2  [frame (0,0): -u^2/30]", lhs.coef(2) == (-(g2 + w1 - R.of(Fr(2, 3)))) % p and lhs.coef(1) == 0)
    log("GI_HAT_EDGE %s: %s   (%.0f s)" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
