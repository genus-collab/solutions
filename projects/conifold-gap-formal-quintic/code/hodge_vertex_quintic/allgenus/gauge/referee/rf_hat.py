# -*- coding: utf-8 -*-
"""rf_hat.py -- REFEREE: the hatted statements of GAUGE_IDENTITY_HAT.md, EXACT over Q and in the FLAT frame (no zeta, no EdgeLeaf):
with Q' = Ghat Psi0 and f = Psi0^-1 Yhat^mod one has E' = Psi0^T Ehat Psi0, Ehat(psi,psi') = [eta^-1 - Ghat(-psi)^T eta^-1 Ghat(-psi')]/(psi+psi')
(uses (Psi0 Psi0^T)^-1 = eta^-1), so
  (E)  Sum_{i,m} phi_{i,m} E'(m,n)_ij = 0 for all j   <=>   Sum_m (-1)^m (Yhat^mod_m)^T Ehat(m, n) = 0      (Yhat^mod_m := [x^-m] Yhat^mod)
  (L)  Sum_m (-1)^m Y_m^T Ehat(m, n) = [psi'^(n+1)] h(psi')^T eta^-1 Ghat(-psi'),  h = [Ghat(x) Y(x)]_(>0), any covector Y
  (T)  Pi(T') = Sum_{d>=2} Sum_b [x^(d-1)] Ghat_0b [x^-d] Yhat^mod_b = u tau^(-3/2)/varpi_0 - t_c
Ghat = canonical gauge (project ur_gauge, certified in rf_gi/rf_frob_run), hatted here.  Completeness by the a-priori degree bound
k <= (5N+4)/4 of [u^N] Ghat (mu(b,k) >= (4k-4)/5).  Controls: no (-1)^m; the STORED normalisation.
Run: py -I -B rf_hat.py al be NZ K      log: rf_hat_<frame>.log
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import systems, quintic_vc, ic_terms, powf, mul, inv, CEXP, EPS
from rf_ur import stored_gauge, p11_zero, normalise
J3 = lambda a: (3 - a) % 5


def mu(b, k): return (4 * k + CEXP[b] - CEXP[(b + k) % 5]) // 5


def main():
    al = Fr(sys.argv[1]); be = Fr(sys.argv[2]); NZ = int(sys.argv[3]); K = int(sys.argv[4])
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "rf_hat_%s.log" % tag), "w", encoding="ascii"); FAIL = []
    def log(s): print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(name, ok, det=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(name)
    t0 = time.time()
    p, prec = stored_gauge(NZ, 3 * K + 24, al, be); pc, ci = normalise(p, p11_zero(p)); S = systems(K + 4, al, be)
    tau = S['tau'][:K + 2]; deg = lambda N: (5 * N + 4) // 4
    def hat(g):
        G = {}
        for b in range(5):
            f = powf(tau[:K], Fr(1 - CEXP[b], 4), K)
            for k in range(NZ + 1):
                m = mu(b, k)
                if m < K: G[(k, (b + k) % 5, b)] = [Fr(0)] * m + mul(f, g[k][b][:K], K)[:K - m]
        return G
    Gc = hat(pc); Gs = hat(p)
    Z = [Fr(0)] * K
    ent = lambda G, k, a, b: G.get((k, a, b), Z) if 0 <= k <= NZ else Z
    # hatted model covector Yhat^mod_m[b] = u tau^(-1/2-m) [t^-m] y_b
    ycoef = [{} for b in range(5)]
    for b in range(5):
        for (n, j, c) in ic_terms(CEXP[b], NZ + 2): ycoef[b][n] = ycoef[b].get(n, 0) + c * EPS[b]
    Ym = []
    for m in range(NZ + 2):
        tp = powf(tau[:K], Fr(-1, 2) - m, K); Ym.append([[Fr(0)] + [ycoef[b].get(m, 0) * x for x in tp[:K - 1]] for b in range(5)])
    log("rf_hat: frame (%s, %s), zhat <= %d, u < %d, exact Q (%.0f s)" % (al, be, NZ, K, time.time() - t0))
    # symplecticity of the hatted canonical gauge: Sum_k (-1)^k Ghat_{M-k}^T eta^-1 Ghat_k = eta^-1 delta_M0
    bad = 0
    for M in range(NZ + 1):
        for c in range(5):
            for d in range(5):
                acc = [Fr(0)] * K
                for k in range(M + 1):
                    for a in range(5):
                        x = mul(ent(Gc, M - k, a, c), ent(Gc, k, J3(a), d), K)
                        acc = [s + (-1) ** k * y / 5 for s, y in zip(acc, x)]
                want = Fr(1, 5) if (M == 0 and (c + d) % 5 == 3) else 0
                if acc[0] != want or any(acc[1:]): bad += 1
    check("hatted canonical gauge symplectic: Ghat(-x)^T eta^-1 Ghat(x) = eta^-1, zhat^0..zhat^%d, mod u^%d" % (NZ, K), bad == 0, bad)
    # hatted (G): positive part of Ghat Yhat^mod vanishes; the rest is Yhat^Q
    NQ = K + (4 * (NZ + 2)) // 5 + 3; YQ = quintic_vc(systems(NQ + 2, al, be)['b'], NQ, NZ + 2)    # y^(n) needed to u^(N + (4n + c_a - 4)/5)
    tested = 0; badp = []; badn = []; okN = [N for N in range(K) if (5 * N - 1) // 4 <= NZ]
    for a in range(5):
        for e in range(-NZ, NZ + 1):
            vec = [Fr(0)] * K
            for k in range(max(e, 0), NZ + 1):
                if k - e > NZ + 1: continue
                b = (a - k) % 5; x = mul(ent(Gc, k, a, b), Ym[k - e][b], K); vec = [s + y for s, y in zip(vec, x)]
            for N in okN:
                if e + (5 * N - 1) // 4 > NZ or (5 * N - 1) // 4 - e > NZ + 1: continue      # all k <= (5N-1)/4, m = k - e <= NZ + 1
                tested += 1
                if e > 0 and vec[N]: badp.append((a, e, N))
                if e <= 0:
                    n = -e; w5 = 4 - CEXP[a] - 4 * n; sh = N - w5 // 5 if w5 % 5 == 0 else None
                    if sh is not None: assert n in YQ and sh < len(YQ[n][a]), (a, n, sh)
                    want = YQ[n][a][sh] if (sh is not None and sh >= 0) else 0
                    if vec[N] != want: badn.append((a, e, N))
    check("hatted (G): no positive power of zhat in Ghat Yhat^mod, and the rest = Yhat^Q (%d complete slots)" % tested, not badp and not badn, (badp + badn)[:4])
    # edge tables in the flat frame
    def F(G, M, N_):
        out = [[[Fr(0)] * K for d in range(5)] for c in range(5)]
        sg = (-1) ** (M + N_)
        for c in range(5):
            for d in range(5):
                acc = [Fr(0)] * K
                for a in range(5):
                    x = mul(ent(G, M, a, c), ent(G, N_, J3(a), d), K); acc = [s - sg * y / 5 for s, y in zip(acc, x)]
                if M == 0 and N_ == 0 and (c + d) % 5 == 3: acc[0] += Fr(1, 5)
                out[c][d] = acc
        return out
    def E(G, cache, m, n):
        if (m, n) not in cache:
            tot = [[[Fr(0)] * K for d in range(5)] for c in range(5)]
            for k in range(n + 1):
                f = F(G, m + 1 + k, n - k)
                tot = [[[s + (-1) ** k * y for s, y in zip(tot[c][d], f[c][d])] for d in range(5)] for c in range(5)]
            cache[(m, n)] = tot
        return cache[(m, n)]
    degF = lambda N: max(deg(n1) + deg(N - n1) for n1 in range(N + 1))
    Ncomp = [N for N in range(1, K) if degF(N - 1) <= NZ]
    for G, name in ((Gc, "canonical"), (Gs, "CONTROL stored normalisation")):
        cache = {}
        for signed in ((True, False) if name == "canonical" else (True,)):
            tested = 0; nz = 0; genuine = 0
            for n in range(NZ):
                for d in range(5):
                    acc = [Fr(0)] * K; terms = [0] * K
                    for m in range(0, NZ - n):
                        Em = E(G, cache, m, n); sg = (-1) ** m if signed else 1
                        for c in range(5):
                            x = mul(Ym[m][c], Em[c][d], K)
                            for N in range(K):
                                if x[N]: acc[N] += sg * x[N]; terms[N] += 1
                    for N in Ncomp:
                        tested += 1; nz += 1 if acc[N] else 0; genuine += 1 if terms[N] >= 2 else 0
            if name == "canonical" and signed:
                check("(E) canonical: Sum_m (-1)^m Ymod_m^T Ehat(m,n) = 0, all n < %d, all columns, u^N for N in %s (%d coefficients, %d sums of >= 2 terms)" % (NZ, Ncomp, tested, genuine), nz == 0, nz)
            else:
                check("(E) %s%s: the contraction does NOT vanish (%d of %d non-zero)" % (name, "" if signed else " CONTROL without (-1)^m", nz, tested), nz > 0)
    # Lemma L with random rational covectors (m <= 3)
    rnd = random.Random(20261113); cache = {}; bad = 0; tested = 0
    for trial in range(2):
        Yr = [[Fr(rnd.randint(-9, 9), rnd.randint(1, 9)) for c in range(5)] for m in range(4)]
        for n in range(NZ - 4):
            for d in range(5):
                lhs = [Fr(0)] * K
                for m in range(4):
                    Em = E(Gc, cache, m, n)
                    for c in range(5): lhs = [s + (-1) ** m * Yr[m][c] * y for s, y in zip(lhs, Em[c][d])]
                rhs = [Fr(0)] * K
                for r in range(1, n + 2):
                    for a in range(5):
                        h = [Fr(0)] * K
                        for m in range(4):
                            b = (a - r - m) % 5; h = [s + Yr[m][b] * y for s, y in zip(h, ent(Gc, r + m, a, b))]
                        x = mul(h, ent(Gc, n + 1 - r, J3(a), d), K); rhs = [s + (-1) ** (n + 1 - r) * y / 5 for s, y in zip(rhs, x)]
                tested += 1; bad += 0 if lhs == rhs else 1
    check("(L) Lemma L in the flat frame, 2 random covectors, %d series identities mod u^%d" % (tested, K), bad == 0, bad)
    # (T)
    Pi = [Fr(0)] * K
    for d in range(2, NZ + 2):
        for b in range(5):
            x = mul(ent(Gc, d - 1, 0, b), Ym[d][b], K); Pi = [s + y for s, y in zip(Pi, x)]
    vhat = [Fr(0)] + mul(powf(tau[:K], Fr(-3, 2), K), inv(S['varpi0'][:K], K), K)[:K - 1]
    Nt = [N for N in range(K) if (5 * (N - 1) + 4) // 4 + 1 <= NZ + 1]
    check("(T) Pi(T') = u tau^(-3/2)/varpi_0 - t_c through u^%d; leading coefficient %s" % (max(Nt), Pi[2]),
          all(Pi[N] == vhat[N] - S['tc'][N] for N in Nt) and Pi[2] == pc[1][4][0] / 2)
    Ps = [Fr(0)] * K
    for d in range(2, NZ + 2):
        for b in range(5):
            x = mul(ent(Gs, d - 1, 0, b), Ym[d][b], K); Ps = [s + y for s, y in zip(Ps, x)]
    dif = [N for N in Nt if Ps[N] != vhat[N] - S['tc'][N]]
    check("(T) CONTROL stored normalisation: Pi(T') differs from vhat - t_c (first at u^%s)" % (min(dif) if dif else None), bool(dif))
    log("RF_HAT %s: %s (%.0f s)" % (tag, "ALL PASS" if not FAIL else "FAIL %s" % FAIL, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
