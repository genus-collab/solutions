# -*- coding: utf-8 -*-
"""rf_direct.py -- REFEREE: independent DIRECT test (analogue of ev_direct / gi_actual_norm, sharing no code with them), mod p.
Actual side: all five colours from rf_core (WKB from the PF operator, QRR of colour i, rows Qh_a^(i)), model rows of colour i (own
recursion), Ghat^act = Qhat R^(tau zhat)^-1 with R^ = (F_b^(i)) inverted as a 5x5 matrix series (no symplecticity/covariance used).
Canonical side: rf_gauge (own Theorem-U solver), P_11(0, z) = 1, hatted; stored normalisation for the control.
Checks: D1 closing relation Qh_0 = (1/C0) A_0 Qh_4 (the PF equation as the step 4 -> 0 of the cyclic system with c = (C1, 1/(u C0^2
C1^2), C1, C0, C0)); D2 covariance of rows and model; D3 model symplecticity R^(-t)^T eta' R^(t) = 5; D4 Lemma P (colour-reduced
pairing) = matrix inverse; D5 grading; D6 Ghat^act = Ghat^can slot by slot; D7 normalisation slots; D8 stored control differs first at
zhat^5; D9 CONNECTION ROUTE (no WKB recursion, no K_q): ph^can_N := C0 [row 0, colour 0 of Ghat^can R^(tau zhat)] is a polynomial of
degree floor(4N/5) and ph^can_N(1) = h_N for every N <= NZ; it also equals the QRR WKB polynomial ph_N.
Run (cwd outside the project): py -I -B rf_direct.py PRIME alpha beta NZ K [NU]    log logs/rf_direct_<P>_<frame>_<NZ>_<K>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
import rf_gauge as G
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def matinv(F, M):
    n = len(M); A = [list(r) + [F.of(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c]); A[c], A[piv] = A[piv], A[c]; iv = F.inv(A[c][c])
        A[c] = [F.red(x * iv) for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]: f = A[r][c]; A[r] = [F.red(x - f * y) for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]


def main():
    a = sys.argv[1:]; P = a[0]; al = Fr(a[1]); be = Fr(a[2]); NZ = int(a[3]); K = int(a[4]); NU = int(a[5]) if len(a) > 5 else NZ + K + 10
    F = C.PF(PR[P]); z = F.zeta5(); fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "rf_direct_%s_%s_%d_%d.log" % (P, fr, NZ, K)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)
    t0 = time.time(); log("rf_direct (referee, independent): p = %s, frame (%s, %s), zhat <= %d, u < %d, gauge precision NU = %d; log logs/%s"
                          % (F.name, al, be, NZ, K, NU, name))
    C0, C1, tc, gam = C.periods(F, K, al, be); iC0 = C.inv(F, C0, K); iC1 = C.inv(F, C1, K)
    k2 = C.mul(F, C.mul(F, C0, C0, K), C.mul(F, C1, C1, K), K); tau = C.inv(F, C.itau(F, K), K)
    KH = [(-1, iC1), (-2, k2), (2, iC1), (1, iC0)]                    # Qh_(a+1) = kh_a A_(e_a) Qh_a
    Qh = [[None] * 5 for _ in range(5)]; Hok = True; close = True; ph0 = None; h0 = None
    for i in range(5):
        zi = pow(z, i, F.p); h = C.qrr(F, NZ, zi, zeta=z); W = C.Wkb(F, zi)
        for n in range(1, NZ + 1): W.step(C.qrr_rule(F, h))
        Hok = Hok and all(F.red(sum(W.ph(n)) - h[n]) == 0 for n in range(1, NZ + 1) if n % 5)
        if i == 0: ph0 = [W.ph(n) for n in range(NZ + 1)]; h0 = h
        Qh[0][i] = [C.mul(F, [F.red(zi * x) for x in W.ph(n)], iC0, K) for n in range(NZ + 1)]
        for j, (e, kh) in enumerate(KH):
            X = Qh[j][i]; Qh[j + 1][i] = [C.mul(F, C.A_at(F, e, zi, X[n], X[n - 1] if n else None, n, K), kh, K) for n in range(NZ + 1)]
        X = Qh[4][i]; close = close and all(C.mul(F, C.A_at(F, 0, zi, X[n], X[n - 1] if n else None, n, K), iC0, K) == Qh[0][i][n]
                                            for n in range(NZ + 1))
    log("   rows of the five colours in %.0f s" % (time.time() - t0))
    check("(H) for all five colours, 5 not | n <= %d (QRR weights zeta^i - zeta^j, -5 zeta^i)" % NZ, Hok)
    check("D1 closing relation Qh_0 = (1/C0) A_0 Qh_4, all colours, n <= %d, u < %d" % (NZ, K), close)
    cov = all(Qh[a_][i][n] == [F.red(pow(z, (i * (a_ + 1 - n)) % 5, F.p) * x) for x in Qh[a_][0][n]]
              for a_ in range(5) for i in range(5) for n in range(NZ + 1))
    Fb = [C.model(F, NZ, pow(z, i, F.p)) for i in range(5)]
    covm = all(Fb[i][b][m] == F.red(pow(z, (i * (b + 1 - m)) % 5, F.p) * Fb[0][b][m]) for i in range(5) for b in range(5) for m in range(NZ + 1))
    check("D2 covariance Qh_(a,n)^(i) = zeta^(i(a+1-n)) Qh_(a,n)^(0) and F_(b,m)^(i) = zeta^(i(b+1-m)) F_(b,m)", cov and covm)
    R = [[[Fb[i][b][m] for i in range(5)] for b in range(5)] for m in range(NZ + 1)]
    sym = all(F.red(sum((-1) ** l * R[l][b][i] * R[m - l][(3 - b) % 5][j] for l in range(m + 1) for b in range(5))
                    - (5 if (m == 0 and i == j) else 0)) == 0 for m in range(NZ + 1) for i in range(5) for j in range(5))
    check("D3 model symplecticity Sum_b R^_b^(i)(-t) R^_(3-b)^(j)(t) = 5 delta_ij, t^0..t^%d" % NZ, sym)
    Ri = [matinv(F, R[0])]
    for m in range(1, NZ + 1):
        S = [[F.red(sum(R[j][r][s] * Ri[m - j][s][c] for j in range(1, m + 1) for s in range(5))) for c in range(5)] for r in range(5)]
        Ri.append([[F.red(-sum(Ri[0][r][s] * S[s][c] for s in range(5))) for c in range(5)] for r in range(5)])
    tp = [[1] + [0] * (K - 1)]
    for m in range(1, NZ + 1): tp.append(C.mul(F, tp[-1], tau, K))
    Gact = []
    for N in range(NZ + 1):
        acc = [[[0] * K for _ in range(5)] for _ in range(5)]
        for n in range(N + 1):
            m = N - n
            for a_ in range(5):
                for c in range(5):
                    y = [F.red(sum(Qh[a_][i][n][j] * Ri[m][i][c] for i in range(5))) for j in range(K)]
                    if any(y): acc[a_][c] = [F.red(s + t) for s, t in zip(acc[a_][c], C.mul(F, y, tp[m], K))]
        Gact.append(acc)
    log("   Ghat^act (matrix inverse, five colours) in %.0f s" % (time.time() - t0))
    lp = all(Gact[N][a_][c] == ([0] * K if (N - a_ + c) % 5 else
             [F.red(s) for s in map(sum, zip(*[[x * (Fb[0][(3 - c) % 5][N - n] if (N - n) % 2 == 0 else -Fb[0][(3 - c) % 5][N - n])
                                                  for x in C.mul(F, Qh[a_][0][n], tp[N - n], K)] for n in range(N + 1)]))])
             for N in range(min(NZ, 30) + 1) for a_ in range(5) for c in range(5))
    check("D4 Lemma P: colour-reduced pairing = matrix-inverse Ghat^act, 25 entries, zhat <= %d, u < %d" % (min(NZ, 30), K), lp)
    check("D5 grading: [zhat^N] Ghat^act_(a,c) = 0 unless a = c + N mod 5", all(not any(Gact[N][a_][c]) for N in range(NZ + 1)
                                                                                for a_ in range(5) for c in range(5) if (a_ - c - N) % 5))
    res = G.gauge(F, NZ, NU, al, be, "can"); can = G.hatted(F, res, K)
    sto = G.hatted(F, G.gauge(F, NZ, NU, al, be, "stored"), K)
    log("   canonical and stored gauges (own Theorem-U solver) in %.0f s" % (time.time() - t0))
    bad = []; dst = []; tested = nz = 0
    for k in range(NZ + 1):
        for b in range(5):
            g = Gact[k][(b + k) % 5][b]
            for j in range(K):
                if j < G.mu(b, k):
                    if g[j]: bad.append((k, b, j))
                    continue
                tested += 1; nz += bool(g[j])
                if g[j] != can[(k, b)][j]: bad.append((k, b, j))
                if g[j] != sto[(k, b)][j]: dst.append((k, b, j))
    check("D6 Ghat^act = Ghat^can in all %d slots (k <= %d, mu(b,k) <= j < %d; %d non-zero), lower slots vanish" % (tested, NZ, K, nz), not bad, bad[:4])
    check("D8 control: stored normalisation differs in %d slots, first at k = %s" % (len(dst), min(dst)[0] if dst else None),
          bool(dst) and min(dst)[0] == 5)
    qs = [q for q in range(1, NZ // 5 + 1) if 4 * q < K]
    check("D7 normalisation slots [u^4q zhat^5q] Ghat^act_11 = Ghat^act_22 = 0 for q = 1..%d (c^act = 1 + O(z^%d))" % (max(qs), 5 * max(qs) + 5),
          all(Gact[5 * q][1][1][4 * q] == 0 and Gact[5 * q][2][2][4 * q] == 0 for q in qs))
    conn = []; degs = []
    for N in range(NZ + 1):
        s = [0] * K
        for k in range(N + 1):
            b = (-k) % 5; c = R[N - k][b][0]
            if c: s = [F.red(x + c * y) for x, y in zip(s, C.mul(F, can[(k, b)], tp[N - k], K))]
        phc = C.mul(F, C0, s, K); d = max([j for j in range(K) if phc[j]] or [-1]); degs.append(d)
        conn.append((N, d == (4 * N // 5 if N >= 3 else 0), F.red(sum(phc) - h0[N]) == 0, phc[:d + 1] == ph0[N] + [0] * (d + 1 - len(ph0[N]))))
    ok = [x for x in conn if 4 * x[0] // 5 + 1 < K]
    check("D9 connection route: ph^can_N polynomial of degree floor(4N/5), ph^can_N(1) = h_N and ph^can_N = ph_N(QRR), N <= %d" % max(x[0] for x in ok),
          all(x[1] and x[2] and x[3] for x in ok), [x[0] for x in ok if not (x[1] and x[2] and x[3])][:5])
    log("RF_DIRECT %s %s NZ=%d K=%d: %s   (%.0f s)" % (F.name, fr, NZ, K, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED,
                                                     time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
