# -*- coding: utf-8 -*-
"""rf_lemmaK.py -- Lemma K in the FLAT frame, exact over Q (a different route from the worker's hatted mod-p check).

Covector components Y_a = eta(s, H^a) of flat sections; the gauge acts by Y^Q = P Y^mod (ur_gauge.gauge, project module, read only;
p_{k,b} = entry (b+k, b) of P_k).  On vectors it acts by R(z) = eta^-1 P(z) eta, so R^-1(w)^T = P(-w) when P(-z)^T eta^-1 P(z) = eta^-1.
Edge V(w,z) = (eta^-1 - P(-w)^T eta^-1 P(-z))/(w+z);  leaf point L(w) = -w R^-1(w) e_0 = -w P(-w)^T e_0;
covector phi(w) = Y^C(v(u), -w), Y^C_b = SGN_b I_(MEXP_b), I_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)! (own transcription).
Pi_m(f) = [w^m] Sum_a phi_a(w) f^a(w).  C_(m,l) := Pi_m of the first slot of [z^l] V.  Reliable region: P known to z^K, u^NU;
[u^N] C_(m,l) needs phi_j for j <= jmax(N) (val_u phi_j >= j - floor(j/5)) and P to z^(m+l+1+jmax).
Modes: canon (P_11(0,z) = 1), stored (ur_gauge as computed), c5 (canon times exp(z^5/7)).  Run: py -I -B rf_lemmaK.py K NU mode
"""
import os, sys, time, random
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory"))
import ur_gauge
MEXP = (0, -1, 3, 2, 1); SGN = (1, -1, -1, 1, 1)


def pmul(a, b, NU):
    out = [Fr(0)] * NU
    for i in range(NU):
        if a[i]:
            for j in range(NU - i):
                if b[j]: out[i + j] += a[i] * b[j]
    return out


def jmax(N):
    j = 0
    while (j + 1) - (j + 1) // 5 <= N: j += 1
    return j


def main():
    K = int(sys.argv[1]); NU = int(sys.argv[2]); mode = sys.argv[3]
    t0 = time.time(); res = ur_gauge.gauge(K, NU + 2 * K + 4)
    ser = lambda s: [Fr(s.coef(n)) for n in range(NU)]
    v = ser(res['v']); tc = ser(res['C'].tc) if hasattr(res['C'], 'tc') else None
    raw = {k: [ser(res['p'][k][b]) for b in range(5)] for k in range(K + 1)}
    c11 = [raw[k][1][0] if k % 5 == 0 else Fr(0) for k in range(K + 1)]
    cz = [Fr(1)] + [Fr(0)] * K
    if mode in ("canon", "c5"):
        for k in range(1, K + 1): cz[k] = -sum(c11[j] * cz[k - j] for j in range(1, k + 1)) / c11[0]
    if mode == "c5":
        e5 = [Fr(0)] * (K + 1)
        for j in range(K // 5 + 1): e5[5 * j] = Fr(1, 7) ** j / factorial(j)
        cz = [sum(cz[i] * e5[k - i] for i in range(k + 1)) for k in range(K + 1)]
    P = {}
    for k in range(K + 1):
        for b in range(5):
            P[(k, (b + k) % 5, b)] = [sum(cz[k - i] * raw[i][b][n] for i in range(k + 1)) for n in range(NU)]
    zero = [Fr(0)] * NU
    E = lambda k, a, b: P.get((k, a, b), zero)
    print("mode %s: stored P_11(0,z) at z^0,5,10,..: %s ; now %s" % (mode, [str(c11[k]) for k in range(0, K + 1, 5)],
          [str(E(k, 1, 1)[0]) for k in range(0, K + 1, 5)]))
    einv = {(a, (3 - a) % 5): Fr(1, 5) * SGN[a] * SGN[(3 - a) % 5] for a in range(5)}
    # symplecticity (orders z^0..z^K)
    sb = None
    for k in range(0, K + 1):
        for a in range(5):
            for b in range(5):
                acc = list(zero)
                for i in range(k + 1):
                    for (c, d), x in einv.items():
                        y = pmul(E(i, c, a), E(k - i, d, b), NU)
                        acc = [s + (-1) ** i * x * t for s, t in zip(acc, y)]
                if k == 0: acc[0] -= einv.get((a, b), 0)
                if any(acc) and sb is None: sb = k
    print("symplecticity through z^%d: %s" % (K, "OK" if sb is None else "FAILS from z^%d" % sb))
    # Y^C coefficients phi: Ycf[j][b] = [z^-j] Y^C_b (u-series); phi_{b,j} = (-1)^j Ycf[j][b]
    J = jmax(NU - 1); vp = {0: [Fr(1)] + [Fr(0)] * (NU - 1)}
    for n in range(1, J + 2): vp[n] = pmul(vp[n - 1], v, NU)
    Ycf = {j: [list(zero) for b in range(5)] for j in range(J + 1)}
    for b in range(5):
        i = 0
        while MEXP[b] + 5 * i + 1 <= J:
            j = MEXP[b] + 5 * i + 1; pw = MEXP[b] + 4 * i + 1
            c = SGN[b] * Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** pw, factorial(pw))
            Ycf[j][b] = [c * x for x in vp[pw]]; i += 1
    # model ODE: z d/dv Y_b = bmod_b Y_(b+1): on monomials c v^pw z^-j: pw c v^(pw-1) z^(1-j) ... checked via the closed form at v-level
    # (G): [z^r] (P Y^C)_a = Sum_j P_(r+j) Y_j, r >= 1, reliable for u^N with r + jmax(N) <= K
    gb = gt = 0
    for r in range(1, K + 1):
        for a in range(5):
            acc = list(zero)
            for j in range(0, J + 1):                    # j = 0 carries Y_1 = -1 (first run omitted it: 129 spurious non-zeros)
                if r + j > K: break
                for b in range(5):
                    if (r + j, a, b) in P: acc = [s + t for s, t in zip(acc, pmul(P[(r + j, a, b)], Ycf[j][b], NU))]
            for N in range(NU):
                if r + jmax(N) <= K: gt += 1; gb += acc[N] != 0
    print("(G): %d coefficients [u^N z^r](P Y^C)_a, r >= 1, in the reliable region; %d non-zero" % (gt, gb))
    # edge: (w+z) V = einv - P(-w)^T einv P(-z);  G(k,l)[a][b] = -(-1)^(k+l) Sum_(c,d) P_k[c][a] einv[c][d] P_l[d][b]  ((k,l) != (0,0))
    def Gm(k, l, a, b):
        acc = list(zero)
        for (c, d), x in einv.items():
            y = pmul(E(k, c, a), E(l, d, b), NU); acc = [s - (-1) ** (k + l) * x * t for s, t in zip(acc, y)]
        if k == 0 and l == 0:
            acc = [s + (einv.get((a, b), 0) if n == 0 else 0) for n, s in enumerate(acc)]
        return acc
    Vc = {}
    def V(k, l, a, b):
        key = (k, l, a, b)
        if key not in Vc:
            acc = list(zero)
            for t in range(l + 1):
                g = Gm(k + 1 + t, l - t, a, b); acc = [s + (-1) ** t * x for s, x in zip(acc, g)]
            Vc[key] = acc
        return Vc[key]
    rnd = random.Random(5); ctrl = {(j, b): Fr(rnd.randrange(1, 97)) for j in range(J + 1) for b in range(5)}
    Mmax = 6; tot = bad = cbad = gen = 0
    for m in range(Mmax + 1):
        for l in range(0, K):
            Ns = [N for N in range(NU) if m + l + 1 + jmax(N) <= K]
            if not Ns: continue
            jm = jmax(max(Ns))
            for b in range(5):
                s = list(zero); sc = list(zero); cnt = [0] * NU
                for j in range(0, jm + 1):
                    for a in range(5):
                        if not any(Ycf[j][a]): continue
                        y = pmul(Ycf[j][a], V(m + j, l, a, b), NU)
                        if any(y):
                            s = [p_ + (-1) ** j * q_ for p_, q_ in zip(s, y)]
                            cnt = [c_ + (q_ != 0) for c_, q_ in zip(cnt, y)]
                        sc = [p_ + ctrl[(j, a)] * q_ for p_, q_ in zip(sc, V(m + j, l, a, b))]
                for N in Ns:
                    tot += 1; bad += s[N] != 0; cbad += sc[N] != 0; gen += cnt[N] >= 2
    print("edge: %d coefficients [u^N] Pi_m([z^l] V) (m <= %d, reliable region); %d non-zero; %d sums of >= 2 terms; control covector non-zero in %d"
          % (tot, Mmax, bad, gen, cbad))
    # leaf: Pi_m(L) = -[w^(m-1)] Sum_a phi_a(w) P(-w)_(0,a) = -[w^(m-1)] (P(-w) phi(w))_0
    for m in range(0, Mmax + 1):
        acc = list(zero); ok_all = True
        for N in range(NU):
            if m + jmax(N) > K: ok_all = False; continue
        for j in range(0, J + 1):
            k = m - 1 + j
            if k < 0 or k > K: continue
            for a in range(5):
                if (k, 0, a) in P:
                    y = pmul(P[(k, 0, a)], Ycf[j][a], NU); acc = [s - (-1) ** k * (-1) ** j * t for s, t in zip(acc, y)]
        rel = [N for N in range(NU) if m + jmax(N) <= K]
        print("leaf Pi_%d(L) mod u^%d: %s%s" % (m, max(rel) + 1, [str(acc[N]) for N in rel][:8],
              ("   -t_c = %s" % [str(-x) for x in tc[:8]]) if (m == 0 and tc) else ""))
    print("rf_lemmaK %s K=%d NU=%d done (%.0f s)" % (mode, K, NU, time.time() - t0))


if __name__ == "__main__":
    main()
