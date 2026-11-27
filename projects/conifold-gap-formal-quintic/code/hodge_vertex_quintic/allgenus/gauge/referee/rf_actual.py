# -*- coding: utf-8 -*-
"""rf_actual.py -- REFEREE: is the gauge of the ACTUAL R-matrix the canonical gauge?  Two routes different from gi_actual_norm.py:
 (A) EXACT over Q: the referee table natural_form/referee/G_fourier_NZ16_K12.json (actual G = Rt M_0(tau zhat)^-1, Fourier components,
     rationals reconstructed at 2^521-1 and checked at two more primes; zhat <= 16, u <= 11, frame (0,0)) against the canonical gauge
     (project ur_gauge, certified by rf_gi / rf_frob_run), hatted and Fourier-transformed HERE (dictionary derived in REFEREE.md s.2).
 (B) mod p: rows natural_form/referee/ref_hat.HatRows (zhat <= NZ, u < K) with M_0 taken from the LG MODEL recursion (LG_MODEL s.3; no
     quintic input), Ghat^act = Qhat R^(tau zhat)^-1, compared slot by slot with the canonical hatted gauge; control: stored normalisation.
Run: py -I -B rf_actual.py al be NZ K     (part A only for frame (0,0))      log: rf_actual_<frame>.log
"""
import os, sys, time, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import powf, mul, CEXP, EPS
from rf_ur import stored_gauge, p11_zero, normalise, PROJ
sys.path.insert(0, os.path.join(PROJ, "natural_form", "referee"))
from ref_hat import HatRows
P89 = (1 << 89) - 1; P61 = (1 << 61) - 1


def mu(b, k): return (4 * k + CEXP[b] - CEXP[(b + k) % 5]) // 5


def lg_rows(NZ):
    """g^a_n (LG_MODEL s.3 closed recursion), exact over Q"""
    g = [[Fr(1)] for a in range(5)]
    for n in range(1, NZ + 1):
        d = [Fr(0)]
        for a in range(4): d.append(d[a] + Fr(4 * n - 3 - CEXP[a], 5) * g[a][n - 1])
        t = -sum((4 * n + 1 - CEXP[a]) * d[a] for a in range(5)) / (20 * n)
        for a in range(5): g[a].append(t + d[a])
    return g


def hatted_can(pc, tau, NZ, K):
    """{(k, b): [Fraction]*K} = [u^0..u^(K-1)] of [zhat^k] Ghat_{b+k, b}"""
    out = {}
    for b in range(5):
        f = powf(tau[:K], Fr(1 - CEXP[b], 4), K)
        for k in range(NZ + 1):
            m = mu(b, k)
            s = mul(f, pc[k][b][:K], K) if m < K else []
            out[(k, b)] = [Fr(0)] * m + s[:K - m] if m < K else [Fr(0)] * K
    return out


def main():
    al = Fr(sys.argv[1]); be = Fr(sys.argv[2]); NZ = int(sys.argv[3]); K = int(sys.argv[4])
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "rf_actual_%s.log" % tag), "w", encoding="ascii"); FAIL = []
    def log(s): print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(name, ok, det=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(name)
    t0 = time.time()
    p, prec = stored_gauge(NZ, 3 * K + 24, al, be); c = p11_zero(p); pc, ci = normalise(p, c)
    assert prec >= K + 1, prec
    tau = [Fr(1)]; i4 = [Fr(4, 5 * n + 4) for n in range(K + 2)]
    for n in range(1, K + 2): tau.append(-sum(i4[j] * tau[n - j] for j in range(1, n + 1)))
    can = hatted_can(pc, tau, NZ, K); sto = hatted_can(p, tau, NZ, K)
    log("rf_actual: frame (%s, %s), zhat <= %d, u < %d; canonical gauge from the periods (precision O(u^%d)), %.0f s" % (al, be, NZ, K, prec, time.time() - t0))
    if al == 0 and be == 0:
        J = json.load(open(os.path.join(PROJ, "natural_form", "referee", "G_fourier_NZ16_K12.json")))
        F = {}
        for key, v in J['G_fourier'].items():
            u_, z_, mm = key.split(" "); m, mp = mm.strip("()").split(",")
            F[(int(u_[1:]), int(z_[1:]), int(m), int(mp))] = Fr(v)
        Kj = J['K']; NZj = J['NZ']
        g0 = [hatted_can(pc, tau, 0, Kj)[(0, a)] for a in range(5)]
        g0i = []
        for a in range(5):                                     # 1/g_a(u) as a power series
            s = g0[a]; o = [1 / s[0]] + [Fr(0)] * (Kj - 1)
            for n in range(1, Kj): o[n] = -sum(s[j] * o[n - j] for j in range(1, n + 1)) / s[0]
            g0i.append(o)
        mine = {}
        for k in range(NZj + 1):
            for b in range(5):
                a = (b + k) % 5; gt = mul(g0i[a], can[(k, b)][:Kj], Kj)
                m_, mp_ = (a + 1) % 5, (-b - 1) % 5
                for n in range(Kj):
                    x = Fr(EPS[a] * EPS[b], 5) * gt[n]
                    if x: mine[(n, k, m_, mp_)] = x
        bad = [key for key in set(F) | set(mine) if F.get(key, 0) != mine.get(key, 0)]
        log("   (A) referee table: %d non-zero Fourier components; canonical gauge: %d; slot u4 z5 (2,3) [= p_{5,1}(0)/5]: table %s, canonical %s"
            % (len(F), len(mine), F.get((4, 5, 2, 3), 0), mine.get((4, 5, 2, 3), 0)))
        check("(A) EXACT: actual G (referee table, Q) = canonical gauge in every Fourier component, zhat <= %d, u <= %d" % (NZj, Kj - 1), not bad, sorted(bad)[:5])
        stmine = {}
        for k in range(NZj + 1):
            for b in range(5):
                a = (b + k) % 5; gt = mul(g0i[a], sto[(k, b)][:Kj], Kj)
                for n in range(Kj):
                    x = Fr(EPS[a] * EPS[b], 5) * gt[n]
                    if x: stmine[(n, k, (a + 1) % 5, (-b - 1) % 5)] = x
        badst = [key for key in set(F) | set(stmine) if F.get(key, 0) != stmine.get(key, 0)]
        check("(A) control: the STORED normalisation differs from the table (%d components), first at %s" % (len(badst), min(badst) if badst else None), len(badst) > 0)
    g = lg_rows(NZ)
    for pr, name in ((P89, "2^89-1"), (P61, "2^61-1")):
        t1 = time.time(); H = HatRows(pr, NZ, K, al, be); z = H.zeta; of = lambda x: x.numerator * pow(x.denominator, pr - 2, pr) % pr
        Rh = [[[of(Fr(EPS[a]) * g[a][n]) * pow(z, (i * (a + 1 - n)) % 5, pr) % pr for i in range(5)] for a in range(5)] for n in range(NZ + 1)]
        check("[%s] model rows (LG recursion) = referee rows at u^0, all a, i, zhat^0..zhat^%d; referee self-checks %s" % (name, NZ, all(H.checks.values())),
              all(H.Q[a][i][n][0] % pr == Rh[n][a][i] for a in range(5) for i in range(5) for n in range(NZ + 1)) and all(H.checks.values()))
        V0 = [[EPS[b] * pow(z, (-i * (b + 1)) % 5, pr) * pow(5, pr - 2, pr) % pr for b in range(5)] for i in range(5)]
        mm = lambda X, Y: [[sum(X[i][j] * Y[j][l] for j in range(5)) % pr for l in range(5)] for i in range(5)]
        V = [V0]
        for n in range(1, NZ + 1):
            acc = [[0] * 5 for _ in range(5)]
            for m in range(1, n + 1):
                t = mm(Rh[m], V[n - m]); acc = [[(acc[i][j] + t[i][j]) % pr for j in range(5)] for i in range(5)]
            V.append([[(-x) % pr for x in row] for row in mm(V0, acc)])
        taup = [[1] + [0] * (K - 1)]; tq = [of(x) for x in tau[:K]]
        for n in range(1, NZ + 1): taup.append([sum(taup[-1][i] * tq[N - i] for i in range(N + 1)) % pr for N in range(K)])
        tested = 0; bad = []; diff_st = []; offgrade = 0
        for k in range(NZ + 1):
            Gk = [[[0] * K for _ in range(5)] for _ in range(5)]
            for n in range(k + 1):
                for a in range(5):
                    for b in range(5):
                        coef = [0] * K
                        for i in range(5):
                            if V[n][i][b]:
                                q = H.Q[a][i][k - n]
                                for N in range(K): coef[N] = (coef[N] + q[N] * V[n][i][b]) % pr
                        for N in range(K): Gk[a][b][N] = (Gk[a][b][N] + sum(coef[j] * taup[n][N - j] for j in range(N + 1))) % pr
            for a in range(5):
                for b in range(5):
                    if (a - b - k) % 5:
                        if any(Gk[a][b]): offgrade += 1
                        continue
                    for N in range(K):
                        tested += 1; x = Gk[a][b][N]
                        if x != of(can[(k, b)][N]): bad.append((k, a, b, N))
                        if x != of(sto[(k, b)][N]): diff_st.append((k, a, b, N))
            if k % 5 == 0 and k:
                log("   [%s] zhat^%d, u^%d: actual diagonal (b = 0..4) %s; canonical %s" % (name, k, 4 * k // 5,
                    [Gk[b][b][4 * k // 5] for b in range(5)] if 4 * k // 5 < K else "-", [str(pc[k][b][0]) for b in range(5)]))
        check("[%s] Ghat^act (M_0 from the model) is graded (%d off-grade entries non-zero) and = canonical in all %d slots zhat^k u^N (k <= %d, N < %d)"
              % (name, offgrade, tested, NZ, K), not bad and offgrade == 0, bad[:4])
        check("[%s] control: the stored normalisation differs in %d slots, first %s" % (name, len(diff_st), min(diff_st) if diff_st else None), len(diff_st) > 0)
        log("   [%s] done in %.0f s" % (name, time.time() - t1))
    log("RF_ACTUAL %s: %s (%.0f s)" % (tag, "ALL PASS" if not FAIL else "FAIL %s" % FAIL, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
