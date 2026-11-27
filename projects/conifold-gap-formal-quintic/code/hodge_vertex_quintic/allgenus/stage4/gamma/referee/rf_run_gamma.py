# -*- coding: utf-8 -*-
"""rf_run_gamma.py -- REFEREE driver: identity (Gamma) c(z;L) = R(z;L), (H), the reduction to free constants, Ginf, the model
coefficients, the single-Gamma lemma, and failing controls.  Own code only (rf_ring, rf_wkb, rf_sp); no worker or project import.
Run (working directory outside the project):  py -I -B rf_run_gamma.py NQ NP [NR]
  NQ: exact order over Q;  NP: order mod p127 = 2^127-1 and p107 = 2^107-1 (primes NOT used by the worker);  NR: order of the
  random-L test (default 60).  Log: logs/rf_run_gamma_<NQ>_<NP>.log (relative to this folder).
"""
import os, sys, time, random
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_ring import FQ, FP, P127, mul, inv, exp0, add, bernoulli_poly, stirling_b, log1, powr
from rf_wkb import wkb, qrr_h, R_at, ph_poly
from rf_sp import fs_series, sp_generic, bexp

P107 = (1 << 107) - 1


def first_bad(F, c, R):
    """first n >= 1 with [z^n](c/R) != 0 (None if c/R == 1)"""
    N = len(c) - 1; D = mul(F, c, inv(F, R[:N + 1], N + 1), N + 1)
    assert D[0] == F.one
    for n in range(1, N + 1):
        if D[n] != F.zero: return n
    return None


def ginf_bern(N):
    lg = [Fr(0)] * (N + 1)
    for l in range(1, N // 5 + 1): lg[5 * l] = Fr((-1) ** (l + 1)) * bernoulli_poly(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1))
    return exp0(FQ(), lg, N + 1)


def ginf_direct(N, a=Fr(1, 5)):
    """5 x [v^n] of (1/v + a - 1/2) log(1 + a v) - a + Sum_k b_k v^(2k-1) (1 + a v)^(1-2k), kept for 5 | n, exponentiated"""
    F = FQ(); n = N + 2
    lg1 = log1(F, [Fr(1), a] + [Fr(0)] * (n - 2), n + 1)                # log(1 + a v)
    ser = [lg1[k + 1] for k in range(n)]                                # (1/v) log(1 + a v)
    ser = add(F, ser, [(a - Fr(1, 2)) * x for x in lg1[:n]]); ser[0] -= a
    for k in range(1, n // 2 + 1):
        m = 2 * k - 1
        for r in range(n - m):                                          # (1 + a v)^(-m) = Sum_r (-1)^r C(m+r-1, r) a^r v^r
            ser[m + r] += stirling_b(k) * (-1) ** r * comb(m + r - 1, r) * a ** r
    assert ser[0] == 0
    lg = [5 * ser[k] if k % 5 == 0 else Fr(0) for k in range(N + 1)]
    return exp0(F, lg, N + 1)


def model_f(N):
    """model F_0: normalised SP of exp(-(x^4/4 + 1/x)/z) dx at x = 1"""
    F = FQ(); nt = 2 * N + 6
    ph = [Fr(0)] * nt                                                   # -( (1+t)^4/4 + 1/(1+t) ) + 5/4
    for k in range(5): ph[k] -= Fr(comb(4, k), 4)
    for k in range(nt): ph[k] -= (-1) ** k
    ph[0] += Fr(5, 4)
    return sp_generic(F, ph, {0: [Fr(1)] + [Fr(0)] * (nt - 1)}, N)


def single_gamma(N, r):
    """normalised SP of Gamma_St(eta) Y^-eta eta^r d eta at eta = Y (eps = 1/Y): phase w log w - w, amplitude
    w^(r - 1/2) exp(Sum_k b_k eps^(2k-1) w^(1-2k)); the lemma predicts Sum_k S(r,k) (-eps)^(r-k)."""
    F = FQ(); nt = 2 * N + 6
    lg = log1(F, [Fr(1), Fr(1)] + [Fr(0)] * (nt - 2), nt)               # log(1+t)
    ph = add(F, mul(F, [Fr(1), Fr(1)], lg, nt), [Fr(0), Fr(-1)])         # (1+t) log(1+t) - t
    terms = {}
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m <= N: terms[m] = [stirling_b(k) * (-1) ** i * comb(m + i - 1, i) for i in range(nt)]
    G = bexp(F, terms, N, nt)
    pre = powr(F, [Fr(1), Fr(1)] + [Fr(0)] * (nt - 2), Fr(2 * r - 1, 2), nt)
    c = sp_generic(F, ph, {m: mul(F, pre, G[m], nt) for m in G}, N)
    S = [[1]]                                                           # Stirling numbers of the 2nd kind
    for a in range(1, r + 1): S.append([0] + [S[a - 1][k - 1] + (k * S[a - 1][k] if k < a else 0) for k in range(1, a + 1)])
    pred = [Fr(0)] * (N + 1)
    for k in range(0, r + 1):
        if r - k <= N and (k > 0 or r == 0): pred[r - k] += S[r][k] * (-1) ** (r - k)
    return c == pred, c[:4]


def main():
    NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 30; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    NR = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "rf_run_gamma_%d_%d.log" % (NQ, NP)), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("rf_run_gamma: NQ = %d (exact over Q), NP = %d (mod 2^127-1 and 2^107-1), NR = %d" % (NQ, NP, NR))
    Q = FQ(); h = qrr_h(Q, NQ); R, info = wkb(Q, NQ, h)
    check("(A1) transport operator: T_0 = 0 and T_1 = (L^6 - L) d/dL (n <= 12)", info['T01'])
    check("(A1) (H) h_auto: R_n(1) = h_n automatically for all 5 not| n <= %d" % NQ, all(info['h_auto'].values()))
    degs = [len(ph_poly(R, n)) - 1 for n in range(NQ + 1)]
    check("(A1) exponents = 4n mod 5 and R_n in Q[L] (min exponent >= 0), deg ph_n <= floor(4n/5)",
          min(info['minexp']) >= 0 and all(d <= 4 * n // 5 for n, d in enumerate(degs)), degs[:16])
    fm = model_f(NQ)
    check("(A2) top coefficients [L^(4n)] R_n = ph_n(0) = model F_0 coefficients (own model SP), n <= %d" % NQ,
          all(R[n].get(4 * n, 0) == fm[n] for n in range(NQ + 1)), fm[:4])
    x1 = info['free'][1]; m5 = Fr(-43397, 807975000)
    check("(A3) h_5 = -257071/8064000000, x_1 = 4537/3937500 = B_6(1/5)/6, x_1 + m5 = 4437977/4039875000",
          h[5] == Fr(-257071, 8064000000) and x1 == Fr(4537, 3937500) == bernoulli_poly(6, Fr(1, 5)) / 6
          and x1 + m5 == Fr(4437977, 4039875000))
    check("(A3) ph_5 = 47735919/128000000 - (203511/250000)u + (548547/1000000)u^2 - (6789/62500)u^3 + (4537/3937500)u^4",
          ph_poly(R, 5) == [Fr(47735919, 128000000), Fr(-203511, 250000), Fr(548547, 1000000), Fr(-6789, 62500), Fr(4537, 3937500)])
    check("(A3) x_2 = -177640111151/31007812500000, x_3 = -284328315069009259339/6226756347656250000000 (sibling EVID values)",
          info['free'][2] == Fr(-177640111151, 31007812500000) and
          (NQ < 15 or info['free'][3] == Fr(-284328315069009259339, 6226756347656250000000)))
    for L in (Fr(2), Fr(1, 2), Fr(7, 3), Fr(5, 4), Fr(2, 3), Fr(-3)):
        c = fs_series(Q, 1 / L, NQ); RL = R_at(Q, R, L)
        bad = [n for n in range(NQ + 1) if c[n] != RL[n]]
        check("(A4) (Gamma) exact: c_n(L) = R_n(L), n <= %d, L = %s" % (NQ, L), not bad, bad[:3])
    log("     c_5(2) = %s" % fs_series(Q, Fr(1, 2), 5)[5])
    gb = ginf_bern(NQ); gd = ginf_direct(NQ)
    check("(A5) Ginf: Bernoulli-polynomial closed form = direct re-expansion of Stirling for Gamma(y + 1/5), z^0..z^%d" % NQ, gb == gd)
    check("(A5) free constants x_m = [z^5m] Ginf, m <= %d (exact)" % (NQ // 5), all(info['free'][m] == gb[5 * m] for m in info['free']))
    for r in range(4):
        ok, c4 = single_gamma(NQ, r)
        check("(A6) single-Gamma lemma, r = %d: Stirling saddle series = (-Y d/dY)^r e^-Y exactly through eps^%d" % (r, NQ), ok, c4)
    c = fs_series(Q, Fr(1, 2), 12, drop_j0=True); fb = first_bad(Q, c, R_at(Q, R, 2)[:13])
    check("(A7) control: drop the j = 0 Stirling term (the factor Gamma((1-sigma)/z)) -> (Gamma) fails", fb is not None, "first at z^%s" % fb)
    c = fs_series(Q, Fr(1, 2), 12, bmod=(2, Fr(1))); fb = first_bad(Q, c, R_at(Q, R, 2)[:13])
    check("(A7) control: b_2 -> b_2 + 1 in every Stirling series -> fails first at z^3", fb == 3, "first at z^%s" % fb)
    log("     (A) done in %.1f s" % (time.time() - t0))
    gbP = ginf_bern(NP) if NP <= 400 else None
    for p in (P127, P107):
        F = FP(p); t1 = time.time(); hP = qrr_h(F, NP); RP, infoP = wkb(F, NP, hP, checks=False)
        check("[%s] (B1) (H) h_auto for all 5 not| n <= %d; min exponent >= 0" % (F.name, NP),
              all(infoP['h_auto'].values()) and min(infoP['minexp']) >= 0)
        if gbP is not None:
            check("[%s] (B4) free constants x_m = [z^5m] Ginf, m <= %d" % (F.name, NP // 5),
                  all(infoP['free'][m] == F.of(gbP[5 * m]) for m in infoP['free']))
        for L in ((Fr(2), Fr(3, 5), Fr(4, 3)) if p == P127 else (Fr(2), Fr(9, 7))):
            fb = first_bad(F, fs_series(F, 1 / L, NP), R_at(F, RP, L))
            check("[%s] (B2) (Gamma): D = c(z;L)/R(z;L) = 1 + O(z^%d), L = %s" % (F.name, NP + 1, L), fb is None, fb)
        log("     [%s] (B) main part %.1f s" % (F.name, time.time() - t1))
    F = FP(P127); RP, infoP = wkb(F, NR, qrr_h(F, NR), checks=False); rng = random.Random(20261114); bad = []
    for _ in range(6):
        L = rng.randrange(2, P127 - 1)
        if first_bad(F, fs_series(F, F.inv(L), NR), R_at(F, RP, L)) is not None: bad.append(L)
    check("[%s] (B3) (Gamma) at 6 random L in F_p (seeded), n <= %d" % (F.name, NR), not bad)
    for (m, d) in ((1, 1), (NR // 5, 3)):
        fr = dict(infoP['free']); fr[m] = F.red(fr[m] + d); Rc, _ = wkb(F, NR, qrr_h(F, NR), free=fr, checks=False)
        fb = first_bad(F, fs_series(F, Fr(1, 2), NR), R_at(F, Rc, 2))
        check("[%s] (B5) control: x_%d -> x_%d + %d makes D - 1 first non-zero at z^%d" % (F.name, m, m, d, 5 * m), fb == 5 * m, fb)
    log("RF_RUN_GAMMA %d %d: %s   (%.0f s)" % (NQ, NP, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
