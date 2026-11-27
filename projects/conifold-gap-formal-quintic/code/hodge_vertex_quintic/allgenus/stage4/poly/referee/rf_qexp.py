# -*- coding: utf-8 -*-
"""rf_qexp.py -- REFEREE: the QRR-normalised colour-0 WKB function computed at LARGE RADIUS (q-expansion), not from (HPF).

Independent route (no use of the worker's hatted recursion).  psi = e^(u_0/z) Sum_n z^n A_n(q), A_n in F_p[[q]], M = zD + L,
D = q d/dq, L = (1 - 3125 q)^(-1/5) (colour 0), and the Picard-Fuchs operator of the formal quintic at H = lambda_0 = 1:
    PF = (M^5 - 1) - q (3125 M^5 + 9375 z M^4 + 10625 z^2 M^3 + 5625 z^3 M^2 + 1370 z^4 M + 120 z^5)      (prod_k (5M + kz)).
z^(n+1) coefficient: T(A_n) + (terms in A_0..A_(n-1)) = 0 with the transport operator (referee's derivation, checked numerically
below):  T(A) = 5 L^4 [ (1 - 3125 q) D A - 625 q A ],  T(L) = 0.  q = 0 is a regular singular point with exponent 0, so A_n is fixed
by A_n(0); the QRR normalisation is A_n(0) = h_n := [z^n] h_{w_0}(z) for EVERY n (no free coefficient on this side).
Hatted form: ph_n(u) = u^((4n+1)/5) A_n = (1 - 3125 q)^((4n+1)/5) A_n.  Claims tested (POLY s.1-s.2, s.4):
  (T)  ph_n is a polynomial of degree <= floor(4n/5): the q-series terminates.  For 5 not| n this IS (H) at order n (otherwise a
       multiple of u^(4n/5), a non-terminating q-series, is present).
  (M)  ph_n(0) equals the model coefficient [t^n] F_0 (Lemma 1).
  (W)  f_m = [u^(4m)] ph_(5m) equal the worker's exact values (read from ../logs/pp_run_norm_Q14.log, data only).
Run: py -I -B rf_qexp.py N       (log: logs/rf_qexp_N<N>.log)
"""
import os, sys, time, re
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import P61, P89, md, Ser, hcoef_mod, model_rows

EP = (9375, 10625, 5625, 1370, 120)            # 5^(5-j) e_j(1..5), j = 1..5


def wkb_q(p, N, Kq=None, hmod=None):
    """A_n (n <= N) as q-series mod q^Kq; returns (S, A, L).  hmod: {n: shift of A_n(0)} (controls)"""
    if Kq is None: Kq = 4 * N // 5 + 12
    S = Ser(p, Kq); L = S.powr([1, (-3125) % p] + [0] * (Kq - 2), Fr(-1, 5))
    h = hcoef_mod(p, N)
    for n, d in (hmod or {}).items(): h[n] = (h[n] + d) % p
    def D(a): return [k * x % p for k, x in enumerate(a)]
    def Mstep(prev_lower, cur):                  # (M^j F)_m from (M^(j-1) F)_(m-1) and (M^(j-1) F)_m
        return S.add(D(prev_lower), S.mul(L, cur))
    zero = S.z()
    L4 = S.mul(S.mul(L, L), S.mul(L, L)); iL4 = S.inv(L4); i5 = pow(5, -1, p)
    # numerical check of the transport operator on a test series
    test = [(7 * k * k + 3) % p for k in range(Kq)]
    def full_T(A):                               # linear part at z^(n+1) of PF(z^n A), computed from the definitions
        Mn = [A]; Mn1 = [zero]
        for j in range(1, 6):
            Mn1.append(S.add(D(Mn[j - 1]), S.mul(L, Mn1[j - 1]))); Mn.append(S.mul(L, Mn[j - 1]))
        r = S.add(Mn1[5], S.shift(Mn1[5], 1), (-3125) % p)
        return S.add(r, S.shift(S.scale(Mn[4], EP[0]), 1), p - 1)
    T_formula = S.scale(S.mul(L4, S.add(S.add(D(test), S.shift(D(test), 1), (-3125) % p), S.shift(test, 1), (-625) % p)), 5)
    assert full_T(test) == T_formula, "transport operator"
    assert not any(full_T(L)), "T(L) != 0"
    MJ = [[] for _ in range(6)]                   # MJ[j][m] = (M^j F)_m, true values
    A = []
    def get(j, m): return MJ[j][m] if m >= 0 else zero
    for n in range(N + 1):
        if n == 0:
            An = L
        else:
            Mn = [zero]; Mn1 = [zero]
            for j in range(1, 6):
                Mn.append(Mstep(get(j - 1, n - 1), Mn[j - 1]))
                Mn1.append(Mstep(Mn[j - 1], Mn1[j - 1]))
            rest = S.add(Mn1[5], S.shift(Mn1[5], 1), (-3125) % p)
            low = S.scale(Mn[4], EP[0])
            for jj, mm in ((2, 1), (3, 2), (4, 3), (5, 4)):
                low = S.add(low, S.scale(get(5 - jj, n - mm), EP[jj - 1]))
            rest = S.add(rest, S.shift(low, 1), p - 1)
            rho = S.scale(S.mul(rest, iL4), (-i5) % p)            # (1-3125q) D A - 625 q A = rho
            assert rho[0] == 0, "obstruction at q^0, n = %d" % n
            An = [h[n]] + [0] * (Kq - 1)
            for k in range(1, Kq):
                An[k] = (rho[k] + (3125 * k - 2500) * An[k - 1]) * pow(k, -1, p) % p
        A.append(An)
        MJ[0].append(An)
        for j in range(1, 6): MJ[j].append(Mstep(get(j - 1, n - 1), MJ[j - 1][n]))
    return S, A, L


def hatted_polys(p, N, Kq=None, hmod=None):
    """ph_n(u) mod p as u-polynomials (lists), plus the termination flags; ph_n via q-series"""
    S, A, L = wkb_q(p, N, Kq, hmod); Kq = S.K
    one_m = [1, (-3125) % p] + [0] * (Kq - 2)
    out = []; term = []
    i3125 = pow(3125, -1, p)
    for n in range(N + 1):
        ph_q = S.mul(S.powr(one_m, Fr(4 * n + 1, 5)), A[n])
        d = 4 * n // 5
        term.append(all(x == 0 for x in ph_q[d + 1:]))
        poly = [0] * (d + 1)                      # q = (1-u)/3125
        for k in range(d + 1):
            c = ph_q[k] * pow(i3125, k, p) % p
            for j in range(k + 1):
                poly[j] = (poly[j] + c * comb(k, j) * (-1) ** j) % p
        while len(poly) > 1 and poly[-1] == 0: poly.pop()
        out.append(poly)
    return out, term


def worker_free_values(path):
    """f_q from the worker's log (data, not code)"""
    vals = {}
    if not os.path.exists(path): return vals
    for line in open(path, encoding='utf-8'):
        m = re.search(r"f_(\d+) \(free coefficient .*\) = (-?\d+)/(\d+)", line)
        if m: vals[int(m.group(1))] = Fr(int(m.group(2)), int(m.group(3)))
    return vals


def main(N):
    t0 = time.time()
    print("rf_qexp: colour-0 WKB at large radius (q-expansion), n <= %d, primes 2^61-1, 2^89-1" % N)
    wf = worker_free_values(os.path.join(HERE, "..", "logs", "pp_run_norm_Q14.log"))
    print("worker's f_q read from ../logs/pp_run_norm_Q14.log: q = %s" % sorted(wf))
    allok = True; res = {}
    for p in (P61, P89):
        polys, term = hatted_polys(p, N)
        F = model_rows(p, N)[0]
        degs = [len(x) - 1 for x in polys]
        ok_T = all(term); ok_deg = all(degs[n] <= 4 * n // 5 for n in range(N + 1))
        ok_M = all(polys[n][0] == F[n] for n in range(N + 1))
        eq = [n for n in range(N + 1) if degs[n] == 4 * n // 5]
        fm = {m: (polys[5 * m][4 * m] if len(polys[5 * m]) > 4 * m else 0) for m in range(1, N // 5 + 1)}
        ok_W = all(fm[m] == md(wf[m], p) for m in fm if m in wf)
        nW = sum(1 for m in fm if m in wf)
        print("p = %d:" % p)
        print("  (T) q-series of ph_n terminate at degree floor(4n/5) for all n <= %d: %s  [(H) for the %d values 5 not| n]"
              % (N, ok_T, sum(1 for n in range(1, N + 1) if n % 5)))
        print("  degrees: deg ph_n <= floor(4n/5): %s; equality fails exactly for n in %s" % (ok_deg, [n for n in range(N + 1) if n not in eq]))
        print("  (M) ph_n(0) = model [t^n] F_0 for n <= %d: %s" % (N, ok_M))
        print("  (W) f_m = worker's exact values (mod p) for %d values m = 1..%d: %s" % (nW, max(fm) if fm else 0, ok_W))
        allok = allok and ok_T and ok_deg and ok_M and ok_W and nW > 0
        res[p] = polys
    # exact spot checks (POLY s.1): ph_5, h_5, f_1
    ph5 = [Fr(47735919, 128000000), Fr(-203511, 250000), Fr(548547, 1000000), Fr(-6789, 62500), Fr(4537, 3937500)]
    s1 = all(res[p][5] == [md(c, p) for c in ph5] for p in res)
    print("ph_5 = 47735919/128000000 - (203511/250000)u + (548547/1000000)u^2 - (6789/62500)u^3 + (4537/3937500)u^4 (POLY s.1), both primes: %s" % s1)
    allok = allok and s1
    # failing controls: a wrong constant at q = 0 must destroy termination (H) at the predicted order
    for hm, first in (({3: 1}, 3), ({5: 1}, 6), ({10: 7}, 11)):
        _, term = hatted_polys(P61, min(N, 15), hmod=hm)
        bad = [n for n, t in enumerate(term) if not t]
        ok = bool(bad) and bad[0] == first
        print("CONTROL A_n(0) shifted by %s: termination first fails at n = %s (predicted %d): %s"
              % (hm, bad[0] if bad else None, first, "fails as predicted" if ok else "UNEXPECTED"))
        allok = allok and ok
    print("RESULT rf_qexp n <= %d: %s  [%.1f s]" % (N, "PASS" if allok else "FAIL", time.time() - t0))
    return res


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
