# -*- coding: utf-8 -*-
"""lg_gm.py -- Gauss-Manin relations, the model connection and the closed recursion for the conifold symbol (exact, over Q).

Conventions (W-form of the model; see LG_MODEL.md s.2):  W(xi) = xi^4/4 + u/xi,  I_k(u, z) := Int exp(-W/z) xi^k dxi,
   critical points xi^5 = u (xi_c = zeta^{-i} s, s = u^{1/5}), W(xi_c) = (5/4) xi_c^4, W'' = 5 xi_c^2, zhat = z/s^4;
   I_k ~ exp(-5/(4 zhat)) (2 pi z/5)^{1/2} s^{k-1} g_k(zhat),   g_k(zhat) := J_k(-4 zhat)   (J_k of lg_sp: lambda = x^4 + u/x, eps = -4 zhat).
Relations (derived in LG_MODEL.md from d(exp(-W/z) xi^m) and d/du under the integral; CHECKED here on the stationary-phase series):
   (R1)  z d/du I_k = -I_{k-1}            <=>   g_{k-1} = g_k - (zhat/5)(k - 1 - 4 theta) g_k,      theta = zhat d/dzhat
   (R2)  I_{m+3} = u I_{m-2} + m z I_{m-1} <=>   g_{k+5} = g_k + (k+2) zhat g_{k+1}                   (k = m - 2)
   (R3)  [u (z d/du)^5 - (k-3) z (z d/du)^4 + 1] I_k = 0
                                          <=>   zhat^5 T (T-1)(T-2)(T-3)(T-k-1) g_k + g_k = 0,  T = -1/zhat + (k-1)/5 - (4/5) theta
Flat frame a = 0..4 <-> k = c_a = (0, -1, 3, 2, 1) (g^a := g_{c_a}):   the MODEL CONNECTION (cyclic in a, g^5 = g^0)
   (MC)  (4/5) zhat theta g^a = g^{a+1} - g^a + (zhat/5)(c_a - 1) g^a     i.e.   g^{a+1}_n - g^a_n = ((4n - 3 - c_a)/5) g^a_{n-1},
   closed by Sum_a (4n + 1 - c_a) g^a_n = 0  (n >= 1), which fixes the common constant of level n:  symbol_rows(N) below.
Run:  py -I lg_gm.py     (prints; writes lg_symbol_rows_Q.json = the rows g^a_n, n <= 40, as exact rationals)
"""
import os, sys, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lg_sp import J_series

CEXP = (0, -1, 3, 2, 1)
NSIGN = (1, -1, -1, 1, 1)


def symbol_rows(N):
    """g[a][n] (a = 0..4, n = 0..N) from the cyclic recursion (MC) alone (no integral, no quintic)"""
    g = [[Fr(1)] for a in range(5)]
    for n in range(1, N + 1):
        d = [Fr(0)]
        for a in range(4): d.append(d[a] + Fr(4 * n - 3 - CEXP[a], 5) * g[a][n - 1])
        # closure of the cycle: d_5 = d_4 + ((4n-3-c_4)/5) g^4_{n-1} must be 0  (it is, by the constraint at level n-1)
        close = d[4] + Fr(4 * n - 3 - CEXP[4], 5) * g[4][n - 1]
        assert close == 0, (n, close)
        t = -sum((4 * n + 1 - CEXP[a]) * d[a] for a in range(5)) / (20 * n)
        for a in range(5): g[a].append(t + d[a])
    return g


# ---------------------------------------------------------------- polynomial helpers (lists of Fractions, variable m)
def padd(A, B):
    n = max(len(A), len(B)); out = [Fr(0)] * n
    for i, x in enumerate(A): out[i] += x
    for i, x in enumerate(B): out[i] += x
    while out and out[-1] == 0: out.pop()
    return out
def pmul(A, B):
    if not A or not B: return []
    out = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B): out[i + j] += x * y
    return out
def peval(A, x): return sum((c * Fr(x) ** i for i, c in enumerate(A)), Fr(0))
def pfmt(A): return " + ".join("(%s) m^%d" % (c, i) for i, c in enumerate(A) if c != 0) or "0"


def scalar_operator(k):
    """zhat^5 T(T-1)(T-2)(T-3)(T-k-1) as {r: P_r(m)}: zhat^m -> Sum_r P_r(m) zhat^{m+5-r}, T zhat^m = -zhat^{m-1} + ((k-1-4m)/5) zhat^m"""
    op = {0: [Fr(1)]}
    for c in (0, 1, 2, 3, k + 1):
        new = {}
        for r, P in op.items():
            new[r + 1] = padd(new.get(r + 1, []), [-x for x in P])
            # diagonal part at the current exponent m - r:  (k - 1 - 4(m - r))/5 - c
            new[r] = padd(new.get(r, []), pmul(P, [Fr(k - 1 + 4 * r, 5) - c, Fr(-4, 5)]))
        op = new
    return op


def main():
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    NE = 10; t0 = time.time()
    KS = list(range(-4, 9))
    g = {k: [x * (-4) ** n for n, x in enumerate(J_series(k, NE))] for k in KS}        # g_k(zhat) = J_k(-4 zhat)
    print("stationary-phase series g_k(zhat) = J_k(-4 zhat), k = -4..8, through zhat^%d (%.1f s)" % (NE, time.time() - t0))
    # (R1)
    bad = []
    for k in KS[1:]:
        for n in range(NE + 1):
            rhs = g[k][n] - (Fr(k - 1 - 4 * (n - 1), 5) * g[k][n - 1] if n >= 1 else 0)
            if rhs != g[k - 1][n]: bad.append((k, n))
    check("(R1) g_{k-1} = g_k - (zhat/5)(k-1-4 theta) g_k, k = -3..8, all orders <= %d" % NE, not bad, bad[:5])
    # (R2)
    bad = [(k, n) for k in KS if k + 5 in g for n in range(NE + 1) if g[k + 5][n] != g[k][n] + ((k + 2) * g[k + 1][n - 1] if n >= 1 else 0)]
    check("(R2) g_{k+5} = g_k + (k+2) zhat g_{k+1}, k = -4..3, all orders <= %d (k = -2: g_3 = g_-2)" % NE, not bad, bad[:5])
    # (R3)
    for k in KS:
        op = scalar_operator(k)
        assert op[5] == [Fr(-1)]
        bad = []
        for n in range(1, NE + 1):       # coefficient of zhat^n of zhat^5 T... g + g, needs g_{n-5+r}, r = 0..5
            tot = g[k][n]
            for r in range(6):
                m = n - 5 + r
                if m >= 0: tot += peval(op[r], m) * g[k][m]
            if tot != 0: bad.append(n)
        if bad: check("(R3) k = %d" % k, False, bad)
    check("(R3) zhat^5 T(T-1)(T-2)(T-3)(T-k-1) g_k + g_k = 0 for k = -4..8, orders 1..%d" % NE, not any(f.startswith("(R3)") for f in FAIL))
    op0 = scalar_operator(0)
    print("   scalar recursion for f = g_0 (identity row): Sum_{r=0}^{4} P_r(n-5+r) f_{n-5+r} = 0 with")
    for r in range(5): print("      P_%d(m) = %s" % (r, pfmt(op0[r])))
    check("(R3) P_4(m) vanishes at m = 0 only among integers 0..200 (f_0 free, f_n determined for n >= 1), k = 0",
          peval(op0[4], 0) == 0 and all(peval(op0[4], m) != 0 for m in range(1, 201)))
    for k in (-1, 1, 2, 3):
        opk = scalar_operator(k)
        check("     same for k = %d" % k, peval(opk[4], 0) == 0 and all(peval(opk[4], m) != 0 for m in range(1, 201)), pfmt(opk[4]))
    # (MC) recursion against the engine
    N = 40; rows = symbol_rows(N)
    bad = [(a, n) for a in range(5) for n in range(NE + 1) if rows[a][n] != g[CEXP[a]][n]]
    check("(MC) cyclic recursion reproduces g_{c_a}, a = 0..4, orders <= %d (50 non-trivial coefficients)" % NE, not bad, bad[:5])
    # symplecticity of the recursion output:  S_r(n) := Sum_{a, m: a+1-m = r mod 5} (-1)^{n-m} g^a_m g^{3-a}_{n-m} = delta_{n0} for every r
    bad = []
    for n in range(N + 1):
        for r in range(5):
            S = sum((-1) ** (n - m) * rows[a][m] * rows[(3 - a) % 5][n - m] for a in range(5) for m in range(n + 1) if (a + 1 - m - r) % 5 == 0)
            if S != (1 if n == 0 else 0): bad.append((n, r))
    check("(SP) symplecticity Sum_ab R^_a^(i)(zhat) eta^ab R^_b^(j)(-zhat) = delta_ij of the recursion output through zhat^%d" % N, not bad, bad[:5])
    # data: both primes
    for tag in ("89", "61"):
        d = json.load(open(os.path.join(HERE, "lg_symbol_%s.json" % tag))); p = d['p']; z = d['zeta']; Rh = d['Rh']
        bad = 0
        for a in range(5):
            for i in range(5):
                for n in range(7):
                    x = rows[a][n]; pred = NSIGN[a] * pow(z, (i * (a + 1 - n)) % 5, p) * (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p
                    if pred != Rh[a][i][n] % p: bad += 1
        check("(DATA) p = 2^%s-1: all 175 symbol residues R^_{a,n}^{(i)} = N_a zeta^{i(a+1-n)} g^a_n, N = (1,-1,-1,1,1), n <= 6" % tag, bad == 0, "%d mismatches" % bad)
    # genus one of the frozen theory from the model alone:  dF_1 = (1/48) Sum_i dlog Delta_i + (1/2) Sum_i (R_1)_ii du_i
    #   Delta_i ~ s^2 (Sum_i dlog Delta_i = 2 du/u),  (R_1)_ii = C_1[i][i]/s^4 = zeta^{-i} c11/s^4,  du_i = d(-(5/4) zeta^i s^4) = -zeta^i s^4 du/u
    c11 = sum(rows[a][1] for a in range(5)) / 5
    f1 = Fr(2, 48) + Fr(1, 2) * 5 * c11 * (-1)
    check("(F1) C_1[0][0] = (1/5) Sum_a g^a_1 = 1/20 and dF_1 = (1/24 - 5 C_1[0][0]/2) du/u = -(1/12) du/u", c11 == Fr(1, 20) and f1 == Fr(-1, 12), (str(c11), str(f1)))
    out = dict(comment="g[a][n] = [zhat^n] of J_{c_a}(-4 zhat), c = (0,-1,3,2,1); symbol R^_a^{(i)}(zhat) = N_a zeta^{i(a+1)} g^a(zeta^{-i} zhat), N = (1,-1,-1,1,1)",
               c=list(CEXP), N=list(NSIGN), g=[[str(x) for x in rows[a]] for a in range(5)])
    json.dump(out, open(os.path.join(HERE, "lg_symbol_rows_Q.json"), "w"), indent=0)
    print("rows g^a through zhat^8:")
    for a in range(5): print("   a = %d (k = %d): %s" % (a, CEXP[a], ", ".join(str(x) for x in rows[a][:9])))
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
