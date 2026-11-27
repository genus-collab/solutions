# -*- coding: utf-8 -*-
"""v4_recursions.py -- ADVERSARIAL checks of the closed recursions of LG_MODEL.md s.3 against the INDEPENDENT stationary-phase series
of v1_sp_indep.J_indep (Morse coordinate + Lagrange inversion), to higher order than the report's own check (zhat^10):
 (S) the five-term scalar recursion for the identity row, typed here from the printed formula of LG_MODEL.md s.3 (not imported):
        4 n f_n = -(1/5)(32m^2+32m+3)|_{m=n-1} f_{n-1} - (4/25)(32m^3+96m^2+73m+12)|_{m=n-2} f_{n-2}
                  - (4/625)(320m^4+1920m^3+3700m^2+2580m+459)|_{m=n-3} f_{n-3} - (1/3125)(4m+1)(4m+6)^2(4m+11)(4m+16)|_{m=n-4} f_{n-4}
     against J_0(-4 zhat) (orders <= N) and against row 0 of the stored lg_symbol_rows_Q.json (orders <= 40);
 (M) the cyclic recursion (MC), re-implemented here from its statement, against J_{c_a}(-4 zhat), a = 0..4, orders <= N; and against the
     stored JSON (orders <= 40);
 (R) the Gauss-Manin relations (R1), (R2) on the independent series, k = -4..7, orders <= N;
 (P) symplecticity of the model rows built from the independent series, orders <= N;
 (Y) the matrix form  theta Y = ((5/4) P/zhat + mu) Y  coefficientwise:  (n - mu_a) g^a_n = (5/4)(g^{a+1}_{n+1} - g^a_{n+1}).
Run:  py -I v4_recursions.py [N]     (default N = 24)
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from v1_sp_indep import J_indep

CEXP = (0, -1, 3, 2, 1)
MU = [Fr(c - 1, 4) for c in CEXP]


def scalar_rec(N):
    f = [Fr(1)]
    for n in range(1, N + 1):
        s = Fr(0)
        m = n - 1; s -= Fr(1, 5) * (32 * m * m + 32 * m + 3) * f[n - 1]
        if n >= 2: m = n - 2; s -= Fr(4, 25) * (32 * m ** 3 + 96 * m ** 2 + 73 * m + 12) * f[n - 2]
        if n >= 3: m = n - 3; s -= Fr(4, 625) * (320 * m ** 4 + 1920 * m ** 3 + 3700 * m ** 2 + 2580 * m + 459) * f[n - 3]
        if n >= 4: m = n - 4; s -= Fr(1, 3125) * (4 * m + 1) * (4 * m + 6) ** 2 * (4 * m + 11) * (4 * m + 16) * f[n - 4]
        f.append(s / (4 * n))
    return f


def cyclic_rec(N):
    """g^{a+1}_n - g^a_n = ((4n-3-c_a)/5) g^a_{n-1} (a in Z/5), Sum_a (4n+1-c_a) g^a_n = 0"""
    g = [[Fr(1)] for _ in range(5)]
    for n in range(1, N + 1):
        d = [Fr(0)]
        for a in range(4): d.append(d[-1] + Fr(4 * n - 3 - CEXP[a], 5) * g[a][n - 1])
        wrap = d[4] + Fr(4 * n - 3 - CEXP[4], 5) * g[4][n - 1]
        if wrap != 0: raise AssertionError("cycle does not close at level %d: %s" % (n, wrap))
        t = -sum((4 * n + 1 - CEXP[a]) * d[a] for a in range(5)) / (20 * n)
        for a in range(5): g[a].append(t + d[a])
    return g


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    t0 = time.time()
    KS = list(range(-4, 8))
    g = {k: [x * (-4) ** n for n, x in enumerate(J_indep(k, N))] for k in KS}
    print("independent series k = -4..7 through zhat^%d: %.1f s" % (N, time.time() - t0))
    dj = json.load(open(os.path.join(os.path.dirname(HERE), "lg", "lg_symbol_rows_Q.json")))
    stored = [[Fr(x) for x in row] for row in dj['g']]
    f = scalar_rec(40)
    check("(S) scalar recursion (typed from the note) == J_0(-4 zhat), orders <= %d" % N, f[:N + 1] == g[0])
    check("(S) scalar recursion == row 0 of the stored JSON, orders <= 40", f == stored[0][:41])
    rows = cyclic_rec(40)
    bad = [(a, n) for a in range(5) for n in range(N + 1) if rows[a][n] != g[CEXP[a]][n]]
    check("(M) cyclic recursion (re-implemented) == J_{c_a}(-4 zhat), a = 0..4, orders <= %d  [Theorem B to this order]" % N, not bad, bad[:5])
    check("(M) cyclic recursion == stored JSON, all rows, orders <= 40", all(rows[a] == stored[a][:41] for a in range(5)))
    bad = [(k, n) for k in KS[1:] for n in range(N + 1)
           if g[k - 1][n] != g[k][n] - (Fr(k - 1 - 4 * (n - 1), 5) * g[k][n - 1] if n else 0)]
    check("(R1) g_{k-1} = g_k - (zhat/5)(k-1-4 theta) g_k, k = -3..7, orders <= %d" % N, not bad, bad[:5])
    bad = [(k, n) for k in KS if k + 5 in g for n in range(N + 1) if g[k + 5][n] != g[k][n] + ((k + 2) * g[k + 1][n - 1] if n else 0)]
    check("(R2) g_{k+5} = g_k + (k+2) zhat g_{k+1}, k = -4..2, orders <= %d" % N, not bad, bad[:5])
    check("(R2) consequence: g_k != g_{k+5} for k = -4, -3, -1, 0, 1, 2 (only g_3 = g_{-2}): 'x^{-a} dx mod 5' is NOT an identity of classes",
          all(g[k] != g[k + 5] for k in (-4, -3, -1, 0, 1, 2)) and g[-2] == g[3])
    bad = []
    for n in range(N + 1):
        for r in range(5):
            S = sum((-1) ** (n - m) * g[CEXP[a]][m] * g[CEXP[(3 - a) % 5]][n - m] for a in range(5) for m in range(n + 1) if (a + 1 - m - r) % 5 == 0)
            if S != (1 if n == 0 else 0): bad.append((n, r))
    check("(P) symplecticity of the model rows (independent series), orders <= %d" % N, not bad, bad[:5])
    bad = [(a, n) for a in range(5) for n in range(N) if (n - MU[a]) * g[CEXP[a]][n] != Fr(5, 4) * (g[CEXP[(a + 1) % 5]][n + 1] - g[CEXP[a]][n + 1])]
    check("(Y) (n - mu_a) g^a_n = (5/4)(g^{a+1}_{n+1} - g^a_{n+1}), mu = (-1/4,-1/2,1/2,1/4,0), orders < %d" % N, not bad, bad[:5])
    print("f_7, f_8 =", f[7], f[8])
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
