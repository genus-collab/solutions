# -*- coding: utf-8 -*-
"""v1_sp_indep.py -- ADVERSARIAL re-derivation of the stationary-phase series of lambda = x^4 + u/x by a DIFFERENT algorithm
(Morse coordinate + Lagrange inversion, no exponential expansion), and comparison with
  (a) the report's engine lg/lg_sp.J_series,   (b) the report's closed recursion lg/lg_gm.symbol_rows,
  (c) the conifold symbol re-derived here from the project pipeline (Rows -> Conifold -> SRows -> Symbol), both primes,
  (d) the top L-coefficients of the WKB series of the Picard-Fuchs operator (fq_rows.WKB), identity row.
Nothing is written to disk (sys.dont_write_bytecode is set before any import).

Method.  x = x_c (1+y), x_c^5 = u/4, lambda/z = phi(y)/eps, phi = (1+y)^4 + 4/(1+y) = 5 + 10 y^2 B(y), B(0) = 1.
Morse coordinate w = y sqrt(B(y)):  lambda/z = (5 + 10 w^2)/eps.  With h = -eps the weight is exp(-10 w^2/h), <w^{2j}> = (2j-1)!! (h/20)^j.
     Int exp(lambda/z) x^k dx  ~  Sum_j <w^{2j}> [w^{2j}] ( (1+y(w))^k dy/dw ),
     [w^{2j}] ((1+y)^k y'(w)) = Res_w (1+y)^k dy / w^{2j+1} = [y^{2j}] (1+y)^k B(y)^{-(2j+1)/2}        (Lagrange),
so   j_{k,n} = (-1)^n (2n-1)!! / 20^n . [y^{2n}] (1+y)^k B(y)^{-(2n+1)/2}       (coefficient of eps^n in J_k).
Run:  py -I v1_sp_indep.py [N]      (default N = 12)
"""
import sys
sys.dont_write_bytecode = True
import os, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
LG = os.path.join(os.path.dirname(HERE), "lg"); PROJ = os.path.dirname(os.path.dirname(HERE))

CEXP = (0, -1, 3, 2, 1); NS = (1, -1, -1, 1, 1)


def mul(A, B, n):
    out = [Fr(0)] * (n + 1)
    for i, x in enumerate(A[:n + 1]):
        if x == 0: continue
        for j, y in enumerate(B[:n + 1 - i]):
            if y != 0: out[i + j] += x * y
    return out


def binom(k, n):
    out = [Fr(1)]
    for j in range(1, n + 1): out.append(out[-1] * (Fr(k) - (j - 1)) / j)
    return out


def power(B, al, n):
    """B(y)^al through y^n for B(0) = 1 and rational al:  P' B = al B' P"""
    P = [Fr(1)] + [Fr(0)] * n
    for m in range(1, n + 1):
        s = Fr(0)
        for j in range(1, m + 1):
            if j < len(B): s += (al * j - (m - j)) * B[j] * P[m - j]
        P[m] = s / m
    return P


def Bser(n):
    """B(y) = (phi(y) - 5)/(10 y^2) through y^n"""
    b4 = binom(4, n + 2); bm = binom(-1, n + 2)
    phi = [b4[m] + 4 * bm[m] for m in range(n + 3)]
    assert phi[0] == 5 and phi[1] == 0 and phi[2] == 10
    return [phi[m + 2] / 10 for m in range(n + 1)]


def dfact(n):
    r = 1
    for k in range(1, 2 * n, 2): r *= k
    return r


def J_indep(k, N):
    """[j_{k,0..N}] coefficients of eps^n"""
    B = Bser(2 * N); bk = binom(k, 2 * N); out = []
    for n in range(N + 1):
        P = power(B, Fr(-(2 * n + 1), 2), 2 * n)
        c = sum(bk[i] * P[2 * n - i] for i in range(2 * n + 1))
        out.append(Fr((-1) ** n * dfact(n), 20 ** n) * c)
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    t0 = time.time()
    KS = list(range(-4, 9))
    J = {k: J_indep(k, N) for k in KS}
    print("independent J_k (Lagrange inversion), k = -4..8, through eps^%d: %.1f s" % (N, time.time() - t0))
    print("   J_0 =", ", ".join(str(x) for x in J[0][:7]))
    # hand value: first correction of a Gaussian integral with Phi''' = 0: A_1 = -(1/8) Phi''''/Phi''^2 = -(1/8) 120/400 = -3/80 in h = -eps
    check("hand value j_{0,1} = 3/80 and j_{0,2} = 9/12800 (task's hand check)", J[0][1] == Fr(3, 80) and J[0][2] == Fr(9, 12800))
    g = {k: [x * (-4) ** n for n, x in enumerate(J[k])] for k in KS}
    check("f = J_0(-4 zhat) = 1 - 3/20 zhat + 9/800 zhat^2 + 11571/80000 zhat^3 (task's identity row)",
          g[0][:4] == [1, Fr(-3, 20), Fr(9, 800), Fr(11571, 80000)], [str(x) for x in g[0][:4]])
    # (a) the report's engine
    sys.path.insert(0, LG)
    from lg_sp import J_series
    bad = [(k, n) for k in KS for n in range(min(N, 10) + 1) if J_series(k, min(N, 10))[n] != J[k][n]]
    check("(a) report engine lg_sp.J_series == independent series, k = -4..8, orders <= %d" % min(N, 10), not bad, bad[:5])
    # (b) the report's closed recursion and its stored JSON
    from lg_gm import symbol_rows
    rows = symbol_rows(N)
    bad = [(a, n) for a in range(5) for n in range(N + 1) if rows[a][n] != g[CEXP[a]][n]]
    check("(b) report recursion (MC) == independent J_{c_a}(-4 zhat), a = 0..4, orders <= %d (THEOREM B to this order)" % N, not bad, bad[:5])
    import json
    dj = json.load(open(os.path.join(LG, "lg_symbol_rows_Q.json")))
    bad = [(a, n) for a in range(5) for n in range(N + 1) if Fr(dj['g'][a][n]) != g[CEXP[a]][n]]
    check("(b) stored lg_symbol_rows_Q.json == independent series, orders <= %d" % N, not bad, bad[:5])
    check("(b) g_3 == g_{-2} identically (orders <= %d)" % N, g[3] == g[-2])
    # (c) the symbol from the project pipeline, re-derived here
    sys.path.insert(0, PROJ)
    from fq_fp import FpRing, P1, P2
    from fq_rows import Rows, WKB
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq_symbol import Symbol
    NZ = min(N, 10)
    for tag, p in (("89", P2), ("61", P1)):
        R = FpRing(p); t0 = time.time()
        rws = Rows(R, Nz=NZ); C = Conifold(R, Nu=10, alpha=Fr(0), beta=Fr(0)); S = SRows(rws, C, prec_abs=34); sym = Symbol(S)
        bad = 0; badguess = {a: 0 for a in range(5)}; badsign = 0
        for a in range(5):
            for i in range(5):
                for n in range(NZ + 1):
                    ph = R.pow(R.zeta, (i * (a + 1 - n)) % 5)
                    if sym.Rh[a][i][n] != NS[a] * ph * R.of(g[CEXP[a]][n]) % p: bad += 1
                    if sym.Rh[a][i][n] != ph * R.of(g[CEXP[a]][n]) % p: badsign += 1
                    if (sym.Rh[a][i][n] - sym.Rh[a][i][0] * R.pow(R.zeta, (-i * n) % 5) * R.of(g[a][n])) % p: badguess[a] += 1
        check("(c) p = 2^%s-1: symbol re-derived (Nz = %d), %d residues == N_a zeta^{i(a+1-n)} [zhat^n] J_{c_a}(-4 zhat)" % (tag, NZ, 25 * (NZ + 1)),
              sym.val_ok and bad == 0, "%d mismatches, %.1f s" % (bad, time.time() - t0))
        check("(c) control: without the signs N the comparison FAILS (mismatches > 0)", badsign > 0, "%d mismatches" % badsign)
        check("(c) control: the labelling 'row a = J_a' fails for a = 1..4 and holds for a = 0", badguess[0] == 0 and all(badguess[a] > 0 for a in range(1, 5)), badguess)
        # which k reproduces the identity row
        unit = [k for k in KS if all(sym.Rh[0][0][n] == R.of(g[k][n]) for n in range(NZ + 1))]
        check("(c) identity row == J_k(-4 zhat) only for k = 0 among k = -4..8", unit == [0], unit)
        # (d) WKB top coefficients, identity row, to order N
        for i in (0, 2):
            w = WKB(R, i, N)
            ok = w.failed is None and all(max(w.phi[n]) == 4 * n and w.phi[n][4 * n] == R.pow(R.zeta, (-i * n) % 5) * R.of(g[0][n]) % p for n in range(N + 1))
            check("(d) p = 2^%s-1, colour %d: deg_L phi_n = 4n and [L^{4n}] phi_n = zeta^{-in} [zhat^n] J_0(-4 zhat), n <= %d" % (tag, i, N), ok)
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
