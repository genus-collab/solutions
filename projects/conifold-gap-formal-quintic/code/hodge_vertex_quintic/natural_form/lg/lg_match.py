# -*- coding: utf-8 -*-
"""lg_match.py -- compare the conifold symbol of the formal quintic (lg_symbol_<tag>.json, dumped by lg_dump_symbol.py from
fq3_par.setup; residues mod p) with the stationary-phase series J_k of the Landau-Ginzburg model lambda = x^4 + u/x (lg_sp.J_series).

Blocks (each run at both primes 2^89-1 and 2^61-1; the J_k are exact rationals, reduced mod p for the comparison):
 (A) AGNOSTIC identification: for every flat index a and every k in -3..6, fit N_a (zhat^0) and c (eps = c zhat, from zhat^1) and test the
     remaining orders zhat^2..zhat^6 of  R^_a^{(0)}(zhat) = N_a J_k(c zhat).  Reports all (a, k) that pass and the fitted (N_a, c).
 (B) the predicted dictionary  R^_a^{(i)}(zhat) = zeta^{i(a+1)} J_{c_a}(-4 zeta^{-i} zhat),  c = (0,-1,3,2,1): all 5 x 5 x 7 residues.
 (C) controls: the identity row against the J_k of the other rank-5 potentials x^p + u/x^(5-p), p = 1, 2, 3 (k = -6..6), and against the
     A_4 / A_5 type x^5 + u x, x^6 + u x ... is not attempted (different covariance); only what is stated is tested.
 (D) Teleman colour form: C_n[j][i] = (1/5) Sum_a zeta^{j(4-a) + i(a+1-n)} j_{c_a,n} (-4)^n against the colour matrices rebuilt from the
     dumped symbol (same formula as fq_symbol.Symbol.colour_matrix).
Run:  py -I lg_match.py        (writes nothing; prints; redirect to a log)
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lg_sp import J_series, J_general

CEXP = (0, -1, 3, 2, 1)
NZ = 6


def fp(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def ratrec(a, p, B=10 ** 12):
    a %= p
    if a == 0: return Fr(0)
    r0, r1 = p, a; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1)


def main():
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    KS = list(range(-3, 7))
    J = {k: J_series(k, NZ) for k in KS}
    print("J_k(eps) coefficients (exact), k = -3..6, through eps^6:")
    for k in KS: print("   J_%d: %s" % (k, ", ".join(str(x) for x in J[k])))
    for tag in ("89", "61"):
        d = json.load(open(os.path.join(HERE, "lg_symbol_%s.json" % tag)))
        p = d['p']; z = d['zeta']; Rh = d['Rh']
        print("=" * 100); print("p = 2^%s-1; dump checks: %s" % (tag, d['checks']))
        # ---------------------------------------------------------------- (A)
        print("(A) agnostic fit R^_a^(0)(zhat) = N_a J_k(c zhat), k in -3..6")
        passes = {}
        for a in range(5):
            row = Rh[a][0]
            for k in KS:
                j = [fp(x, p) for x in J[k]]
                N = row[0]
                if j[1] == 0: continue
                c = row[1] * pow(N * j[1] % p, p - 2, p) % p
                ok = [(N * pow(c, n, p) * j[n] - row[n]) % p == 0 for n in range(NZ + 1)]
                if all(ok):
                    passes.setdefault(a, []).append((k, ratrec(N, p), ratrec(c, p)))
            print("   a = %d: row (rational reconstruction where possible) %s" % (a, [str(ratrec(x, p)) for x in row]))
            print("          matching (k, N_a, c): %s" % [(k, str(N), str(c)) for (k, N, c) in passes.get(a, [])])
        # first run (see RUNLOG.md): the naive expectation N_a = 1 for all a FAILED for a = 1, 2 (N_1 = N_2 = -1); the signs
        # N = (1,-1,-1,1,1) are those of the chain X_{a+1} = (C1 kappa_a)^{-1} (zD + zeta^i L) X_a with D -> -d/du, C1 -> -1 (LG_MODEL.md s.4)
        exp = {0: 0, 1: -1, 2: 3, 3: 2, 4: 1}; NS = (1, -1, -1, 1, 1)
        for a in range(5):
            good = [(k, N, c) for (k, N, c) in passes.get(a, []) if N == NS[a] and c == -4]
            allowed = {exp[a], -2} if exp[a] == 3 else {exp[a]}
            check("(A) row a = %d equals N_a J_k(-4 zhat), N_a = %d, exactly for k in %s (expected c_a = %d; J_3 = J_-2 identically)" % (a, NS[a], [g[0] for g in good], exp[a]),
                  set(g[0] for g in good) == allowed and len(passes.get(a, [])) == len(good))
        # how many orders are a genuine test: for the matched k the orders 2..6 (5 coefficients per row), 25 in all
        # ---------------------------------------------------------------- (B)
        bad = 0
        for a in range(5):
            j = [fp(x, p) for x in J[CEXP[a]]]
            for i in range(5):
                for n in range(NZ + 1):
                    pred = NS[a] * pow(z, (i * (a + 1 - n)) % 5, p) * pow(p - 4, n, p) * j[n] % p
                    if pred != Rh[a][i][n] % p: bad += 1
        check("(B) all 175 residues R^_{a,n}^{(i)} = N_a zeta^{i(a+1-n)} (-4)^n j_{c_a,n}, N = (1,-1,-1,1,1)", bad == 0, "%d mismatches" % bad)
        # the task's first guess H^a <-> x^a dx (rows J_a): which rows would pass?
        guess = [a for a in range(5) if all((Rh[a][0][0] * pow(p - 4, n, p) * fp(J[a][n], p) - Rh[a][0][n]) % p == 0 for n in range(NZ + 1))]
        check("(B) NEGATIVE: the guess H^a <-> x^a dx (row a = N_a J_a) holds only for a = 0 (the correct labelling is k = -a mod 5)", guess == [0], guess)
        # which unit: the identity row against every k (primitive form x^m dx <-> identity row J_m)
        unit = [k for k in KS if all((pow(p - 4, n, p) * fp(J[k][n], p) - Rh[0][0][n]) % p == 0 for n in range(NZ + 1))]
        check("(B) identity row f = J_k(-4 zhat) only for k = 0 (primitive form dx; k = -1 would be dx/x, k = 1 x dx)", unit == [0], unit)
        first_diff = {k: next((n for n in range(NZ + 1) if (pow(p - 4, n, p) * fp(J[k][n], p) - Rh[0][0][n]) % p), None) for k in KS}
        print("      first zhat-order at which J_k(-4 zhat) differs from f:", first_diff)
        # ---------------------------------------------------------------- (C)
        print("(C) controls: f against N J_k(c zhat) for lambda = x^P + u/x^Q, P + Q = 5 (c fitted at zhat^1)")
        for (P, Q) in ((4, 1), (3, 2), (2, 3), (1, 4)):
            hits = []
            for k in range(-6, 7):
                lc, a2, ser = J_general(k, NZ, power=P, pole=Q)
                j = [fp(x, p) for x in ser]
                if j[1] == 0: continue
                c = Rh[0][0][1] * pow(j[1], p - 2, p) % p
                nok = next((n for n in range(NZ + 1) if (pow(c, n, p) * j[n] - Rh[0][0][n]) % p), None)
                if nok is None: hits.append((k, str(ratrec(c, p))))
            print("   (P, Q) = (%d, %d): matching (k, c) = %s" % (P, Q, hits))
            if (P, Q) == (4, 1): check("(C) (4,1): only k = 0, c = -4", hits == [(0, "-4")], hits)
            elif (P, Q) == (1, 4): check("(C) (1,4): only k = -2, c = -1 (the SAME model under x -> 1/x: x + u/x^4 with x^-2 dx = -d(1/x))", hits == [(-2, "-1")], hits)
            else: check("(C) (%d,%d): no k in -6..6 matches" % (P, Q), hits == [], hits)
        # ---------------------------------------------------------------- (D)
        eta = [(a, b, e) for (a, b, e) in d['eta_inv']]
        bad = 0; inv5 = pow(5, p - 2, p)
        for n in range(NZ + 1):
            for jj in range(5):
                for i in range(5):
                    C = sum(Rh[b][jj][0] * e * Rh[a][i][n] for (a, b, e) in eta) % p
                    pred = sum(pow(z, (jj * (4 - a) + i * (a + 1 - n)) % 5, p) * fp(J[CEXP[a]][n], p) for a in range(5)) * pow(p - 4, n, p) * inv5 % p
                    if C != pred: bad += 1
        check("(D) colour matrices C_0..C_6 (175 entries) = (1/5) Sum_a zeta^{j(4-a)+i(a+1-n)} (-4)^n j_{c_a,n}", bad == 0, "%d mismatches" % bad)
        check("(D) eta^{ab} of the dump is (1/5) delta_{a+b = 3 mod 5}", sorted((a, b) for (a, b, e) in eta) == sorted((a, (3 - a) % 5) for a in range(5)) and all(e == inv5 for (a, b, e) in eta))
    print("=" * 100)
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
