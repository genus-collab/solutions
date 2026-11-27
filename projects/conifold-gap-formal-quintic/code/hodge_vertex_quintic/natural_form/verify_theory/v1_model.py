# -*- coding: utf-8 -*-
"""v1_model.py -- adversarial re-derivation for UNIFORM_REDUCTION.md s.1-2 (own code, no import from natural_form/theory).

(1) Stationary phase, done DIRECTLY on the five oscillatory integrals I_m = Int x^m exp(-(x^4/4 + v/x)/z) dx, m = 0, -1, 3, 2, 1
    (and m = -2 for the alternative representative v dx/x^2 of the H^2 row), at v = 1, critical point x = 1 (colour 0).
    Algorithm (different from ur_model_check.stationary_phase): y = h w, h = sqrt(z); integrand = exp(-(5/2) w^2) exp(S(h)),
    S = Sum_j s_j h^j, s_j = L_{j+2} w^(j+2) + m (-1)^(j+1) w^j / j; exponential by the recurrence n E_n = Sum j s_j E_{n-j};
    Gaussian moments <w^(2r)> = (2r-1)!!/5^r.  No ODE is used.
(2) WKB solution of the SCALAR equation (v th^5 + 3 z th^4 + 1) Y = 0, Y = e^{U/z} sigma^-1 F(t) (own operator calculus on monomials
    sigma^k t^n), compared with (1), m = 0.
(3) Comparison with the project's conifold symbol sym.Rh[b][i][n] (all rows, colours, n <= 6, both primes):
        Rh[b][i][n] = zeta^(i(b+1-n)) sign_b [z^n] (stationary phase of x^(m_b)),  m = (0,-1,3,2,1), sign = (1,-1,-1,1,1).
Run: py -I -B v1_model.py
"""
import os, sys
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
OUT = open(os.path.join(HERE, "v1_model.log"), "w", encoding="utf-8")


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def padd(a, b, c=1):
    for k, v in b.items():
        a[k] = a.get(k, 0) + c * v
    return a


def pmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] = out.get(i + j, 0) + x * y
    return out


def sp(m, N):
    """[z^0..z^N] of Int (1+y)^m exp((lam(1+y) - lam(1))/z) dy / sqrt(2 pi z/5), lam = -(x^4/4 + 1/x)"""
    J = 2 * N
    L = lambda k: Fr(-comb(4, k), 4) - Fr((-1) ** k)
    assert L(1) == 0 and L(2) == Fr(-5, 2)
    s = [None] + [{j + 2: L(j + 2), j: Fr(m * (-1) ** (j + 1), j)} for j in range(1, J + 1)]
    E = [{0: Fr(1)}]
    for n in range(1, J + 1):
        acc = {}
        for j in range(1, n + 1): padd(acc, pmul(s[j], E[n - j]), Fr(j, n))
        E.append(acc)
    out = []
    for n in range(N + 1):
        tot = Fr(0)
        for deg, c in E[2 * n].items():
            if deg % 2: continue
            r = deg // 2; dd = 1
            for q in range(1, deg, 2): dd *= q
            tot += c * dd * Fr(1, 5) ** r
        out.append(tot)
    return out


def scalar_wkb(N):
    """F with Y = e^{U/z} sigma^-1 F(t) solving (v th^5 + 3 z th^4 + 1) Y = 0, colour 0, F(0) = 1"""
    def th(f):
        out = {}
        for (k, n), c in f.items():
            out[(k - 1, n)] = out.get((k - 1, n), 0) - c
            out[(k - 1, n + 1)] = out.get((k - 1, n + 1), 0) + c * Fr(k - 4 * n, 5)
        return out

    def op(f):
        g = f
        for _ in range(4): g = th(g)
        g4 = g; g5 = th(g4)
        res = dict(f)
        for (k, n), c in g5.items(): res[(k + 5, n)] = res.get((k + 5, n), 0) + c
        for (k, n), c in g4.items(): res[(k + 4, n + 1)] = res.get((k + 4, n + 1), 0) + 3 * c
        return res

    F = [Fr(1)]
    for n in range(1, N + 1):
        base = op({(-1, j): F[j] for j in range(n)})
        unit = op({(-1, n): Fr(1)})
        assert unit.get((-1, n), 0) == 0
        piv = unit[(-1, n + 1)]
        F.append(-base.get((-1, n + 1), 0) / piv)
    res = op({(-1, j): F[j] for j in range(N + 1)})
    assert all(v == 0 for (k, n), v in res.items() if n <= N + 1), "scalar equation residual"
    assert all(k == -1 for (k, n) in res), "weights"
    return F


def main():
    N = 10; ok_all = True
    ms = (0, -1, 3, 2, 1); sg = (1, -1, -1, 1, 1)
    Fsp = [[sg[b] * x for x in sp(ms[b], N)] for b in range(5)]
    alt = [-x for x in sp(-2, N)]
    Fw = scalar_wkb(N)
    log("stationary phase, identity row (m = 0): " + ", ".join(str(x) for x in Fsp[0][:8]))
    quoted = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000),
              Fr(52989974037, 25600000000), Fr(-7749624763881, 512000000000)]
    for name, ok in (("stationary phase (m = 0) = F_0 quoted in the note through t^7", Fsp[0][:8] == quoted),
                     ("scalar-equation WKB = stationary phase (m = 0) through t^%d" % N, Fw == Fsp[0]),
                     ("H^2 row: -I_3 and -v I_{-2} have the same expansion through t^%d" % N, alt == Fsp[2]),
                     ("F_b(0) = (1,-1,-1,1,1)", [Fsp[b][0] for b in range(5)] == [1, -1, -1, 1, 1])):
        log(("PASS " if ok else "FAIL ") + name); ok_all = ok_all and ok
    q = {1: [Fr(1), Fr(1, 20), Fr(-111, 800), Fr(13191, 80000)], 2: [Fr(1), Fr(9, 20), Fr(-63, 800), Fr(-9009, 80000)],
         3: [Fr(1), Fr(1, 20), Fr(81, 800), Fr(-16569, 80000)], 4: [Fr(1), Fr(-3, 20), Fr(21, 160), Fr(-5229, 80000)]}
    ok = all([Fsp[b][n] / Fsp[b][0] for n in range(4)] == q[b] for b in q)
    log(("PASS " if ok else "FAIL ") + "rows 1..4 from the stationary phase of x^(-1), x^3, x^2, x^1 = the note's F_b/F_b(0) through t^3"); ok_all = ok_all and ok
    from fq_fp import P1, P2
    from fq3_par import setup
    for prime, tag in ((P2, "2^89-1"), (P1, "2^61-1")):
        st = setup((prime, "0", "0", 34, 6)); R = st['R']; p = R.p; sym = st['sym']
        bad = [(b, i, n) for b in range(5) for i in range(5) for n in range(7)
               if R.of(Fsp[b][n]) * R.z(i * (b + 1 - n)) % p != sym.Rh[b][i][n] % p]
        log(("PASS " if not bad else "FAIL ") + "p = %s: symbol Rh[b][i][n] = zeta^(i(b+1-n)) x stationary phase of sign_b x^(m_b) dx, all 175 entries %s" % (tag, bad[:5]))
        ok_all = ok_all and not bad
    log("V1_MODEL " + ("ALL PASS" if ok_all else "FAILURES")); OUT.close()


if __name__ == "__main__":
    main()
