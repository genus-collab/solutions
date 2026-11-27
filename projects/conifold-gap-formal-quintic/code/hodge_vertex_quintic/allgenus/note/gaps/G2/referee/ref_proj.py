# -*- coding: utf-8 -*-
"""ref_proj.py -- referee check (key gaps/G2/referee) of lem:gau-comm (c)-(e) and of the two projection facts of
lem:gau-scal (a) (G1-G3), by a route independent of GAPS_G2 and of its scripts:
  * GENERAL weights (random polynomials in u, no quintic, no model, no pole needed: the facts are identities of rational
    functions in the weights and their first derivatives), evaluated at a random rational point u0 with DUAL NUMBERS
    (value + derivative), so that d/du is exact and no series truncation occurs;
  * pr computed by honest linear algebra at u0: X = Sum_j c_j Y^j + [At, W] solved for (c, W) (30 unknowns, 25 equations);
    nothing of the closed formula (lambda_m, calB) is used for pr.
Checks: commutant = span(Y^0..Y^4) (rank of ad At = 20); the decomposition is direct; closed formula (lem:gau-comm (e)) =
linear-algebra pr on random X; pr((Y^j)') = -(j/5)(Bprod'/Bprod) Y^j; pr(dlogTheta Y^j) = (tr dlogTheta/5) Y^j;
lem:gau-scal (a) on random c = Sum c_j Y^j.  Failing controls: factor j/4; non-diagonal matrix in place of dlogTheta;
Frobenius-orthogonal projection in place of pr; (EP) broken (Theta does not close).
Usage: py -I -B ref_proj.py [npoints].  Stdlib only; imports nothing from the project.
"""
import sys, random
from fractions import Fraction as Fr

FAIL = []


def check(name, ok):
    print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
    if not ok:
        FAIL.append(name)


class Du(object):
    """dual number a + b eps, eps^2 = 0: value and exact first derivative at u0"""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a, self.b = Fr(a), Fr(b)

    def __add__(self, o):
        o = o if isinstance(o, Du) else Du(o)
        return Du(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __neg__(self):
        return Du(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Du) else Du(o)))

    def __mul__(self, o):
        o = o if isinstance(o, Du) else Du(o)
        return Du(self.a * o.a, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def inv(self):
        assert self.a != 0
        return Du(1 / self.a, -self.b / (self.a * self.a))

    def __truediv__(self, o):
        o = o if isinstance(o, Du) else Du(o)
        return self * o.inv()


def poly_at(cs, u0):
    """random polynomial sum cs[i] u^i at u0, as a dual number"""
    return Du(sum(c * u0 ** i for i, c in enumerate(cs)), sum(i * c * u0 ** (i - 1) for i, c in enumerate(cs) if i))


def zero():
    return [[Du(0) for _ in range(5)] for _ in range(5)]


def mm(X, Y):
    return [[sum((X[i][k] * Y[k][j] for k in range(5)), Du(0)) for j in range(5)] for i in range(5)]


def val(X):
    return [[X[i][j].a for j in range(5)] for i in range(5)]


def der(X):
    return [[X[i][j].b for j in range(5)] for i in range(5)]


def solve(rows, rhs):
    """exact solution of a consistent linear system over Q (any solution); returns (solution, rank)"""
    n = len(rows[0]); M = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]; M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    assert all(all(x == 0 for x in M[i][:n]) <= (M[i][n] == 0) for i in range(len(M))), "inconsistent"
    sol = [Fr(0)] * n
    for i, c in enumerate(piv):
        sol[c] = M[i][n]
    return sol, r


def pr_linalg(At, Ypw, X):
    """coefficients c_j with X - Sum c_j Y^j in im(ad At), by solving 25 equations in 30 unknowns at the point"""
    cols = []
    for j in range(5):
        cols.append([Ypw[j][r][s] for r in range(5) for s in range(5)])
    for a in range(5):
        for b in range(5):
            E = [[Fr(int(r == a and s == b)) for s in range(5)] for r in range(5)]
            C = [[sum(At[r][k] * E[k][s] - E[r][k] * At[k][s] for k in range(5)) for s in range(5)] for r in range(5)]
            cols.append([C[r][s] for r in range(5) for s in range(5)])
    rows = [[cols[c][e] for c in range(30)] for e in range(25)]
    sol, rk = solve(rows, [X[r][s] for r in range(5) for s in range(5)])
    return sol[:5], rk
