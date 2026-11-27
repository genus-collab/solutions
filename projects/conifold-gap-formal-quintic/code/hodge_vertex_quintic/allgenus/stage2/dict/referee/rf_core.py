"""Referee of Stage 2 / D(iii): independent core.  Imports NOTHING from ../dict_*.py or ../../../critic/.
Pure stdlib; exact over Q (P = None) or modulo a prime P.  Run the scripts with  py -I -B.

Model curve x = W = p^4/4 + v/p, y = p.  Weights p:1, v:5, x:4, q = v/p:4, log p:0.
A homogeneous element of weight w is a pair (A0, A1) of dicts {j: c} (only j = w mod 5 occur):
    Sum_j A0[j] v^((w-j)/5) p^j  +  log p . Sum_j A1[j] v^((w-j)/5) p^j          (v-exponent implicit).
Two leg constructions (b = 0..4 <-> 1, H, H^2, H^3, H^4):
  legs_iter  : Phi_0 = p, log p, v/p, p^3, p^2 ;  Phi_(k+1) = -Int Phi_k dx  (termwise in p, constant 0).
  legs_closed: generating function.  G = Sum_k z^k Phi_k solves dG/dp = dPhi_0/dp - z G dx/dp, so
               e^(zx) G = Sum_n z^n (A_n + C_n),  A_n := Int x^n dPhi_0 / n!  (constant 0),  C_n p-constants;
               hence Phi_k = Sum_j (-x)^j/j! A_(k-j)  modulo Q[v][x].   (Different algorithm from legs_iter.)
"""
from fractions import Fraction
from math import comb, factorial

W0 = (1, 0, 4, 3, 2)                                        # weight of Phi_0^b
PHI0 = (({1: 1}, {}), ({}, {0: 1}), ({-1: 1}, {}), ({3: 1}, {}), ({2: 1}, {}))
DPHI0 = ({0: 1}, {-1: 1}, {-2: -1}, {2: 3}, {1: 2})          # dPhi_0^b/dp (weight W0[b]-1); b=2: -v p^-2
XP = {3: 1, -2: -1}                                         # dx/dp = p^3 - v p^-2  (weight 3)
P61, P89 = 2 ** 61 - 1, 2 ** 89 - 1


class Ar:
    """Coefficient arithmetic: Q (P None, Fractions) or Z/P."""

    def __init__(self, P=None):
        self.P = P

    def q(self, a, b=1):
        if self.P is None:
            return Fraction(a, b)
        return (a % self.P) * pow(b % self.P, self.P - 2, self.P) % self.P

    def inv(self, a):
        return 1 / Fraction(a) if self.P is None else pow(a % self.P, self.P - 2, self.P)

    def nz(self, c):
        return c != 0 if self.P is None else c % self.P != 0

    def clean(self, A):
        if self.P is None:
            return {k: c for k, c in A.items() if c != 0}
        return {k: c % self.P for k, c in A.items() if c % self.P}

    def lift(self, A):
        return self.clean({k: self.q(c) for k, c in A.items()})


def acc(D, k, c):
    D[k] = D.get(k, 0) + c


def lin(ar, *terms):
    """Sum of s * A over (s, A) pairs."""
    C = {}
    for s, A in terms:
        for k, c in A.items():
            acc(C, k, s * c)
    return ar.clean(C)


def mul(ar, A, B):
    C = {}
    for j1, c1 in A.items():
        for j2, c2 in B.items():
            acc(C, j1 + j2, c1 * c2)
    return ar.clean(C)


_XP = {}


def xpow(ar, n):
    """x^n = Sum_t C(n,t) 4^(t-n) v^t p^(4n-5t)  (weight 4n); lowest key -n has coefficient 1."""
    key = (ar.P, n)
    if key not in _XP:
        _XP[key] = ar.clean({4 * n - 5 * t: ar.q(comb(n, t), 4 ** (n - t)) for t in range(n + 1)})
    return _XP[key]


def integ(ar, A0, A1):
    """Termwise p-antiderivative of A0 + log p A1 (constant 0): p^-1 -> log p; p^a log p -> p^(a+1) log p/(a+1) - p^(a+1)/(a+1)^2."""
    B0, B1 = {}, {}
    for a, c in A0.items():
        if a == -1:
            acc(B1, 0, c)
        else:
            acc(B0, a + 1, c * ar.q(1, a + 1))
    for a, c in A1.items():
        if a == -1:
            raise ValueError("(log p)^2 term: p^-1 log p in the integrand")
        acc(B1, a + 1, c * ar.q(1, a + 1))
        acc(B0, a + 1, -c * ar.q(1, (a + 1) ** 2))
    return ar.clean(B0), ar.clean(B1)


def I_step(ar, A0, A1):
    """Descendant step F -> -Int F dx (weight +4)."""
    xp = ar.lift(XP)
    B0, B1 = integ(ar, mul(ar, A0, xp), mul(ar, A1, xp))
    return lin(ar, (-1, B0)), lin(ar, (-1, B1))


def legs_iter(ar, K, phi0=None):
    phi0 = phi0 or PHI0
    T = {}
    for b in range(5):
        A0, A1 = ar.lift(phi0[b][0]), ar.lift(phi0[b][1])
        for k in range(K + 1):
            T[(b, k)] = (A0, A1)
            if k < K:
                A0, A1 = I_step(ar, A0, A1)
    return T


def legs_closed(ar, K):
    T = {}
    for b in range(5):
        d = ar.lift(DPHI0[b])
        A = [integ(ar, lin(ar, (ar.q(1, factorial(n)), mul(ar, xpow(ar, n), d))), {}) for n in range(K + 1)]
        for k in range(K + 1):
            S0, S1 = {}, {}
            for j in range(k + 1):
                f, X = ar.q((-1) ** j, factorial(j)), xpow(ar, j)
                S0 = lin(ar, (1, S0), (f, mul(ar, X, A[k - j][0])))
                S1 = lin(ar, (1, S1), (f, mul(ar, X, A[k - j][1])))
            T[(b, k)] = (S0, S1)
    return T


def wt(b, k):
    return W0[b] + 4 * k


# ---------------- classes and normal forms ----------------
def in_Qpq(A0, w):
    """Support rule of Lemma 1(c): every term v^m p^j has m >= 0 and j >= -m, i.e. -w <= 4j and j <= w."""
    return all(-w <= 4 * j and j <= w for j in A0)


def x_coords(ar, A1, w):
    """Write A1 (weight w) as Sum_n a_n v^((w-4n)/5) x^n by top-down elimination.  Returns ({n: a_n}, ok) where
    ok = (A1 is a polynomial in x with v-exponents >= 0); on failure returns (partial, False)."""
    G, res = dict(A1), {}
    while G:
        jt = max(G)
        if jt < 0 or jt % 4 or jt > w:
            return ar.clean(res), False
        n = jt // 4
        a = G[jt] * ar.q(4 ** n)
        res[n] = a
        G = lin(ar, (1, G), (-a, xpow(ar, n)))
    return ar.clean(res), True                  # clean() also reduces mod P


def nf_vx(ar, A0, w):
    """Bottom-up elimination of negative p-powers by multiples of v^((w-4n)/5) x^n.  Returns (rem, negv):
    rem = remainder without its p^0 term (keys j >= 1; j > w means a negative v-exponent), negv = True iff some
    multiplier had a negative v-exponent (4n > w), i.e. a Q(v)[x]-, not Q[v][x]-, subtraction was needed."""
    G, negv = dict(A0), False
    while True:
        neg = [j for j in G if j < 0]
        if not neg:
            break
        jb = min(neg)
        n = -jb
        negv = negv or 4 * n > w
        G = lin(ar, (1, G), (-G[jb], xpow(ar, n)))
    return {j: c for j, c in G.items() if j != 0}, negv


def in_L(ar, A0, w):
    """Is the log-free A0 (weight w >= 0) in L = Q[v][p] + Q(v)[x]?  (remainder must have v-exponents >= 0)."""
    rem, _ = nf_vx(ar, A0, w)
    return all(j <= w for j in rem)


def in_M(ar, A0, A1, w):
    """A0 + log p A1 in M = Q[p,q] (+) log p Q[v,x]?"""
    return in_Qpq(A0, w) and x_coords(ar, A1, w)[1]


def shiftv(T, b, k, m):
    """leg of v^m tau_k(H^b): same dicts (keys are p-exponents; v-exponent implicit in the weight)."""
    return T[(b, k)]


def triples(w):
    """(b, k, m) with W0[b] + 4k + 5m = w, m, k >= 0."""
    out = []
    for b in range(5):
        for k in range(w // 4 + 1):
            r = w - W0[b] - 4 * k
            if r >= 0 and r % 5 == 0:
                out.append((b, k, r // 5))
    return out


# ---------------- linear algebra (rows = dicts col -> coefficient) ----------------
def rank(ar, rows):
    piv = {}                                    # col -> reduced row with leading col
    r = 0
    for row in rows:
        row = ar.clean(dict(row))
        while row:
            col = min(row)
            if col not in piv:
                piv[col] = row
                r += 1
                break
            pr = piv[col]
            f = row[col] * ar.inv(pr[col])
            row = lin(ar, (1, row), (-f, pr))
    return r


def left_kernel(ar, rows):
    """Basis of {y : Sum_i y_i rows_i = 0}, generic elimination (no structure used)."""
    piv = {}
    kern = []
    for i, row in enumerate(rows):
        row, comb_ = ar.clean(dict(row)), {i: ar.q(1)}
        while row:
            col = min(row)
            if col not in piv:
                piv[col] = (row, comb_)
                break
            pr, pc = piv[col]
            f = row[col] * ar.inv(pr[col])
            row, comb_ = lin(ar, (1, row), (-f, pr)), lin(ar, (1, comb_), (-f, pc))
        if not row:
            kern.append(comb_)
    return kern
