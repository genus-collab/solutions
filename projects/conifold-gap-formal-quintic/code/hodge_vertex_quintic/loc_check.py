#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
loc_check.py -- step 3 of PLAN_2026-11-03.md: INDEPENDENT control by direct
Graber-Pandharipande virtual localisation for the formal quintic of
Lho-Pandharipande (arXiv:1702.06096, sec. 0.2).

Self-contained: own exact arithmetic (Fractions; Q(zeta_5) as Q[x]/(Phi_5)),
own Witten-Kontsevich numbers (string/dilaton/DVV), own Mumford-GRR recursion
for ch_k of the Hodge bundle, own graph enumeration.  Imports nothing from the
project.  Run:  py -I loc_check.py   (pure Python, ~1 min).

======================================================================
CONVENTIONS AND SIGN DECISIONS (all re-derived here; the graph-sum code
of step 3 is compared against the numbers produced by this file)
======================================================================
Torus.  T = (C^*)^5 acts on C^5 = H^0(P^4,O(1))^vee so that the coordinate
FUNCTION x_i (a section of O(1)) has weight lambda_i.  Hence
  O(1)|_{p_i} = lambda_i   (x_i trivialises O(1) near p_i),  H|_{p_i} = lambda_i,
  O(5)|_{p_i} = 5 lambda_i,
  T_{p_i}P^4 = Hom(l_i, C^5/l_i) has weights  lambda_i - lambda_j  (j != i),
  e(T_{p_i}) = e_i := prod_{j != i} (lambda_i - lambda_j)
  (consistent with prod_j (H - lambda_j) = 0 in H_T^*(P^4)).
The T-line L_ij = p_i p_j has tangent weight lambda_i - lambda_j at p_i; on
its degree-d cover C_e the tangent weight at the point over p_i is
  omega_(e,i) = (lambda_i - lambda_j)/d        (the FLAG weight).

Twist.  N_{g,d} = Int_{[Mbar_g(P^4,d)]^vir} e(R pi_* f^* O(5)), Euler class of
the K-theory class R^0 - R^1 in the NUMERATOR.  More generally the code takes
a list of line bundles O(a) (x) C_eps on P^r with a sign s = +1 (e(R pi_*) in
the numerator, the quintic case) or s = -1 (e(-R pi_*), i.e. e(R^1 - R^0),
the local-CY case used as a control: X = O(-1)+O(-1) -> P^1).

Weights of the cohomology of a line bundle on a degree-d cover C_e of L_ij
with fibre weights w_0 at p_i and w_0 + m*omega at p_j (m = deg on C_e):
  m >= 0 : H^0 weights  w_0 - k*omega,  k = 0..m   (sections z^k e_0; the
           coordinate z has weight -omega),  H^1 = 0;
  m <  0 : H^0 = 0,  H^1 weights  w_0 + k*omega,  k = 1..|m|-1  (Cech).
For O(5), d = 1: {(5-k) lambda_i + k lambda_j, k = 0..5}  (as in the task).
For f^*TP^4 = O(2) + O(1)^3 on L_ij: O(2d): (d-k)omega, k = 0..2d, the ZERO
weight (k = d) being removed (it is the infinitesimal reparametrisation of the
cover, cancelled against the weight-0 part of Ext^0(Omega_C(D),O), i.e. the
C^* of the edge); O(d)_m (m != i,j): lambda_i - lambda_m - k*omega, k = 0..d.

Graber-Pandharipande.  [M]^vir = sum_F [F]^vir / e(N_F^vir), with
  T^vir = Ext^1(Omega_C(D),O) - Ext^0(Omega_C(D),O) + H^0(f^*TX) - H^1(f^*TX)
and the normalisation sequence for a fixed-locus curve C = (contracted
components C_v, genus g_v, at p_{i(v)}) U (edge covers C_e), glued at nodes.
Per piece, with E = Hodge bundle of C_v and Lambda^vee_h(w) := e(E^vee (x) w)
= sum_k (-1)^k lambda_k w^{h-k}:
  contracted vertex:  H^0 - H^1 (f^*TX) = T_{p_i} - E^vee (x) T_{p_i}
                      -> factor prod_{j != i} Lambda^vee(lambda_i - lambda_j) / e_i
                      twist O(5):  5lambda_i - E^vee (x) 5lambda_i
                      -> factor [5lambda_i / Lambda^vee(5lambda_i)]^{s}
  node at p_i:        -T_{p_i} and -(5lambda_i) in the K-classes
                      -> factor e_i * (5lambda_i)^{-s}
  node smoothing (Ext^1):  edge-vertex node: weight omega_F - psi_F -> 1/(omega_F - psi_F)
                           edge-edge node (2-valent genus-0 vertex): 1/(omega_1 + omega_2)
  Ext^0 moving part:  1-valent genus-0 vertex (the edge point over p_i is not
                      a node, the vector field d/dz at it has weight omega_F)
                      -> factor omega_F in the numerator
  edge:               1/d (cover automorphisms, part of 1/|Aut|)
                      x [e(H^0 - H^1)(f^*O(5))]^{s} / e(H^0(f^*TX))_{moving}.
Checked on Mbar_{0,0}(P^1,1) = pt: edge 1/(omega*(-omega)) x vertices omega*(-omega) = 1.

Mumford identity used to move the numerator twist into Hodge form:
  1/Lambda^vee_h(w) = (-1)^h Lambda^vee_h(-w) / w^{2h}
(equivalent to c(E)c(E^vee) = 1, i.e. ch_{2k}(E) = 0 for k >= 1; it is
verified EXACTLY in this file at the level of polynomials in the ch_k).
With Hodge_v := prod_a Lambda^vee(a)/a^h over a weight multiset v, a stable
vertex (h, n flags) of colour i becomes
  V = W_i^{h+n-1} Int_{Mbar_{h,n}} Hodge_{wt_i} prod_F 1/(omega_F - psi_F),
  wt_i = {lambda_i - lambda_j (j != i)} U {-5 lambda_i}   (numerator twist
         enters with a MINUS sign; a denominator twist would enter as +w),
  W_i  = e_i / (5 lambda_i)   (= e(T_{p_i}) x e(denominator-twist fibres)
         / e(numerator-twist fibres);  = W_i^0 of the plan; = zeta^{3i} at
         the zeta point with lambda = 1),
with the unstable conventions  (0,1): omega_F,  (0,2): W_i/(omega_1+omega_2)
(these are exactly the GP factors derived above: no node for a 1-valent
genus-0 vertex, one node for a 2-valent one).
Derivation of the exponent: prefactor e_i^{n-1} (nodes and T_{p_i})
x (5lambda_i)^{1-n} (nodes and the H^0 term of the twist) x (-1)^h
(5lambda_i)^{-2h} (Mumford) x (e_i (-5lambda_i))^h (Hodge normalisation)
= e_i^{h+n-1} (5lambda_i)^{1-n-h}.

CONSTANT MAPS (d = 0), verified by hand from the obstruction theory: the fixed
locus at p_i is Mbar_g x {p_i}; T^vir = T Mbar_g (fixed) + [T_{p_i} - E^vee (x)
T_{p_i}] (moving, no edges, no nodes), E|_F = 5lambda_i - E^vee (x) 5lambda_i,
so the contribution is
  Int_{Mbar_g} [5lambda_i/Lambda^vee(5lambda_i)] [prod_{j!=i} Lambda^vee(lambda_i-lambda_j)/e_i]
  = (5lambda_i/e_i) (-1)^g (5lambda_i)^{-2g} Int Lambda^vee(-5lambda_i) prod_j Lambda^vee(lambda_i-lambda_j)
  = (5lambda_i/e_i) (-1)^g (5lambda_i)^{-2g} e_i^g (-5lambda_i)^g Int Hodge_{wt_i}
  = (e_i/(5lambda_i))^{g-1} Int_{Mbar_g} Hodge_{wt_i}  =  (W_i^0)^{g-1} Int Hodge_{wt_i}.
Hence F_g(0) = sum_i (W_i^0)^{g-1} Int_{Mbar_g} Hodge_{wt_i}.  At the zeta
point p_k(wt_i) = zeta^{-ik} p_k(wt_0) so Int Hodge_{wt_i} = zeta^{-3i(g-1)}
Int Hodge_{wt_0} and F_g(0) = 5 Int Hodge_{wt_0}, a rational number.

Hodge integrals.  prod_a Lambda^vee_h(a) is expanded as a polynomial in the
ch_k(E) (k odd; ch_{2k} := 0) via c(E) = exp(sum_{k odd} (k-1)! ch_k t^k),
c_k(E^vee) = (-1)^k c_k(E); each monomial prod psi^a prod ch_k is integrated by
Mumford's GRR formula on Mbar_{g,n} (Faber-Pandharipande form)
  ch_{2l-1} = B_{2l}/(2l)! [ kappa_{2l-1} - sum_i psi_i^{2l-1}
              + 1/2 sum_{p+q = 2l-2} (-1)^p ( iota_*(psi_.^p psi_*^q) over the
                irreducible boundary and over ORDERED splittings (g1,I|g2,J) ) ],
kappa_j = pi_*(psi_{n+1}^{j+1}) (Arbarello-Cornalba), with the string, dilaton
and DVV recursions for the psi-numbers.  Bernoulli numbers: B_2 = 1/6, B_4 = -1/30.

Zeta point.  lambda_i = zeta^i (lambda = 1 by degree-0 homogeneity, which is
also checked by rescaling a generic sample), computed in Q(zeta) =
Q[x]/(x^4+x^3+x^2+x+1) exactly; inverses via the Galois norm.  Every N_{g,d}
is symmetric in the lambda_i, hence Galois-invariant, hence in Q: checked.

Weight choices (must avoid lambda_i = lambda_j, lambda_i = 0 and, for d = 2,
lambda_i + lambda_j = 2 lambda_k; any accidental zero raises ZeroDivisionError):
  sample A = (1, 3, 7, 12, 20),  sample B = (2/3, -5, 4, 15/2, 9),  3*A,  zeta.
"""
import sys, os, math, itertools, time
from fractions import Fraction as Fr
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
_LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    sys.stdout.flush()
    _LOG.append(s)


def prod(xs, one=1):
    r = one
    for x in xs:
        r = r * x
    return r


# ----------------------------------------------------------------------
# Bernoulli numbers (B_1 = -1/2 convention; only even ones are used)
# ----------------------------------------------------------------------
def bernoulli(nmax):
    B = [Fr(1)]
    for m in range(1, nmax + 1):
        B.append(-sum(math.comb(m + 1, k) * B[k] for k in range(m)) / Fr(m + 1))
    return B


BERN = bernoulli(24)


def dfact(n):
    """double factorial n!! for odd n >= -1 ((-1)!! = 1)."""
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def splits(t):
    """all ORDERED splittings (I, J) of the tuple t (by position)."""
    t = tuple(t)
    n = len(t)
    for mask in range(1 << n):
        yield (tuple(t[i] for i in range(n) if mask >> i & 1),
               tuple(t[i] for i in range(n) if not mask >> i & 1))


# ----------------------------------------------------------------------
# Witten-Kontsevich numbers <tau_{d_1} ... tau_{d_n}>_g
# ----------------------------------------------------------------------
@lru_cache(maxsize=None)
def wk(g, d):
    """d: sorted tuple (ascending).  String/dilaton/DVV; base <tau_1>_1 = 1/24."""
    n = len(d)
    if g < 0 or 2 * g - 2 + n <= 0:
        return Fr(0)
    if sum(d) != 3 * g - 3 + n:
        return Fr(0)
    if g == 0:
        r = Fr(math.factorial(n - 3))
        for x in d:
            r /= math.factorial(x)
        return r
    if g == 1 and n == 1:
        return Fr(1, 24)
    if d[0] == 0:                                   # string equation
        rest = d[1:]
        s = Fr(0)
        for j in range(len(rest)):
            if rest[j] > 0:
                s += wk(g, tuple(sorted(rest[:j] + (rest[j] - 1,) + rest[j + 1:])))
        return s
    if d[0] == 1:                                   # dilaton equation
        return (2 * g - 2 + n - 1) * wk(g, d[1:])
    # DVV on the largest insertion tau_{k+1}
    k = d[-1] - 1
    rest = d[:-1]
    s = Fr(0)
    for j in range(len(rest)):
        s += Fr(dfact(2 * k + 2 * rest[j] + 1), dfact(2 * rest[j] - 1)) * \
            wk(g, tuple(sorted(rest[:j] + (rest[j] + k,) + rest[j + 1:])))
    for r in range(k):
        q = k - 1 - r
        c = Fr(dfact(2 * r + 1) * dfact(2 * q + 1), 2)
        t = wk(g - 1, tuple(sorted(rest + (r, q))))
        for I, J in splits(rest):
            for g1 in range(g + 1):
                x = wk(g1, tuple(sorted(I + (r,))))
                if x:
                    t += x * wk(g - g1, tuple(sorted(J + (q,))))
        s += c * t
    return s / dfact(2 * k + 3)


# ----------------------------------------------------------------------
# Int_{Mbar_{g,n}} prod psi_i^{a_i} prod_{k in K} ch_k(E)   (Mumford GRR)
# ----------------------------------------------------------------------
@lru_cache(maxsize=None)
def chint(g, a, K):
    """a: sorted tuple of psi exponents; K: sorted tuple of (odd) ch indices."""
    n = len(a)
    if g < 0 or 2 * g - 2 + n <= 0:
        return Fr(0)
    if sum(a) + sum(K) != 3 * g - 3 + n:
        return Fr(0)
    if any(k % 2 == 0 for k in K):
        return Fr(0)
    if not K:
        return wk(g, a)
    if g == 0:
        return Fr(0)
    k = K[-1]
    KK = K[:-1]
    l = (k + 1) // 2
    tot = chint(g, tuple(sorted(a + (k + 1,))), KK)          # kappa_k term
    for i in range(n):                                        # -sum psi_i^k
        tot -= chint(g, tuple(sorted(a[:i] + (a[i] + k,) + a[i + 1:])), KK)
    for p in range(k):                                        # boundary terms
        q = k - 1 - p
        b = chint(g - 1, tuple(sorted(a + (p, q))), KK)      # irreducible
        for I, J in splits(a):                                # reducible, ordered
            for K1, K2 in splits(KK):
                for g1 in range(g + 1):
                    x = chint(g1, tuple(sorted(I + (p,))), tuple(sorted(K1)))
                    if x:
                        b += x * chint(g - g1, tuple(sorted(J + (q,))), tuple(sorted(K2)))
        tot += Fr((-1) ** p, 2) * b
    return tot * BERN[2 * l] / math.factorial(2 * l)


# ----------------------------------------------------------------------
# Q(zeta_5) = Q[x]/(x^4+x^3+x^2+x+1)
# ----------------------------------------------------------------------
class Qz:
    """element c0 + c1 x + c2 x^2 + c3 x^3, x = zeta_5; coefficients Fractions."""
    __slots__ = ("c",)

    def __init__(self, c):
        if isinstance(c, (int, Fr)):
            c = (Fr(c), Fr(0), Fr(0), Fr(0))
        self.c = tuple(Fr(x) for x in c)
        assert len(self.c) == 4

    @staticmethod
    def _co(x):
        if isinstance(x, Qz):
            return x
        if isinstance(x, (int, Fr)):
            return Qz(x)
        return NotImplemented

    def __add__(self, o):
        o = Qz._co(o)
        if o is NotImplemented:
            return o
        return Qz(tuple(a + b for a, b in zip(self.c, o.c)))
    __radd__ = __add__

    def __neg__(self):
        return Qz(tuple(-a for a in self.c))

    def __sub__(self, o):
        o = Qz._co(o)
        if o is NotImplemented:
            return o
        return self + (-o)

    def __rsub__(self, o):
        return Qz._co(o) - self

    def __mul__(self, o):
        o = Qz._co(o)
        if o is NotImplemented:
            return o
        r = [Fr(0)] * 5
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    if b:
                        r[(i + j) % 5] += a * b           # x^5 = 1
        c4 = r[4]                                          # x^4 = -(1+x+x^2+x^3)
        return Qz((r[0] - c4, r[1] - c4, r[2] - c4, r[3] - c4))
    __rmul__ = __mul__

    def galois(self, k):
        """sigma_k: x -> x^k (k = 1..4)."""
        r = [Fr(0)] * 5
        for m, a in enumerate(self.c):
            r[(k * m) % 5] += a
        c4 = r[4]
        return Qz((r[0] - c4, r[1] - c4, r[2] - c4, r[3] - c4))

    def inverse(self):
        conj = self.galois(2) * self.galois(3) * self.galois(4)
        nrm = self * conj
        assert nrm.c[1] == 0 and nrm.c[2] == 0 and nrm.c[3] == 0, "norm not rational?"
        if nrm.c[0] == 0:
            raise ZeroDivisionError("division by zero in Q(zeta)")
        return conj * (1 / nrm.c[0])

    def __truediv__(self, o):
        o = Qz._co(o)
        if o is NotImplemented:
            return o
        return self * o.inverse()

    def __rtruediv__(self, o):
        return Qz._co(o) * self.inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse() ** (-n)
        r = Qz(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def __eq__(self, o):
        o = Qz._co(o)
        if o is NotImplemented:
            return False
        return self.c == o.c

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(self.c)

    def __bool__(self):
        return any(self.c)

    def is_rational(self):
        return self.c[1] == 0 and self.c[2] == 0 and self.c[3] == 0

    def __repr__(self):
        if self.is_rational():
            return str(self.c[0])
        return "(" + " + ".join(f"({a})z^{m}" if m else f"({a})" for m, a in enumerate(self.c) if a) + ")"


ZETA = Qz((0, 1, 0, 0))


# ----------------------------------------------------------------------
# polynomials in ch_1, ch_3, ch_5, ... with coefficients in a field
# ----------------------------------------------------------------------
class ChPoly:
    """dict {sorted tuple of odd k (multiset of ch_k factors): coefficient}."""

    def __init__(self, d=None):
        self.d = {}
        if d:
            for k, v in d.items():
                if v != 0:
                    self.d[k] = v

    @staticmethod
    def const(c):
        return ChPoly({(): c})

    def __add__(self, o):
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = r.get(k, 0) + v
        return ChPoly(r)

    def scale(self, c):
        return ChPoly({k: v * c for k, v in self.d.items()})

    def mul(self, o, maxdeg=None):
        r = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                k = tuple(sorted(k1 + k2))
                if maxdeg is not None and sum(k) > maxdeg:
                    continue
                r[k] = r.get(k, 0) + v1 * v2
        return ChPoly(r)

    def __eq__(self, o):
        return self.d == o.d

    def integrate(self, g, psis):
        """Int_{Mbar_{g,n}} (self) prod psi_i^{psis_i}."""
        psis = tuple(sorted(psis))
        tot = 0
        for K, v in self.d.items():
            x = chint(g, psis, K)
            if x:
                tot = tot + v * x
        return tot


def odd_partitions(n, mx=None):
    """partitions of n into odd parts, non-increasing."""
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    p = min(n, mx)
    if p % 2 == 0:
        p -= 1
    while p >= 1:
        for rest in odd_partitions(n - p, p):
            yield (p,) + rest
        p -= 2


@lru_cache(maxsize=None)
def chern_E(k):
    """c_k(E) as ChPoly: coefficient of t^k in exp(sum_{p odd} (p-1)! ch_p t^p)."""
    r = {}
    for part in odd_partitions(k):
        c = Fr(1)
        mult = {}
        for p in part:
            mult[p] = mult.get(p, 0) + 1
            c *= math.factorial(p - 1)
        for m in mult.values():
            c /= math.factorial(m)
        key = tuple(sorted(part))
        r[key] = r.get(key, 0) + c
    return ChPoly(r)


def lambda_vee(h, w):
    """Lambda^vee_h(w) = sum_k (-1)^k lambda_k w^{h-k} = sum_k c_k(E^vee) w^{h-k}."""
    P = ChPoly()
    for k in range(h + 1):
        P = P + chern_E(k).scale(((-1) ** k) * (w ** (h - k)))
    return P


def lambda_monomial(h, ks):
    """prod_k lambda_k as a ChPoly (for validation)."""
    P = ChPoly.const(Fr(1))
    for k in ks:
        P = P.mul(chern_E(k))
    return P


def hodge_poly(h, wts, maxdeg):
    """Hodge_v = prod_a Lambda^vee_h(a) / a^h, truncated above degree maxdeg."""
    P = ChPoly.const(Fr(1))
    for a in wts:
        P = P.mul(lambda_vee(h, a), maxdeg)
    return P.scale(1 / prod([a ** h for a in wts], Fr(1)))


# ----------------------------------------------------------------------
# Graber-Pandharipande localisation on P^r with line-bundle twists
# ----------------------------------------------------------------------
class Loc:
    """lam: weights lambda_0..lambda_r (Fractions or Qz).
    twists: list of (a, eps, s): line bundle O(a) (x) C_eps on P^r, fibre weight
    a*lambda_i + eps at p_i; s = +1: e(R pi_*) in the numerator; s = -1: e(-R pi_*)."""

    def __init__(self, lam, twists):
        self.lam = list(lam)
        self.r = len(lam) - 1
        self.tw = list(twists)
        R = range(self.r + 1)
        self.eT = [prod([lam[i] - lam[j] for j in R if j != i]) for i in R]
        self.fw = [[a * lam[i] + eps for (a, eps, s) in twists] for i in R]
        self.W = []
        self.hw = []
        for i in R:
            Wi = self.eT[i]
            hw = [lam[i] - lam[j] for j in R if j != i]
            for (a, eps, s), w in zip(twists, self.fw[i]):
                if s > 0:
                    Wi = Wi / w
                    hw.append(-w)
                else:
                    Wi = Wi * w
                    hw.append(w)
            self.W.append(Wi)
            self.hw.append(hw)
        self._hp = {}
        self._hint = {}
        self._vert = {}

    # Hodge data at p_i
    def hodge(self, i, h, nmax=3):
        key = (i, h)
        if key not in self._hp:
            self._hp[key] = hodge_poly(h, self.hw[i], 3 * h - 3 + nmax)
        return self._hp[key]

    def hint(self, i, h, ks):
        """Int_{Mbar_{h,n}} Hodge_{wt_i} prod psi^{k}."""
        key = (i, h, tuple(sorted(ks)))
        if key not in self._hint:
            self._hint[key] = self.hodge(i, h).integrate(h, ks)
        return self._hint[key]

    def vertex(self, i, h, oms):
        """vertex factor of colour i, genus h, flag weights oms (unstable conventions included)."""
        n = len(oms)
        key = (i, h, tuple(oms))
        if key in self._vert:
            return self._vert[key]
        if h == 0 and n == 1:
            r = oms[0]
        elif h == 0 and n == 2:
            r = self.W[i] / (oms[0] + oms[1])
        else:
            dim = 3 * h - 3 + n
            tot = 0
            for ks in itertools.product(range(dim + 1), repeat=n):
                if sum(ks) > dim:
                    continue
                hv = self.hint(i, h, ks)
                if hv == 0:
                    continue
                c = hv
                for o, k in zip(oms, ks):
                    c = c / o ** (k + 1)
                tot = tot + c
            r = tot * self.W[i] ** (h + n - 1)
        self._vert[key] = r
        return r

    def edge(self, i, j, d):
        """degree-d cover of L_ij. Returns (factor, omega at p_i); omega at p_j is -omega."""
        lam = self.lam
        om = (lam[i] - lam[j]) / d
        f = Fr(1, d)
        for k in range(2 * d + 1):                     # O(2d): (d-k) omega, zero removed
            if k != d:
                f = f / ((d - k) * om)
        for m in range(self.r + 1):                    # O(d)_m
            if m in (i, j):
                continue
            for k in range(d + 1):
                f = f / (lam[i] - lam[m] - k * om)
        for (a, eps, s) in self.tw:
            w0 = a * lam[i] + eps
            m = a * d
            if m >= 0:
                val = prod([w0 - k * om for k in range(m + 1)], Fr(1))
            else:
                val = 1 / prod([w0 + k * om for k in range(1, -m)], Fr(1))
            f = f * val if s > 0 else f / val
        return f, om

    def constant_maps(self, g):
        assert g >= 2
        return sum((self.W[i] ** (g - 1) * self.hint(i, g, ()) for i in range(self.r + 1)), 0)

    def gw(self, g, d):
        """N_{g,d}: sum over connected fixed-locus graphs with total edge degree d."""
        if d == 0:
            return self.constant_maps(g)
        tot = 0
        for nv in range(2, d + 2):
            for edges in edge_multisets(nv, d):
                if not connected(nv, edges):
                    continue
                h1 = len(edges) - nv + 1
                if h1 > g:
                    continue
                wgt = Fr(1, math.factorial(nv))
                cnt = {}
                for e in edges:
                    cnt[e] = cnt.get(e, 0) + 1
                for m in cnt.values():
                    wgt /= math.factorial(m)
                for gens in compositions(g - h1, nv):
                    for cols in itertools.product(range(self.r + 1), repeat=nv):
                        if any(cols[u] == cols[v] for u, v, _ in edges):
                            continue
                        tot = tot + self.config(cols, gens, edges) * wgt
        return tot

    def config(self, cols, gens, edges):
        nv = len(cols)
        oms = [[] for _ in range(nv)]
        val = Fr(1)
        for u, v, dg in edges:
            f, om = self.edge(cols[u], cols[v], dg)
            val = val * f
            oms[u].append(om)
            oms[v].append(-om)
        for x in range(nv):
            val = val * self.vertex(cols[x], gens[x], oms[x])
        return val


def partitions_upto(rem):
    def parts(n, mx):
        if n == 0:
            yield ()
            return
        for p in range(min(n, mx), 0, -1):
            for r in parts(n - p, p):
                yield (p,) + r
    for n in range(rem + 1):
        yield from parts(n, n)


def edge_multisets(nv, d):
    """edge multisets (u<v, degree) on labelled vertices 0..nv-1 with total degree d."""
    pairs = [(u, v) for u in range(nv) for v in range(u + 1, nv)]

    def rec(idx, rem):
        if rem == 0:
            yield ()
            return
        if idx == len(pairs):
            return
        for part in partitions_upto(rem):
            for r in rec(idx + 1, rem - sum(part)):
                yield tuple((pairs[idx][0], pairs[idx][1], q) for q in part) + r
    return list(rec(0, d))


def connected(nv, edges):
    par = list(range(nv))

    def f(x):
        while par[x] != x:
            x = par[x]
        return x
    for u, v, _ in edges:
        par[f(u)] = f(v)
    return len({f(i) for i in range(nv)}) == 1


def compositions(n, k):
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for r in compositions(n - i, k - 1):
            yield (i,) + r


# ----------------------------------------------------------------------
# validation of the intersection-number engine
# ----------------------------------------------------------------------
def check(name, got, exp):
    ok = (got == exp)
    say(f"  [{'OK' if ok else 'FAIL'}] {name}: got {got}, expected {exp}")
    return ok


def validate_engine():
    say("== Engine validation (Witten-Kontsevich, Mumford ch_k, lambda classes) ==")
    ok = True
    ok &= check("<tau_0^3>_0", wk(0, (0, 0, 0)), Fr(1))
    ok &= check("<tau_1>_1", wk(1, (1,)), Fr(1, 24))
    ok &= check("<tau_4>_2 = Int_{Mbar_{2,1}} psi^4", wk(2, (4,)), Fr(1, 1152))
    ok &= check("<tau_2 tau_3>_2", wk(2, (2, 3)), Fr(29, 5760))
    ok &= check("<tau_2^3>_2", wk(2, (2, 2, 2)), Fr(7, 240))
    ok &= check("<tau_7>_3", wk(3, (7,)), Fr(1, 82944))
    ok &= check("<tau_2 tau_6>_3", wk(3, (2, 6)), Fr(77, 414720))
    ok &= check("<tau_3 tau_5>_3", wk(3, (3, 5)), Fr(503, 1451520))
    ok &= check("<tau_4 tau_4>_3", wk(3, (4, 4)), Fr(607, 1451520))
    ok &= check("lambda_1 on Mbar_{1,1}", lambda_monomial(1, (1,)).integrate(1, (0,)), Fr(1, 24))
    ok &= check("Int_{Mbar_2} lambda_1^3", lambda_monomial(2, (1, 1, 1)).integrate(2, ()), Fr(1, 2880))
    ok &= check("Int_{Mbar_2} lambda_1 lambda_2", lambda_monomial(2, (1, 2)).integrate(2, ()), Fr(1, 5760))
    ok &= check("Int_{Mbar_{2,1}} lambda_1 lambda_2 psi (= 2 Int lambda_1 lambda_2, dilaton)",
                lambda_monomial(2, (1, 2)).integrate(2, (1,)), Fr(1, 2880))
    ok &= check("Int_{Mbar_{2,1}} lambda_2 psi^2 (lambda_g conjecture 7/5760)",
                lambda_monomial(2, (2,)).integrate(2, (2,)), Fr(7, 5760))
    ok &= check("Int_{Mbar_3} lambda_1 lambda_2 lambda_3", lambda_monomial(3, (1, 2, 3)).integrate(3, ()),
                Fr(1, 1451520))
    ok &= check("Int_{Mbar_3} lambda_2^3", lambda_monomial(3, (2, 2, 2)).integrate(3, ()), Fr(1, 725760))
    ok &= check("Int_{Mbar_3} lambda_1^6", lambda_monomial(3, (1,) * 6).integrate(3, ()), Fr(1, 90720))
    ok &= check("Int_{Mbar_{3,1}} lambda_3 psi^4 (lambda_g conjecture 31/967680)",
                lambda_monomial(3, (3,)).integrate(3, (4,)), Fr(31, 967680))
    ok &= check("Int_{Mbar_{1,1}} Lambda^vee(w) = -1/24 (w=7)", lambda_vee(1, Fr(7)).integrate(1, (0,)), Fr(-1, 24))
    # Mumford identity Lambda^vee(w) Lambda^vee(-w) = (-1)^h w^{2h}.  As polynomials in the ch_k it
    # holds only modulo c_k(E) = 0 for k > h, so it is tested under the integral, against psi-monomials
    # and against another Lambda^vee, on Mbar_{h,n}.
    for h in (1, 2, 3):
        for w in (Fr(3), Fr(-5, 7)):
            P = lambda_vee(h, w).mul(lambda_vee(h, -w))
            Q = lambda_vee(h, Fr(11, 2))
            for n in (1, 2):
                dim = 3 * h - 3 + n
                for ks in itertools.product(range(dim + 1), repeat=n):
                    if sum(ks) > dim or sum(ks) < dim - h:
                        continue
                    lhs = P.mul(Q).integrate(h, ks)
                    rhs = Q.scale(((-1) ** h) * w ** (2 * h)).integrate(h, ks)
                    ok &= check(f"Int Lambda^vee_{h}(w)Lambda^vee_{h}(-w) Lambda^vee(11/2) psi^{ks} "
                                f"= (-1)^h w^(2h) Int Lambda^vee(11/2) psi^{ks}, w={w}", lhs, rhs)
    # exp form of Hodge_v vs product form (3 weights, g = 2, 3), again under the integral
    for h in (2, 3):
        wts = [Fr(2), Fr(-3), Fr(5, 7)]
        P1 = hodge_poly(h, wts, 3 * h)
        P2 = ChPoly()                    # exp(-sum_{k odd} (k-1)! p_k ch_k)
        for part in itertools.chain.from_iterable(odd_partitions(n) for n in range(3 * h + 1)):
            c = Fr(1)
            mult = {}
            for q in part:
                mult[q] = mult.get(q, 0) + 1
                c *= -math.factorial(q - 1) * sum(Fr(1) / x ** q for x in wts)
            for m in mult.values():
                c /= math.factorial(m)
            P2 = P2 + ChPoly({tuple(sorted(part)): c})
        for ks in ((), (0,), (1,), (3 * h - 2,), (1, 2), (0, 3)):
            ok &= check(f"Int Hodge_v psi^{ks} product form == exp(-sum (k-1)! p_k ch_k) form, h={h}",
                        P1.integrate(h, ks), P2.integrate(h, ks))
    say("  engine validation:", "ALL OK" if ok else "FAILURES")
    return ok


# ----------------------------------------------------------------------
# control: local P^1 = O(-1)+O(-1) -> P^1 with a CY torus (known GV answers)
# ----------------------------------------------------------------------
def control_local_P1():
    say("== Control: X = O(-1)+O(-1) -> P^1, CY torus (eps_1+eps_2 = lambda_0+lambda_1), "
        "e(-R pi_* f^*N) twist ==")
    say("   expected from GV (F = sum_k (1/k)(2 sin(k g_s/2))^-2 Q^k): N_{0,1}=1, N_{1,1}=1/12, "
        "N_{2,1}=1/240, N_{0,2}=1/8, N_{1,2}=1/24, N_{2,2}=1/120; constant maps "
        "(-1)^g chi/2 Int lambda_{g-1}^3 with chi=2: N_{2,0}=1/2880, N_{3,0}=-1/725760")
    exp = {(0, 1): Fr(1), (1, 1): Fr(1, 12), (2, 1): Fr(1, 240), (0, 2): Fr(1, 8), (1, 2): Fr(1, 24),
           (2, 2): Fr(1, 120), (2, 0): Fr(1, 2880), (3, 0): Fr(-1, 725760)}
    ok = True
    for lam, e1 in (((Fr(2), Fr(-3, 7)), Fr(5, 3)), ((Fr(1, 3), Fr(11, 2)), Fr(-7))):
        e2 = lam[0] + lam[1] - e1
        L = Loc(lam, [(-1, e1, -1), (-1, e2, -1)])
        say(f"   weights lambda={tuple(map(str, lam))}, eps=({e1},{e2})")
        for (g, d), v in exp.items():
            ok &= check(f"local P1 N_{{{g},{d}}}", L.gw(g, d), v)
    say("  local P^1 control:", "ALL OK" if ok else "FAILURES")
    return ok


# ----------------------------------------------------------------------
# the formal quintic
# ----------------------------------------------------------------------
def quintic_table():
    say("== Formal quintic: T=(C^*)^5 on P^4, e(R pi_* f^*O(5)) in the numerator ==")
    A = [Fr(1), Fr(3), Fr(7), Fr(12), Fr(20)]
    B = [Fr(2, 3), Fr(-5), Fr(4), Fr(15, 2), Fr(9)]
    A3 = [3 * x for x in A]
    Z = [ZETA ** i for i in range(5)]
    samples = [("A=(1,3,7,12,20)", A), ("B=(2/3,-5,4,15/2,9)", B), ("3A (homogeneity)", A3),
               ("zeta: lambda_i = zeta^i", Z)]
    items = [(0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
    results = {}
    ok = True
    for name, lam in samples:
        t0 = time.time()
        L = Loc(lam, [(5, 0, +1)])
        say(f"-- sample {name}")
        say("   e_i = e(T_{p_i}):", [str(x) for x in L.eT])
        say("   W_i^0 = e_i/(5 lambda_i):", [str(x) for x in L.W])
        say("   Hodge weights wt_0:", [str(x) for x in L.hw[0]])
        for (g, d) in items:
            v = L.gw(g, d)
            results[(name, g, d)] = v
            say(f"   N_{{{g},{d}}} = {v}")
        for g in (2, 3):
            v = L.constant_maps(g)
            results[(name, "F0", g)] = v
            say(f"   F_{g}(0) = sum_i (W_i^0)^(g-1) Int_Mbar_{g} Hodge_wt_i = {v}")
        say(f"   ({time.time() - t0:.1f}s)")
    say("-- checks")
    for name, _ in samples:
        ok &= check(f"N_{{0,1}} = 2875 [{name}]", results[(name, 0, 1)], Fr(2875))
        ok &= check(f"N_{{0,2}} = 4876875/8 [{name}]", results[(name, 0, 2)], Fr(4876875, 8))
    for key in [(1, 1), (1, 2), (2, 1), (2, 2), ("F0", 2), ("F0", 3)]:
        vA = results[("A=(1,3,7,12,20)",) + key]
        vB = results[("B=(2/3,-5,4,15/2,9)",) + key]
        vA3 = results[("3A (homogeneity)",) + key]
        lab = f"N_{{{key[0]},{key[1]}}}" if key[0] != "F0" else f"F_{key[1]}(0)"
        ok &= check(f"{lab} degree-0 homogeneity (A vs 3A)", vA3 == vA, True)
        say(f"  [info] {lab} lambda-independent between samples A and B: {vA == vB}")
    say("  zeta point: rationality (LP remark; Galois invariance of the symmetric sum):")
    for key in items + [("F0", 2), ("F0", 3)]:
        v = results[("zeta: lambda_i = zeta^i",) + key]
        lab = f"N_{{{key[0]},{key[1]}}}" if key[0] != "F0" else f"F_{key[1]}(0)"
        ok &= check(f"{lab} at zeta in Q", v.is_rational(), True)
    # true-quintic constant maps for comparison
    l13 = lambda_monomial(2, (1, 1, 1)).integrate(2, ())
    l23 = lambda_monomial(3, (2, 2, 2)).integrate(3, ())
    say(f"  [info] true quintic constant maps (-1)^g chi/2 Int lambda_(g-1)^3, chi=-200: "
        f"g=2: {-100 * l13},  g=3: {100 * l23}")
    # formal zeta-point F_g(0) = 5 Int Hodge_{wt_0}
    Lz = Loc(Z, [(5, 0, +1)])
    for g in (2, 3):
        say(f"  [info] zeta point: 5*Int_Mbar_{g} Hodge_wt_0 = {5 * Lz.hint(0, g, ())}, "
            f"F_{g}(0) = {results[('zeta: lambda_i = zeta^i', 'F0', g)]}")
    # d = 1 breakdown by genus split (h at the i-end, g-h at the j-end), summed over the 10 lines,
    # and the single-line values i=0, j=1 (for the graph-by-graph comparison with the graph sum)
    say("-- d = 1 breakdown: edge p_i p_j (degree-1 cover), contracted genera (h_i, h_j); "
        "'all lines' = sum over ordered pairs i != j with weight 1/2")
    for name, lam in samples:
        L = Loc(lam, [(5, 0, +1)])
        for g in (1, 2):
            parts = []
            for hi in range(g + 1):
                hj = g - hi
                tot = 0
                for i in range(5):
                    for j in range(5):
                        if i != j:
                            tot = tot + L.config((i, j), (hi, hj), ((0, 1, 1),)) * Fr(1, 2)
                v01 = L.config((0, 1), (hi, hj), ((0, 1, 1),))
                parts.append(f"(h_i,h_j)=({hi},{hj}): all lines {tot}; line p_0p_1 {v01}")
            say(f"   [{name}] g={g}: " + " | ".join(parts))
    say("== SUMMARY TABLE (formal quintic) ==")
    hdr = f"{'quantity':<10}" + "".join(f"{n:<46}" for n, _ in samples)
    say(hdr)
    for key in items + [("F0", 2), ("F0", 3)]:
        lab = f"N_{{{key[0]},{key[1]}}}" if key[0] != "F0" else f"F_{key[1]}(0)"
        say(f"{lab:<10}" + "".join(f"{str(results[(n,) + key]):<46}" for n, _ in samples))
    say("  quintic checks:", "ALL OK" if ok else "FAILURES")
    return ok


def main():
    t0 = time.time()
    say(f"loc_check.py  --  {time.strftime('%Y-%m-%d %H:%M:%S')}  Python {sys.version.split()[0]}")
    ok1 = validate_engine()
    ok2 = control_local_P1()
    ok3 = quintic_table()
    say(f"total time {time.time() - t0:.1f}s; overall:", "ALL OK" if (ok1 and ok2 and ok3) else "FAILURES")
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    with open(os.path.join(HERE, "logs", "loc_check.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LOG) + "\n")


if __name__ == "__main__":
    main()
