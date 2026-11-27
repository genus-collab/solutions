# -*- coding: utf-8 -*-
"""fq_phi.py -- constant, symplectic, Z/5-covariant scalar factorisation of the formal-quintic R-matrix
(PLAN_scalar_2026-11-04.md):  R = R' phi_i(z),  phi_i(z) = phi_0(zeta^-i z),  phi_0(z) = exp(a1 z + a2 z^3)
(a3 and higher never enter at genus 2), with the vertex CohFT phi_i . 1 (PPZ action of the scalar on the trivial theory).

Shortcut for the vertex classes (no PPZ graph sum needed; Mumford/Teleman form):
    phi . 1 = exp( Sum_k (2k)!/B_{2k} a_k ch_{2k-1}(E) ) = exp(12 a1 ch_1 - 720 a2 ch_3 + ...),
so  I_phi(g; d) = Int_{Mbar_{g,n}} (phi.1) prod psi^{d_i} = Sum_{odd partitions K of D} hi(g,d,K) prod_k s_k^{m_k}/m_k!
with s_1 = 12 a1, s_3 = -720 a2 (fq_hodge.VertexTable with these s_k), D = 3g-3+n-Sum d; I_phi is homogeneous of a-weight
D (weight 1 for a1, 3 for a2).  At colour i the phases are s_k^{(i)} = zeta^{-ik} s_k.
Hodge points: a Hodge class Hodge_w = h_w . 1 corresponds to a1 = -(1/12) Sum_a 1/w_a, a2 = (1/360) Sum_a 1/w_a^3
(fq_hodge.hv conventions), e.g. w = {1}: (-1/12, 1/360) = frame F4; w = wt_0: (-3/20, -7/2500) = true class F0.
Rows: R' = R phi_i^-1, coefficients  zeta^{-in} c_n  with c_n = [z^n] exp(-a1 z - a2 z^3)  (c_n has a-weight n).

API:  phi_inv_coeffs(R, a1, a2, i, Nz) -> list of z-coefficients of phi_i(z)^-1;
      PhiVertexTable(R, a1, a2, i=0)   -> VertexTable subclass (.val(g, d)), asserts k <= 3;
      I_phi_poly(g, d)                 -> {(e1, e2): Fraction}  the polynomial Sum c a1^e1 a2^e2 over Q;
      hodge_point(R, w)                -> (a1, a2) in R for the weight multiset w;
      phi_frame(rows, a1, a2, name)    -> gsum_q.Frame with rows R phi_i^-1 and vertex tables phi_i . 1.
Self-tests:  py -I fq_phi.py
"""
import math, os, sys, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_hodge import VertexTable, hi, partitions_odd, dmax, series_exp, hv, wt, wk


def phi_inv_coeffs(R, a1, a2, i, Nz):
    """[z^0..z^Nz] of phi_i(z)^-1 = exp(-a1 zeta^-i z - a2 zeta^-3i z^3)"""
    a1 = R.of(a1); a2 = R.of(a2)
    ex = [R.zero] * (Nz + 1)
    if Nz >= 1: ex[1] = R.neg(R.mul(a1, R.z(-i)))
    if Nz >= 3: ex[3] = R.neg(R.mul(a2, R.z(-3 * i)))
    return series_exp(R, ex, Nz)


class PhiVertexTable(VertexTable):
    """Int_{Mbar_{g,n}} (phi_i . 1) prod psi^d (mod p): VertexTable with s_1 = 12 a1 zeta^-i, s_3 = -720 a2 zeta^-3i."""
    def __init__(self, R, a1, a2, i=0):
        VertexTable.__init__(self, R, None)
        self.w = ('phi', i)                       # not None: the general (Hodge-polynomial) branch of val()
        a1 = R.of(a1); a2 = R.of(a2)
        self.a1 = a1; self.a2 = a2; self.colour = i
        s1 = R.mul(R.of(12), R.mul(a1, R.z(-i)))
        s3 = R.mul(R.of(-720), R.mul(a2, R.z(-3 * i)))
        def s(k):
            if k == 1: return s1
            if k == 3: return s3
            raise AssertionError("PhiVertexTable: s_%d requested (only s_1, s_3 are defined: genus <= 2)" % k)
        self.s = s


def I_phi_poly(g, d):
    """{(e1, e2): Fraction}: I_phi(g; d) as a polynomial in (a1, a2) over Q (s_1 = 12 a1, s_3 = -720 a2)"""
    d = tuple(sorted(d)); n = len(d); D = 3 * g - 3 + n - sum(d)
    out = {}
    if D < 0 or 2 * g - 2 + n <= 0 or D > dmax(g): return out
    for K in partitions_odd(D, max(2 * g - 1, 1)):
        x = hi(g, d, K)
        if x == 0: continue
        m1 = K.count(1); m3 = K.count(3)
        assert m1 + m3 == len(K), "odd part > 3 in a partition (genus > 2?): %s" % (K,)
        c = x * Fr(12) ** m1 / math.factorial(m1) * Fr(-720) ** m3 / math.factorial(m3)
        key = (m1, m3)
        out[key] = out.get(key, Fr(0)) + c
    return {k: v for k, v in out.items() if v != 0}


def poly_fmt(P):
    if not P: return "0"
    parts = []
    for (e1, e2), c in sorted(P.items()):
        mon = "".join([("a1^%d" % e1 if e1 > 1 else "a1" if e1 == 1 else ""), ("a2^%d" % e2 if e2 > 1 else "a2" if e2 == 1 else "")])
        parts.append("(%s)%s" % (c, " " + mon if mon else ""))
    return " + ".join(parts)


def hodge_point(R, w):
    """(a1, a2) of the Hodge class of the weight multiset w: a1 = -(1/12) Sum 1/w, a2 = (1/360) Sum 1/w^3"""
    p1 = R.zero; p3 = R.zero
    for x in w:
        xi = R.inv(R.of(x))
        p1 = R.add(p1, xi); p3 = R.add(p3, R.pow(xi, 3))
    return R.mul(R.of(Fr(-1, 12)), p1), R.mul(R.of(Fr(1, 360)), p3)


def phi_frame(rows, a1, a2, name):
    """gsum_q.Frame: rows R phi_i^-1 (colour-dependent coefficients zeta^-in c_n), vertex classes phi_i . 1"""
    from gsum_q import Frame
    R = rows.R
    hfac = {i: phi_inv_coeffs(R, a1, a2, i, rows.Nz) for i in range(5)}
    F = Frame(rows, name, hfac, None)
    F.vt = {i: PhiVertexTable(R, a1, a2, i) for i in range(5)}
    F.phi_point = (R.of(a1), R.of(a2))
    return F


# ------------------------------------------------------------------ self-tests
def _selftest():
    from fq_fp import FpRing, P1, P2
    R = FpRing(P2)
    # (1) polynomials over Q
    assert I_phi_poly(1, (0,)) == {(1, 0): Fr(1, 2)}, I_phi_poly(1, (0,))          # I(1;0) = a1/2
    assert I_phi_poly(1, (1,)) == {(0, 0): Fr(1, 24)}                              # I(1;1) = 1/24
    I2 = I_phi_poly(2, ())
    assert I2 == {(3, 0): Fr(1, 10), (0, 1): Fr(1, 48)}, I2                         # I(2;empty) = a1^3/10 + a2/48
    assert I_phi_poly(2, (4,)) == {(0, 0): Fr(1, 1152)}                            # I(2;4) = 1/1152
    print("I(1;0) =", poly_fmt(I_phi_poly(1, (0,))), "; I(2;empty) =", poly_fmt(I2), "; I(2;4) =", poly_fmt(I_phi_poly(2, (4,))))
    # homogeneity of a-weight D
    for g in (1, 2):
        for n in range(0, 4):
            for d in set(tuple(sorted(t)) for t in __import__('itertools').product(range(4), repeat=n)):
                P = I_phi_poly(g, d); D = 3 * g - 3 + n - sum(d)
                assert all(e1 + 3 * e2 == D for (e1, e2) in P), (g, d, P)
    # values at the two Hodge points
    ev = lambda P, a1, a2: sum((c * a1 ** e1 * a2 ** e2 for (e1, e2), c in P.items()), Fr(0))
    assert ev(I2, Fr(-1, 12), Fr(1, 360)) == 0                                     # Int_{M_2} Lambda^vee(1) = 0
    assert ev(I2, Fr(-3, 20), Fr(-7, 2500)) == Fr(-19, 48000)                      # true class: F_2(0)/5
    # (2) the two points in F_p: w = {1} and w = wt_0 (p_1 = 9/5, p_3 = -126/125)
    lams = [R.z(k) for k in range(5)]
    assert hodge_point(R, [1]) == (R.of(Fr(-1, 12)), R.of(Fr(1, 360)))
    w0 = wt(R, lams, 0)
    p1 = sum(R.inv(x) for x in w0) % R.p; p3 = sum(R.pow(R.inv(x), 3) for x in w0) % R.p
    assert p1 == R.of(Fr(9, 5)) and p3 == R.of(Fr(-126, 125)), (R.fmt(p1), R.fmt(p3))
    assert hodge_point(R, w0) == (R.of(Fr(-3, 20)), R.of(Fr(-7, 2500)))
    for i in range(5):
        wi = wt(R, lams, i)
        a1, a2 = hodge_point(R, wi)
        # wt_i = zeta^i wt_0, so a_k^{(i)} = zeta^{-i(2k-1)} a_k: exactly phi_i(z) = phi_0(zeta^-i z)
        assert a1 == R.mul(R.of(Fr(-3, 20)), R.z(-i)) and a2 == R.mul(R.of(Fr(-7, 2500)), R.z(-3 * i))
    # (3) PhiVertexTable vs VertexTable at the Hodge points and at a random two-weight multiset, all colours
    rng = random.Random(20261104)
    tests = [([1], 'w={1}'), (w0, 'w=wt_0')]
    for _ in range(3):
        tests.append(([rng.randrange(1, R.p), rng.randrange(1, R.p)], 'random pair'))
    ds = [(g, d) for g in (0, 1, 2) for n in range(0, 4) for d in set(tuple(sorted(t)) for t in __import__('itertools').product(range(5), repeat=n))]
    for w, label in tests:
        a1, a2 = hodge_point(R, w)
        for i in range(5):
            wi = [R.mul(R.of(x), R.z(i)) for x in w]
            VT = VertexTable(R, wi); PT = PhiVertexTable(R, a1, a2, i)
            for (g, d) in ds:
                assert VT.val(g, d) == PT.val(g, d), (label, i, g, d)
    print("PhiVertexTable = VertexTable at", [t[1] for t in tests], "for all colours and (g, d) with g <= 2, n <= 3")
    # a1 = a2 = 0 gives the point class
    VT0 = VertexTable(R, None); PT0 = PhiVertexTable(R, 0, 0, 2)
    for (g, d) in ds: assert VT0.val(g, d) == PT0.val(g, d), (g, d)
    # (4) phi_inv_coeffs vs hv(..., -1) up to z^3 at the Hodge points, all colours (z^4.. differ: a3 absent)
    Nz = 6
    for i in range(5):
        c = phi_inv_coeffs(R, Fr(-1, 12), Fr(1, 360), i, Nz); h = hv(R, [lams[i]], Nz, -1)
        assert c[:4] == h[:4], (i, c, h)
        a1, a2 = hodge_point(R, w0)
        c = phi_inv_coeffs(R, a1, a2, i, Nz); h = hv(R, wt(R, lams, i), Nz, -1)
        assert c[:4] == h[:4], (i, c, h)
        c0 = phi_inv_coeffs(R, 7, 11, 0, Nz); ci = phi_inv_coeffs(R, 7, 11, i, Nz)
        assert all(ci[n] == R.mul(R.z(-i * n), c0[n]) for n in range(Nz + 1))      # covariance zeta^-in c_n
    c = phi_inv_coeffs(R, 7, 11, 0, 3)
    assert c == [1, R.of(-7), R.of(Fr(49, 2)), R.of(Fr(-343, 6) - 11)], c          # c_3 = -a1^3/6 - a2
    # (5) same at the second prime (spot)
    R1 = FpRing(P1)
    a1, a2 = hodge_point(R1, wt(R1, [R1.z(k) for k in range(5)], 0))
    assert (a1, a2) == (R1.of(Fr(-3, 20)), R1.of(Fr(-7, 2500)))
    assert PhiVertexTable(R1, a1, a2, 3).val(2, ()) == VertexTable(R1, wt(R1, [R1.z(k) for k in range(5)], 3)).val(2, ())
    print("fq_phi self-tests ok")


if __name__ == "__main__":
    _selftest()
