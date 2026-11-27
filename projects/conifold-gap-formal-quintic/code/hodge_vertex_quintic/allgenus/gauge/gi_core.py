# -*- coding: utf-8 -*-
"""gi_core.py -- common tools for the gauge identity (G) of GAUGE_IDENTITY.md (exact, over Q; stdlib only).

Objects (flat frame, index a in Z/5, columns = vectors Y = (Y_0..Y_4), th = z d/du):
  quintic system   th Y_a = b_a(u) Y_{a+1}          (ur_gauge.coefficients; conifold frame (alpha, beta))
  model system     z d/dv Y_b = bmod_b Y_{b+1},  bmod = (1, -1/v, 1, -1, -1), pulled back by v(u) = u tau^(-5/4)
  regular gauge    P(u, z) = Sum_k z^k P_k(u),  P_k = Sum_b p_{k,b}(u) E_{b+k,b}     (ur_gauge.gauge, read-only import)
  vanishing solution of the model (residue at the pole x = 0):
      Y^C = (I_0, -I_{-1}, -I_3, I_2, I_1),   I_m = Sum_{i>=0} (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!
  vanishing solution of the quintic: the unique solution in C((z))[[u]]^5 with value -e_1 at u = 0, built by iterated integrals
      with zero constants:  y^(0) = -e_1,  y^(n)_a = Int_0^u b_a y^(n-1)_{a+1}   (coefficient of z^-n).
A "z-Laurent vector" is a dict {m: [Ser]*5} (coefficient of z^m, Ser in u over Q); zero entries may be absent.
Nothing here writes a file.  No existing project file is modified (ur_gauge.Conifold is rebound IN MEMORY to select the frame).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
THEORY = os.path.join(PROJ, "natural_form", "theory")
sys.path.insert(0, PROJ); sys.path.insert(0, THEORY)
from fq_fp import QRing
from fq_series import Ser
import fq_conifold
import ur_gauge
Q = QRing()
MEXP = (0, -1, 3, 2, 1)          # Y_b = SGN_b I_{MEXP_b}
SGN = (1, -1, -1, 1, 1)


def gauge_frame(K, NU, alpha=0, beta=0, log=None):
    """ur_gauge.gauge in the conifold frame (alpha, beta): varpi_0 = f_0 + alpha g + beta f_2."""
    a = Fr(alpha); b = Fr(beta)
    orig = ur_gauge.Conifold
    ur_gauge.Conifold = lambda R, Nu=40: fq_conifold.Conifold(R, Nu=Nu, alpha=a, beta=b)
    try:
        res = ur_gauge.gauge(K, NU, log)
    finally:
        ur_gauge.Conifold = orig
    assert res['C'].alpha == a and res['C'].beta == b
    return res


def fact(n):
    out = 1
    for k in range(2, n + 1): out *= k
    return out


def ic_terms(m, nmax):
    """I^C_m as a list of (n, j, coefficient): coefficient * v^j * z^(-n), n = i + j, j = m + 4i + 1 >= 0, n <= nmax"""
    out = []
    i = 0
    while True:
        j = m + 4 * i + 1
        if j >= 0:
            n = i + j
            if n > nmax: break
            c = Fr((-1) ** i, 4 ** i * fact(i)) * Fr((-1) ** j, fact(j))
            out.append((n, j, c))
        i += 1
        if i > nmax + 2: break
    return out


def model_vanishing(v, NU, nmax):
    """Y^C(v(u), z) as a z-Laurent vector {-n: [Ser]*5}, truncated at O(u^NU), all z^-n with n <= nmax.
    (the term v^j z^-n has j >= 4n/5, so u^N needs n <= 5N/4 only)"""
    vp = [Ser.const(Q, 1, NU)]
    jmax = max(j for b in range(5) for (_, j, _) in ic_terms(MEXP[b], nmax)) if nmax >= 0 else 0
    vt = v.truncate(NU)
    for j in range(1, jmax + 1): vp.append(vp[-1].mul(vt).truncate(NU))
    out = {}
    for b in range(5):
        for (n, j, c) in ic_terms(MEXP[b], nmax):
            if j >= NU: continue
            vec = out.setdefault(-n, [Ser.zero(Q, NU) for _ in range(5)])
            vec[b] = vec[b].add(vp[j].scale(c * SGN[b]))
    return out


def integrate0(f):
    """Int_0^u f du for a Ser f with valuation >= 0 (no u^-1 term allowed)"""
    if not f.c: return Ser.zero(Q, f.prec + 1)
    assert f.v >= 0, "integrand has a pole: valuation %d" % f.v
    return Ser(Q, f.v + 1, [x / (f.v + i + 1) for i, x in enumerate(f.c)])


def quintic_vanishing(bq, NU, nmax):
    """Y^C_Q as a z-Laurent vector {-n: [Ser]*5}: y^(0) = -e_1, y^(n)_a = Int_0 b_a y^(n-1)_{a+1}, n <= nmax.  Asserts that no
    logarithm arises (the integrand b_1 y_2 is a power series)."""
    y = [Ser.zero(Q, NU) for _ in range(5)]; y[1] = Ser.const(Q, -1, NU)
    out = {0: y}
    for n in range(1, nmax + 1):
        new = []
        for a in range(5):
            src = y[(a + 1) % 5]
            if not src.c: new.append(Ser.zero(Q, NU)); continue
            new.append(integrate0(bq[a].mul(src)).truncate(NU))
        y = new; out[-n] = y
        if all(not x.c for x in y): break
    return out


def apply_gauge(p, Y, mmin, mmax, NU):
    """(P Y)[z^m] for mmin <= m <= mmax as a z-Laurent vector; uses every P_k, k <= K = len(p) - 1.  The caller is responsible for
    completeness (see complete())."""
    K = len(p) - 1; out = {}
    for m in range(mmin, mmax + 1):
        vec = [Ser.zero(Q, NU) for _ in range(5)]; nz = False
        for k in range(K + 1):
            src = Y.get(m - k)
            if src is None: continue
            for b in range(5):
                if not src[b].c: continue
                a = (b + k) % 5
                if not p[k][b].c: continue
                vec[a] = vec[a].add(p[k][b].mul(src[b]).truncate(NU)); nz = True
        if nz: out[m] = vec
    return out


def mul_scalar_z(c, W, mmin, mmax, NU):
    """(c(z) W)[z^m], c = {k: Fraction} a polynomial in z (k >= 0)"""
    out = {}
    for m in range(mmin, mmax + 1):
        vec = [Ser.zero(Q, NU) for _ in range(5)]; nz = False
        for k, ck in c.items():
            src = W.get(m - k)
            if src is None or ck == 0: continue
            for a in range(5):
                if src[a].c: vec[a] = vec[a].add(src[a].scale(ck)); nz = True
        if nz: out[m] = vec
    return out


def p11_at_zero(p):
    """P_11(0, z) = Sum_{k = 0 mod 5} p_{k,1}(0) z^k as {k: Fraction}"""
    return {k: p[k][1].coef(0) for k in range(0, len(p), 5)}


def inv_poly(c, K):
    """1/c(z) mod z^(K+1) for c = {k: Fraction}, c[0] != 0"""
    out = {0: 1 / c[0]}
    for k in range(1, K + 1):
        s = Fr(0)
        for j in range(1, k + 1):
            if j in c and (k - j) in out: s += c[j] * out[k - j]
        if s != 0: out[k] = -s / c[0]
    return out


def complete(m, n, K):
    """is the coefficient [z^m u^n] of P Y^C (or of c(z) P Y^C, c a series in z) determined by P_k, k <= K?
    The terms are P_k [u^n1] * Y^C[z^(m-k)] [u^n2], n1 + n2 = n, and Y^C[z^-n'] = O(u^ceil(4 n'/5)): n' <= floor(5n/4)."""
    return m + (5 * n) // 4 <= K


def coef_or_none(s, n):
    return s.coef(n) if n < s.prec else None


def normalise(p, cinv, NU=None):
    """the gauge P c(z): pc[k][b] = Sum_j cinv[j] p[k-j][b]   (cinv = {j: Fraction}, a series in z^5 for the canonical choice)"""
    out = []
    for k in range(len(p)):
        row = []
        for b in range(5):
            acc = None
            for j, cj in cinv.items():
                if j > k or cj == 0: continue
                t = p[k - j][b].scale(cj)
                acc = t if acc is None else acc.add(t)
            if acc is None: acc = Ser.zero(Q, p[k][b].prec)
            row.append(acc if NU is None else acc.truncate(NU))
        out.append(row)
    return out
