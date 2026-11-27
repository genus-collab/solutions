# -*- coding: utf-8 -*-
"""fq_ourec.py -- exact reconstruction helpers for genus2_ou.py: Q (fq_fp.ratrec) and Q(zeta) (Galois conjugates, fq_scalar_lib.qzeta_coords)
values of the O(u) data, and the cross-prime record.

Galois action sigma_m (zeta -> zeta^m) = colour relabelling i -> m i (CONIFOLD_SYMBOL s.1).  On a linear functional J on sp (a J row,
coordinates v = 0..39 in the fq_mfit convention): sigma_m(J_v) = sgn J_{pi_m(v)} with (pi_m, sgn) = fq_mfit.var_perm(m) (the sign for
nu_2 when the pair flips) -- the degree-1 case of fq_leadfit.conj_coeffs3.  On a functional restricted to the CIRCULANT class with the
basis E_{k,d} of fq_leadsolve.class_circulant (c_k(d) = 1, partner c_k(-d) = s_k zeta^{kd}, s_2 = -1, s_1 = s_3 = 1; d in {0,1,2}):
   sigma_m(l(E_{k,0})) = l(E_{k,0});   sigma_m(l(E_{k,d})) = l(E_{k,d'}) if md = d' mod 5 with d' in {1,2},
                                                             = s_k zeta^{-k d'} l(E_{k,d'}) if md = -d' mod 5 with d' in {1,2}
(E_{k,d}^sigma has c(md) = 1, c(-md) = s_k zeta^{kmd}, which is E_{k,d'} resp. s_k zeta^{-kd'} E_{k,d'}).  Functionals restricted to the
covariant-scalar class or to the diagonal class in Fourier modes are Galois-INVARIANT (rational).
Rec: per prime, named residues; at the reporting prime also their Q / Q(zeta) reconstructions; cross(recA, recB) maps every reconstruction
of A to F_p of B and compares with B's residue (the only reliable use of 2^61-1, whose ratrec bound is ~1e7).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import var_perm
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt

CIRC_NAMES = ["c1_0", "c1_1", "c1_2", "c2_1", "c2_2", "c3_0", "c3_1", "c3_2"]
CIRC_KD = [(1, 0), (1, 1), (1, 2), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]


def ratbig(R, a, B=10 ** 13):
    r = R.ratrec(a, B)
    return r


def fmt_rat(R, a, B=None):
    r = R.ratrec(a) if B is None else R.ratrec(a, B)
    return str(r) if r is not None else "[%d mod p]" % (a % R.p)


def jrow_conjugates(row, p):
    """[sigma_1(row), .., sigma_4(row)] with sigma_m(J_v) = sgn J_{pi_m(v)}"""
    out = []
    for m in range(1, 5):
        vp = var_perm(m); out.append([row[vp[v][0]] * vp[v][1] % p for v in range(40)])
    return out


def jrow_qzeta(R, row):
    """Q(zeta) coordinates (F_p 4-tuples) of the 40 entries of a J row"""
    conj = jrow_conjugates(row, R.p)
    return [qzeta_coords(R, [conj[m - 1][v] for m in range(1, 5)]) for v in range(40)]


def circ_conjugates(R, lam):
    """lam: the 8 restricted coefficients l(E_{k,d}) in CIRC_KD order; returns [sigma_1(lam), .., sigma_4(lam)]"""
    p = R.p; out = []
    for m in range(1, 5):
        row = []
        for (k, d) in CIRC_KD:
            if d == 0: row.append(lam[CIRC_KD.index((k, 0))]); continue
            md = (m * d) % 5
            if md in (1, 2): row.append(lam[CIRC_KD.index((k, md))])
            else:
                dp = (-md) % 5; s = (p - 1) if k == 2 else 1
                row.append(s * R.z(-k * dp) * lam[CIRC_KD.index((k, dp))] % p)
        out.append(row)
    return out


def circ_qzeta(R, lam):
    conj = circ_conjugates(R, lam)
    return [qzeta_coords(R, [conj[m - 1][j] for m in range(1, 5)]) for j in range(8)]


def qz_eval(Rb, cs):
    """Q(zeta) coordinates (Fractions) -> F_p element of the ring Rb"""
    return sum(Rb.of(c) * Rb.z(k) for k, c in enumerate(cs)) % Rb.p


class Rec(object):
    """named residues at one prime, with reconstructions (Fraction, or tuple of 4 Fractions for Q(zeta), or None)"""
    def __init__(self, R):
        self.R = R; self.res = {}; self.rec = {}; self.kind = {}

    def rat(self, name, x, B=None):
        """rational reconstruction at the default bound; at the large prime a 1e13 fallback (spurious rate ~1e-1, so these
        reconstructions count only once confirmed by the cross-prime block); None if both fail"""
        R = self.R; x %= R.p; self.res[name] = x; self.kind[name] = 'Q'
        r = R.ratrec(x) if B is None else R.ratrec(x, B)
        if r is None and B is None and R.p > 2 ** 80: r = R.ratrec(x, 10 ** 13); self.kind[name] = 'Qbig' if r is not None else 'Q'
        self.rec[name] = r
        return r

    def qz(self, name, conjugates):
        """conjugates: [sigma_1(x), .., sigma_4(x)]; stores x = sigma_1(x)"""
        R = self.R; self.res[name] = conjugates[0] % R.p; self.kind[name] = 'Qz'
        cs = qzeta_ratrec(R, qzeta_coords(R, conjugates)); self.rec[name] = cs
        return cs

    def fmt(self, name):
        r = self.rec.get(name)
        if r is None: return "[%d mod p]" % self.res[name]
        if self.kind[name] in ('Q', 'Qbig'): return str(r)
        return qzfmt(r)

    def names(self): return list(self.res)


def qzfmt(cs):
    if all(c == 0 for c in cs[1:]): return str(cs[0])
    parts = []
    for k, c in enumerate(cs):
        if c == 0: continue
        parts.append(str(c) if k == 0 else "(%s) z%s" % (c, "" if k == 1 else "^%d" % k))
    return " + ".join(parts) if parts else "0"


def cross(recA, recB):
    """compare A's reconstructions with B's residues: (n compared, mismatches, names not reconstructed at A, names missing at B)"""
    RB = recB.R; n = 0; bad = []; norec = []; missing = []
    for name in recA.names():
        if name not in recB.res: missing.append(name); continue
        r = recA.rec.get(name)
        if r is None: norec.append(name); continue
        n += 1
        val = RB.of(r) if recA.kind[name] in ('Q', 'Qbig') else qz_eval(RB, r)
        if val != recB.res[name] % RB.p: bad.append(name)
    return n, bad, norec, missing


def fmt_qz_fp(R, cs):
    """F_p Q(zeta) coordinates -> string"""
    return qzeta_fmt(R, cs)
