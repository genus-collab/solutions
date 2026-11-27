# -*- coding: utf-8 -*-
"""fq3_rec.py -- Q / Q(zeta) reconstruction helpers for the genus-3 programme and the cross-prime comparison of two state files.

Galois action sigma_m (zeta -> zeta^m) = colour relabelling i -> m i (CONIFOLD_SYMBOL s.1) on the 75 sp coordinates of fq3_m:
sp_var_perm(m) -> list over v of (pi_m(v), sign): coordinate (k, (i, j)) -> (k, sorted(m i, m j)) with sign -1 on the antisymmetric
(even k) blocks when the order of the pair flips (the degree-1 case of fq_mfit.var_perm, extended to the six blocks).  For a linear
functional J on sp (a Jacobian row): sigma_m(J_v) = sign . J_{pi_m(v)}; fq_scalar_lib.qzeta_coords then gives the Q(zeta) coordinates
(c_0..c_3 in the basis 1, zeta, zeta^2, zeta^3) of every entry, reconstructed with fq_fp.ratrec.
Rec(R): fmt(a) (Q with the default bound, 1e13 fallback at 2^89-1, then r*zeta^k, then the residue), jrow_qzeta(row) -> list of
(c_0..c_3) Fractions or None per entry.
cross_states(dA, dB, pA, pB): every residue in state A (a JSON tree of ints / lists / dicts) that reconstructs to a rational r at pA is
reduced mod pB and compared with the residue at the same position in state B; returns (compared, mismatches, not_reconstructed, missing, big) where big lists the
1e13-fallback reconstructions with their confirmation status at pB.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import PAIRS1, PAIRS2
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt
from fq3_m import BLOCKS, NSP


def sp_var_perm(m):
    out = []
    for k, pairs, off in BLOCKS:
        for (i, j) in pairs:
            a, b = (m * i) % 5, (m * j) % 5
            if k % 2 == 1: out.append((off + pairs.index(tuple(sorted((a, b)))), 1))
            elif a < b: out.append((off + pairs.index((a, b)), 1))
            else: out.append((off + pairs.index((b, a)), -1))
    return out


class Rec(object):
    def __init__(self, R): self.R = R; self.perms = {m: sp_var_perm(m) for m in range(1, 5)}
    def rat(self, a):
        R = self.R; a %= R.p
        r = R.ratrec(a)
        if r is None and R.p > 10 ** 20: r = R.ratrec(a, 10 ** 13)
        return r
    def fmt(self, a):
        """Q at the default bound, then r*zeta^k at the default bound, then the 1e13 fallback (marked ~, confirm cross-prime), else residue"""
        R = self.R; a %= R.p; r = R.ratrec(a)
        if r is not None: return str(r)
        for k in range(1, 5):
            r = R.ratrec(a * R.z(-k) % R.p)
            if r is not None: return "%s*zeta^%d" % (r, k)
        r = self.rat(a)
        if r is not None: return "~%s" % r
        return "[%d mod p]" % a
    def jrow_conj(self, row):
        p = self.R.p
        return [[row[self.perms[m][v][0]] * self.perms[m][v][1] % p for v in range(NSP)] for m in range(1, 5)]
    def jrow_qzeta(self, row):
        conj = self.jrow_conj(row); out = []
        for v in range(NSP):
            cs = qzeta_coords(self.R, [conj[m - 1][v] for m in range(1, 5)])
            out.append(qzeta_ratrec(self.R, cs))
        return out
    def jrow_fmt(self, row, names):
        q = self.jrow_qzeta(row)
        return "; ".join("%s: %s" % (names[v], qzeta_fmt(self.R, qzeta_coords(self.R, [c[v] for c in self.jrow_conj(row)]))) for v in range(NSP) if row[v])


def _walk(x, path=()):
    if isinstance(x, dict):
        for k, v in x.items(): yield from _walk(v, path + (str(k),))
    elif isinstance(x, list):
        for i, v in enumerate(x): yield from _walk(v, path + (i,))
    elif isinstance(x, int) and not isinstance(x, bool): yield path, x


def _lookup(d, path):
    for k in path:
        if isinstance(d, dict):
            if k not in d: return None
            d = d[k]
        elif isinstance(d, list):
            if not isinstance(k, int) or k >= len(d): return None
            d = d[k]
        else: return None
    return d


def cross_states(dA, dB, RA, pB, skip_prefixes=('time_',)):
    """compare the reconstructible residues of state A with the residues of state B"""
    rec = Rec(RA); compared = 0; mism = []; norec = []; missing = []; big = []
    for path, a in _walk(dA):
        if any(str(path[0]).startswith(s) for s in skip_prefixes): continue
        if abs(a) < 10 ** 6: continue                       # small ints (ranks, dims, graph indices): structural, compared separately
        r = RA.ratrec(a)                                    # default bound (1e10): spurious rate ~ B^2/p ~ 1e-7
        rbig = rec.rat(a) if r is None else None            # 1e13 fallback: reported separately (confirmed or not by the other prime)
        b = _lookup(dB, path)
        if b is None: missing.append(path); continue
        if r is None:
            if rbig is not None:
                ok = (rbig.numerator * pow(rbig.denominator, pB - 2, pB) - b) % pB == 0
                big.append((path, str(rbig), ok))
            norec.append(path); continue
        compared += 1
        if (r.numerator * pow(r.denominator, pB - 2, pB) - b) % pB != 0: mism.append((path, str(r), b))
    return compared, mism, norec, missing, big


if __name__ == "__main__":
    from fq_fp import FpRing, P2
    import random
    R = FpRing(P2); p = R.p; rng = random.Random(5)
    # sigma_m acts as a permutation with signs; sigma_1 = id; sigma_m sigma_n = sigma_{mn}
    for m in range(1, 5):
        vp = sp_var_perm(m); assert sorted(x[0] for x in vp) == list(range(NSP))
    assert all(x == (v, 1) for v, x in enumerate(sp_var_perm(1)))
    for m in range(1, 5):
        for n in range(1, 5):
            vm = sp_var_perm(m); vn = sp_var_perm(n); vmn = sp_var_perm((m * n) % 5)
            assert all((vm[vn[v][0]][0], vm[vn[v][0]][1] * vn[v][1]) == vmn[v] for v in range(NSP))
    # a covariant functional (rational on the cov class) reconstructs with rational Q(zeta) coordinates
    rec = Rec(R); row = [0] * NSP
    from fq3_m import class_cov
    A, names = class_cov(R)
    for v in range(NSP):
        for t, c in A[v].items(): row[v] = (row[v] + c * (t + 1)) % p          # functional with J(b-direction) = Sum (t+1) zeta^{-ki}...
    q = rec.jrow_qzeta(row); assert all(x is not None for x in q)
    print("sp_var_perm group law ok; Q(zeta) reconstruction of a covariant row ok:", [qzeta_fmt(R, qzeta_coords(R, [c[v] for c in rec.jrow_conj(row)])) for v in (0, 5, 9, 12, 14)])
    print("fq3_rec ok")
