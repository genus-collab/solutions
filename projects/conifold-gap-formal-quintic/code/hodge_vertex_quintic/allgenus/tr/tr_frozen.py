# -*- coding: utf-8 -*-
"""tr_frozen.py -- the project's frozen theory M_0.omega^ for comparison with TR (read-only use of project files and modules).

M_0 = R~ (Teleman colour-colour form, zhat-series), C_n[j][i] = (1/5) Sum_a zeta^(j(4-a) + i(a+1-n)) g^a_n  (LG_MODEL.md s.2 (D)),
g^a_n = rows of natural_form/lg/lg_symbol_rows_Q.json (n <= 40, exact).  M_0^-1(zhat) = M_0(-zhat)^T.
Cross-check with the project's own symbol: natural_form/lg/lg_symbol_<89|61>.json (Rh[a][i][n], n <= 6, the dump of fq3_par.setup):
C_n[j][i] = Sum_ab Rh[b][j][0] eta^{ba} Rh[a][i][n].
Frozen correlators: fq3_mvertex.MVertex6(R, MConst(M_0)) (project code, mod p, psi-deficit <= 6), in the colour basis et_i;
flat insertions H^a = N_a Sum_i zeta^(i(a+1)) et_i, N = (1,-1,-1,1,1).  The zeta of fq_fp.FpRing is used throughout, so residues
are directly comparable with tr_field.make_fp(p, zeta=R.zeta).
"""
import os, sys, json, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
if HERE not in sys.path: sys.path.insert(0, HERE)
ROWS = os.path.join(PROJ, "natural_form", "lg", "lg_symbol_rows_Q.json")
NSIGN = (1, -1, -1, 1, 1)


def rows_Q(nmax):
    d = json.load(open(ROWS))
    return [[Fr(x) for x in d['g'][a][:nmax + 1]] for a in range(5)]


def M0_coeffs(F, nmax):
    """[C_0, .., C_nmax] with C_n[j][i] in the field F (tr_field)"""
    g = rows_Q(nmax); out = []
    for n in range(nmax + 1):
        M = [[F.zero] * 5 for _ in range(5)]
        for j in range(5):
            for i in range(5):
                s = F.zero
                for a in range(5): s = s + F.z(j * (4 - a) + i * (a + 1 - n)) * g[a][n]
                M[j][i] = s * Fr(1, 5)
        out.append(M)
    return out


def dump_check(tag):
    """C_n from the rows (mod p, the dump's zeta) == C_n from the project's dumped symbol, n <= 6: returns (ok, number of entries)"""
    d = json.load(open(os.path.join(PROJ, "natural_form", "lg", "lg_symbol_%s.json" % tag)))
    p = d['p']; z = d['zeta']; Rh = d['Rh']; Nz = d['Nz']
    ei = {(a, b): e for a, b, e in d['eta_inv']}
    from tr_field import make_fp
    F = make_fp(p, zeta=z); C = M0_coeffs(F, Nz); bad = 0; cnt = 0
    for n in range(Nz + 1):
        for j in range(5):
            for i in range(5):
                v = 0
                for (a, b), e in ei.items(): v += Rh[b][j][0] * e * Rh[a][i][n]
                cnt += 1
                if (v - C[n][j][i].v) % p: bad += 1
    return bad == 0, cnt


class Frozen(object):
    """frozen theory mod p via the project's MVertex6 (psi-deficit <= kmax = 6)"""
    def __init__(self, p):
        if PROJ not in sys.path: sys.path.insert(0, PROJ)
        from fq_fp import FpRing
        from fq3_mvertex import MVertex6, MConst
        from tr_field import make_fp
        self.R = FpRing(p); self.p = p; R = self.R
        self.F = make_fp(p, zeta=R.zeta); C = M0_coeffs(self.F, 6)
        M = [[[C[k][j][i].v for i in range(5)] for j in range(5)] for k in range(7)]
        Minv = [[[(C[k][j][i].v * (1 if k % 2 == 0 else p - 1)) % p for j in range(5)] for i in range(5)] for k in range(7)]
        self.M0 = MConst(R, M, Minv, check=True)
        self.mv = MVertex6(R, self.M0)

    def colour(self, g, d):
        """{(i_1..i_n): <prod tau_{d_l}(et_{i_l})>_g} (F_p ints)"""
        T = self.mv.tensor(g, tuple(d)); n = len(d)
        return {idx: T[sum(i * 5 ** l for l, i in enumerate(idx))] for idx in itertools.product(range(5), repeat=n)}

    def flat(self, g, d, avec):
        """<prod tau_{d_l}(H^{a_l})>_g (F_p int)"""
        R = self.R; p = self.p; T = self.mv.tensor(g, tuple(d)); n = len(d); tot = 0
        for idx in itertools.product(range(5), repeat=n):
            w = 1
            for l, i in enumerate(idx): w = w * NSIGN[avec[l]] * R.z(i * (avec[l] + 1)) % p
            tot += w * T[sum(i * 5 ** l for l, i in enumerate(idx))]
        return tot % p
