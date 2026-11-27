# -*- coding: utf-8 -*-
"""m_1 in the left-placed Fourier coordinates: m_1 = x hat c^D_1(3) + y hat c^D_2(4)?  (edge-dressed genus-1 one-point function);
   and the other stored circulant solutions (circ_b_solutions, circ_b_extra) in the left-placed Fourier coordinates"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
from rlib import *
from fq3_rec import Rec
from fq_mfit import fp_solve
from fq3_quad import rank_fp
for tag in ("61", "89"):
    S, d, cs = load(tag); p = S.p; rec = Rec(S.R); M = d['M']; Mi = d['Minv']
    rowsF = []
    for j in range(15):
        e = [0] * 15; e[j] = 1; F, _ = fourier(S, Ad(M, Mi, ser_from_circ(S, e), p)); rowsF.append((F[(1, 3)], F[(2, 4)]))
    sol = fp_solve([[rowsF[j][0], rowsF[j][1]] for j in range(15)], S.m1, p)
    x, y = sol[0] if sol[0] is not None else (None, None)
    print("p = 2^%s-1: m_1 = x hat c^D_1(3) + y hat c^D_2(4) solvable on all 15 coordinates: %s; y/x = %s" % (tag, sol[0] is not None, rec.fmt(y * pow(x, p - 2, p) % p) if x else None))
    for key in ('circ_b_solutions', 'circ_b_extra'):
        for v in cs.get(key, []):
            if isinstance(v, dict): v = v.get('v') or v.get('point')
            if not isinstance(v, list) or len(v) != 15: print("   %s: unexpected entry %s" % (key, str(v)[:80])); continue
            res = [k for k, r in S.residuals(v).items() if r]
            F, _ = fourier(S, Ad(M, Mi, ser_from_circ(S, v), p))
            pat = {k: "".join("0" if F[(k, m)] == 0 else "*" for m in range(5)) for k in range(1, 7)}
            print("   %s: c1 = %s; failing conditions %s; left-placed Fourier pattern (m = 0..4 per order): %s; hat c_1(3) = %s, hat c^D_2(4) = %s" % (key, [rec.fmt(t) for t in v[:3]], res, pat, rec.fmt(F[(1, 3)]), rec.fmt(F[(2, 4)])))
