"""is the rescaling tau(u) forced by regularity?  rows of the O(u), O(u^2), O(u^3) linear conditions applied to the rescaling direction M_0^-1 zhat M_0'"""
import sys, json, os
import os
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, HERE)
from fq_fp import P2, P1
from fq3_par import setup
from fq3_rec import Rec
from fq_symbol import mat_series_mul
from fq3_m import KMAX, sp_params, in_sp, NSP
from fq3_quad import QuadForm, IDX_2EDGE

for prime, tag in ((P1, "61"), (P2, "89")):
    st = setup((prime, "0", "0", 34, 6)); R = st['R']; p = R.p; rec = Rec(R); M0 = st['M0']
    dM = [[[k * x % p for x in r] for r in M0.M[k]] for k in range(KMAX + 1)]
    N = mat_series_mul(M0.Minv, dM, p, KMAX); assert in_sp(N, p); d = sp_params(N, p)          # rescaling direction (right placed)
    g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % tag)))
    us = json.load(open(os.path.join(HERE, "g3_nat_%s_nat.json" % tag))); n1 = us['n1']
    dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    rows = {"genus-2 graph 2": g3['J2']['2'], "genus-2 graph 4": g3['J2']['4'], "genus-3 graph 2": g3['J']['2'], "genus-3 graph 5": g3['J']['5']}
    print("p = 2^%s-1: linear rows applied to the rescaling direction M_0^-1 zhat M_0':" % tag)
    for nm, r in rows.items(): print("    J[%s] . d = %s" % (nm, rec.fmt(dot(r, d))))
    Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
    for i in IDX_2EDGE: print("    grad Q_%d(N_1^nat) . d = %s" % (i, rec.fmt(dot(Q[i].gradient(n1), d))))
    # also: the pure-M_0 values and the tau_1-dependence: g1 + a J.d = 0 ?
    from fractions import Fraction as Fr
    g1 = {"genus-2 graph 2": Fr(1, 750), "genus-2 graph 4": Fr(-13, 8000)}
    for nm in g1: print("    %s: [u^1] with M_0 alone = %s; -g1/(J.d) = %s" % (nm, g1[nm], rec.fmt((-R.of(g1[nm])) * pow(dot(rows[nm], d), p - 2, p) % p)))
    for nm, key in (("genus-3 graph 2", '2'), ("genus-3 graph 5", '5')):
        g = g3['A_M0'][key][1]; print("    %s: [u^1] with M_0 alone = %s; -g1/(J.d) = %s" % (nm, rec.fmt(g), rec.fmt((-g) * pow(dot(rows[nm], d), p - 2, p) % p)))
