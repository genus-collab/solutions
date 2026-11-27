# -*- coding: utf-8 -*-
"""ref13_stored.py -- is Delta* related to the stored O(u^2) circulant solutions of g3_quad_<p>.json (rational_solution_n, circ_solution_n)?"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *
from fq_mfit import fp_solve
from fq3_rec import Rec
for tag in ("61", "89"):
    S = load(tag); p = S.p; I = IRing(p); rec = Rec(S.R)
    qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag))); d = json.load(open(os.path.join(HERE, "ref9_newsol_%s.json" % tag))); dl = d['delta']
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    for key in ('rational_solution_n', 'circ_solution_n'):
        n = qs[key]; diff = [(n[k] - S.nt[k]) % p for k in range(len(n))]
        sol = fp_solve([[S.cols[j][k] for j in range(15)] for k in range(len(n))], diff, p)
        if sol[0] is None: print(tag, key, "not in the circulant class"); continue
        v = sol[0]; res = S.residuals(v); bad = [k for k, x in res.items() if x]
        print("p=2^%s-1 %s: circulant coords low part %s" % (tag, key, [rec.fmt(x) for x in v[:8]]))
        print("     non-zero residuals: %s; rank with Delta* on coords 0..9: %d; rank with (A, B) on coords 0..9: %d" % (bad, rank_field(I, [v[:10], dl[:10]]), rank_field(I, [v[:10], A[:10], B[:10]])))
