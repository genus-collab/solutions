# -*- coding: utf-8 -*-
"""FROZEN theory (rows R^_0, u = 0): the edge-graph sums G^0_Gamma(eps) for the constant M_eps = exp(eps zhat D_1) M_0 (R' = Psi exp(-eps zhat D_1)),
as polynomials in eps (degree <= 6: every zhat carries at most one eps), for D_1 in Pi.  Which eps-coefficients vanish?
usage: py -I -B t9_frozen_eps.py <tag> <beta|gamma|mix|alpha> <graphs> <outfile>"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, time
from evlib import *
from rlib import load, from_fourier, zeros, mul, scale
from fq_symbol import mat_series_exp
from fq3_mvertex import MConst
from fq3_msum import frozen_sums6
from fq3_circsys import interp_simple
tag, which, gl, outf = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], sys.argv[4]
S, d, cs = load(tag); p = S.p; R = S.R
F = {'beta': {(1, 2): 1, (1, 4): 1}, 'gamma': {(1, 0): 1, (1, 1): 1}, 'alpha': {(1, 3): 1},
     'mix': {(1, 2): 1234567891234567, (1, 4): 1234567891234567, (1, 0): 987654321987654321, (1, 1): 987654321987654321}}[which]
D1 = from_fourier(S, F)[1]; st = setup(frame_of(p, (0, 0), 34, 4)); fr = st['frozen']
nodes = list(range(0, 9)); vals = {i: [] for i in gl}; t = time.time()
for eps in nodes:
    D = zeros(); D[1] = [[x * eps % p for x in r] for r in D1]
    E = mat_series_exp(D, R, 6); Ei = mat_series_exp(scale(D, p - 1, p), R, 6)
    Me = MConst(R, mul(E, d['M'], p), mul(d['Minv'], Ei, p), check=True)
    res = frozen_sums6(fr, Me, [st['graphs'][i] for i in gl])
    for i in gl: vals[i].append(res[i])
out = {}
for i in gl:
    c = interp_simple(nodes, vals[i], p); out[i] = c
    print("   %s graph %2d: non-zero eps-orders %s (9 nodes: eps^7, eps^8 must vanish: %s)" % (which, i, [k for k in range(9) if c[k]], c[7] == 0 and c[8] == 0), flush=True)
json.dump(dict(which=which, coeffs={str(i): v for i, v in out.items()}), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), outf), "w"))
print("%s p = 2^%s-1: %.0f s" % (which, tag, time.time() - t))
