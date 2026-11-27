# -*- coding: utf-8 -*-
"""PREDICTION from the explanation of Pi: in the FULL class (no covariance) every Delta = M0^-1 (zhat D_1) M0 with D_1 symmetric and
D_1 phi_3 = 0 (10 parameters: the propagator Sum_j D_1 e_j (x) e_j lies in Sym^2 W, W = span(1, H^2, H^3, H^4)) satisfies ALL N_1-conditions.
Direct evaluations (frozen vertex tensors, exact for the leading coefficients), plus the edge tables, plus controls."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, random, time
from evlib import *
from rlib import Ad, zeros, fourier, from_fourier, ser_from_circ, load
from fq3_m import sp_params, in_sp
from fq3_quad import QuadForm, IDX_2EDGE
from fq3_cubic import IDX_3EDGE
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
S, d, cs = load(tag); p = S.p; R = S.R; z = R.z; rng = random.Random(99); M = d['M']; Mi = d['Minv']; nt = S.nt
qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag))); g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag)))
Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
rows = dict(L2_0=g3['J2']['2'], L2_1=g3['J2']['4'], J3_0=g3['J']['2'], J3_1=g3['J']['5'])
dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
inv5 = pow(5, p - 2, p); phi = lambda m: [z(m * j) for j in range(5)]
def sym_random():
    X = [[0] * 5 for _ in range(5)]
    for i in range(5):
        for j in range(i, 5): X[i][j] = X[j][i] = rng.randrange(p)
    return X
def mm(A_, B_): return [[sum(A_[i][k] * B_[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
def T(A_): return [[A_[j][i] for j in range(5)] for i in range(5)]
Pi_ = [[((1 if i == j else 0) - phi(3)[i] * phi(2)[j] * inv5) % p for j in range(5)] for i in range(5)]      # Pi phi_3 = 0
def D_in_W(): return mm(mm(T(Pi_), sym_random()), Pi_)
def D_Hrow():
    """propagator H (x) x + x (x) H, x in W random (no (H, H) term): D = phi_2 x^T + x phi_2^T, x.phi_3 = 0"""
    x = [rng.randrange(p) for _ in range(5)]; c = sum(x[j] * phi(3)[j] for j in range(5)) * inv5 % p; x = [(x[j] - c * phi(2)[j]) % p for j in range(5)]
    return [[(phi(2)[i] * x[j] + x[i] * phi(2)[j]) % p for j in range(5)] for i in range(5)]
def n_of(D1):
    D = zeros(); D[1] = D1; assert in_sp(D, p)
    return sp_params(Ad(Mi, M, D, p), p)
ALL = [2, 5] + IDX_2EDGE + IDX_3EDGE; EOF_ = {2: 1, 5: 1}; EOF_.update({i: 2 for i in IDX_2EDGE}); EOF_.update({i: 3 for i in IDX_3EDGE})
def report(name, D1, direct=True):
    dl = n_of(D1); n1 = [(nt[k] + dl[k]) % p for k in range(NSP)]
    kills = all(sum(D1[i][j] * phi(3)[j] for j in range(5)) % p == 0 for i in range(5))
    lin = [k for k, r in rows.items() if dot(r, dl)]; qu = [i for i in IDX_2EDGE if Q[i].eval(n1)]
    msg = "%s: D_1 phi_3 = 0: %s; non-zero linear rows %s; non-zero stored quadrics %s" % (name, kills, lin, qu)
    if direct:
        t = time.time(); v, _ = evaluate(tag, [n1], ALL)
        nz = [(i, k) for i in ALL for k in range(EOF_[i] + 1) if v[i][k]]
        msg += "; direct evaluation, non-zero [u^k], k <= E, (graph, k): %s (%.0f s)" % (nz, time.time() - t)
    print(msg, flush=True)
    return n1
for t in range(2): n1 = report("full class, random D_1 in Sym^2 W (#%d)" % t, D_in_W())
# edge and leaf tables at the last point: which (m, n) carry a u^1 term?
from fq3_steps import frame_of
st = setup(frame_of(p, (0, 0), 34, 4)); GS = GraphSummer6(st['H'], MU(R, st['M0'], [sp_mats(n1, p)], 4, 'right'), st['mv0'])
nzE = sorted((m, n) for m in range(6) for n in range(6 - m) if any(GS.EL.E(m, n)[i][j].coef(1) for i in range(5) for j in range(5)))
nz0 = sorted((m, n) for m in range(6) for n in range(6 - m) if any(GS.EL.E(m, n)[i][j].coef(0) for i in range(5) for j in range(5)))
nzT = [dd for dd in range(2, 8) if any(GS.EL.T(dd)[i].coef(1) for i in range(5))]
print("edge tables E'(m, n) at this point: non-zero at u^0: %s; non-zero at u^1: %s (expected only (0, 0): a primary-primary propagator); leaves T'_d non-zero at u^1: d in %s" % (nz0, nzE, nzT))
report("control: random symmetric D_1 (no constraint)", sym_random())
report("control: propagator H (x) x + x (x) H, x in W (no (H, H) term, but H appears)", D_Hrow())
report("circulant Pi point A + 2B (stored)", None) if False else None
v15 = [(a + 2 * b) % p for a, b in zip(cs['Pi_A'], cs['Pi_B'])]; n1 = S.to_sp(v15); t = time.time(); v, _ = evaluate(tag, [n1], ALL)
print("stored circulant Pi point A + 2B: direct evaluation, non-zero [u^k], k <= E: %s (%.0f s)" % ([(i, k) for i in ALL for k in range(EOF_[i] + 1) if v[i][k]], time.time() - t))
print("   its [u^(E+1)] coefficients (frozen vertices: NOT the true values, shown only as non-vanishing of higher terms): %s" % {i: bool(v[i][EOF_[i] + 1]) for i in ALL if EOF_[i] < 3})
