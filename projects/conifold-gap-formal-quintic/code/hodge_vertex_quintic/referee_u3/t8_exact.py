# -*- coding: utf-8 -*-
"""EXACT left dressing of the tautological factorisation: M(zhat; u) = exp(u zhat D_1) R~(zhat; u) mod (zhat^7, u^4), R' = Psi(u) exp(-u zhat D_1).
Prediction: for the string direction D_beta (propagator 1 (x) H^2 + H^2 (x) 1; D_beta^2 = 0, D_beta 1 = 0) EVERY edge graph vanishes identically in u
(string equation of the vertex CohFT), i.e. this direction of Pi is a trivial (gauge) direction.  D_gamma ((H^3, H^4) propagator): tested, no prediction.
usage: py -I -B t8_exact.py <tag> <beta|gamma|mix> <NN> <graphs> <outfile>"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, time
from evlib import *
from rlib import load, from_fourier, zeros
from fq_series import Ser
from fq_sser import mat_mul, mat_T
from fq3_m import Z5, KMAX
tag, which, NN, gl, outf = sys.argv[1], sys.argv[2], int(sys.argv[3]), [int(x) for x in sys.argv[4].split(",")], sys.argv[5]
S, d, cs = load(tag); p = S.p; R = S.R; K = 4
F = {'beta': {(1, 2): 1, (1, 4): 1}, 'gamma': {(1, 0): 1, (1, 1): 1}, 'mix': {(1, 2): 1234567891234567, (1, 4): 1234567891234567, (1, 0): 987654321987654321, (1, 1): 987654321987654321}}[which]
D1 = from_fourier(S, F)[1]
st = setup(frame_of(p, (0, 0), 34, K)); H = st['H']
from fq3_taut import colour_series_u
C = [colour_series_u(H, k) for k in range(K)]                    # C[k][n] = [u^k] of the colour matrices of R~
# bivariate product exp(u zhat D_1) . C(u): keys (n, j) = (zhat-order, u-order)
Mb = {}
pw = [[1 if i == j else 0 for j in range(5)] for i in range(5)]; fact = 1
for j in range(K):
    if j > 0: pw = mat_mul(pw, D1, p); fact = fact * j % p
    c = pow(fact, p - 2, p)
    for k in range(K - j):
        for n in range(KMAX + 1 - j):
            P = mat_mul(pw, C[k][n], p); key = (n + j, j + k); cur = Mb.get(key, Z5())
            Mb[key] = [[(cur[a][b] + c * P[a][b]) % p for b in range(5)] for a in range(5)]
class MExact(MU):
    def __init__(self):
        self.R = R; self.p = p; self.K = K; self.kmax = KMAX; self.placement = 'exact'; self.M0 = st['M0']; self.is_identity = False
        self.M = [[[[Mb.get((n, j), Z5())[a][b] for j in range(K)] for b in range(5)] for a in range(5)] for n in range(KMAX + 1)]
        sg = lambda n: 1 if n % 2 == 0 else p - 1
        self.Minv = [[[[Mb.get((n, j), Z5())[b][a] * sg(n) % p for j in range(K)] for b in range(5)] for a in range(5)] for n in range(KMAX + 1)]
        self.M_ser = [[[Ser(R, 0, list(self.M[n][a][b])) for b in range(5)] for a in range(5)] for n in range(KMAX + 1)]
        self.Minv_ser = [[[Ser(R, 0, list(self.Minv[n][a][b])) for b in range(5)] for a in range(5)] for n in range(KMAX + 1)]
M = MExact(); assert M.check(), "M not symplectic / inverse wrong / M(0) != M_0"
# row 0 of Psi(u) proportional to the unit phi_1 at every u-order?
unit_ok = all((H.Q[0][i][0].coef(k) - H.Q[0][0][0].coef(k) * R.z(i)) % p == 0 for k in range(K) for i in range(5))
t = time.time(); vt = UVertex(R, M, K, NN, degree_check=False); GS = GraphSummer6(H, M, vt)
nzE = sorted((m, n, k) for m in range(6) for n in range(6 - m) for k in range(K) if any(GS.EL.E(m, n)[i][j].coef(k) for i in range(5) for j in range(5)))
nzT = sorted((dd, k) for dd in range(2, 8) for k in range(K) if any(GS.EL.T(dd)[i].coef(k) for i in range(5)))
out = {}
for i in gl:
    s = GS.graph(st['graphs'][i]); out[i] = [s.coef(k) if k < s.prec else None for k in range(K)]
json.dump(dict(which=which, NN=NN, maxdeg=vt.maxdeg, vals={str(i): v for i, v in out.items()}, nzE=nzE, nzT=nzT, unit_ok=unit_ok), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), outf), "w"))
print("%s p=2^%s-1 NN=%d maxdeg=%d (%.0f s): row 0 of Psi(u) prop. to the unit: %s; non-zero edge tables (m, n, u-order): %s; non-zero leaves (d, u-order): %s" % (which, tag, NN, vt.maxdeg, time.time() - t, unit_ok, nzE, nzT))
for i in gl: print("   graph %2d: non-zero u-orders among u^0..u^3: %s" % (i, [k for k in range(K) if out[i][k]]))
