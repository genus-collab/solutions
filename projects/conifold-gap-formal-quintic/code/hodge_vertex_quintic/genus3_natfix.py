# -*- coding: utf-8 -*-
"""genus3_natfix.py -- which left dressing A(u, zhat) = exp(u^2 a_2 + u^3 a_3) of the natural-variable model M_nat = M_0(tau(u) zhat) is needed for
all genus-3 edge graphs to be regular?  (companion of genus3_natural.py; linear algebra on stored pass tables, no graph sums)

Input: g3_nat_<p>_<base>.json (pass table A3 with M = A_base M_nat, and its n1, n2, n3).  Since M_0^-1 exp(u^k a_k) M = exp(u^k M_0^-1 a_k M_0 + O(u^{k+1}))
exp(Sum u^j N_j), a further left dressing by u^k a_k shifts N_k by delta_k = M_0^-1 a_k M_0 and nothing of lower order.
  --order 2:  residuals r = ([u^2] of the one-edge graphs 2, 5; [u^3] of the five two-edge graphs); equations  J . delta_2 = -r (2),
              grad Q_Gamma(N_1) . delta_2 = -r (5), solved on several classes of a_2:
                 z1hh  a_2 = zhat x E_33                                   (psi-free H (x) H propagator, 1 parameter)
                 z1    a_2 = zhat (x_33 E_33 + x_24 E_24 + x_01 E_01)       (all covariant psi-free propagators, 3 parameters)
                 low   a_2 = any covariant element of orders 1, 2, 3        (8 parameters)
                 circ  a_2 = any covariant element (15 parameters)
              E_mm' = the symmetric bivector with Fourier components (m, m') and (m', m) equal to 25.
  --order 3:  residuals r = [u^3] of the one-edge graphs 2, 5; equations J . delta_3 = -r, same classes.
Writes g3_nat_dress_<name>_<p>.json (the accumulated dressing: base dressing + the new term) for the chosen class (--cls), to be run with
genus3_natural.py --variant <name>.
Run: py -I genus3_natfix.py --prime 89|61 --base nat --order 2 --cls z1 --name d2
"""
import os, sys, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq_symbol import mat_series_mul
from fq3_steps import frame_of
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE
from fq3_par import setup
from fq3_m import NSP, sp_params, in_sp, class_columns, sp_mats, Z5, KMAX
from fq3_rec import Rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--base", default="nat"); ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--cls", default="z1"); ap.add_argument("--name", default="")
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R); z = R.zeta
    LOG = open(os.path.join(HERE, "logs", "genus3_natfix_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 110); log("genus3_natfix: base %s, order u^%d, p = 2^%s-1" % (a.base, a.order, a.prime))
    us = json.load(open(os.path.join(HERE, "g3_nat_%s_%s.json" % (a.prime, a.base)))); A3 = {int(k): v for k, v in us['A3'].items()}; n1 = us['n1']
    qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % a.prime))); g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % a.prime)))
    Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
    J3 = [g3['J']['2'], g3['J']['5']]
    st = setup(frame_of(p, (0, 0), 34, 6)); M0 = st['M0']
    inv25 = pow(25, p - 2, p)
    def E(m, mp):
        """symmetric (m != mp: both components) bivector with Fourier components 25 at (m, mp), (mp, m)"""
        Xh = {(m, mp): 25, (mp, m): 25}
        return [[inv25 * sum(pow(z, (-mm * i) % 5, p) * v * pow(z, (-mq * j) % 5, p) for (mm, mq), v in Xh.items()) % p for j in range(5)] for i in range(5)]
    def left_dir(k, X):
        """right-placed sp coordinates of M_0^-1 (zhat^k X) M_0"""
        S = [Z5() for _ in range(KMAX + 1)]; S[k] = X
        N = mat_series_mul(mat_series_mul(M0.Minv, S, p, KMAX), M0.M, p, KMAX); assert in_sp(N, p)
        return sp_params(N, p)
    circ, cnames = class_columns(R, 'circ')
    def left_circ(j):
        N = mat_series_mul(mat_series_mul(M0.Minv, sp_mats(circ[j], p), p, KMAX), M0.M, p, KMAX); return sp_params(N, p)
    classes = {'z1hh': ([left_dir(1, E(3, 3))], ["x33"], [(1, E(3, 3))]),
               'z1': ([left_dir(1, E(3, 3)), left_dir(1, E(2, 4)), left_dir(1, E(0, 1))], ["x33", "x24", "x01"], [(1, E(3, 3)), (1, E(2, 4)), (1, E(0, 1))]),
               'low': ([left_circ(j) for j in range(8)], cnames[:8], None), 'circ': ([left_circ(j) for j in range(15)], cnames, None)}
    dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    if a.order == 2:
        ok_low = all(A3[i][k] == 0 for i in (2, 5) for k in range(2)) and all(A3[i][k] == 0 for i in IDX_2EDGE for k in range(3))
        log("   base pass: one-edge graphs O(u^2) and two-edge graphs O(u^3): %s" % ok_low)
        rows = [J3[0], J3[1]] + [Q[i].gradient(n1) for i in IDX_2EDGE]; res = [A3[2][2], A3[5][2]] + [A3[i][3] for i in IDX_2EDGE]
        names = ["[u^2] graph 2", "[u^2] graph 5"] + ["[u^3] graph %d" % i for i in IDX_2EDGE]
    else:
        ok_low = all(A3[i][k] == 0 for i in (2, 5) for k in range(3))
        log("   base pass: one-edge graphs O(u^3): %s" % ok_low)
        rows = [J3[0], J3[1]]; res = [A3[2][3], A3[5][3]]; names = ["[u^3] graph 2", "[u^3] graph 5"]
    log("   residuals with the base M: %s" % ", ".join("%s = %s" % (nm, rec.fmt(r)) for nm, r in zip(names, res)))
    out = {}
    for cls in ('z1hh', 'z1', 'low', 'circ'):
        cols, cn, gens = classes[cls]
        Jr = [[dot(r, c) for c in cols] for r in rows]; rhs = [(-r) % p for r in res]
        x, ns = fp_solve(Jr, rhs, p)[:2]; rk = rank_fp(Jr, p)
        log("   class %-5s (%2d parameters): rank %d, %s" % (cls, len(cols), rk, ("CONSISTENT, %d free; particular solution %s" % (len(cols) - rk, ["%s: %s" % (cn[j], rec.fmt(x[j])) for j in range(len(cols)) if x[j]])) if x is not None else "INCONSISTENT"))
        out[cls] = (x, cols, gens)
    if a.name:
        x, cols, gens = out[a.cls]; assert x is not None and gens is not None, "class has no solution or no explicit generators"
        term = [Z5() for _ in range(KMAX + 1)]
        for c, (k, X) in zip(x, gens): term[k] = [[(term[k][i][j] + c * X[i][j]) % p for j in range(5)] for i in range(5)]
        prev = []
        if a.base != 'nat': prev = json.load(open(os.path.join(HERE, "g3_nat_dress_%s_%s.json" % (a.base, a.prime))))['a']
        acc = [prev[n] if n < len(prev) else [Z5() for _ in range(KMAX + 1)] for n in range(a.order)]
        acc[a.order - 1] = [[[(acc[a.order - 1][k][i][j] + term[k][i][j]) % p for j in range(5)] for i in range(5)] for k in range(KMAX + 1)]
        desc = "base %s + u^%d term of class %s: %s" % (a.base, a.order, a.cls, {nm: rec.fmt(c) for nm, c in zip(classes[a.cls][1], x)})
        json.dump(dict(a=acc, desc=desc, coeffs={nm: c for nm, c in zip(classes[a.cls][1], x)}), open(os.path.join(HERE, "g3_nat_dress_%s_%s.json" % (a.name, a.prime)), "w"))
        log("   written g3_nat_dress_%s_%s.json: %s" % (a.name, a.prime, desc))
    LOG.close()


if __name__ == "__main__":
    main()
