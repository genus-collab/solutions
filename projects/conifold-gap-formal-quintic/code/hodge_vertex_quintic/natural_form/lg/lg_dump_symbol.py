# -*- coding: utf-8 -*-
"""lg_dump_symbol.py -- dump the conifold symbol rows and the hatted rows (u-series) of the formal quintic at one prime
to a JSON file in THIS folder (read-only use of the project modules; nothing in the project folder is written).

Output lg_symbol_<tag>.json:
   p, zeta, eta_inv ([(a, b, e)]), Rh[a][i][n] (symbol residues), HQ[a][i][n] = list of the u-coefficients [u^0..u^{K-1}] of the hatted
   row Qhat_{a,n}^{(i)}, xhat (u-coefficients of x/s), checks.
Run:  py -I lg_dump_symbol.py 89|61
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(os.path.dirname(HERE)); sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_par import setup


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "89"
    p = P2 if tag == "89" else P1; t0 = time.time()
    st = setup((p, "0", "0", 34, 6)); R = st['R']; sym = st['sym']; H = st['H']
    K = st['K']
    def coefs(s): return [(s.coef(k) % p) if k < s.prec else None for k in range(K)]
    out = dict(p=p, zeta=R.zeta, K=K, Nz=sym.Nz,
               eta_inv=[[a, b, e % p] for (a, b), e in sym.eta_inv.items()],
               Rh=[[[sym.Rh[a][i][n] % p for n in range(sym.Nz + 1)] for i in range(5)] for a in range(5)],
               HQ=[[[coefs(H.Q[a][i][n]) for n in range(H.Nz + 1)] for i in range(5)] for a in range(5)],
               xhat=coefs(H.x),
               checks=dict(val_ok=sym.val_ok, hat_val_ok=H.val_ok, symbol_check=H.symbol_check(sym),
                           symplectic_failures=len(sym.symplectic_check()),
                           cov_a_plus_1_minus_n=sym.covariance()['phase_a_plus_1_minus_n']),
               time=time.time() - t0)
    path = os.path.join(HERE, "lg_symbol_%s.json" % tag)
    json.dump(out, open(path, "w"))
    print("wrote", path, "checks:", out['checks'], "time %.1f s" % out['time'])


if __name__ == "__main__":
    main()
