# -*- coding: utf-8 -*-
"""lg_rows_highorder.py -- ALL five rows of the conifold symbol beyond the stored order zhat^6: the project pipeline
(fq_rows.Rows -> fq_conifold.Conifold -> fq_sser.SRows -> fq_symbol.Symbol; read-only imports, alpha = beta = 0) is re-run in memory with a
larger z-order Nz and the symbol is compared, residue by residue, with the closed LG recursion (lg_gm.symbol_rows):
        R^_{a,n}^{(i)} = N_a zeta^{i(a+1-n)} g^a_n,   N = (1,-1,-1,1,1),   g^a = J_{c_a}(-4 zhat),  c = (0,-1,3,2,1).
Also the symplecticity and covariance checks of the Symbol object at that order.  Nothing is written to the project folder.
Run:  py -I lg_rows_highorder.py [Nz] [alpha beta]     (default Nz = 12; optional conifold-frame parameters as fractions, e.g. 1/3 -2/7)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_conifold import Conifold
from fq_sser import SRows
from fq_symbol import Symbol
from lg_gm import symbol_rows, CEXP, NSIGN


def main():
    NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    al = Fr(sys.argv[2]) if len(sys.argv) > 3 else Fr(0); be = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0)
    g = symbol_rows(NZ); allok = True
    for tag, p in (("89", P2), ("61", P1)):
        R = FpRing(p); t0 = time.time()
        rows = Rows(R, Nz=NZ); C = Conifold(R, Nu=10, alpha=al, beta=be); S = SRows(rows, C, prec_abs=34); sym = Symbol(S)
        bad = [(a, i, n) for a in range(5) for i in range(5) for n in range(NZ + 1)
               if sym.Rh[a][i][n] != NSIGN[a] * R.pow(R.zeta, (i * (a + 1 - n)) % 5) * R.of(g[a][n]) % p]
        sp = len(sym.symplectic_check()); cov = sym.covariance()['phase_a_plus_1_minus_n']
        ok = sym.val_ok and not bad and sp == 0 and cov
        print("%s p = 2^%s-1, (alpha, beta) = (%s, %s), Nz = %d: valuations ok %s; %d residues compared, %d mismatches %s; symplectic failures %d; covariance %s  (%.1f s)"
              % ("PASS" if ok else "FAIL", tag, al, be, NZ, sym.val_ok, 25 * (NZ + 1), len(bad), bad[:4], sp, cov, time.time() - t0), flush=True)
        allok = allok and ok
    print("ALL PASS" if allok else "FAILURES")


if __name__ == "__main__":
    main()
