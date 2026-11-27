# -*- coding: utf-8 -*-
"""ms_regress3.py -- FULL cross-check (symbol to zhat^KZ, default 9, so every tensor up to (3,3) is reachable) of the flat-basis engine against the project's MVertex6 (independent code), mod 2^89-1:
EVERY one-point tensor at genus 1..4 and EVERY two-point tensor at genus 1..3 (and every three-point tensor at genus 2) whose psi-deficit
3h - 3 + n - Sum d is <= 6 (MVertex6 needs the symbol to zhat^deficit), all flat index tuples.  Same construction as ms_regress part B.
Run: py -I -B ms_regress3.py [KZ] [g,n ...]   (built from ms_regress2.py by sed; same logic)
"""
import itertools, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr
import ms_engine as me


def main():
    sys.path.insert(1, mr.PROJ)
    from fractions import Fraction as Fr
    from fq_fp import FpRing, P2
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq_symbol import Symbol
    from fq3_mvertex import MVertex6
    R = FpRing(P2); p = R.p; KZ = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    rows = Rows(R, Nz=KZ); Cf = Conifold(R, Nu=10, alpha=Fr(0), beta=Fr(0)); S = SRows(rows, Cf, prec_abs=34)
    sym = Symbol(S); assert sym.val_ok
    Cm = sym.colour_matrix()
    Minv = [[[Cm[k][j][i] * (1 if k % 2 == 0 else p - 1) % p for j in range(5)] for i in range(5)] for k in range(KZ + 1)]
    class _M(object): pass
    M0 = _M(); M0.Minv = Minv
    mv = MVertex6(R, M0, kmax=KZ)
    Psi = [[mr.N_SIGN[a] * R.z(i * (a + 1)) % p for i in range(5)] for a in range(5)]
    jobs = []
    for g, n in [tuple(map(int, x.split(","))) for x in sys.argv[2:]] or ((3, 3), (3, 2), (4, 1), (2, 3), (1, 4)):
        if not os.path.exists(me._cache_path(g, n)): print("   skip (%d,%d): no amplitude" % (g, n)); continue
        D = 3 * g - 3 + n
        for d in itertools.combinations_with_replacement(range(D + 1), n):
            if sum(d) > D or D - sum(d) > KZ: continue
            jobs.append((g, d))
    nt = bad = 0; t0 = time.time()
    for g, d in jobs:
        n = len(d); T = mv.tensor(g, d)
        for idx in itertools.product(range(5), repeat=n):
            val = 0
            for col in itertools.product(range(5), repeat=n):
                x = T[sum(c * 5 ** l for l, c in enumerate(col))]
                if not x: continue
                for l in range(n): x = x * Psi[idx[l]][col[l]] % p
                val = (val + x) % p
            mine = me.corr(g, tuple(zip(idx, d)))
            nt += 1
            if R.of(mine) != val: bad += 1; print("   mismatch", g, d, idx, mine)
    print(("PASS" if bad == 0 and nt else "FAIL") + " %d flat correlators on %d tensors (genus 1-4, n <= 3, psi-deficit <= 6) = MVertex6 mod 2^89-1 (%.0f s)" % (nt, len(jobs), time.time() - t0))


if __name__ == "__main__":
    main()
