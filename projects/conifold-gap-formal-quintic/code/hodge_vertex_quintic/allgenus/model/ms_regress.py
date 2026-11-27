# -*- coding: utf-8 -*-
"""ms_regress.py -- regression of the flat-basis engine ms_engine against the project's own code and numbers on record.

(A) numbers on record: F_2 = 1/1200, F_3 = 1/25200 (F_4 = 1/180000 if the (4,0) amplitude is cached), <tau_1(1)>_1 = 5/24,
    <H>_1; the primary H-sector  <h^n>_g = d_t^n (c_g t^(2-2g)) at t = -1 with h = h^1 = -H (UNIFORM_REDUCTION s.6.3:
    there "H" is the Fourier class without the sign N_1), i.e. <H^n>_g = d_v^n (c_g v^(2-2g)) at v = 1.
(B) the project's vertex tensors: fq3_mvertex.MVertex6 on the symbol M_0 produced by the project pipeline
    (fq_rows -> fq_conifold -> fq_sser -> fq_symbol, as in fq4_split.Split.__init__), prime 2^89-1, canonical basis e~_i;
    flat components H^a = Sum_i N_a zeta^(i(a+1)) e~_i.  Compared with ms_engine.corr reduced mod p for every tensor
    (g, d) listed (psi-deficit <= 6) and every flat index tuple.
(C) string and dilaton equations, symmetry under permutation of the insertions (sample).
Everything read-only on the project; run:  py -I -B ms_regress.py [--noB]
"""
import itertools, os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr
import ms_engine as me

FAIL = []
def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok: FAIL.append(name)


def cg(g):
    from fq_fp import bernoulli
    return (-1) ** (g - 1) * Fr(bernoulli(2 * g)) / (2 * g * (2 * g - 2)) / 5 ** (g - 1)


def part_A():
    check("(A) F_2 = 1/1200", me.corr(2, ()) == Fr(1, 1200) == cg(2))
    check("(A) F_3 = 1/25200", me.corr(3, ()) == Fr(1, 25200) == cg(3))
    if os.path.exists(me._cache_path(4, 0)):
        check("(A) F_4 = 1/180000", me.corr(4, ()) == Fr(1, 180000) == cg(4))
    check("(A) <tau_1(1)>_1 = 5/24", me.corr(1, ((0, 1),)) == Fr(5, 24))
    check("(A) <H>_1 = -1/12  (F_1 = -(1/12) log v)", me.corr(1, ((1, 0),)) == Fr(-1, 12))
    check("(A) <H,H,H>_0 = -5  (F_0''' = -5/v)", me.corr(0, ((1, 0),) * 3) == -5)
    for g, nmax in ((2, 3), (3, 3), (4, 1)):
        for n in range(1, nmax + 1):
            if not os.path.exists(me._cache_path(g, n)): continue
            exp_ = cg(g)
            for j in range(n): exp_ *= (2 - 2 * g - j)
            check("(A) <H^%d>_%d = d_v^%d (c_g v^(2-2g)) at v = 1 = %s" % (n, g, n, exp_), me.corr(g, ((1, 0),) * n) == exp_)
            check("(A) <h^%d>_%d = d_t^%d (c_g t^(2-2g)) at t = -1" % (n, g, n), me.corr_h(g, ((1, 0),) * n) == (-1) ** n * exp_)
    for n, val in ((1, Fr(-1, 12)), (2, Fr(1, 12)), (3, Fr(-1, 6)), (4, Fr(1, 2))):
        if os.path.exists(me._cache_path(1, n)):
            check("(A) <H^%d>_1 = d_v^%d (-(1/12) log v) = %s" % (n, n, val), me.corr(1, ((1, 0),) * n) == val)


def part_B():
    sys.path.insert(1, mr.PROJ)
    from fq_fp import FpRing, P2
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq_symbol import Symbol
    from fq3_mvertex import MVertex6
    t0 = time.time(); R = FpRing(P2); p = R.p; KZ = 6; K = 2
    rows = Rows(R, Nz=KZ); Cf = Conifold(R, Nu=K + 4, alpha=Fr(0), beta=Fr(0)); S = SRows(rows, Cf, prec_abs=5 * K + 4)
    sym = Symbol(S); assert sym.val_ok
    Cm = sym.colour_matrix()
    Minv = [[[Cm[k][j][i] * (1 if k % 2 == 0 else p - 1) % p for j in range(5)] for i in range(5)] for k in range(KZ + 1)]
    class _M(object): pass
    M0 = _M(); M0.Minv = Minv
    mv = MVertex6(R, M0, kmax=KZ)
    print("   project symbol and MVertex6 built (%.1f s)" % (time.time() - t0), flush=True)
    # the project's colour matrices against the Q rows of the LG model: C_n[j][i] = (1/5) Sum_a zeta^(j(4-a)+i(a+1-n)) g^a_n
    ok = True
    for n in range(KZ + 1):
        for j in range(5):
            for i in range(5):
                x = sum(R.z(j * (4 - a) + i * (a + 1 - n)) * R.of(mr.g(a, n)) for a in range(5)) * R.inv(5) % p
                ok = ok and x == Cm[n][j][i] % p
    check("(B) project colour matrices C_0..C_6 = Fourier transform of the Q rows (2^89-1)", ok)
    Psi = [[mr.N_SIGN[a] * R.z(i * (a + 1)) % p for i in range(5)] for a in range(5)]
    tests = [(0, (0, 0, 0)), (0, (0, 0, 0, 0)), (0, (1, 0, 0, 0)), (1, (0,)), (1, (1,)), (1, (0, 0)), (1, (1, 0)), (1, (2, 0)), (1, (0, 0, 0)), (1, (1, 1, 0)),
             (2, ()), (2, (0,)), (2, (1,)), (2, (3,)), (2, (0, 0)), (2, (1, 0)), (2, (2, 1)), (2, (0, 0, 0)), (2, (1, 0, 0)), (2, (2, 2, 1)),
             (3, ()), (3, (1,)), (3, (2,)), (3, (1, 1)), (3, (2, 1)), (3, (3, 0)), (3, (1, 1, 1)), (3, (2, 1, 0))]
    nt = 0; bad = 0; t0 = time.time()
    for g, d in tests:
        n = len(d)
        if not os.path.exists(me._cache_path(g, n)): continue
        T = mv.tensor(g, d)
        for idx in itertools.product(range(5), repeat=n):
            val = 0
            for col in itertools.product(range(5), repeat=n):
                x = T[sum(c * 5 ** l for l, c in enumerate(col))]
                if not x: continue
                for l in range(n): x = x * Psi[idx[l]][col[l]] % p
                val = (val + x) % p
            mine = R.of(me.corr(g, tuple(zip(idx, d))))
            nt += 1
            if mine != val: bad += 1; print("   mismatch", g, d, idx, R.fmt(val), me.corr(g, tuple(zip(idx, d))))
    check("(B) %d flat correlators (%d tensors, psi-deficit <= 6, genus 0..3) = MVertex6 on the project symbol, mod 2^89-1 (%.0f s)" % (nt, len(tests), time.time() - t0), bad == 0 and nt > 0)


def part_C():
    unit = (0, 0); dil = (0, 1); cnt = 0; ok_s = ok_d = ok_p = True
    samples = [(1, ((1, 0),)), (1, ((0, 2),)), (1, ((4, 2),)), (2, ((4, 2),)), (2, ((3, 3),)), (2, ((1, 0), (0, 2))), (2, ((2, 4),)),
               (2, ((3, 1), (4, 1))), (3, ((1, 5),)), (3, ((0, 6),)), (2, ((1, 1), (0, 1))), (1, ((2, 0), (0, 1))), (2, ((2, 2), (4, 2))),
               (0, ((1, 0), (1, 0), (1, 0))), (0, ((2, 1), (0, 0), (0, 0))), (1, ((3, 0), (4, 0))), (2, ((3, 0), (4, 0)))]
    for g, ins in samples:
        n = len(ins)
        if not os.path.exists(me._cache_path(g, n + 1)): continue
        cnt += 1
        lhs = me.corr(g, ins + (unit,)); rhs = Fr(0)
        for j in range(n):
            if ins[j][1] > 0:
                rhs += me.corr(g, ins[:j] + ((ins[j][0], ins[j][1] - 1),) + ins[j + 1:])
        ok_s = ok_s and lhs == rhs
        ok_d = ok_d and me.corr(g, ins + (dil,)) == (2 * g - 2 + n) * me.corr(g, ins)
        for perm in itertools.permutations(range(n)):
            ok_p = ok_p and me.corr(g, tuple(ins[i] for i in perm)) == me.corr(g, ins)
    check("(C) string equation on %d samples" % cnt, ok_s and cnt > 0)
    check("(C) dilaton equation on %d samples" % cnt, ok_d and cnt > 0)
    check("(C) symmetry under permutations of the insertions", ok_p)


if __name__ == "__main__":
    print("rows: symplectic %s, edge symmetry %s, edge g.f. %s" % (mr.symplectic_ok(20), mr.edge_symmetric_ok(), mr.edge_gf_ok()))
    part_A(); part_C()
    if "--noB" not in sys.argv: part_B()
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)
