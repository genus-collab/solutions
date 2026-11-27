# -*- coding: utf-8 -*-
"""genus3_quadan.py -- analysis of the O(u^2) quadrics Q_Gamma (g3_quad_<p>.json) on the O(u) solution families:
(1) K_4 = kernel of the four O(u) equations in sp (75 coordinates; dim 71) and its circulant part (dim 11); the affine O(u) families are
    N_1^taut + K_4 (full) and N_1^taut + K_4^circ; (2) the restricted forms q_Gamma(t) = Q_Gamma(N_1^taut + B t): homogeneous (cone with
    vertex at the tautological point), their ranks, the dimension of their span; (3) Hilbert function of the ideal (q_3, q_6, q_7, q_8, q_17)
    in degrees 3 (and 4 for the circulant class) against a complete intersection of five quadrics: if equal, the common zero locus X in
    P(K_4) has dimension exactly dim K_4 - 1 - 5 over the algebraic closure (65 full, 5 circulant): non-tautological complex solutions exist;
(4) sectors: for every subset S of the zhat-order blocks of the circulant class, K_4 intersected with the span of S; where that has
    dimension <= 3 the common zeros of the five forms are determined exactly (gcd / resultants); the diagonal sector of the full class.
Run:  py -I genus3_quadan.py --prime 89|61   -> appends to logs/genus3_quad_<p>.log, writes keys into g3_quad_<p>.json.
"""
import os, sys, json, argparse, itertools, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_steps import State
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE, sqrt_fp
from fq3_m import NSP, class_columns, BLOCKS, SP_NAMES
from fq3_rec import Rec
import fq3_quadlib as L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); a = ap.parse_args()
    p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R)
    qs = State(os.path.join(HERE, "g3_quad_%s.json" % a.prime)); g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % a.prime)))
    LOG = open(os.path.join(HERE, "logs", "genus3_quad_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    Q = {i: QuadForm(p, qs.d['Q'][str(i)]['c0'], qs.d['Q'][str(i)]['b'], qs.d['Q'][str(i)]['A']) for i in IDX_2EDGE}
    nt = qs.d['n_taut']
    log("-" * 110); log("(A) analysis of the O(u^2) quadrics, p = 2^%s-1" % a.prime)
    rows = [g3['J2']['2'], g3['J2']['4'], g3['J']['2'], g3['J']['5']]
    rhs = [(-R.of(Fr(1, 750))) % p, (-R.of(Fr(-13, 8000))) % p, (-g3['A_M0']['2'][1]) % p, (-g3['A_M0']['5'][1]) % p]
    check("N_1^taut satisfies the four O(u) equations", all(sum(r[v] * nt[v] for v in range(NSP)) % p == b for r, b in zip(rows, rhs)))
    K4 = L.kernel_basis(rows, NSP, p); log("   dim K_4 (full class) = %d" % len(K4)); check("dim K_4 = 71", len(K4) == 71)
    qf, bad = L.restricted_hom(Q, nt, K4)
    check("full class: every restricted form is homogeneous in t (cone with vertex N_1^taut): no constant or linear terms", bad == 0)
    ranks_full = {i: rank_fp(qf[i].A, p) for i in IDX_2EDGE}; log("   ranks of the five forms on K_4 (71 x 71): %s" % ranks_full)
    rng = random.Random(5); cf = {i: rng.randrange(p) for i in IDX_2EDGE}; comb = [[sum(cf[i] * qf[i].A[r][c] for i in IDX_2EDGE) % p for c in range(71)] for r in range(71)]
    log("   rank of a random linear combination: %d; dimension of the span of the five forms: %d" % (rank_fp(comb, p), L.span_dim_forms([qf[i] for i in IDX_2EDGE], p)))
    qs.put('full_ranks', ranks_full)
    r3, nm3 = L.ideal_rank_degree([qf[i] for i in IDX_2EDGE], 3, p); e3 = nm3 - L.expected_ci(71, 3, 5)
    check("full class: Hilbert function in degree 3 equals the complete-intersection value (ideal rank %d = 5 x 71: no linear syzygies; NOT a proof of the dimension)" % r3, r3 == e3, (r3, e3))
    # circulant class
    cols, names = class_columns(R, 'circ'); rows_c = [[sum(r[v] * col[v] for v in range(NSP)) % p for col in cols] for r in rows]
    Kc = L.kernel_basis(rows_c, len(cols), p); Bc = L.compose_cols(cols, Kc, p)
    log("   dim K_4^circ = %d (rank of the 4 x 15 restricted system %d)" % (len(Kc), rank_fp(rows_c, p))); check("dim K_4^circ = 11", len(Kc) == 11)
    qc, badc = L.restricted_hom(Q, nt, Bc)
    check("circulant class: restricted forms homogeneous", badc == 0)
    ranks_c = {i: rank_fp(qc[i].A, p) for i in IDX_2EDGE}; log("   ranks of the five forms on K_4^circ (11 x 11): %s" % ranks_c)
    cf = {i: rng.randrange(p) for i in IDX_2EDGE}; comb = [[sum(cf[i] * qc[i].A[r][c] for i in IDX_2EDGE) % p for c in range(11)] for r in range(11)]
    log("   rank of a random combination: %d; span dimension: %d" % (rank_fp(comb, p), L.span_dim_forms([qc[i] for i in IDX_2EDGE], p)))
    qs.put('circ_ranks', ranks_c); qs.put('circ_forms', {str(i): qc[i].A for i in IDX_2EDGE}); qs.put('circ_kernel_coords', Kc)
    hf = {}
    for d in (3, 4):
        r, nm = L.ideal_rank_degree([qc[i] for i in IDX_2EDGE], d, p); hf[d] = (nm - r, L.expected_ci(11, d, 5))
        check("circulant class: Hilbert function H(%d) = %d equals the complete-intersection value %d (consistency, not a dimension proof)" % (d, nm - r, L.expected_ci(11, d, 5)), nm - r == L.expected_ci(11, d, 5))
    qs.put('circ_hilbert', hf)
    log("   => the Hilbert function agrees with a complete intersection in degrees 3 and 4 (no linear syzygies, only the Koszul quadratic ones); dim X^circ >= 5 holds in any case and the exact dimension of the main component follows from the low/high decomposition (genus3_quadsolve2)")
    # sectors of the circulant class
    log("(B) sectors: K_4 intersected with the span of subsets of zhat-order blocks of the circulant class; exact common zeros where dim <= 3")
    blocks = L.circ_block_indices(names); found = []
    for r_ in range(1, 7):
        for S in itertools.combinations(sorted(blocks), r_):
            idxs = [j for k in S for j in blocks[k]]
            sub = [cols[j] for j in idxs]
            rows_s = [[sum(r[v] * col[v] for v in range(NSP)) % p for col in sub] for r in rows]
            Ks = L.kernel_basis(rows_s, len(sub), p)
            if not Ks: continue
            Bs = L.compose_cols(sub, Ks, p); qs_, badS = L.restricted_hom(Q, nt, Bs)
            dim = len(Ks)
            if dim <= 3:
                pts, how = L.common_zeros_small([qs_[i] for i in IDX_2EDGE], p)
                msg = "blocks %s: dim %d, common zeros %s (%s)" % (S, dim, pts, how)
                if pts and pts != []: found.append((S, dim, str(pts)))
            else:
                rk = [rank_fp(qs_[i].A, p) for i in IDX_2EDGE]; msg = "blocks %s: dim %d, form ranks %s (not solved: dim > 3)" % (S, dim, rk)
            log("      " + msg)
    qs.put('circ_sector_solutions', found)
    log("   block sectors of dimension <= 3 with common zeros (all of them are pure-high directions, on which q_HH = 0, or their mixtures): %s" % found)
    check("the common zeros found in block sectors of dimension <= 3 all involve the high blocks 4, 5, 6 (no low-order-only sector solution)", all(any(k >= 4 for k in f[0]) for f in found), found)
    # the diagonal sector of the full class (colour-dependent scalars, not Z/5-covariant)
    dcols, dnames = L.diag_columns(R); rows_d = [[sum(r[v] * col[v] for v in range(NSP)) % p for col in dcols] for r in rows]
    Kd = L.kernel_basis(rows_d, 15, p); Bd = L.compose_cols(dcols, Kd, p); qd, badd = L.restricted_hom(Q, nt, Bd)
    log("(C) diagonal sector (15 colour-dependent scalar coordinates, odd orders): dim K_4^diag = %d; form ranks %s; span %d" % (len(Kd), [rank_fp(qd[i].A, p) for i in IDX_2EDGE], L.span_dim_forms([qd[i] for i in IDX_2EDGE], p)))
    r3d, nm3d = L.ideal_rank_degree([qd[i] for i in IDX_2EDGE], 3, p)
    log("   diagonal sector Hilbert function H(3) = %d (complete intersection: %d)" % (nm3d - r3d, L.expected_ci(len(Kd), 3, 5)))
    qs.put('diag_dim', len(Kd)); qs.put('diag_hilbert3', (nm3d - r3d, L.expected_ci(len(Kd), 3, 5)))
    # per-block sub-sectors of the diagonal sector (k = 1, 3, 5 separately, pairs)
    dblocks = {}
    for j, nm in enumerate(dnames): dblocks.setdefault(int(nm[1]), []).append(j)
    for r_ in (1, 2):
        for S in itertools.combinations(sorted(dblocks), r_):
            idxs = [j for k in S for j in dblocks[k]]; sub = [dcols[j] for j in idxs]
            rows_s = [[sum(r[v] * col[v] for v in range(NSP)) % p for col in sub] for r in rows]; Ks = L.kernel_basis(rows_s, len(sub), p)
            if not Ks: log("      diagonal blocks %s: K_4 intersection trivial" % (S,)); continue
            Bs = L.compose_cols(sub, Ks, p); qs_, _ = L.restricted_hom(Q, nt, Bs)
            if len(Ks) <= 3: log("      diagonal blocks %s: dim %d, common zeros %s" % (S, len(Ks), L.common_zeros_small([qs_[i] for i in IDX_2EDGE], p)))
            else: log("      diagonal blocks %s: dim %d, ranks %s" % (S, len(Ks), [rank_fp(qs_[i].A, p) for i in IDX_2EDGE]))
    log("QUADAN ALL CHECKS PASSED" if not FAILED else "QUADAN FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
