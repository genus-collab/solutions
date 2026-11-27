# -*- coding: utf-8 -*-
"""genus3_quadsolve2.py -- the residual low-order hypersurface of the circulant class and explicit O(u^2) solutions (genus 3).

genus3_quadsolve showed: for every Delta_L in L_6 (low-order perturbations compatible with the genus-2 O(u) equations) the 5 x 5 system
C(Delta_L) h = -q(Delta_L) for the high part Delta_H (h in the 5-dimensional J3-kernel H_5 of the zhat^{4,5,6} coordinates) has rank 4, with
a Delta_L-DEPENDENT left kernel lambda(Delta_L) (degree 4 in Delta_L, by 4 x 4 cofactors).  So the two-edge O(u^2) conditions in the
circulant class are equivalent to  F(Delta_L) := lambda(Delta_L) . q(Delta_L) = 0,  a homogeneous polynomial of degree <= 6 on L_6 (a
hypersurface Y of dimension 4 in P^5), each point of which carries an affine line of admissible Delta_H.  This script: computes F along random
lines (interpolation, degree check), finds F_p-points of Y by univariate root finding, extends one to a full solution N_1 = N_1^taut +
Delta_L + Delta_H (all five Q_Gamma = 0, all four O(u) equations), stores it for the independent pass verification (--verify), and tests
the analogous decomposition on the FULL class (low = zhat^{1,2,3}, 40 coordinates, genus-2 kernel of dim 38; high = zhat^{4,5,6}, 35
coordinates, genus-3 kernel of dim 33): there the 5 x 33 cross-term system is expected to have rank 5, so every Delta_L extends.
Run: py -I genus3_quadsolve2.py --prime 89|61 [--verify]
"""
import os, sys, json, argparse, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq_hat import vander_inv, interp_fp
from fq3_steps import State, frame_of, hat_pass, IDX_BY_E, GD
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE, nullspace_fp
from fq3_m import NSP, class_columns, sp_mats, BLOCKS
from fq3_rec import Rec
import fq3_quadlib as L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--verify", action="store_true"); ap.add_argument("--nproc", type=int, default=24)
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R)
    qs = State(os.path.join(HERE, "g3_quad_%s.json" % a.prime)); g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % a.prime)))
    LOG = open(os.path.join(HERE, "logs", "genus3_quad_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    Q = {i: QuadForm(p, qs.d['Q'][str(i)]['c0'], qs.d['Q'][str(i)]['b'], qs.d['Q'][str(i)]['A']) for i in IDX_2EDGE}; nt = qs.d['n_taut']
    J2 = [g3['J2']['2'], g3['J2']['4']]; J3 = [g3['J']['2'], g3['J']['5']]
    rhs4 = [(-R.of(Fr(1, 750))) % p, (-R.of(Fr(-13, 8000))) % p, (-g3['A_M0']['2'][1]) % p, (-g3['A_M0']['5'][1]) % p]
    dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    log("-" * 110); log("(E) residual hypersurface of the circulant class and explicit O(u^2) solutions, p = 2^%s-1" % a.prime)
    cols, names = class_columns(R, 'circ'); Lidx = [j for j, nm in enumerate(names) if int(nm[1]) <= 3]; Hidx = [j for j, nm in enumerate(names) if int(nm[1]) >= 4]
    Lcols = [cols[j] for j in Lidx]; Hcols = [cols[j] for j in Hidx]
    L6 = L.compose_cols(Lcols, L.kernel_basis([[dot(r, c) for c in Lcols] for r in J2], len(Lcols), p), p)
    H5 = L.compose_cols(Hcols, L.kernel_basis([[dot(r, c) for c in Hcols] for r in J3], len(Hcols), p), p)
    check("dim L_6 = 6 and dim H_5 = 5", len(L6) == 6 and len(H5) == 5)
    M3H = [[dot(r, c) for c in Hcols] for r in J3]
    def base(dL):
        """nt + dL + h_part with J3 . (dL + h_part) = 0 (h_part in the high coordinates; particular solution)"""
        xH = fp_solve(M3H, [(-dot(r, dL)) % p for r in J3], p)[0]
        return [(nt[v] + dL[v] + sum(Hcols[j][v] * xH[j] for j in range(7))) % p for v in range(NSP)]
    def Cq(dL):
        x = base(dL)
        C = [[sum(Q[i].gradient(x)[v] * h[v] for v in range(NSP)) % p for h in H5] for i in IDX_2EDGE]
        return C, [Q[i].eval(x) for i in IDX_2EDGE]
    rng = random.Random(3)
    def lam_cof(C, c0):
        """left-kernel vector by cofactors: lambda_g = (-1)^g det(C without row g and column c0)"""
        out = []
        for g in range(5):
            M = [[C[r][c] for c in range(5) if c != c0] for r in range(5) if r != g]
            out.append((-1) ** g * L.det_fp(M, p) % p)
        return out
    # choose the column c0 for the cofactor formula: one that keeps rank 4 at a random point
    dL0 = (lambda cf: [sum(cf[k] * L6[k][v] for k in range(6)) % p for v in range(NSP)])([rng.randrange(p) for _ in range(6)]); C0, q0 = Cq(dL0)
    c0 = next(c for c in range(5) if rank_fp([[C0[r][cc] for cc in range(5) if cc != c] for r in range(5)], p) == 4)
    lam0 = lam_cof(C0, c0); check("cofactor left-kernel vector annihilates C at a random Delta_L (rank 4 confirmed)", any(lam0) and all(sum(lam0[g] * C0[g][c] for g in range(5)) % p == 0 for c in range(5)))
    def F(dL):
        C, q = Cq(dL); lam = lam_cof(C, c0); return sum(lam[g] * q[g] for g in range(5)) % p
    # F along random lines: degree (interpolation at 9 nodes: degree <= 8 resolved; expect <= 6), roots, points of Y
    degs = []; pts = []
    for trial in range(4):
        A_ = (lambda cf: [sum(cf[k] * L6[k][v] for k in range(6)) % p for v in range(NSP)])([rng.randrange(p) for _ in range(6)]); B_ = (lambda cf: [sum(cf[k] * L6[k][v] for k in range(6)) % p for v in range(NSP)])([rng.randrange(p) for _ in range(6)])
        nodes = list(range(9)); Vinv = vander_inv(nodes, p)
        vals = [F([(A_[v] + s * B_[v]) % p for v in range(NSP)]) for s in nodes]; poly = interp_fp(vals, Vinv, p)
        deg = max([k for k, c in enumerate(poly) if c] or [0]); degs.append(deg)
        roots = L.roots_fp(poly, p, rng)
        log("   random line #%d in L_6: F restricted has degree %d, %d root(s) in F_p" % (trial, deg, len(roots)))
        for s in roots: pts.append([(A_[v] + s * B_[v]) % p for v in range(NSP)])
    check("F has degree <= 6 along every random line (degree 4 cofactors x degree 2 quadrics)", all(d <= 6 for d in degs), degs)
    qs.put('F_degrees_on_lines', degs)
    check("F_p-points of the residual hypersurface Y found on random lines", len(pts) > 0, len(pts))
    log("   (F has %d F_p-roots on %d random lines; a degree-6 polynomial has on average one root in F_p)" % (len(pts), 4))
    # extend the first point to a full solution
    sol_n = None
    for dL in pts:
        C, q = Cq(dL); sol = fp_solve(C, [(-x) % p for x in q], p)
        if sol[0] is None: continue
        x0 = base(dL); n = [(x0[v] + sum(sol[0][k] * H5[k][v] for k in range(5))) % p for v in range(NSP)]
        okQ = all(Q[i].eval(n) == 0 for i in IDX_2EDGE); okL = [dot(r, n) for r in J2 + J3] == rhs4
        log("   point of Y: C h = -q consistent (rank %d, %d-dim fibre); Q_Gamma(N_1) all zero: %s; O(u) equations: %s" % (rank_fp(C, p), len(sol[1]), okQ, okL))
        if okQ and okL and sol_n is None: sol_n = n
    check("an explicit non-tautological circulant N_1 with all five two-edge quadrics zero and the four O(u) equations satisfied was constructed", sol_n is not None)
    if sol_n is not None:
        qs.put('circ_solution_n', sol_n)
        cc = fp_solve([[cols[j][v] for j in range(15)] for v in range(NSP)], [(sol_n[v] - nt[v]) % p for v in range(NSP)], p)[0]
        log("   circulant coordinates of Delta = N_1 - N_1^taut (%s): %s" % (names, [rec.fmt(x) for x in cc]))
    # full class: low/high decomposition
    Lf = [v for k, pairs, off in BLOCKS if k <= 3 for v in range(off, off + len(pairs))]; Hf = [v for k, pairs, off in BLOCKS if k >= 4 for v in range(off, off + len(pairs))]
    def unit(v): e = [0] * NSP; e[v] = 1; return e
    L38 = L.compose_cols([unit(v) for v in Lf], L.kernel_basis([[r[v] for v in Lf] for r in J2], len(Lf), p), p)
    H33 = L.compose_cols([unit(v) for v in Hf], L.kernel_basis([[r[v] for v in Hf] for r in J3], len(Hf), p), p)
    check("full class: dim of the genus-2 kernel in the low part = 38, of the genus-3 kernel in the high part = 33", len(L38) == 38 and len(H33) == 33)
    qHf = [Q[i].restrict([0] * NSP, [unit(v) for v in Hf]) for i in IDX_2EDGE]
    check("full class: q_HH = 0 on the 35 high coordinates", all(not any(x for r in q.A for x in r) for q in qHf))
    rks = []; nf = None
    M3Hf = [[r[v] for v in Hf] for r in J3]
    for trial in range(3):
        dLf = (lambda cf: [sum(cf[k] * L38[k][v] for k in range(38)) % p for v in range(NSP)])([rng.randrange(p) for _ in range(38)])
        xHf = fp_solve(M3Hf, [(-dot(r, dLf)) % p for r in J3], p)[0]
        xf = [(nt[v] + dLf[v]) % p for v in range(NSP)]
        for j, v in enumerate(Hf): xf[v] = (xf[v] + xHf[j]) % p
        Cf = [[sum(Q[i].gradient(xf)[v] * h[v] for v in range(NSP)) % p for h in H33] for i in IDX_2EDGE]; rks.append(rank_fp(Cf, p))
        sol = fp_solve(Cf, [(-Q[i].eval(xf)) % p for i in IDX_2EDGE], p)
        if sol[0] is not None and nf is None:
            nf = [(xf[v] + sum(sol[0][k] * H33[k][v] for k in range(33))) % p for v in range(NSP)]
    check("full class: the 5 x 33 cross-term system has rank 5 for random Delta_L (every low-order perturbation extends to an O(u^2) solution: X^full is rational over Q(zeta))", all(r == 5 for r in rks), rks)
    if nf is not None:
        check("full class: the extended random N_1 satisfies all five quadrics and the four O(u) equations", all(Q[i].eval(nf) == 0 for i in IDX_2EDGE) and [dot(r, nf) for r in J2 + J3] == rhs4)
        qs.put('full_solution_n', nf)
    if a.verify and sol_n is not None:
        frame = frame_of(p, (0, 0), 34, 6); idxs = IDX_BY_E[0] + IDX_BY_E[1] + IDX_BY_E[2] + IDX_BY_E[3]
        for tag, n in (("circ", sol_n), ("full", nf)):
            if n is None: continue
            t = time.time(); A, info = hat_pass(frame, ('M0exp', [sp_mats(n, p)], 'right'), idxs, a.nproc, log, 21, label="verify " + tag)
            log("   verification pass (%s solution; %d graphs, full hatted theory, 21-node vertex interpolation): %.0f s" % (tag, len(idxs), time.time() - t))
            for i in idxs: log("      graph %2d (E=%d): %s" % (i, GD[i]['E'], [rec.fmt(x) for x in A[i][:5]]))
            check("VERIFY (%s, independent full pass): one-edge graphs vanish through u^1 and all five two-edge graphs through u^2 (regular at O(u^2))" % tag, all(A[i][k] == 0 for i in (2, 5) for k in range(2)) and all(A[i][k] == 0 for i in IDX_2EDGE for k in range(3)))
            log("   three-edge graphs [u^3] (cubic conditions, not imposed): %s" % {i: rec.fmt(A[i][3]) for i in IDX_BY_E[3]})
            qs.put('verify_table_' + tag, {str(i): A[i] for i in A})
    log("QUADSOLVE2 ALL CHECKS PASSED" if not FAILED else "QUADSOLVE2 FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
