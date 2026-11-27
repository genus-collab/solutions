# -*- coding: utf-8 -*-
"""genus3_quadsolve3.py -- the structure of the residual condition found by the referees, and an explicit Q(zeta)-rational, Z/5-covariant,
non-tautological O(u^2) solution (genus 3, circulant class).

Facts (both primes): on K_4^circ the rank-2 forms factor, q_6 = c_1(0) . ell(Delta) (the loop at the trivalent genus-0 vertex of graph 6 sees
only the scalar zhat^1 coordinate) and q_17 = m_1(Delta_L) . m_2(Delta) (chain g1-g1-g1), with m_1 a linear form in the LOW circulant
coordinates only.  The residual sextic of genus3_quadsolve2 is F = c_1(0) . m_1 . F_4: the factor c_1(0) is spurious (there the 5 x 5 system
drops to rank 3 and is generically inconsistent), m_1 is genuine: on the hyperplane {m_1 = 0} of L_6 the equation of graph 17 holds for
every Delta_H and the other four have rank 4 on the four visible zhat^{4,5} directions, so EVERY point of {m_1 = 0} extends (uniquely up to the
invisible zhat^6 direction) to an O(u^2) solution.  m_1 is defined over Q(zeta) (its coefficients, in a normalisation, are rationals times
powers of zeta), so Q(zeta)-rational solutions form a 5-parameter family.  This script constructs m_1 from q_17, checks the facts, builds the
solution with the small integer free coordinates (c_2(1), c_2(2), c_3(0), c_3(1), c_3(2)) = (1, 2, -1, 3, 1) (the remaining low coordinates
c_1(0), c_1(1), c_1(2) fixed by the two genus-2 rows and m_1 = 0; zhat^6 part set to zero), reconstructs its 15 circulant coordinates in Q(zeta)
where possible, stores it (key 'rational_solution_n') and optionally verifies it by an independent full pass (--verify).
Run:  py -I genus3_quadsolve3.py --prime 89|61 [--verify]
"""
import os, sys, json, argparse, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq3_steps import State, frame_of, hat_pass, IDX_BY_E, GD
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE, sqrt_fp
from fq3_m import NSP, class_columns, sp_mats
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
    log("-" * 110); log("(F) factorisation of q_6, q_17 and the Q(zeta)-rational hyperplane {m_1 = 0}; explicit rational solution, p = 2^%s-1" % a.prime)
    cols, names = class_columns(R, 'circ')
    # circulant coordinates -> 75-vector, and the four rows on the 15 circulant coordinates
    def circ_to_sp(c): return [sum(cols[j][v] * c[j] for j in range(15)) % p for v in range(NSP)]
    rows15 = [[dot(r, cols[j]) for j in range(15)] for r in J2 + J3]
    # K_4^circ: kernel of the four rows on the 15 circulant coordinates (11-dim); Delta = circ_to_sp(Kc . t)
    Kc = L.kernel_basis(rows15, 15, p); check("dim K_4^circ = 11", len(Kc) == 11)
    Bc = [circ_to_sp(k) for k in Kc]                                   # 75-vectors
    def hom11(i):
        q = Q[i].restrict(nt, Bc); assert q.c0 == 0 and not any(q.b); return q.A
    # pure-high sub-space of K_4^circ (Delta_L = 0): t with Kc . t supported on the high coordinates
    hi = [j for j, nm in enumerate(names) if int(nm[1]) >= 4]; lo = [j for j, nm in enumerate(names) if int(nm[1]) <= 3]
    Mlow = [[Kc[t][j] for t in range(11)] for j in lo]; Tpure = L.kernel_basis(Mlow, 11, p); check("pure-high part of K_4^circ has dimension 5", len(Tpure) == 5)
    def factor_rank2(A, n):
        rows_nz = [r for r in A if any(r)]; basis = []
        for r in rows_nz:
            if rank_fp(basis + [r], p) > len(basis): basis.append(r)
        assert len(basis) == 2, len(basis)
        u, v = basis; eqs = []; rhs = []
        for i in range(n):
            for j in range(n):
                eqs.append([u[i] * u[j] % p, (u[i] * v[j] + v[i] * u[j]) % p, v[i] * v[j] % p]); rhs.append(A[i][j])
        alpha, beta, gamma = fp_solve(eqs, rhs, p)[0]
        D = (beta * beta - alpha * gamma) % p; r = sqrt_fp(D, p); assert r is not None, "discriminant not a square"
        if alpha:
            t1 = (-beta + r) * pow(alpha, p - 2, p) % p; t2 = (-beta - r) * pow(alpha, p - 2, p) % p
            return [(u[k] - t1 * v[k]) % p for k in range(n)], [(u[k] - t2 * v[k]) % p for k in range(n)], alpha
        return list(v), [(2 * beta * u[k] + gamma * v[k]) % p for k in range(n)], 1
    def check_factor(A, m1, m2, scale, n):
        rng = random.Random(9); ok = True
        for _ in range(3):
            c = [rng.randrange(p) for _ in range(n)]
            lhs = sum(c[i] * A[i][j] * c[j] for i in range(n) for j in range(n)) % p
            ok = ok and lhs == scale * sum(m1[k] * c[k] for k in range(n)) * sum(m2[k] * c[k] for k in range(n)) % p
        return ok
    m1_17 = None
    for g in (6, 17):
        A = hom11(g); check("q_%d has rank 2 on K_4^circ" % g, rank_fp(A, p) == 2)
        m1, m2, sc = factor_rank2(A, 11); check("q_%d = a . m . m' with linear forms over F_p on K_4^circ" % g, check_factor(A, m1, m2, sc, 11))
        lowonly = [m for m in (m1, m2) if all(sum(m[t] * w[t] for t in range(11)) % p == 0 for w in Tpure)]
        log("   q_%d: number of factors vanishing on the pure-high directions (i.e. depending on Delta_L only): %d" % (g, len(lowonly)))
        check("q_%d has exactly one factor depending on Delta_L only" % g, len(lowonly) == 1)
        if g == 17: m1_17 = lowonly[0]
        else:
            # the Delta_L-only factor of q_6 is c_1(0): it must be proportional to the functional t -> (Kc . t)[c1_0]
            f6 = lowonly[0]; c10 = [Kc[t][0] for t in range(11)]
            check("the Delta_L-only factor of q_6 is c_1(0)", rank_fp([f6, c10], p) == 1)
    # m_1 as a functional on the 15 circulant coordinates: extend from K_4^circ (choose the representative vanishing on the pivot complement)
    # solve: find row vector m on 15 coords with m . Kc[t] = m1_17[t] for all t and m supported on the low coordinates (8 unknowns, 11 equations)
    eqs = [[Kc[t][j] for j in lo] for t in range(11)]; msol = fp_solve(eqs, m1_17, p)
    check("m_1 extends to a functional of the 8 low circulant coordinates (consistent)", msol[0] is not None)
    m1n = [0] * 15
    for j, val in zip(lo, msol[0]): m1n[j] = val
    piv = next(k for k in range(15) if m1n[k]); inv = pow(m1n[piv], p - 2, p); m1n = [x * inv % p for x in m1n]
    log("   m_1 (graph 17, Delta_L-only factor) normalised at %s: %s" % (names[piv], ["%s: %s" % (names[k], rec.fmt(m1n[k])) for k in range(15) if m1n[k]]))
    qs.put('m1_circ_coords', m1n)
    # explicit solution: t in K_4^circ with m_1 = 0 and the five free low coordinates fixed; then the 4 visible equations
    free = {'c2_1': 1, 'c2_2': 2, 'c3_0': -1, 'c3_1': 3, 'c3_2': 1}
    lin = [[sum(m1n[j] * Kc[t][j] for j in range(15)) % p for t in range(11)]] + [[Kc[t][names.index(nm)] for t in range(11)] for nm in free]
    rhs_l = [0] + [val % p for val in free.values()]
    t0, tns, rk = fp_solve(lin, rhs_l, p); check("m_1 = 0 and the five free coordinates cut a 5-dimensional affine family in K_4^circ", t0 is not None and len(tns) == 5)
    check("the five remaining directions are pure-high", all(all(sum(Kc[t][j] * w[t] for t in range(11)) % p == 0 for j in lo) for w in tns))
    x0 = [(nt[v] + sum(Bc[t][v] * t0[t] for t in range(11))) % p for v in range(NSP)]
    dirs = [[sum(Bc[t][v] * w[t] for t in range(11)) % p for v in range(NSP)] for w in tns]
    rowsC = []; rhsC = []
    for i in IDX_2EDGE:
        g = Q[i].gradient(x0); rowsC.append([sum(g[v] * d[v] for v in range(NSP)) % p for d in dirs]); rhsC.append((-Q[i].eval(x0)) % p)
    check("on {m_1 = 0}: the equation of graph 17 holds at every point of the fibre (row and rhs zero)", not any(rowsC[IDX_2EDGE.index(17)]) and rhsC[IDX_2EDGE.index(17)] == 0)
    solC = fp_solve(rowsC, rhsC, p); check("the other four equations have rank 4 on the five fibre directions and are consistent (solution = line, the invisible zhat^6 direction)", solC[0] is not None and len(solC[1]) == 1 and rank_fp(rowsC, p) == 4)
    # fix the invisible direction by c6_2 = 0
    sline = solC[0]; w6 = solC[1][0]
    def coords(sv):
        tt = [(t0[t] + sum(sv[k] * tns[k][t] for k in range(5))) % p for t in range(11)]
        return [sum(Kc[t][j] * tt[t] for t in range(11)) % p for j in range(15)]
    c_a = coords(sline); c_b = coords([(sline[k] + w6[k]) % p for k in range(5)]); j62 = names.index('c6_2')
    if (c_b[j62] - c_a[j62]) % p:
        lam = (-c_a[j62]) * pow((c_b[j62] - c_a[j62]) % p, p - 2, p) % p; sfin = [(sline[k] + lam * w6[k]) % p for k in range(5)]
    else: sfin = sline
    cfin = coords(sfin); n = [(nt[v] + circ_to_sp(cfin)[v]) % p for v in range(NSP)]
    check("the explicit N_1 = N_1^taut + Delta satisfies all five two-edge quadrics and the four O(u) equations", all(Q[i].eval(n) == 0 for i in IDX_2EDGE) and [dot(r, n) for r in J2 + J3] == rhs4)
    qs.put('rational_solution_n', n); qs.put('rational_solution_circ', cfin)
    log("   circulant coordinates of Delta (Q(zeta) where reconstructible; free (c2_1, c2_2, c3_0, c3_1, c3_2) = (1, 2, -1, 3, 1), c6_2 = 0): %s" % ["%s: %s" % (names[k], rec.fmt(cfin[k])) for k in range(15)])
    if a.verify:
        frame = frame_of(p, (0, 0), 34, 6); idxs = IDX_BY_E[0] + IDX_BY_E[1] + IDX_BY_E[2] + IDX_BY_E[3]
        t = time.time(); A, info = hat_pass(frame, ('M0exp', [sp_mats(n, p)], 'right'), idxs, a.nproc, log, 21, label="verify rational")
        log("   verification pass (%d graphs, full hatted theory, 21-node vertex interpolation): %.0f s" % (len(idxs), time.time() - t))
        for i in idxs: log("      graph %2d (E=%d): %s" % (i, GD[i]['E'], [rec.fmt(x) for x in A[i][:5]]))
        check("VERIFY (independent full pass): one-edge graphs vanish through u^1 and all five two-edge graphs through u^2 for the explicit Q(zeta)-rational circulant N_1", all(A[i][k] == 0 for i in (2, 5) for k in range(2)) and all(A[i][k] == 0 for i in IDX_2EDGE for k in range(3)))
        qs.put('verify_table_rational', {str(i): A[i] for i in A})
    log("QUADSOLVE3 ALL CHECKS PASSED" if not FAILED else "QUADSOLVE3 FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
