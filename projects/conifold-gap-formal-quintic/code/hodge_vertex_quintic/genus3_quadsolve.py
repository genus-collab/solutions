# -*- coding: utf-8 -*-
"""genus3_quadsolve.py -- first attack on the O(u^2) conditions in the circulant class by the low/high decomposition (genus 3).

Dimension count: a two-edge graph at genus 3 has Sum of half-edge psi-exponents <= 4, so it cannot carry two edge factors of psi-order
>= 4; hence Q_Gamma(Delta) = q_LL(Delta_L) + 2 B_LH(Delta_L, Delta_H) with Delta_L the zhat^{1,2,3} part and Delta_H the zhat^{4,5,6} part of
Delta = N_1 - N_1^taut (q_HH = 0, verified).  The genus-2 O(u) equations constrain Delta_L only (6-dimensional space L_6 inside the 8 low
circulant coordinates), the genus-3 O(u) equations are solved by Delta_H in a 5-dimensional affine space.  For fixed Delta_L the five O(u^2)
conditions are affine-linear in Delta_H: a 5 x 5 system.  RESULT (both primes): this system has rank 4 for every Delta_L, so the mechanism
"every Delta_L extends" FAILS by one condition; the KEY check of this script therefore reports FAIL by design, and the residual condition
(one hypersurface of degree <= 6 on L_6) is analysed in genus3_quadsolve2.py.  Run:  py -I genus3_quadsolve.py --prime 89|61
"""
import os, sys, json, argparse, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq3_steps import State
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE
from fq3_m import NSP, class_columns
from fq3_rec import Rec
import fq3_quadlib as L


def lincomb(vecs, coef, p):
    return [sum(c * v[k] for c, v in zip(coef, vecs)) % p for k in range(NSP)]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89")
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
    log("-" * 110); log("(D) low/high decomposition of the circulant class at O(u^2), p = 2^%s-1" % a.prime)
    cols, names = class_columns(R, 'circ'); Lidx = [j for j, nm in enumerate(names) if int(nm[1]) <= 3]; Hidx = [j for j, nm in enumerate(names) if int(nm[1]) >= 4]
    J2 = [g3['J2']['2'], g3['J2']['4']]; J3 = [g3['J']['2'], g3['J']['5']]
    dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    check("genus-2 rows vanish on the high blocks (zhat^4, zhat^5, zhat^6) of the circulant class", all(not any(dot(r, cols[j]) for r in J2) for j in Hidx))
    Hcols = [cols[j] for j in Hidx]; Lcols = [cols[j] for j in Lidx]
    qH = [Q[i].restrict([0] * NSP, Hcols) for i in IDX_2EDGE]
    check("the quadratic forms vanish identically on the 7 high circulant coordinates (q_HH = 0): a two-edge graph cannot carry two edges of psi-order >= 4", all(not any(x for r in q.A for x in r) for q in qH))
    M2 = [[dot(r, c) for c in Lcols] for r in J2]; L6 = L.compose_cols(Lcols, L.kernel_basis(M2, len(Lcols), p), p)
    log("   dim L_6 = %d (genus-2 O(u) equations on the 8 low circulant coordinates, rank %d)" % (len(L6), rank_fp(M2, p))); check("dim L_6 = 6", len(L6) == 6)
    M3H = [[dot(r, c) for c in Hcols] for r in J3]
    check("genus-3 rows have rank 2 on the 7 high coordinates (Delta_H can always absorb them)", rank_fp(M3H, p) == 2)
    def solve_high(dL):
        xH, nsH, rk = fp_solve(M3H, [(-dot(r, dL)) % p for r in J3], p)
        base = [sum(Hcols[j][v] * xH[j] for j in range(7)) % p for v in range(NSP)]
        dirs = [[sum(Hcols[j][v] * n[j] for j in range(7)) % p for v in range(NSP)] for n in nsH]
        x0 = [(nt[v] + dL[v] + base[v]) % p for v in range(NSP)]
        rows = []; rhs = []
        for i in IDX_2EDGE:
            g = Q[i].gradient(x0); rows.append([sum(g[v] * d[v] for v in range(NSP)) % p for d in dirs]); rhs.append((-Q[i].eval(x0)) % p)
        sol = fp_solve(rows, rhs, p); rk5 = rank_fp(rows, p)
        if sol[0] is None: return None, rk5, False
        return [(base[v] + sum(sol[0][k] * dirs[k][v] for k in range(len(dirs)))) % p for v in range(NSP)], rk5, True
    rng = random.Random(11); ranks = []
    for trial in range(5):
        dL = lincomb(L6, [rng.randrange(p) for _ in range(6)], p)
        check("random Delta_L #%d lies in L_6 (genus-2 rows vanish)" % trial, not any(dot(r, dL) for r in J2))
        dH, rk5, ok = solve_high(dL); ranks.append(rk5)
        log("   random Delta_L #%d: 5 x 5 system rank %d, %s" % (trial, rk5, "solved" if ok else "inconsistent"))
    for k in range(6):
        dH, rk5, ok = solve_high(L6[k]); ranks.append(rk5); log("   basis vector %d of L_6: rank %d, %s" % (k, rk5, "solved" if ok else "inconsistent"))
    check("KEY: for generic Delta_L in L_6 the 5 x 5 system for Delta_H is INVERTIBLE (then every low-order perturbation would extend to an O(u^2) solution)", all(r == 5 for r in ranks), "ranks %s: by design this FAILS when the rank is 4 (one residual condition on Delta_L, see genus3_quadsolve2)" % ranks)
    qs.put('quadsolve_ranks', ranks)
    log("QUADSOLVE ALL CHECKS PASSED" if not FAILED else "QUADSOLVE FAILURES (expected for the KEY line): %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
