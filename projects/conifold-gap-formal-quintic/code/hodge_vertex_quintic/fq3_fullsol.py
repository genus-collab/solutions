# -*- coding: utf-8 -*-
"""fq3_fullsol.py -- ABANDONED ATTEMPT (kept as a record; inconsistent at stage F3, rank 11 of 12, both primes; see OU_EXTENSION_genus3.md
s.6d, where full-class solutions are obtained instead as left dressings by constant propagators avoiding H).
Original purpose: an explicit N_1 in the FULL class (no Z/5 covariance) with non-zero low-order part Delta_L satisfying every N_1-condition
of genus 2 and 3: the four O(u) rows, the five two-edge quadrics and the nine three-edge cubics (genus3_u3.py --cls full).

Structure used (weights: an order-k coordinate has weight k - 1; a three-edge graph sees weight <= 3, a two-edge graph weight <= 4):
  * graph 32 (genus-0 trivalent centre, three genus-1 leaves) sees the orders 1, 2 only and is a sum of five cubes of linear forms
    (the genus-0 three-point tensor is diagonal in the colour basis): c_32 = Sum_i tau_i l_i(Delta_1, Delta_2)^3.  The span W of the l_i is
    found from gradients of c_32 at random points; on S = ker W (20 of the 25 coordinates of orders 1, 2) c_32 vanishes.
  * graphs 4, 18, 20 do not see order 4 and are affine-linear in order 3: with the two genus-2 rows, 5 linear equations for Delta_3.
  * graphs 9, 10, 11, 12, 19 are affine-linear in order 4; with the two genus-3 rows and the five quadrics (affine-linear in the orders
    4, 5, 6): 12 linear equations for the 35 high coordinates.
full_n1(...) returns the 75-vector n1 = n_taut + Delta; every stage is verified by direct evaluation.
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_solve
from fq3_m import NSP, BLOCKS
from fq3_quad import rank_fp, IDX_2EDGE
from fq3_cubic import job_cubic, IDX_3EDGE
from fq3_par import run_pool
import fq3_quadlib as L

G_LOW = [4, 18, 20]            # no order-4 dependence, affine-linear in order 3
G_HIGH = [9, 10, 11, 12, 19]   # affine-linear in order 4


def small_ints(n): return [((7 * k + 3) % 11) - 5 for k in range(n)]


def full_n1(R, Q, nt, J2, J3, frame4, nproc, log, check, us):
    p = R.p; dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    o = {k: list(range(off, off + len(pairs))) for k, pairs, off in BLOCKS}
    o12 = o[1] + o[2]; o3 = o[3]; o4 = o[4]; Hf = o[4] + o[5] + o[6]
    def at(delta): return [(nt[v] + delta[v]) % p for v in range(NSP)]
    def plus(delta, v, s=1):
        d = list(delta); d[v] = (d[v] + s) % p; return d
    def run(jobs, n, label):
        t = time.time(); res = dict(run_pool(job_cubic, jobs, min(n, len(jobs)), log, label=label)); log("   %s: %d evaluations in %.0f s" % (label, len(jobs), time.time() - t))
        return res
    # ---------------- stage F1: the span of the linear forms of graph 32
    rng = random.Random(41); pts = []
    for m in range(7):
        a_ = [0] * NSP
        for v in o12: a_[v] = rng.randrange(p)
        pts.append(a_)
    jobs = []
    for m, a_ in enumerate(pts):
        jobs.append((frame4, at(a_), [32], ('a', m)))
        for v in o12:
            for s in (1, 2, 3): jobs.append((frame4, at(plus(a_, v, s)), [32], ('g', m, v, s)))
    res = run(jobs, min(nproc, 10), "F1 gradients of c_32")
    inv6 = pow(6, p - 2, p); grads = []
    for m in range(7):
        f0 = res[('a', m)][32][3]; g = []
        for v in o12:
            d1, d2, d3 = [(res[('g', m, v, s)][32][3] - f0) % p for s in (1, 2, 3)]
            g.append((18 * d1 - 9 * d2 + 2 * d3) * inv6 % p)
        grads.append(g)
    r5 = rank_fp(grads[:5], p); r7 = rank_fp(grads, p)
    check("F1: the gradients of c_32 at 5 and at 7 random points span the same 5-dimensional space (c_32 is a cubic in five linear forms of the orders 1, 2)", r5 == 5 and r7 == 5, (r5, r7))
    Sb = L.kernel_basis(grads, len(o12), p); check("F1: S = common kernel of the five linear forms has dimension 20", len(Sb) == 20)
    cf = small_ints(20); d12 = [0] * NSP
    for j, v in enumerate(o12): d12[v] = sum(cf[k] * Sb[k][j] for k in range(20)) % p
    log("   (Delta_1, Delta_2) = Sum_k a_k S_k, a = %s: %d non-zero coordinates of order 1, %d of order 2" % (cf, sum(1 for v in o[1] if d12[v]), sum(1 for v in o[2] if d12[v])))
    # ---------------- stage F2: Delta_3 from the genus-2 rows and graphs 4, 18, 20
    idx2 = G_LOW + [32]
    jobs = [(frame4, at(d12), idx2, ('c0',))] + [(frame4, at(plus(d12, v)), idx2, ('c', v)) for v in o3]
    res = run(jobs, nproc, "F2 order-3 rows")
    check("F2: c_32 = 0 on S, independently of the order-3 coordinates", all(res[t][32][3] == 0 for t in res))
    rows = [[r[v] for v in o3] for r in J2]; rhs = [(-dot(r, d12)) % p for r in J2]
    for i in G_LOW:
        c0 = res[('c0',)][i][3]; rows.append([(res[('c', v)][i][3] - c0) % p for v in o3]); rhs.append((-c0) % p)
    rk = rank_fp(rows, p); x, ns, _ = fp_solve(rows, rhs, p)
    log("   F2: 2 genus-2 rows + cubics %s on the 15 order-3 coordinates: rank %d; c_Gamma(Delta_1, Delta_2) = %s" % (G_LOW, rk, {i: ("0" if res[('c0',)][i][3] == 0 else "non-zero") for i in G_LOW}))
    check("F2: the equations for Delta_3 are consistent", x is not None)
    for k, vec in enumerate(ns[:6]):
        w = (k + 1) * (-1) ** k; x = [(x[j] + w * vec[j]) % p for j in range(15)]
    dL = list(d12)
    for j, v in enumerate(o3): dL[v] = x[j]
    check("F2: Delta_L satisfies the genus-2 O(u) rows", all(dot(r, dL) == 0 for r in J2))
    # ---------------- stage F3: the orders 4, 5, 6
    x0 = at(dL)
    jobs = [(frame4, x0, G_HIGH, ('c0',))] + [(frame4, at(plus(dL, v)), G_HIGH, ('c', v)) for v in o4]
    res = run(jobs, min(nproc, 6), "F3 order-4 rows")
    rows = []; rhs = []
    for r in J3: rows.append([r[v] for v in Hf]); rhs.append((-dot(r, dL)) % p)
    for i in IDX_2EDGE:
        g = Q[i].gradient(x0); rows.append([g[v] for v in Hf]); rhs.append((-Q[i].eval(x0)) % p)
    for i in G_HIGH:
        c0 = res[('c0',)][i][3]; rows.append([((res[('c', v)][i][3] - c0) % p) if v in o4 else 0 for v in Hf]); rhs.append((-c0) % p)
    rk = rank_fp(rows, p); rkc = rank_fp(rows[7:], p); h, ns, _ = fp_solve(rows, rhs, p)
    log("   F3: 2 genus-3 rows + 5 quadrics + cubics %s on the 35 high coordinates: rank %d (cubic rows alone: %d)" % (G_HIGH, rk, rkc))
    check("F3: the equations for the orders 4, 5, 6 are consistent", h is not None)
    n1 = list(x0)
    for j, v in enumerate(Hf): n1[v] = (n1[v] + h[j]) % p
    # ---------------- verification of all nine cubics at the final point
    res = run([(frame4, n1, IDX_3EDGE, ('final',))], 1, "final cubic check")
    check("FULL N_1: all nine three-edge cubics vanish at the constructed point (direct frozen-vertex evaluation)", all(res[('final',)][i][3] == 0 for i in IDX_3EDGE))
    us.put('dL', dL); us.put('full_info', dict(rank_F3=rk, rank_F3_cubic=rkc, low_nonzero=[sum(1 for v in o[k] if dL[v]) for k in (1, 2, 3)]))
    log("   Delta_L: non-zero coordinates by order (1, 2, 3): %s of (15, 10, 15)" % [sum(1 for v in o[k] if dL[v]) for k in (1, 2, 3)])
    return n1
