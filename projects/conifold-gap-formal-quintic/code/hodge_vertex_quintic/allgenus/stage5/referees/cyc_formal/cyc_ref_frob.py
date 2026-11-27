# -*- coding: utf-8 -*-
"""cyc_ref_frob.py -- referee cyc_formal: Frobenius basis of PF4 at u = 0 and the pairing B (stdlib only).

  rec_polys()       c_s(n), s = -3..1, with PF4 u^n = Sum_s c_s(n) u^(n+s) (symbolic in n; cyc_ref_jet.sm_apply)
  solve(N, a012, r) the power series a with PF4 a = -r (r: dict exponent -> coefficient), a_0..a_2 given; the conditions
                    at u^-3, u^-2, u^-1 are returned (not assumed) as 'obstruction'
  logpart(f)        r(f) = PF4(f l) - l PF4(f), the l-free part (an element of u^-3 Q[[u]])
  basis(N)          f_0, f_1, f_2 (f_k = u^k + O(u^3)), the 3x3 obstruction matrix of r on span(f_0, f_1, f_2), its kernel,
                    gamma, g = f_1 + gamma f_2, h (h_0 = h_1 = h_2 = 0), y_3 = g l + h
  pairing(P, y, z)  B(y, z) = Sum_ab D^a y P_ab D^b z with the EXACT jet matrix P_ab of cyc_ref_jet.concomitant_jet,
                    returned as (value, log_free, constant) with the precision used
The method (symbolic recursion; B from the jet matrix) differs from critic/cr5_cyc.py (residual recursion; nested formula).
"""
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U
from cyc_ref_ops import Op, REl
from cyc_ref_jet import pf4_from_large_radius, sm_apply, pn_eval, jet_eval, jets_of, concomitant_jet, det_inv

PF4 = pf4_from_large_radius()
CS = sm_apply(PF4)
PJ = concomitant_jet(PF4.A)


def solve(N, a012, r=None):
    r = r or {}
    a = {0: Fr(a012[0]), 1: Fr(a012[1]), 2: Fr(a012[2])}
    ob = []
    for m in (-3, -2, -1):
        ob.append(sum(pn_eval(CS.get(s, []), m - s) * a.get(m - s, 0) for s in range(-3, 2) if m - s >= 0) + r.get(m, 0))
    for m in range(0, N - 3):
        acc = sum(pn_eval(CS.get(s, []), m - s) * a.get(m - s, 0) for s in range(-2, 2) if m - s >= 0) + r.get(m, 0)
        a[m + 3] = -acc / pn_eval(CS[-3], m + 3)
    return Ser(0, [a[n] for n in range(N)], N), ob


def logpart(f):
    out = PF4(REl({1: f}))
    return out.F[0], out.F.get(1)


def kernel(M):
    """kernel of a 3x3 matrix over Q (rows = conditions), by elimination; returns (rank, basis list)"""
    A = [list(map(Fr, r)) for r in M]; piv = []; row = 0
    for col in range(3):
        p = next((i for i in range(row, len(A)) if A[i][col] != 0), None)
        if p is None:
            continue
        A[row], A[p] = A[p], A[row]; pv = A[row][col]; A[row] = [x / pv for x in A[row]]
        for i in range(len(A)):
            if i != row and A[i][col] != 0:
                fc = A[i][col]; A[i] = [x - fc * y for x, y in zip(A[i], A[row])]
        piv.append(col); row += 1
    free = [c for c in range(3) if c not in piv]
    ker = []
    for fcol in free:
        vec = [Fr(0)] * 3; vec[fcol] = Fr(1)
        for i, pc in enumerate(piv):
            vec[pc] = -A[i][fcol]
        ker.append(vec)
    return len(piv), ker


def basis(N):
    f = [solve(N, [int(k == i) for k in range(3)])[0] for i in range(3)]
    r = [logpart(fi)[0] for fi in f]
    Ob = [[r[i].coef(m) for i in range(3)] for m in (-3, -2, -1)]
    rank, ker = kernel(Ob)
    out = {'f': f, 'r': r, 'Ob': Ob, 'rank': rank, 'ker': ker}
    if rank == 2 and ker[0][1] != 0:
        k = [x / ker[0][1] for x in ker[0]]
        gam = k[2]; out['gamma'] = gam; out['kvec'] = k
        g = f[1] + f[2].sc(gam)
        rg = r[1] + r[2].sc(gam)
        assert rg.N >= N - 3, ('precision of r(g)', rg.N, N)
        h, ob = solve(N, [0, 0, 0], {e: rg.coef(e) for e in range(-3, rg.N)})
        out.update(g=g, h=h, h_ob=ob, y3=REl({1: g, 0: h}))
        out['Y'] = [REl.of(f[0]), REl.of(f[1]), REl.of(f[2]), out['y3']]
    return out


def pairing(y, z, P=PJ):
    b = jet_eval(P, jets_of(y), jets_of(z))
    prec = min(s.N for s in b.F.values())
    log_free = all(s.zero_below(prec) for k, s in b.F.items() if k)
    s0 = b.F.get(0, Ser(0, []))
    const = all(x == 0 for i, x in enumerate(s0.c) if s0.v + i != 0 and s0.v + i < prec)
    return (s0.coef(0) if prec > 0 else None), log_free, const, prec


def bmat(Y, P=PJ):
    B = [[None] * len(Y) for _ in Y]; ok_lf, ok_c, prec = True, True, 10 ** 9
    for i, y in enumerate(Y):
        for j, z in enumerate(Y):
            val, lf, c, p = pairing(y, z, P)
            B[i][j] = val; ok_lf &= lf; ok_c &= c; prec = min(prec, p)
    return B, ok_lf, ok_c, prec


def frame(bas, alpha, beta, mu=1, nu=1, param='project'):
    """y_0, y_1 and t: project parametrisation y_0 = mu (f_0 + alpha g + beta f_2) (fq_conifold, pb_core), y_1 = nu g;
    param = 'critic': y_0 = mu (f_0 + alpha f_1 + beta f_2) (cr5_cyc_run)."""
    f, g = bas['f'], bas['g']
    mid = g if param == 'project' else f[1]
    y0 = (f[0] + mid.sc(alpha) + f[2].sc(beta)).sc(mu)
    y1 = g.sc(nu)
    t = y1 * y0.inv()
    return y0, y1, t


def combo(bas, coeffs):
    """the element Sum c_i Y_i of S (basis f_0, f_1, f_2, y_3)"""
    out = None
    for c, y in zip(coeffs, bas['Y']):
        if c:
            t = y.sc(Fr(c)); out = t if out is None else out + t
    return out
