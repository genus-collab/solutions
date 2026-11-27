# -*- coding: utf-8 -*-
"""fq_diag_solve.py -- the polar systems of the diagonal colour-dependent scalar frame (genus2_diag.py step 7).

cf = {(idx, col): {mono: Ser}} the Cf images of the coefficient polynomials.  Unknowns: the ten a1_i (indices 0..4) and a2_i (5..9).
  solve_full(cf, keys, R, log, check)     {all negative powers of s of all 110 coloured edge graphs = 0}: the s^-10, s^-5 coefficients
                                          are a-independent constants, non-zero -> inconsistent; reported per graph.
  solve_weak(cf, keys, R, log, check)     the a-dependent polar coefficients only: linear forms in the a1_i (s^-6, s^-1) -> Gaussian
                                          elimination; quadratic forms (s^-2) restricted to the kernel; exact solution set.
  solve_useries(cf, keys, R, log, check)  all orders (to the series precision): which frames leave EVERY coloured term a u-series
                                          (weight classes 4, 3, 2 are distinct, so each a-weight part must vanish separately).
"""
import itertools
from fq_diag_lib import local_vars, local_monos, mono_weight, mono_global, gmono_name, fp_rref, fp_nullspace, GNAMES
from fq_scalar_lib import ser_coef


def _rows_linear(cf, keys, exps, R):
    """[(key, e, row)] with row[c] = coefficient of a1_c (a-weight-1 local monomials) at s^e"""
    out = []
    for key in keys:
        vs = local_vars(key[1])
        for e in exps:
            row = [0] * 10
            for mono in local_monos(key[1]):
                if mono_weight(mono, vs) != 1: continue
                g = mono_global(mono, vs); c = g.index(1)
                row[c] = (row[c] + ser_coef(cf[key][mono], e)) % R.p
            out.append((key, e, row))
    return out


def _quads(cf, keys, exps, R):
    """[(key, e, {g: coef})] a-weight-2 (quadratic in the a1_i) parts at s^e, g = global exponent 10-tuple"""
    out = []
    for key in keys:
        vs = local_vars(key[1])
        for e in exps:
            Q = {}
            for mono in local_monos(key[1]):
                if mono_weight(mono, vs) != 2: continue
                x = ser_coef(cf[key][mono], e)
                if x: g = mono_global(mono, vs); Q[g] = (Q.get(g, 0) + x) % R.p
            if Q: out.append((key, e, Q))
    return out


def _rank(rows, p):
    return len(fp_rref([r for r in rows], p)[0]) if rows else 0


def _fmt_form(row, R, names=GNAMES):
    parts = ["(%s) %s" % (R.fmt(x), n) for x, n in zip(row, names) if x]
    return " + ".join(parts) if parts else "0"


def solve_full(cf, keys, R, log, check):
    edge = [k for k in keys if k[0] != 1]
    const = lambda key: next(m for m in local_monos(key[1]) if sum(m) == 0)
    nz10 = [k for k in edge if ser_coef(cf[k][const(k)], -10)]
    nz5 = [k for k in edge if ser_coef(cf[k][const(k)], -5)]
    log("   full system over the %d coloured edge graphs: the s^-10 and s^-5 coefficients are the a-independent constants (point-frame term);" % len(edge))
    log("   non-zero s^-10 coefficient for %d of %d, non-zero s^-5 coefficient for %d of %d" % (len(nz10), len(edge), len(nz5), len(edge)))
    log("   RESULT: the full system is %s (no choice of the ten constants changes an a-independent coefficient)"
        % ("INCONSISTENT: non-zero a-independent s^-10 coefficient for every coloured edge graph" if len(nz10) == len(edge)
           else "inconsistent for %d edge graphs; %d edge graphs have a vanishing s^-10 constant" % (len(nz10), len(edge) - len(nz10))))
    return len(nz10), len(nz5)


def _substitute_quads(quads, basis, p):
    """restrict the quadratic forms to a1 = Sum_j t_j basis[j]: forms in the t_j as {(j, l): coef}, j <= l"""
    d = len(basis); out = []
    for key, e, Q in quads:
        F = {}
        for g, c in Q.items():
            idx = [i for i in range(10) for _ in range(g[i])]           # the two variable indices
            i1, i2 = idx[0], idx[1]
            for j in range(d):
                for l in range(d):
                    v = c * basis[j][i1] % p * basis[l][i2] % p
                    if v:
                        kk = (min(j, l), max(j, l)); F[kk] = (F.get(kk, 0) + v) % p
        F = {k: v for k, v in F.items() if v}
        if F: out.append((key, e, F))
    return out


def solve_weak(cf, keys, R, log, check, label="weak system"):
    p = R.p
    edge = [k for k in keys if k[0] != 1]
    lin = _rows_linear(cf, edge, (-6, -1), R)
    quads = _quads(cf, edge, (-2,), R)
    nzlin = [(k, e, r) for (k, e, r) in lin if any(r)]
    log("   %s: %d linear forms in (a1_0..a1_4) from s^-6 and s^-1 of the %d coloured edge graphs (%d non-zero), %d quadratic forms from s^-2"
        % (label, len(lin), len(edge), len(nzlin), len(quads)))
    rows = [r for (_, _, r) in nzlin]
    A, piv = fp_rref(rows, p)
    log("   Gaussian elimination of the linear part: rank %d (of 10 unknowns); pivot variables %s" % (len(A), [GNAMES[c] for c in piv]))
    for r in A: log("      %s = 0" % _fmt_form(r, R))
    # structure of the linear part
    for sub, name in [((-6,), "s^-6 forms only"), ((-1,), "s^-1 forms only")]:
        rk = _rank([r for (_, e, r) in nzlin if e in sub], p); log("      rank of the %s: %d" % (name, rk))
    mono = [k for k in edge if len(set(k[1])) == 1]; mixed = [k for k in edge if len(set(k[1])) > 1]
    log("      rank from the %d monochrome edge graphs alone: %d; from the %d mixed colourings alone: %d"
        % (len(mono), _rank([r for (k, e, r) in nzlin if k in mono], p), len(mixed), _rank([r for (k, e, r) in nzlin if k in mixed], p)))
    for idx in sorted(set(k[0] for k in edge)):
        log("      graph %d alone (all its colourings): rank %d" % (idx, _rank([r for (k, e, r) in nzlin if k[0] == idx], p)))
    basis = fp_nullspace(rows, 10, p) if rows else [[int(i == j) for i in range(10)] for j in range(10)]
    a2free = [GNAMES[c] for c in range(5, 10) if c not in piv]
    kernel_a1 = [b for b in basis if any(b[:5])]
    log("   kernel of the linear part: dimension %d; the a2_i are %s (no polar coefficient involves them: a-weight 3 has valuation >= 2)"
        % (len(basis), "all free" if len(a2free) == 5 else "constrained: %s" % a2free))
    if not kernel_a1:
        log("   => linear part forces a1_0 = a1_1 = a1_2 = a1_3 = a1_4 = 0; the quadratic forms then vanish identically")
        sol = "a1_i = 0 for all i, a2_i arbitrary (solution set of dimension 5)"
    else:
        log("   kernel vectors with a1-components: %d; substituting a1 = Sum t_j v_j into the %d quadratic forms" % (len(kernel_a1), len(quads)))
        F = _substitute_quads(quads, kernel_a1, p)
        d = len(kernel_a1); tm = [(j, l) for j in range(d) for l in range(j, d)]
        rows2 = [[f.get(k, 0) for k in tm] for (_, _, f) in F]
        rk2 = _rank(rows2, p)
        log("   quadratic forms in t (%d non-zero), rank in the %d quadratic monomials: %d" % (len(F), len(tm), rk2))
        if rk2 == len(tm): sol = "a1_i = 0 for all i (the quadratics span all t_j t_l, so t = 0), a2_i arbitrary"; log("   => " + sol)
        else:
            sol = "a1 in a non-trivial variety inside a kernel of dimension %d (quadratic rank %d < %d): NOT resolved here" % (d, rk2, len(tm))
            for key, e, f in F[:10]: log("      %s s^%d: %s" % (str(key), e, " + ".join("(%s) t%d t%d" % (R.fmt(v), j, l) for (j, l), v in sorted(f.items()))))
    log("   => exact solution set of the %s: %s" % (label, sol))
    return dict(rank=len(A), pivots=piv, kernel_dim=len(basis), sol=sol, nzlin=len(nzlin), nquad=len(quads))


def solve_useries(cf, keys, R, log, check):
    """all orders to the precision: a-weight-1 series coefficients (class 4) give linear forms in a1; given a1 = 0 the a-weight-3
    a2-part (p, q series, class 2) gives linear forms in a2.  Reported for all 115 coloured graphs and for the 110 edge graphs."""
    p = R.p
    for keyset, name in [(keys, "all 115 coloured graphs"), ([k for k in keys if k[0] != 1], "the 110 coloured edge graphs")]:
        rows1 = []; rows3 = []
        for key in keyset:
            vs = local_vars(key[1])
            for mono in local_monos(key[1]):
                W = mono_weight(mono, vs); S = cf[key][mono]
                if W == 1:
                    c = mono_global(mono, vs).index(1)
                    for i, x in enumerate(S.c):
                        if x: r = [0] * 10; r[c] = x; rows1.append(r)
                elif W == 3 and sum(mono_global(mono, vs)[5:]):
                    c = mono_global(mono, vs).index(1)
                    for i, x in enumerate(S.c):
                        if x: r = [0] * 10; r[c] = x; rows3.append(r)
        A1, piv1 = fp_rref(rows1, p); A3, piv3 = fp_rref(rows3, p)
        log("   u-series system over %s: a-weight-1 series coefficients: %d non-zero, rank %d, forced zero: %s" % (name, len(rows1), len(A1), [GNAMES[c] for c in piv1]))
        log("      given a1 = 0, the a2-part (a-weight 3, p and q series): %d non-zero coefficients, rank %d, forced zero: %s"
            % (len(rows3), len(A3), [GNAMES[c] for c in piv3]))
        if name.startswith("all"):
            log("   RESULT: every coloured term is a u-series %s" % ("only in the point frame (all ten a_k^(i) = 0)" if len(A1) == 5 and len(A3) == 5
                else "on a set of dimension %d (a1 rank %d, a2 rank %d)" % (10 - len(A1) - len(A3), len(A1), len(A3))))
        else:
            log("      (the a2 constraints over the edge graphs alone: %s)" % ("a2_i = 0 for all i" if len(A3) == 5 else "rank %d" % len(A3)))
    # the a-weight-2 and a1-cubic parts are not needed once a1 = 0; report their non-vanishing counts for completeness
    n2 = sum(1 for key in keys for mono in local_monos(key[1]) if mono_weight(mono, local_vars(key[1])) == 2 and cf[key][mono].c)
    log("   (a-weight-2 series non-zero for %d local monomials over the 115 coloured graphs; irrelevant once a1 = 0)" % n2)
