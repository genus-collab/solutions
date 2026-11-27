# -*- coding: utf-8 -*-
"""fq_scalar_report.py -- tables and solves for genus2_scalar.py (PLAN_scalar_2026-11-04.md steps 4-6).

Input everywhere:  cf = {(idx, col): [S0..S4]}  the Cf images (Laurent Ser in s) of T0..T4 of every coloured genus-2 graph,
reps = the 23 colour-shift orbit representatives (col[0] = 0), R = FpRing, log = line logger.
  coef_table(cf, reps, R, log)   the 23-row table (T0_-10, T0_-5, T1_-6, T1_-1, T2_-2): residues/rationals, Q(zeta)
                                 coordinates via the Galois conjugates col -> m col, per-graph colour sums, totals row.
                                 Returns comparable data (rationals) for the cross-prime check.
  solve_full(cf, R, log)         full system {all negative-power coefficients of the 110 coloured edge graphs = 0}.
  solve_weak(cf, reps, R, log)   weak system a1 T1_-6 = a1^2 T2_-2 = a1 T1_-1 = 0 and the u-series statement.
  solve_qdep(cf, reps, R, log)   OPTIONAL q-dependent leading system T0_-10 + al T1_-6 + al^2 T2_-2 + al^3 T3_2 + be T4_2 = 0.
"""
from fractions import Fraction as Fr
from fq_scalar_lib import (ser_coef, qzeta_coords, qzeta_ratrec, qzeta_fmt, galois_colouring, p_gcd, p_roots, p_sub, p_scale,
                           p_trim, p_eval, p_fmt, p_monic)

COLS = [(0, -10), (0, -5), (1, -6), (1, -1), (2, -2)]          # (k, s-exponent)
COLNAMES = ["T0_-10", "T0_-5", "T1_-6", "T1_-1", "T2_-2"]


def entry(cf, key, k, e):
    return ser_coef(cf[key][k], e)


def qz(cf, key, k, e, R):
    """Q(zeta) coordinates of the entry via its Galois conjugates (colourings m * col)"""
    idx, col = key
    xs = [entry(cf, (idx, galois_colouring(col, m)), k, e) for m in range(1, 5)]
    return qzeta_coords(R, xs)


def coef_table(cf, reps, R, log):
    p = R.p
    log("   %-6s %-8s | %s" % ("graph", "colour", " | ".join("%-28s" % n for n in COLNAMES)))
    log("   (each entry: residue-based R.fmt  [Q(zeta) coordinates c0 + c1 z + c2 z^2 + c3 z^3 from the Galois conjugates])")
    data = {}; nrec_fail = 0; ngal_fail = 0
    for key in reps:
        idx, col = key
        cells = []; row = []
        for (k, e) in COLS:
            x = entry(cf, key, k, e)
            cs = qz(cf, key, k, e, R)
            rs = qzeta_ratrec(R, cs)
            if rs is None: nrec_fail += 1
            # consistency: the coordinates must reproduce x itself (m = 1)
            back = sum(R.mul(c, R.z(kk)) for kk, c in enumerate(cs)) % p
            if back != x: ngal_fail += 1
            cells.append("%s [%s]" % (R.fmt(x), qzeta_fmt(R, cs)))
            row.append(rs)
        data[key] = tuple(row)
        log("   %-6d %-8s | %s" % (idx, str(col), " | ".join(cells)))
    log("   Q(zeta) reconstruction failures: %d of %d entries; Galois back-substitution failures: %d" % (nrec_fail, 5 * len(reps), ngal_fail))
    # per-graph colour sums (all colourings) and totals
    log("   per-graph colour sums (all 5^V colourings), and the totals row:")
    gsums = {}; tot = [0] * 5
    for key in cf:
        idx = key[0]
        acc = gsums.setdefault(idx, [0] * 5)
        for j, (k, e) in enumerate(COLS):
            acc[j] = (acc[j] + entry(cf, key, k, e)) % p
    for idx in sorted(gsums):
        for j in range(5): tot[j] = (tot[j] + gsums[idx][j]) % p
        log("   graph %-3d sum   | %s" % (idx, " | ".join("%-28s" % R.fmt(x) for x in gsums[idx])))
    log("   %-16s | %s" % ("TOTAL", " | ".join("%-28s" % R.fmt(x) for x in tot)))
    data['gsums'] = {idx: tuple(R.ratrec(x) for x in v) for idx, v in gsums.items()}
    data['tot'] = tuple(R.ratrec(x) for x in tot)
    return data, tot, nrec_fail, ngal_fail


def solve_full(cf, R, log):
    """all negative-power coefficients of the 110 coloured edge graphs = 0: inconsistent at s^-10 (a-independent)"""
    edge = [key for key in cf if key[0] != 1]
    nz = [key for key in edge if not R.iszero(entry(cf, key, 0, -10))]
    log("   full system: %d coloured edge graphs; the s^-10 coefficient is T0_-10 (a-independent); non-zero for %d of %d"
        % (len(edge), len(nz), len(edge)))
    ok = len(nz) == len(edge)
    log("   => the full system is %s (a-independent leading coefficients listed in the table above)" % ("INCONSISTENT" if nz else "consistent at s^-10"))
    return ok, len(nz)


def solve_weak(cf, reps, R, log):
    """a1 T1_-6 = a1^2 T2_-2 = a1 T1_-1 = 0 for all coloured edge graphs; and the u-series statement"""
    edge = [key for key in cf if key[0] != 1]
    cnt = {}
    for (k, e), name in zip(COLS[2:], COLNAMES[2:]):
        cnt[name] = sum(1 for key in edge if not R.iszero(entry(cf, key, k, e)))
    log("   weak system (a-dependent polar terms only), over the %d coloured edge graphs: non-vanishing coefficients: %s"
        % (len(edge), ", ".join("%s: %d" % (n, c) for n, c in cnt.items())))
    any_nz = any(cnt.values())
    if any_nz: log("   => solution set of the weak system: a1 = 0, a2 arbitrary (a2 never enters a polar term)")
    else: log("   => solution set of the weak system: all (a1, a2) (no a-dependent polar coefficient)")
    # u-series statement: exponents of T1, T2 are = 4, 3 mod 5; T3, T4 = 2 mod 5
    ser_nz = {}
    for k in (1, 2, 3, 4):
        ser_nz[k] = sum(1 for key in cf if cf[key][k].c)
    log("   u-series statement: non-zero Cf series among the 115 coloured graphs: T1: %d, T2: %d, T3: %d, T4: %d"
        % (ser_nz[1], ser_nz[2], ser_nz[3], ser_nz[4]))
    conds = []
    if ser_nz[1] or ser_nz[2]: conds.append("a1 = 0")
    if ser_nz[4]: conds.append("a2 = 0 (given a1 = 0)" if (ser_nz[1] or ser_nz[2]) else "a1^3 T3 + a2 T4 = 0")
    log("   => the constant covariant scalar frames in which EVERY coloured term is a u-series: %s"
        % ("only the point frame (a1, a2) = (0, 0)" if conds == ["a1 = 0", "a2 = 0 (given a1 = 0)"] else " and ".join(conds) if conds else "all"))
    return cnt, ser_nz


def solve_qdep(cf, reps, R, log):
    """leading system over the 23 representatives: T0_-10 + al T1_-6 + al^2 T2_-2 + al^3 T3_2 + be T4_2 = 0.
    Cubic in al, linear in be: eliminate be with one equation (T4_2 != 0), gcd of the remaining cubics over F_p."""
    p = R.p
    eqs = {}
    for key in reps:
        A = entry(cf, key, 0, -10); B = entry(cf, key, 1, -6); C = entry(cf, key, 2, -2)
        D = entry(cf, key, 3, 2); E = entry(cf, key, 4, 2)
        eqs[key] = ([A, B, C, D], E)                  # cubic in al, coefficient of be
    log("   leading system over %d representatives: coefficients (T0_-10, T1_-6, T2_-2, T3_2, T4_2):" % len(reps))
    for key in reps:
        c, E = eqs[key]
        log("      %-14s %s" % (str(key), ", ".join(R.fmt(x) for x in c + [E])))
    pivot = next((key for key in reps if not R.iszero(eqs[key][1])), None)
    if pivot is None:
        log("   no equation involves al_2 (all T4_2 = 0): al_2 is free; gcd of all cubics taken directly")
        polys = {key: p_trim(list(eqs[key][0])) for key in reps}
    else:
        cP, EP = eqs[pivot]; Einv = pow(EP, p - 2, p)
        log("   eliminating al_2 with equation %s: al_2 = -(%s)/T4_2" % (str(pivot), p_fmt(cP, R, "al_1")))
        polys = {}
        for key in reps:
            if key == pivot: continue
            c, E = eqs[key]
            polys[key] = p_sub(p_trim(list(c)), p_scale(cP, E * Einv % p, p), p)
    degs = {key: len(f) - 1 for key, f in polys.items()}
    log("   resulting polynomials in al_1: degrees %s" % sorted(set(degs.values())))
    g = None
    for key, f in polys.items():
        if not f: continue
        g = f if g is None else p_gcd(g, f, p)
    if g is None:
        log("   all resulting polynomials vanish identically: every al_1 works (degenerate)"); return dict(deg=None, roots=[])
    g = p_monic(g, p)
    log("   gcd over F_p: degree %d: %s" % (len(g) - 1, p_fmt(g, R, "al_1")))
    out = dict(deg=len(g) - 1, roots=[], roots_fmt=[])
    if len(g) <= 1:
        log("   => NO common root: no covariant scalar frame (constant or q-dependent, a_k in L^4k-... Q[[u]]) kills the leading s^-10 pole of all 23 representatives at once")
        return out
    roots = p_roots(g, p, R)
    log("   F_p-roots of the gcd: %d" % len(roots))
    for al in roots:
        be = None
        if pivot is not None:
            be = (-p_eval(cP, al, p)) * Einv % p
        bad = [key for key in reps if (p_eval(eqs[key][0], al, p) + (eqs[key][1] * be if be is not None else 0)) % p]
        ratl = R.fmt(al); ratb = R.fmt(be) if be is not None else "free"
        log("      al_1 = %s, al_2 = %s; all 23 equations satisfied: %s; Q(zeta)-rational (r or r*zeta^k): %s"
            % (ratl, ratb, not bad, not ratl.startswith('[')))
        out['roots'].append((al, be)); out['roots_fmt'].append((ratl, ratb))
    if not roots: log("   => the gcd has no root in F_p (so none in Q(zeta)): no common root")
    return out
