# -*- coding: utf-8 -*-
"""fq_diag_report.py -- Cf stage, valuation checks, tables and the comparison with the covariant scan, for genus2_diag.py.

Input: T = {(idx, col): {mono: NEl}} (local monomials of fq_diag_lib), cf = {(idx, col): {mono: Ser}} their Cf images.
  cf_all(T, ring, imu0)                     Cf of every coefficient (gsum_q.to_sseries)
  exact_classes(T, keys, log, check)        exact exponent classes {-m mod 5} and max L-degree per a-weight (no truncation)
  val_checks(cf, keys, R, log, check)       series valuation bounds and classes per a-weight
  zero_pattern(T, keys, log)                which local monomials are non-zero, per orbit representative
  polar_entries(key, cf)                    [(label, mono, exponent)] of all polar coefficients of a coloured term
  table_reps(cf, reps, R, log, check)       the 23-row table with Q(zeta) coordinates via the Galois conjugates
  table_all(cf, keys, R, log)               compact residue table of the polar coefficients of all 115 coloured graphs
  parse_prev_log(path, p)                   the '(6) leading system' block of logs/genus2_scalar.log for prime p
  compare_prev(cf, reps, R, prev, log, check)  covariant specialisation a_k^(i) = zeta^-i(2k-1) a_k reproduces the previous scan
"""
import os, re
from fractions import Fraction as Fr
from fq_diag_lib import (local_vars, local_monos, mono_weight, mono_global, local_name, gmono_name, gmono_eval, covariant_point)
from fq_scalar_lib import ser_coef, ser_classes, ser_val, qzeta_coords, qzeta_ratrec, qzeta_fmt, galois_colouring
from gsum_q import to_sseries, ser_fmt

VAL_MIN = {0: -10, 1: -6, 2: -2, 3: 2}; CLASSES = {0: {0}, 1: {4}, 2: {3}, 3: {2}}; MAXL = {0: 10, 1: 6, 2: 2, 3: -2}
POLAR = {0: (-10, -5), 1: (-6, -1), 2: (-2,)}          # polar s-exponents per a-weight


def cf_all(T, ring, imu0):
    return {key: {mono: to_sseries(el, ring, imu0) for mono, el in T[key].items()} for key in T}


def exact_classes(T, keys, log, check):
    """exact (all orders): exponent classes of Cf(T_mono) are {-m mod 5 : m an L-power of the numerator}; val >= -max L-degree"""
    check("all coefficient NEls are Lam-free (no log q)", all(not el.has('Lam') for key in keys for el in T[key].values()))
    for W in range(4):
        cls = set(); maxL = None
        for key in keys:
            vs = local_vars(key[1])
            for mono, el in T[key].items():
                if mono_weight(mono, vs) != W or el.iszero(): continue
                cls |= set((-m) % 5 for (m, *_r) in el.num)
                hi = el.Ldeg()[1]; maxL = hi if maxL is None else max(maxL, hi)
        check("a-weight %d: exact exponent classes {-m mod 5 : m an L-power} = %s" % (W, sorted(CLASSES[W])), cls == CLASSES[W], sorted(cls))
        check("a-weight %d: max L-degree over all coefficients = %d (so val Cf >= %d exactly)" % (W, MAXL[W], -MAXL[W]), maxL == MAXL[W], maxL)


def val_checks(cf, keys, R, log, check):
    for W in range(4):
        vals = []; nz = 0; ntot = 0; bad = []; badc = []
        for key in keys:
            vs = local_vars(key[1])
            for mono, S in cf[key].items():
                if mono_weight(mono, vs) != W: continue
                ntot += 1; v = ser_val(S)
                if v is None: continue
                nz += 1; vals.append(v)
                if v < VAL_MIN[W]: bad.append((key, local_name(mono, vs), v))
                if not ser_classes(S, R) <= CLASSES[W]: badc.append((key, local_name(mono, vs)))
        log("   a-weight %d: %d coefficient series, %d non-zero, min valuation %s" % (W, ntot, nz, min(vals) if vals else None))
        if W == 0: check("a-weight 0 (constant part): valuation exactly -10 for all 115 coloured graphs", all(v == -10 for v in vals) and nz == 115, (nz, bad[:3]))
        else: check("a-weight %d: valuation >= %d for all non-zero coefficient series" % (W, VAL_MIN[W]), not bad, bad[:3])
        check("a-weight %d: exponent classes mod 5 in %s (to the series precision)" % (W, sorted(CLASSES[W])), not badc, badc[:3])


def zero_pattern(T, keys, log):
    log("   non-zero local monomials per coloured graph (x, y = a1 of the first/second colour, p, q = a2):")
    pat = {}
    for key in keys:
        vs = local_vars(key[1])
        pat[key] = tuple(local_name(m, vs) for m in local_monos(key[1]) if not T[key][m].iszero())
    groups = {}
    for key, pt in pat.items(): groups.setdefault(pt, []).append(key)
    for pt, ks in sorted(groups.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        byidx = {}
        for (idx, col) in ks: byidx.setdefault(idx, []).append(col)
        log("      %-44s %3d coloured graphs: %s" % ("{%s}" % ", ".join(pt), len(ks),
            "; ".join("graph %d: %s" % (i, "all 5 monochrome" if len(c) == 5 and all(len(set(x)) == 1 for x in c) else
                                        "all 20 mixed" if len(c) == 20 else "all 25" if len(c) == 25 else str(c)) for i, c in sorted(byidx.items()))))
    return pat


def polar_entries(key, cf):
    vs = local_vars(key[1]); out = []
    for mono in local_monos(key[1]):
        W = mono_weight(mono, vs)
        if W not in POLAR: continue
        for e in POLAR[W]: out.append((local_name(mono, vs), mono, e))
    return out


def table_reps(cf, reps, R, log, check):
    """per orbit representative: every polar coefficient as residue/rational [Q(zeta) coordinates]; Galois back-substitution"""
    p = R.p; data = {}; nrec = 0; ngal = 0
    log("   polar coefficients per orbit representative (col[0] = 0); entry = coefficient of (monomial) s^e as R.fmt [Q(zeta) coords c0 + c1 z + c2 z^2 + c3 z^3]")
    for key in reps:
        idx, col = key; cells = []
        for name, mono, e in polar_entries(key, cf):
            x = ser_coef(cf[key][mono], e)
            xs = [ser_coef(cf[(idx, galois_colouring(col, m))][mono], e) for m in range(1, 5)]
            cs = qzeta_coords(R, xs); rs = qzeta_ratrec(R, cs)
            if rs is None: nrec += 1
            if sum(R.mul(c, R.z(k)) for k, c in enumerate(cs)) % p != x: ngal += 1
            data[(key, mono, e)] = (x, cs, rs)
            if x or any(cs): cells.append("(%s) s^%d: %s [%s]" % (name, e, R.fmt(x), qzeta_fmt(R, cs)))
        log("   %-2d %-7s | %s" % (idx, str(col), " | ".join(cells) if cells else "no non-zero polar coefficient"))
    check("Galois back-substitution of the Q(zeta) coordinates reproduces every table entry", ngal == 0, ngal)
    import math
    log("   Q(zeta) coordinate reconstruction failures: %d of %d entries (ratrec bound %d)" % (nrec, len(data), min(math.isqrt(p // 2) // 100, R.RATREC_BOUND)))
    return data, nrec


def table_all(cf, keys, R, log):
    log("   all coloured graphs, polar coefficients (R.fmt), the a-independent ones first: s^-10, s^-5; then (x) s^-6, (x) s^-1, (y) s^-6, (y) s^-1, (x^2) s^-2, (x y) s^-2, (y^2) s^-2")
    for key in keys:
        cells = ["%s" % R.fmt(ser_coef(cf[key][mono], e)) for name, mono, e in polar_entries(key, cf)]
        log("   %-2d %-7s | %s" % (key[0], str(key[1]), " | ".join(cells)))


def leading_w3(cf, reps, R, log):
    log("   leading a-weight-3 coefficients at s^2 (x^3, x^2 y, x y^2, y^3, p, q; these are never polar) per representative:")
    for key in reps:
        vs = local_vars(key[1]); cells = []
        for mono in local_monos(key[1]):
            if mono_weight(mono, vs) == 3:
                S = cf[key][mono]
                cells.append("(%s) %s" % (local_name(mono, vs), R.fmt(ser_coef(S, 2)) if S.c else "0 series"))
        log("   %-2d %-7s | %s" % (key[0], str(key[1]), " | ".join(cells)))


# ------------------------------------------------------------------ comparison with the covariant scan (genus2_scalar.log)
def _parse_fmt(R, s):
    s = s.strip()
    m = re.match(r"^\[(\d+) mod p\]$", s)
    if m: return int(m.group(1)) % R.p
    m = re.match(r"^(-?\d+(?:/\d+)?)\*zeta\^(\d)$", s)
    if m: return R.mul(R.of(Fr(m.group(1))), R.z(int(m.group(2))))
    return R.of(Fr(s))


def parse_prev_log(path, p, R):
    """{(idx, col): [T0_-10, T1_-6, T2_-2, T3_2, T4_2]} from the '(6)' block of the prime-p section of genus2_scalar.log"""
    if not os.path.exists(path): return None
    out = {}; inblock = False; inprime = False
    for line in open(path, encoding="utf-8"):
        if line.startswith("prime p = "): inprime = (int(line.split("=")[1]) == p); inblock = False; continue
        if not inprime: continue
        if "leading system over 23 representatives" in line: inblock = True; continue
        if inblock:
            m = re.match(r"^\s+\((\d+), \(([\d, ]+)\)\)\s+(.*)$", line)
            if not m: inblock = False; continue
            col = tuple(int(x) for x in m.group(2).split(",") if x.strip())
            out[(int(m.group(1)), col)] = [_parse_fmt(R, x) for x in m.group(3).split(", ")]
    return out


def covariant_combination(cf, key, R, phase_pt, W, a2part):
    """Sum over local monomials of a-weight W (a1-only if not a2part, else p, q) of phase(mono) * Cf(coef): the series T_W of the
    covariant scan, phase(mono) = mono evaluated at the covariant point (a1, a2) = (1, 1)"""
    vs = local_vars(key[1]); acc = None
    for mono in local_monos(key[1]):
        if mono_weight(mono, vs) != W: continue
        g = mono_global(mono, vs)
        if (sum(g[5:]) > 0) != a2part: continue
        t = cf[key][mono].scale(gmono_eval(R, g, phase_pt))
        acc = t if acc is None else acc.add(t)
    return acc


def compare_prev(cf, reps, R, prev, log, check):
    """the covariant specialisation reproduces (T0_-10, T1_-6, T2_-2, T3_2, T4_2) of genus2_scalar.log for all 23 representatives"""
    if prev is None: log("   (previous log not found: comparison skipped)"); return
    ph = covariant_point(R, 1, 1)
    spec = [(0, False, -10), (1, False, -6), (2, False, -2), (3, False, 2), (3, True, 2)]
    bad = []; n = 0
    for key in reps:
        if key not in prev: bad.append((key, 'missing')); continue
        for j, (W, a2part, e) in enumerate(spec):
            S = covariant_combination(cf, key, R, ph, W, a2part)
            x = ser_coef(S, e) if S is not None else 0
            n += 1
            if x != prev[key][j]: bad.append((key, j, R.fmt(x), R.fmt(prev[key][j])))
    check("covariant specialisation a_k^(i) = zeta^-i(2k-1) a_k reproduces the %d leading coefficients (T0_-10, T1_-6, T2_-2, T3_2, T4_2) of genus2_scalar.log" % n,
          not bad, bad[:4])
