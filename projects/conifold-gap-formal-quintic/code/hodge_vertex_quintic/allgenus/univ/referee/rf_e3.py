# -*- coding: utf-8 -*-
"""rf_e3.py -- REFEREE (univ): E3 re-derived independently.  Quintic rows from the referee's own route (rf_pf, rf_gauge canonical,
rf_rows; equal to the project rows, rf_rows test), dressings D from rf_e3core (Cayley, own iteration, properties checked on D),
graph sums with the module-disjoint engine (rf_indep).
Cases:  A<seed>   random flat-regular symplectic X (k >= 1), column-1 orbits solved: annihilating, not a gauge
        K<seed>   the same with u-dependent k = 0 diagonal terms
        B<seed>   A<seed> plus (1/7) u^2 zhat H(x)H added to X after solving (violates Delta_1 at u^2)   [control]
        R<seed>   flat-regular symplectic X with RANDOM column-1 orbits (annihilation not imposed)            [control]
        N<seed>, P<seed>   as L with k = 0 terms; P: random X with random column-1 orbits (NOT in the kernel)   [control]
        L<seed>   FIRST ORDER: X in the kernel of the LINEARISED annihilation map, D_t = Cayley(t X), t = 0..4; the polar
                  coefficients are polynomials of degree <= 2g-3 in t (every X is O(u)); prints their t^1 coefficients
Prefix x (e.g. xA1): quintic base, but the column-1 orbits are solved against the MODEL covector (wrong-covector control).
Prefix m (e.g. mA7): PURE MODEL base (P = 1: the model x^4 + v/x at v(u), no quintic input), model vanishing covector.
Run: py -I -B rf_e3.py CASE GENUS [max_outer|full] [nproc]   -> rf_e3.log
"""
import os, sys, json, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_e4 import make_rows, fmtr, SCR, CG
from rf_e3core import random_X, solve_ann, cayley, checks, set_orbit, madd, ident
from rf_rows import dress, red


def yhat(g, K, p, model=False):
    from rf_pf import system
    from rf_gchk import yc_system
    S = system(30); KZ = 3 * g - 3; ys = yc_system(S['bt'] if model else S['b'], KZ + 2, K + KZ + 4)
    out = []
    for n in range(KZ + 2):
        a, y = ys[n]; nu = n - n // 5
        assert all(y.coef(j) == 0 for j in range(nu)) and y.coef(nu) != 0, "valuation of y_%d" % n
        out.append([red(y.coef(nu + j), p) for j in range(K)])
    return out


def build_D(case, g, K, p, yh, log):
    KZ = 3 * g - 3; kind = case[0]; seed = int(case[1:]); rng = random.Random(1000 + seed)
    X = random_X(rng, KZ, K, p, with_k0=(kind in "KN"), skip_col1=(kind not in "RP"))
    if kind == "R": return X, cayley(X, KZ, p, K)
    if kind == "P": return X, None
    lin = kind in "LN"
    X, D, it, skip = solve_ann(X, yh, KZ, p, K, linear=lin)
    log("   column-1 orbits solved in %d iterations (Delta_m not imposed for m = %s)" % (it, skip))
    if kind == "B" and K > 2:                                       # the u^2 term is outside the genus-2 window
        cur = X.get((2, 1, 1), [0] * K); add = [0] * K; add[2] = red(Fr(1, 7), p)
        set_orbit(X, 1, 1, [(x + y) % p for x, y in zip(cur, add)], p)
    return X, (None if lin else cayley(X, KZ, p, K))


def main():
    case, g = sys.argv[1], int(sys.argv[2]); mo = None if (len(sys.argv) <= 3 or sys.argv[3] == "full") else int(sys.argv[3])
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 1; K = 2 * g - 2; KZ = 3 * g - 3
    LOG = open(os.path.join(HERE, "rf_e3.log"), "a", encoding="ascii")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    from rf_indep import run, polar_in_T
    model = case.startswith("m"); wrong = case.startswith("x"); case_ = case
    if model or wrong: case = case[1:]
    t0 = time.time(); base, T1q, p = (model_rows(g, K) if model else make_rows("q", "canon", g, K)); d = json.load(open(base)); yh = yhat(g, K, p, model)
    PR = os.environ.get("RF_PRIME", "89")
    log("=" * 100); log(("rf_e3: case %s (%s), genus %d, max_outer %s, p = 2^" + PR + "-1   [%s]") % (case_, "PURE MODEL base P = 1, model vanishing covector" if model else "quintic base", g, mo, time.strftime("%Y-%m-%dT%H:%M:%S")))
    X, D = build_D(case, g, K, p, yhat(g, K, p, True) if wrong else yh, log)    # x: built against the MODEL covector
    ts = [0, 1, 2, 3, 4] if case[0] in "LNP" else [None]
    polys = {}
    for t in ts:
        if t is None: Dt = D
        else: Dt = cayley({k: [x * t % p for x in v] for k, v in X.items()}, KZ, p, K) if t else ident(K)
        sd, irr, ann, d0, Tu = checks(Dt, yh, KZ, p, K)
        if t in (None, 1):
            log("   D: %d entries; symplectic: %s; flat-regular: %s; annihilating mod u^%d for m = 1..%d: %s%s" % (len(Dt), not sd, not irr, K, KZ, not ann, "" if not ann else " (violated at m = %s)" % sorted(ann)))
            log("   delta_0 (z^0 coefficient / e_1) = %s  -> %s" % ([fmtr(x, p) for x in d0], "not constant: no z-independent target connection" if any(d0[1:]) else "constant"))
        Q = dress(d['Q'], Dt, K, KZ, p); path = os.path.join(SCR, "rows_e3_%s_g%d_t%s_%s.json" % (case_, g, t, PR))
        json.dump(dict(prime=PR, p=p, zeta=d['zeta'], K=K, NZ=KZ, Q=Q), open(path, "w"))
        tag = "e3_%s_g%d_t%s_%s_%s" % (case_, g, t, "full" if mo is None else "E%d" % mo, PR)
        tab, L, complete = run(path, g, mo, nproc, 1600.0, os.path.join(SCR, "parts_%s.json" % tag), log)
        if not complete: log("   INCOMPLETE (rerun to resume)"); return
        shown = [i for i in range(len(L)) if mo is None or L[i]['E'] <= mo]; v0 = [i for i in shown if L[i]['E'] == 0][0]
        if t is not None:
            polys[t] = {i: tab.get(i, [0] * K) for i in shown}; continue
        bad = {i: tab[i] for i in shown if L[i]['E'] > 0 and i in tab and any(tab[i])}
        if not bad: log("   RESULT: all %d outer graphs with edges%s are O(u^%d): REGULAR" % (len(shown) - 1, "" if mo is None else " (<= %d edges)" % mo, K))
        else:
            first = {}
            for i, c in bad.items(): first.setdefault((L[i]['E'], next(n for n in range(K) if c[n])), []).append(i)
            log("   RESULT: %d of %d outer graphs with edges are SINGULAR; (E, first order): %s" % (len(bad), len(shown) - 1, {k: v for k, v in sorted(first.items())}))
            for i in sorted(bad)[:6]: log("      outer graph %d (E=%d, gv=%s): %s" % (i, L[i]['E'], list(L[i]['gv']), [fmtr(x, p) for x in bad[i]]))
        c = polar_in_T(tab.get(v0, [0] * K), Tu, g, p, K)
        okv = c[0] == red(CG[g], p) and not any(c[1:])
        log("   vertex graph in T_eff = -u Delta_-1: %s -> %s" % (" + ".join("(%s) T^%d" % (fmtr(c[j], p), j + 2 - 2 * g) for j in range(K)), "c_g T_eff^(2-2g) + 0" if okv else "NOT c_g T_eff^(2-2g) + 0"))
    if case[0] in "LNP":
        # Lagrange interpolation in t over the points 0..4 (exact mod p); report the t^1 and top coefficients
        pts = ts; nz1 = {}; top = {}
        for i in polys[0]:
            for n in range(K):
                vals = [polys[t][i][n] for t in pts]
                co = interp(pts, vals, p)
                if co[1]: nz1.setdefault(i, []).append(n)
                if co[4]: top.setdefault(i, []).append(n)
        e_nz = {i: v for i, v in nz1.items() if L[i]['E'] > 0}
        log("   FIRST ORDER (t^1) on the edge graphs: %s" % ("ZERO for all %d outer graphs with edges" % (len(polys[0]) - 1) if not e_nz else "NON-ZERO on %s" % e_nz))
        log("   consistency: t^4 coefficients vanish everywhere: %s" % (not top))
    log("   (%.0f s)" % (time.time() - t0))


def model_rows(g, K):
    """the pure model at v(u) in the quintic's hatted frame: P = 1, Ghat = diag(tau^((1-c_b)/4))"""
    from rf_series import tau_series, LS
    from rf_rows import hatted_gauge, rows, prime_zeta
    KZ = 3 * g - 3; tau = tau_series(K + KZ + 8); PR = os.environ.get("RF_PRIME", "89"); p, z = prime_zeta(PR)
    pg = [[(LS.const(1, K + KZ + 8) if k == 0 else LS.zero(K + KZ + 8)) for b in range(5)] for k in range(KZ + 1)]
    Gh = hatted_gauge(pg, tau, K, KZ, p); Q = rows(Gh, tau, K, KZ, p, z)
    path = os.path.join(SCR, "rows_model_g%d_%s.json" % (g, PR))
    json.dump(dict(prime=PR, p=p, zeta=z, K=K, NZ=KZ, Q=Q), open(path, "w"))
    return path, None, p


def interp(xs, ys, p):
    """coefficients of the interpolating polynomial (degree < len(xs)) mod p"""
    n = len(xs); co = [0] * n
    for i in range(n):
        num = [1]; den = 1
        for j in range(n):
            if j == i: continue
            num = [((num[k - 1] if k > 0 else 0) - xs[j] * (num[k] if k < len(num) else 0)) % p for k in range(len(num) + 1)]
            den = den * (xs[i] - xs[j]) % p
        f = ys[i] * pow(den, p - 2, p) % p
        for k in range(n): co[k] = (co[k] + f * num[k]) % p
    return co


if __name__ == "__main__":
    main()
