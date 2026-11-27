# -*- coding: utf-8 -*-
"""ref7_cands.py -- task (3): examine EVERY root alpha of g = gcd(U_23, U_24) (ref6_plane.py) over its root field K = F_p[alpha]/(f):
  1. b: gcd over K[b] of R_k(alpha, b) = Res_tau(D_1, D_k), k = 2, 3, 4 (interpolated over K); the factor (b - b_L)^2 of the line of Pi removed;
  2. tau: gcd over K[tau] of ALL 35 minors of the seven-row system at (a, b);
  3. (c4_1, c4_2, w) from the seven-row linear system; then the full high system (all conditions) and all residuals over K.
Run: py -I ref7_cands.py --prime 61|89 [--seed 1]"""
import sys, os, random, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *
from ref_kpoly import *
from ref6_plane import setup


def eval2K(K, c, apow, bpow):
    """evaluate a dict in (a, b) (3-tuples, t exponent ignored) at K-elements via power tables"""
    s = K.zero
    for (i, j, _), v in c.items(): s = K.add(s, K.scal(v, K.mul(apow[i], bpow[j])))
    return s


def powers(K, x, n):
    out = [K.one]
    for _ in range(n): out.append(K.mul(out[-1], x))
    return out


def examine(st, Y, f, log):
    p = st['p']; K = KRing(f, p); al = K.gen(); dt = st['dt']; mu = st['mu']; P = st['P']; D = st['D']; t0 = time.time(); out = dict(deg=K.d)
    apow = powers(K, al, 10)
    # D_k(alpha, b, tau) over K: table [k][tau-degree][b-degree]
    tab = []
    for fk in st['DD']:
        cs = mcoeffs(fk, 2); cs = cs + [{}] * (dt + 1 - len(cs)); rowk = []
        for c in cs:
            db = max((e[1] for e in c), default=-1); u = [K.zero] * (db + 1)
            for (i, j, _), v in c.items(): u[j] = K.add(u[j], K.scal(v, apow[i]))
            rowk.append(u)
        tab.append(rowk)
    nb = 42 + 1 + 3; vals = [[], [], []]
    for b in range(1, nb + 1):
        ft = [[keval_int(K, u, b) for u in rowk] for rowk in tab]
        for k in range(3): vals[k].append(kres(K, ft[0], ft[k + 1], dt, dt))
    Rk = [knewton(K, 1, v, p) for v in vals]; assert all(len(r) - 1 <= 42 for r in Rk), [len(r) - 1 for r in Rk]
    G = kgcdl(K, Rk); l0 = K.add(K.const(st['Lc'][0]), K.scal(st['Lc'][1], al)); bL = K.scal((-pow(st['Lc'][2], p - 2, p)) % p, l0)
    mL, rest = kmult_root(K, G, bL)
    log("   factor of degree %d: b-degrees of R_2..R_4 %s; gcd in b: degree %d = (b - b_L)^%d x rest of degree %d   (%.0f s)" % (K.d, [len(r) - 1 for r in Rk], len(G) - 1, mL, len(rest) - 1, time.time() - t0))
    out.update(mL=mL, rest=len(rest) - 1); cands = []
    if len(rest) == 2: cands.append(("off L", K.sub(K.zero, K.mul(rest[0], K.inv(rest[1])))))
    elif len(rest) > 2: log("      rest of degree %d in b: NOT examined (would need a tower of extensions)" % (len(rest) - 1)); out['unexamined'] = True
    if mL > 2: cands.append(("on L (higher multiplicity %d)" % mL, bL))
    out['points'] = []
    for where, b in cands:
        a = K.add(al, K.scal(mu, b)); apw = powers(K, a, 10); bpw = powers(K, b, 10); T = []; a4v = eval2K(K, D['a4'], apw, bpw)
        nz = 0
        for m in st['M7'].values():
            if not m: continue
            ft = ktrim(K, [eval2K(K, c, apw, bpw) for c in mcoeffs(m, 2)])
            if ft: T = kgcd(K, T, ft) if T else kmonic(K, ft); nz += 1
        log("      candidate %s: a4 %s 0; tau-gcd of all 35 minors: degree %s (%d non-zero minors at this (a, b))   (%.0f s)" % (where, "=" if K.iszero(a4v) else "!=", len(T) - 1 if T else "ALL ZERO", nz, time.time() - t0))
        rec = dict(where=where, taudeg=(len(T) - 1 if T else None), a4zero=K.iszero(a4v), sols=[])
        taus = []
        if T and len(T) == 2: taus.append(K.sub(K.zero, T[0]))
        elif T and len(T) == 3:
            half = K.scal(pow(2, p - 2, p), T[1]); disc = K.sub(K.mul(half, half), T[0])
            if K.iszero(disc): taus.append(K.sub(K.zero, half)); log("         the tau-gcd is a perfect square (double root)")
            else: log("         tau-gcd of degree 2 with two distinct roots: NOT examined (tower)"); rec['unexamined'] = True
        elif T and len(T) > 3: log("         tau-gcd of degree %d: NOT examined" % (len(T) - 1)); rec['unexamined'] = True
        for tau in taus:
            tpw = powers(K, tau, 6)
            def ev3(c):
                s = K.zero
                for (i, j, k), v in c.items(): s = K.add(s, K.scal(v, K.mul(K.mul(apw[i], bpw[j]), tpw[k])))
                return s
            rows = [[ev3(e) for e in r] for r in D['D7']]
            sol, rk, free = solve_field(K, [r[:3] for r in rows], [K.sub(K.zero, r[3]) for r in rows])
            if sol is None: log("         seven-row system INCONSISTENT at the tau-root (rank %d)" % rk); rec['sols'].append('7-row inconsistent'); continue
            w = sol[2]; low = [K.one, a, b, tau, w]
            v, hrk, nfree, bad = full_lift(Y, P, low, K, p)
            # which of the remaining conditions fail?  solve without each of q_6, q_7, q_8, J3_0, J3_1 in turn
            hr, hh, names = high_system(P, low, K, p); drop = {}
            for nm in ('q_6', 'q_7', 'q_8', 'J3_0', 'J3_1'):
                keep = [i for i, n_ in enumerate(names) if n_ != nm]; s2 = solve_field(K, [hr[i] for i in keep], [hh[i] for i in keep])[0]; drop[nm] = s2 is not None
            keep = [i for i, n_ in enumerate(names) if n_.startswith('c_') or n_ == 'q_3']; s3 = solve_field(K, [hr[i] for i in keep], [hh[i] for i in keep])[0]
            log("         seven-row system: rank %d, free %d -> w, c4 found; cubics + q_3 consistent: %s; FULL high system: rank %d, %s; consistent after dropping one condition: %s   (%.0f s)" % (rk, len(free), s3 is not None, hrk, "ALL CONDITIONS HOLD: A SOLUTION" if not bad else "fails: %s" % bad, drop, time.time() - t0))
            rec['sols'].append('SOLUTION' if not bad else 'not a solution')
        out['points'].append(rec)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--seed", type=int, default=1); a = ap.parse_args()
    st = setup(a.prime, a.seed); p = st['p']; Y = Sys(st['S']); d = json.load(open(os.path.join(HERE, "ref6_plane_%s_seed%d.json" % (a.prime, a.seed))))
    assert d['mu'] == st['mu'] and d['lam'] == st['lam'], "setup mismatch"
    print("=" * 100); print("ref7_cands, p = 2^%s-1, seed %d: g of degree %d, factors of degrees %s" % (a.prime, a.seed, len(d['g']) - 1, fmt_fac(d['factors'])), flush=True)
    res = []
    for f in d['factors']: res.append(examine(st, Y, f, lambda s: print(s, flush=True)))
    json.dump(res, open(os.path.join(HERE, "ref7_cands_%s_seed%d.json" % (a.prime, a.seed)), "w"))
    nsol = sum(r['deg'] for r in res for pt in r['points'] if 'SOLUTION' in pt['sols'])
    print("REF7 summary: candidates (alpha-roots, with field degrees) %d; solution points of the full system found: %d; unexamined: %s" % (sum(r['deg'] for r in res), nsol, [r['deg'] for r in res if r.get('unexamined') or any(pt.get('unexamined') for pt in r['points'])]))


if __name__ == "__main__":
    main()
