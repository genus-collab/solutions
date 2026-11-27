# -*- coding: utf-8 -*-
"""pb_mech.py NU -- Stage-5 'pb' audit: a DIRECT test of the mechanism of Theorem P(b), steps 1-4 and 6, mod 2^61-1.

Step 1 (CFU): Y = X y(X^5;u), y = 1 + sum u^k Y_k(T), V(u), with V/y + T y^4/4 = What(T;u) := sum_m (u^(m-1)/m)(1 + T^m/(5m-1)),
  i.e. What(X) = Xhat-phase times X; solved order by order, (T-1) Y_k = R_k(T) - V_k; checked V = tau^(-5/4).
Step 2: beta(Y) dY = B(X(Y)) X'(Y) dY, B the hatted amplitude, log B = -(1/2) log(1-uX^5)
  + sum_k b_k (uXw)^(2k-1) [5^(1-2k) - 5 sum_P C(2k-2+5P,5P) u^P X^(5P)]; reduction Y^n == (n-3) w Y^(n-4) + V Y^(n-5), n >= 4,
  down to pi_m, m = -1..3.
Step 4 (SP side, with the poly referee's bookkeeping): eps_b pi_(c_b) = C0 sum_(k = -b mod 5) w^k u^mu(b,k) p_(k,b)(u), the
  canonical gauge of pb_gauge (built independently from the periods).  Grading: pi_(c_b) contains w^k with k = -b mod 5 only.
Step 6 (residue side): sum_b pi_(c_b) s^(-c_b) I^C_(c_b)(v(u), z) = V (pb_res.residue_V), I^C_m the model residues of GI Lemma V.
Controls: the stored normalisation of the gauge must fail step 4 (first at w^5 u^4); a reduction with V -> 1 must fail.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C, pb_gauge as G, pb_sp as S, pb_res as R
from pb_core import LS, red, finv
from fractions import Fraction as Fr
from math import comb, factorial

LOG = []
def log(s=''):
    print(s); LOG.append(s)


def md_mul(a, b, NU):
    r = {}
    for (i1, j1, k1), c1 in a.items():
        for (i2, j2, k2), c2 in b.items():
            if i1 + i2 <= NU:
                key = (i1 + i2, j1 + j2, k1 + k2); r[key] = (r.get(key, 0) + c1 * c2) % C.MOD
    return {k: v for k, v in r.items() if v}


def md_add(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) + s * v) % C.MOD
    return {k: v for k, v in r.items() if v}


def md_exp(g, NU):
    out = {(0, 0, 0): 1}; term = {(0, 0, 0): 1}
    for n in range(1, NU + 1):
        term = md_mul(term, g, NU); inv = finv(n)
        term = {k: v * inv % C.MOD for k, v in term.items()}
        out = md_add(out, term)
    return out


def cfu(NU):
    """y[i] (T-polynomials), V[i]: the CFU map to O(u^(NU+1))."""
    y = {(0, 0, 0): 1}; V = [1] + [0] * NU
    for k in range(1, NU + 1):
        # [u^k] (V/y + T y^4/4) with Y_k = V_k = 0; use the variable slot 2 for T
        yi = md_inv(y, k)
        t1 = md_mul({(i, 0, 0): V[i] for i in range(k + 1) if V[i]}, yi, k)
        y4 = md_mul(md_mul(y, y, k), md_mul(y, y, k), k)
        Nk = {}
        for (i, j, t), c in t1.items():
            if i == k: Nk[t] = (Nk.get(t, 0) + c) % C.MOD
        for (i, j, t), c in y4.items():
            if i == k: Nk[t + 1] = (Nk.get(t + 1, 0) + c * finv(4)) % C.MOD
        Rk = {0: finv(k + 1), k + 1: finv((k + 1) * (5 * k + 4))}
        for t, c in Nk.items():
            Rk[t] = (Rk.get(t, 0) - c) % C.MOD
        Vk = sum(Rk.values()) % C.MOD
        V[k] = Vk
        Rk[0] = (Rk.get(0, 0) - Vk) % C.MOD
        # divide sum_t Rk[t] T^t by (T - 1)
        top = max(Rk); q = [0] * top; carry = 0
        for t in range(top, 0, -1):
            carry = (carry + Rk.get(t, 0)) % C.MOD; q[t - 1] = carry
        assert (carry + Rk.get(0, 0)) % C.MOD == 0
        for t, c in enumerate(q):
            if c: y[(k, 0, t)] = c
    return y, V


def md_inv(a, NU):
    """1/a for a = 1 + O(u) (keys (i, j, t))."""
    d = md_add(a, {(0, 0, 0): 1}, -1)
    out = {(0, 0, 0): 1}; term = {(0, 0, 0): 1}
    for n in range(1, NU + 1):
        term = md_mul(term, d, NU); term = {k: (-v) % C.MOD for k, v in term.items()}
        out = md_add(out, term)
    return out


def main():
    NU = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    p = S.P61; C.set_field(p); t0 = time.time()
    log('pb_mech.py %d   (mod 2^61-1; u-orders <= %d)' % (NU, NU))
    y, V = cfu(NU)
    Vt = C.model_v(NU + 2)[0]
    okV = all((V[i] - Vt.coef(i)) % p == 0 for i in range(NU + 1))
    log('  step 1 (CFU): Y_k polynomial in T solved for k <= %d; V = tau^(-5/4) through u^%d: %s; Y_1 = (1+T)/18: %s'
        % (NU, NU, 'PASS' if okV else 'FAIL', 'PASS' if y.get((1, 0, 0)) == finv(18) == y.get((1, 0, 1)) else 'FAIL'))
    # inverse: x(S) with X = Y x(Y^5): iterate x <- 1/y(S x^5)
    x = {(0, 0, 0): 1}
    for it in range(NU + 1):
        x5 = md_mul(md_mul(x, x, NU), md_mul(md_mul(x, x, NU), x, NU), NU)
        Sx5 = {(i, j, t + 1): c for (i, j, t), c in x5.items()}
        comp = {}
        maxdeg = max(t for (_, _, t) in y)
        pw = [{(0, 0, 0): 1}]
        for d in range(1, maxdeg + 1):
            pw.append(md_mul(pw[-1], Sx5, NU))
        for (i, j, t), c in y.items():
            comp = md_add(comp, {(i + a, b, s): c * v % p for (a, b, s), v in pw[t].items() if i + a <= NU})
        x = md_inv(comp, NU)
    # trivariate (u, w, Y): X = Y x(Y^5), X' = sum (1+5s) x_(i,s) Y^(5s)
    X = {(i, 0, 1 + 5 * s): c for (i, j, s), c in x.items()}
    Xp = {(i, 0, 5 * s): c * (1 + 5 * s) % p for (i, j, s), c in x.items()}
    X5 = md_mul(md_mul(md_mul(X, X, NU), md_mul(X, X, NU), NU), X, NU)
    U5 = {(i + 1, j, k): c for (i, j, k), c in X5.items() if i + 1 <= NU}
    bern = C.bernoulli(2 * NU + 4)
    logB = {}; pw = {(0, 0, 0): 1}
    for m in range(1, NU + 1):
        pw = md_mul(pw, U5, NU); logB = md_add(logB, {k: v * finv(2 * m) % p for k, v in pw.items()})
    U5p = [{(0, 0, 0): 1}]
    for P in range(1, NU + 1):
        U5p.append(md_mul(U5p[-1], U5, NU))
    for k in range(1, (NU + 1) // 2 + 1):
        r = 2 * k - 1
        Xr = {(0, 0, 0): 1}
        for _ in range(r):
            Xr = md_mul(Xr, X, NU)
        bk = red(bern[2 * k] / (2 * k * (2 * k - 1)))
        brack = {(0, 0, 0): red(Fr(1, 5 ** r))}
        for P in range(0, NU + 1):
            brack = md_add(brack, {kk: v * (-5 * comb(2 * k - 2 + 5 * P, 5 * P)) % p for kk, v in U5p[P].items()})
        term = md_mul(Xr, brack, NU - r)
        logB = md_add(logB, {(i + r, j + r, kk): v * bk % p for (i, j, kk), v in term.items() if i + r <= NU})
    B = md_exp(logB, NU)
    beta = md_mul(B, Xp, NU)
    log('  step 2: beta(Y) built (%d terms), %.1f s' % (len(beta), time.time() - t0))
    results = {}
    for tag, Vred in (('true V', V), ('control V -> 1', [1] + [0] * NU)):
        pi = reduce_forms(beta, Vred, NU)
        results[tag] = pi
    pi = results['true V']
    grad = all((j - G.CW[b]) % 5 == 0 for b in range(5) for (i, j) in pi[G.CW[b]])
    log('  grading: pi_(c_b) contains only w^k with k = -b mod 5: %s' % ('PASS' if grad else 'FAIL'))
    # gauge side
    K = C.Conifold(NU + 40)
    b_, bt_ = G.systems(K, NU + 34)
    jmax = max(j for m in pi for (i, j) in pi[m]) + 1
    ok4 = {}
    for norm in ('can', 'stored'):
        Pg = G.gauge(b_, bt_, jmax, norm)
        C0 = K.varpi0
        bad = None
        for bb in range(5):
            m = G.CW[bb]
            gs = {}
            for k in range(jmax + 1):
                if (-k) % 5 != bb: continue
                term = C0.mul(Pg[k][bb]).shift(G.mu0(k))
                for i in range(NU + 1):
                    c = term.coef(i)
                    if c: gs[(i, k)] = c
            lhs = {kk: v * G.EPS[bb] % p for kk, v in pi[m].items()}
            keys = sorted(set(lhs) | set(gs))
            for kk in keys:
                if (lhs.get(kk, 0) - gs.get(kk, 0)) % p:
                    bad = min(bad, (kk[1], kk[0], bb)) if bad else (kk[1], kk[0], bb)
        ok4[norm] = bad
    log('  step 4: eps_b pi_(c_b) = C0 x hatted first row of the CANONICAL gauge, all b, u^0..u^%d, w^0..w^%d: %s'
        % (NU, jmax, 'PASS' if ok4['can'] is None else 'FAIL first (w,u,b) = %s' % (ok4['can'],)))
    log('  control [stored normalisation]: first mismatch (w, u, b) = %s (predicted (5, 4, 0)): %s'
        % (ok4['stored'], 'PASS (fails as it must)' if ok4['stored'] == (5, 4, 0) else 'CONTROL DID NOT BEHAVE'))
    # step 6: residue side
    vser = C.model_v(NU + 2)[1]
    Vres = {}
    for bb in range(5):
        m = G.CW[bb]
        IC = icm(m, vser, NU)
        for (i, j), c in pi[m].items():
            sh = (4 * j + m) // 5
            for (e, f), d in IC.items():
                key = (i - sh + e, j + f)
                if 1 <= key[0] <= NU:
                    Vres[key] = (Vres.get(key, 0) + c * d) % p
    C.set_field(0); Vex = R.residue_V(NU, C.bernoulli(2 * NU + 4)); C.set_field(p)
    Vex = {k: red(v) for k, v in Vex.items()}
    # the key (u^E, z^b) receives contributions from pi at u-order E + (4b-1)/5 only (PB_AUDIT s.3): compare where that is <= NU
    inrange = lambda k: 1 <= k[0] <= NU and (k[1] - 4) % 5 == 0 and k[0] + (4 * k[1] - 1) // 5 <= NU
    stray = [k for k in set(Vres) | set(Vex) if 1 <= k[0] <= NU and (k[1] - 4) % 5]
    keys = sorted(set(k for k in set(Vres) | set(Vex) if inrange(k)))
    allk = sorted(set((E, bz) for E in range(1, NU + 1) for bz in range(-(5 * NU) // 4 - 6, 5 * NU) if inrange((E, bz))))
    bad6 = [k for k in allk if (Vres.get(k, 0) - Vex.get(k, 0)) % p]
    npos = len([k for k in allk if k[1] > 0])
    log('  step 6: sum_b pi_(c_b) s^(-c_b) I^C_(c_b)(v(u),z) = V on all %d admissible keys (u^1..u^%d; %d of them with z^b, b > 0,'
        ' where both sides must vanish; %d non-zero keys): %s' % (len(allk), NU, npos, len(keys),
        'PASS' if not bad6 else 'FAIL first %s' % (min(bad6),)))
    log('    z-powers outside 4 mod 5 present: %s' % ('none (PASS)' if not stray else 'FAIL %s' % stray[:3]))
    bad6 = bad6 + stray
    pic = results['control V -> 1']
    diff = any((pic[m].get(k, 0) - pi[m].get(k, 0)) % p for m in pi for k in set(pi[m]) | set(pic[m]))
    log('  control [reduction with V -> 1]: coefficients pi_m change: %s' % ('PASS (fails as it must)' if diff else 'CONTROL DID NOT BEHAVE'))
    ok = okV and grad and ok4['can'] is None and ok4['stored'] == (5, 4, 0) and not bad6 and diff
    log('total %.1f s' % (time.time() - t0))
    log('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    with open(os.path.join(HERE, 'logs', 'pb_mech_%d.log' % NU), 'w') as f:
        f.write('\n'.join(LOG) + '\n')
    print('log written to logs/pb_mech_%d.log' % NU)


def reduce_forms(beta, Vs, NU):
    p = C.MOD
    poly = {}
    for (i, j, k), c in beta.items():
        poly.setdefault(k, {})[(i, j)] = c
    for n in range(max(poly), 3, -1):
        if n not in poly: continue
        for (i, j), c in poly.pop(n).items():
            tgt = poly.setdefault(n - 4, {}); tgt[(i, j + 1)] = (tgt.get((i, j + 1), 0) + (n - 3) * c) % p
            tgt = poly.setdefault(n - 5, {})
            for l, vl in enumerate(Vs):
                if vl and i + l <= NU:
                    tgt[(i + l, j)] = (tgt.get((i + l, j), 0) + vl * c) % p
    assert all(m >= -1 for m in poly), 'reduction left the range -1..3'
    return {m: {k: v for k, v in poly.get(m, {}).items() if v} for m in range(-1, 4)}


def icm(m, vser, NU):
    """I^C_m(v(u), z) = sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!  as {(u-power, z-power): c} mod p."""
    out = {}; i = 0
    while True:
        e = m + 4 * i + 1
        if e > NU: break
        if e >= 0:
            vp = LS.const(1, NU + 1)
            for _ in range(e): vp = vp.mul(vser)
            coef = red(Fr((-1) ** (i + e), 4 ** i * factorial(i) * factorial(e)))
            for n in range(NU + 1):
                c = vp.coef(n)
                if c: out[(n, -(i + e))] = (out.get((n, -(i + e)), 0) + coef * c) % C.MOD
        i += 1
    return out


if __name__ == '__main__':
    main()
