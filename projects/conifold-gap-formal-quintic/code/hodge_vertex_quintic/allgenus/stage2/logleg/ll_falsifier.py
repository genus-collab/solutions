# -*- coding: utf-8 -*-
"""ll_falsifier.py -- falsifier (b) of the Stage-2 route: TR correlators of x = W = p^4/4 + u/p, y = p paired with ALL L-legs (p^m, m >= 1)
and with the log leg -log p, in every weight that can be singular.   py -I -B ll_falsifier.py <out.log> <field> <g,n> [<g,n> ..]
field: 89 | 61 | both (2^89-1 and 2^61-1, cross-checked by rational reconstruction) | exact (Q(zeta5)).
Homogeneity (TR_MODEL s.5): <p^m_1 .. p^m_a, (-log p)^k>_h = value(u=1) * u^e,  e = 2 - 2h - n + Sum m/5,  n = a + k.
ML_K (via D(iii), D(ii)) predicts: (P1) every pairing with a >= 1 and e < 0 vanishes; (P2) a = 0: <(-log p)^k>_h = F_h^(k)(1) (Gap Theorem,
TR normalisation, ll_core.Fh_deriv); (P3) Z/5 selection: non-zero only if Sum m = 0 mod 5.  Controls (must be SINGULAR non-zero):
legs p^-1, p^-2, p^-5 (single valued, outside L); and functions of x (x, x^2) must pair to zero (must be exactly zero).
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from ll_core import QzField, make_fp, P89, P61, setup, Legs, contract_last, pair_table, compositions, Fh_deriv, at_infinity, rat


class Log(object):
    def __init__(self, path):
        self.f = open(path, 'w'); self.t0 = time.time(); self.fails = 0
    def __call__(self, s=''):
        print(s, flush=True); self.f.write(s + '\n'); self.f.flush()
    def check(self, ok, what):
        self('  [%s] %s' % ('PASS' if ok else 'FAIL', what))
        if not ok: self.fails += 1
    def close(self):
        self('== %s (%.1f s, %d failures)' % ('ALL CHECKS PASSED' if not self.fails else 'SOME CHECKS FAILED', time.time() - self.t0, self.fails))
        self.f.close()


def nsing(h, n, a):
    """number of ordered L-weight tuples (m_i >= 1) with Sum m = 0 mod 5 and e < 0"""
    lim = 5 * (2 * h - 2 + n)
    return sum(sum(1 for _ in compositions(s, a)) for s in range(5, lim, 5))


def run_field(F, cases, log, store):
    kmx = max(6 * g - 5 + 2 * n for g, n in cases); N = 2 * kmx + 6
    cur, tr = setup(F, N); c = cur.c; L = Legs(F, c)
    for (h, n) in cases:
        t0 = time.time(); W = tr.omega(h, n); tw = time.time() - t0
        log('%s  omega_{%d,%d}: %d e-basis coefficients (%.1f s, incl. lower omegas)' % (F.name[:22], h, n, len(W), tw))
        lim = 5 * (2 * h - 2 + n)
        for k in range(n + 1):
            a = n - k; t1 = time.time()
            T = pair_table(F, c, W, n, k, lim, L)
            if a == 0:
                val = T.get((), F.zero); pred = Fh_deriv(h, k) if (h >= 2 or (h == 1 and k >= 1) or k >= 3) else None
                if pred is None: continue
                log.check((val - F.of(pred)).iszero(),
                          '(%d,%d) a=0 k=%d: <(-log p)^%d>_%d = %s = F_%d^(%d)(1) = %s   [frozen: %s]' % (
                              h, n, k, k, h, rat(F, val), h, k, pred, pred * Fr(-5) ** (1 - h)))
                store[(h, n, k, ())] = val
                continue
            sel = sum(1 for m, x in T.items() if sum(m) % 5 and not x.iszero())
            sing = [(m, x) for m, x in T.items() if sum(m) % 5 == 0 and sum(m) < lim and not x.iszero()]
            e0 = [(m, x) for m, x in T.items() if sum(m) == lim and not x.iszero()]
            for m, x in T.items():
                if sum(m) % 5 == 0: store[(h, n, k, m)] = x
            log.check(sel == 0 and not sing,
                      '(%d,%d) a=%d L-legs, k=%d log legs: %d singular weight tuples tested (e<0), non-zero: %d; selection violations: %d; '
                      'e=0 non-zero: %d (%.1f s)' % (h, n, a, k, nsing(h, n, a), len(sing), sel, len(e0), time.time() - t1))
            if sing: log('     NON-ZERO SINGULAR (ML_K violated): %s' % [(m, rat(F, x)) for m, x in sing[:6]])
            if e0 and len(e0) <= 4: log('     e=0 values: %s' % [(m, rat(F, x)) for m, x in sorted(e0)])
        controls(F, c, W, h, n, L, log)


def controls(F, c, W, h, n, L, log):
    # functions of x pair to zero (one slot x or x^2, other slots log legs):
    xl = lambda jk: L.mono(4)(jk) * F.of(Fr(1, 4)) + L.mono(-1)(jk)
    x2 = lambda jk: L.mono(8)(jk) * F.of(Fr(1, 16)) + L.mono(3)(jk) * F.of(Fr(1, 2)) + L.mono(-2)(jk)
    z1 = contract_last(W, [xl] + [L.logleg] * (n - 1)); z2 = contract_last(W, [x2] + [L.logleg] * (n - 1))
    zz = at_infinity(F, c, contract_last(W, [xl]), n - 1, 5 * (2 * h - 2 + n) + 10) if n >= 2 else {}
    log.check(not z1 and not z2 and not zz, '(%d,%d) control: legs x, x^2 (functions of x) pair to exactly zero (with log legs, and x with all '
              'p^m, Sum m <= %d)' % (h, n, 5 * (2 * h - 2 + n) + 10))
    # off-L single-valued legs: p^-5 (+ log legs), p^-1 (+ one p^m, m = 1 mod 5, + log legs), p^-2: must give singular non-zero values
    out = []
    for w in ((-5, -10, -15) if n == 1 else (-5, -1, -2)):
        if n == 1:
            T = contract_last(W, [L.mono(w)]); vals = {(w,): T.get((), F.zero)}
        else:
            T = contract_last(W, [L.mono(w)] + [L.logleg] * (n - 2))
            vals = {(w, m[0]): x for m, x in at_infinity(F, c, T, 1, 5 * (2 * h - 2 + n) + 5).items() if (m[0] + w) % 5 == 0}
        nz = [(m, x) for m, x in vals.items() if not x.iszero() and 2 - 2 * h - n + Fr(sum(m), 5) < 0]
        out.append((w, len(nz), sorted(nz)[:2]))
    ok = sum(cnt for _, cnt, _ in out) > 0
    log.check(ok, '(%d,%d) control (must be singular, at least one leg): off-L legs give non-zero e<0 values: %s' % (
        h, n, '; '.join('p^%d: %d [%s]' % (w, cnt, ', '.join('%s=%s' % (m, rat(F, x)) for m, x in ex)) for w, cnt, ex in out)))


def parse(args):
    return [tuple(int(x) for x in a.split(',')) for a in args]


if __name__ == '__main__':
    out, fld = sys.argv[1], sys.argv[2]; cases = parse(sys.argv[3:])
    log = Log(out); log('ll_falsifier.py %s  (allgenus/stage2/logleg, %s)' % (' '.join(sys.argv[2:]), time.strftime('%Y-%m-%d %H:%M')))
    fields = {'89': [make_fp(P89)], '61': [make_fp(P61)], 'both': [make_fp(P89), make_fp(P61)], 'exact': [QzField()]}[fld]
    stores = []
    for F in fields:
        st = {}; run_field(F, cases, log, st); stores.append((F, st))
    if len(stores) == 2:
        (F1, s1), (F2, s2) = stores; bad = 0; nrec = 0
        for key in set(s1) | set(s2):
            x1 = s1.get(key, F1.zero); x2 = s2.get(key, F2.zero)
            if x1.iszero() != x2.iszero(): bad += 1; continue
            if x1.iszero(): continue
            r = F1.ratrec(x1)
            if r is None: continue
            nrec += 1
            if not (F2.of(r) - x2).iszero(): bad += 1
        log.check(bad == 0, 'two primes: identical zero patterns on %d stored pairings (e <= 0); %d non-zero values rationally reconstructed '
                  'at 2^89-1 and confirmed at 2^61-1' % (len(set(s1) | set(s2)), nrec))
    log.close()
