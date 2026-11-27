"""g5_t2.py -- Stage 5 (general), test T2 (+ its T4 controls).

O(5)-twisted P^4 with DISTINCT rational weights lambda_0..lambda_4 (generic; sets with Sum lambda_j != 0 and = 0).
Summand Gamma(1+5s) prod_j Gamma(lambda_j/z - s) (-q)^s, all Gamma factors Stirling-expanded.  For every colour i
(critical point sigma_c -> lambda_i as q -> 0) and n <= N:  (R) c_n(lambda_i + eps) has no negative power of eps;
(V) [eps^0] c_n = [z^n] h_(w_i), w_i = {lambda_i - lambda_j}_(j != i) U {-5 lambda_i}.
Also: (Z) the z-power bookkeeping of the chart argument (Euler substitution x -> x/z in every factor): the exponent of
z is -Sum_f (c_f sigma + d_f)/z - Sum_f e_f; its sigma-coefficient -Sum_f c_f/z vanishes (5 - 5 = 0), leaving the
sigma-independent factor z^(-1 - Sum_j lambda_j/z); (I) the summand reproduces the twisted I-function (residues).
Controls: one weight perturbed in h only (lambda_k -> lambda_k + 1/11, k != i); twist sign flipped in h (+5 lambda_i).
Usage: py -I -B g5_t2.py N [NC]      Log: logs/g5_t2_<N>.log
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D, g5_theory as T

N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
NC = int(sys.argv[2]) if len(sys.argv) > 2 else min(N, 8)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_t2_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []
WSETS = [('W1', [Fr(1), Fr(3), Fr(-2), Fr(5, 2), Fr(-7, 3)]),
         ('W2', [Fr(2), Fr(-1), Fr(3), Fr(-5), Fr(1)]),
         ('W3', [Fr(1), Fr(1, 2), Fr(1, 3), Fr(1, 4), Fr(1, 5)]),
         ('W4', [Fr(7), Fr(-3), Fr(4, 3), Fr(-11, 2), Fr(2, 5)])]


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def main():
    out('g5_t2: O(5)-twisted P^4, distinct rational weights, every colour, z^0..z^%d, exact over Q' % N)
    for name, lam in WSETS:
        fac = T.base(lam) + [T.convex(5)]
        csum = sum(f['c'] for f in fac); dsum = sum(f['d'] for f in fac); esum = sum(f['e'] for f in fac)
        check('(Z) %s: Sum lambda_j = %s; z-exponent of the Euler substitutions = (%s) sigma/z + (%s)/z + (%s): '
              'sigma-independent' % (name, dsum, -csum, -dsum, -esum), csum == 0)
        zt = Fr(1, 7)
        rc = T.residue_check(fac, 0, zt, 3, lambda d: T.poch(1 + 5 * lam[0] / zt, 5 * d) * zt ** (5 * d))
        check('(I) %s: summand residues at s = lambda_0/z + d equal the I-function coefficients, d <= 3, z = 1/7' % name,
              all(r == 1 for _, r in rc), [str(r) for _, r in rc])
        for i in range(5):
            t0 = time.time()
            first, bad, vals, pred = T.run_case(fac, i, N)
            w = T.weights(fac, i)
            check('%s colour %d (lambda_i = %s, w_i = {%s}): (R) no pole and (V) value = h_(w_i), n <= %d'
                  % (name, i, lam[i], ', '.join(str(x) for x in w), N), first is None and not bad,
                  '' if (first is None and not bad) else 'first pole z^%s, value mismatches %s' % (first, bad[:4]))
            out('      c_1 = %s, c_2 = %s, c_5 = %s   (%.1f s)' % (vals[1], vals[2], str(vals[5])[:40] if N >= 5 else '-',
                                                             time.time() - t0))
    out('Controls (T4), W1, every colour, z^0..z^%d:' % NC)
    lam = WSETS[0][1]; fac = T.base(lam) + [T.convex(5)]
    for i in range(5):
        k = (i + 1) % 5
        w = T.weights(fac, i)
        wp = [x + (Fr(-1, 11) if j == (k if k < i else k - 1) else 0) for j, x in enumerate(w)]
        first, bad, _, _ = T.run_case(fac, i, NC, pred_weights=wp)
        check('control colour %d: weight lambda_%d -> lambda_%d + 1/11 in h only: value fails (first at z^%s), no pole'
              % (i, k, k, bad[0] if bad else None), first is None and bool(bad))
        wq = w[:-1] + [-w[-1]]
        first, bad, _, _ = T.run_case(fac, i, NC, pred_weights=wq)
        check('control colour %d: twist weight +5 lambda_i in h: value fails (first at z^%s)' % (i, bad[0] if bad else None),
              bool(bad))
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main(); out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
