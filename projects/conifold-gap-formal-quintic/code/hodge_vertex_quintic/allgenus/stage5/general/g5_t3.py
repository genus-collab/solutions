"""g5_t3.py -- Stage 5 (general), test T3: other theories (+ T4-type controls).

Each theory is a list of Gamma factors (g5_theory), derived from its twisted I-function (GENERAL.md s.2); the summand
is checked against the I-function by residues (I).  For every colour i (sigma_c -> lambda_i as q -> 0):
(R) c_n(lambda_i + eps) has no negative power of eps, n <= N;  (V) [eps^0] c_n = [z^n] h_(w_i) with the PREDICTED
weights w_i = {-(c_f lambda_i + d_f)}_(f != i): tangent weights lambda_i - lambda_j, a convex O(l) (Euler class,
c_1 = l H + mu) contributes -(l lambda_i + mu) (minus the fibre weight), a concave O(-l) (inverse Euler class,
c_1 = -l H + mu) contributes mu - l lambda_i (plus the fibre weight = a tangent weight of the total space).
Control per theory: the twist weight with the opposite sign in h (or, untwisted, one tangent weight perturbed).
(Z) prints the z-exponent sigma-coefficient of the Euler substitutions, -Sum_f c_f, which is cancelled exactly by
Q^s = (z^(Sum c_f) q)^s; with the naive (-q)^s it would remain (non-CY cases).
Usage: py -I -B g5_t3.py N      Log: logs/g5_t3_<N>.log
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D, g5_theory as T

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
NC = min(N, 6)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_t3_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def conv_mod(l, mu):
    return lambda lam, z, d: T.poch(1 + (l * lam + mu) / z, l * d) * z ** (l * d)


def conc_mod(l, mu):
    def f(lam, z, d):
        r = Fr(1)
        for k in range(l * d):
            r *= (-l * lam + mu - k * z)
        return r
    return f


THEORIES = [
    ('P^1, weights (1, -2)', T.base([1, -2]), []),
    ('P^1, weights (3, 1/2)', T.base([3, Fr(1, 2)]), []),
    ('P^2, weights (1, 3, -2)', T.base([1, 3, -2]), []),
    ('P^2, weights (2, -1, 1/3)', T.base([2, -1, Fr(1, 3)]), []),
    ('local P^2 = K_P2, weights (1, 3, -2), lift c_1 = -3H', T.base([1, 3, -2]) + [T.concave(3)], [conc_mod(3, 0)]),
    ('local P^2, weights (2, -1, 1/3), lift c_1 = -3H + 1/2', T.base([2, -1, Fr(1, 3)]) + [T.concave(3, Fr(1, 2))],
     [conc_mod(3, Fr(1, 2))]),
    ('P^1 with O(-1)+O(-1), weights (1, -2), lifts -H, -H', T.base([1, -2]) + [T.concave(1), T.concave(1)],
     [conc_mod(1, 0), conc_mod(1, 0)]),
    ('P^1 with O(-1)+O(-1), weights (3, 1/2), lifts -H + 1/3, -H - 2/5', T.base([3, Fr(1, 2)]) +
     [T.concave(1, Fr(1, 3)), T.concave(1, Fr(-2, 5))], [conc_mod(1, Fr(1, 3)), conc_mod(1, Fr(-2, 5))]),
]
LAM5 = [Fr(1), Fr(3), Fr(-2), Fr(5, 2), Fr(-7, 3)]
for k in (3, 4, 6):
    THEORIES.append(('P^4 with O(%d), weights W1' % k, T.base(LAM5) + [T.convex(k)], [conv_mod(k, 0)]))
THEORIES.append(('P^4 with O(4), weights W1, lift c_1 = 4H + 2/3', T.base(LAM5) + [T.convex(4, Fr(2, 3))],
                 [conv_mod(4, Fr(2, 3))]))


def main():
    out('g5_t3: other theories, every colour, z^0..z^%d, exact over Q' % N)
    for name, fac, mods in THEORIES:
        nb = sum(1 for f in fac if f['kind'] == 'base'); csum = sum(f['c'] for f in fac)
        out('-- %s: factors %s' % (name, ', '.join('Gamma(%s%s s %s %s/z)' % ('1 + ' if f['e'] else '', f['c'],
                                                                         '+', f['d']) for f in fac)))
        out('   (Z) sigma-coefficient of the z-exponent of the Euler substitutions: %s; Q = %s z^(%s) q cancels it'
            % (-csum, '+' if (nb + sum(int(f['c']) for f in fac if f['kind'].startswith('concave'))) % 2 == 0
               else '-', csum))
        zt = Fr(1, 7)
        for i in range(nb):
            lam = fac[i]['d']

            def imod(d, lam=lam):
                r = Fr(1)
                for m in mods:
                    r *= m(lam, zt, d)
                return r
            rc = T.residue_check(fac, i, zt, 3, imod)
            okI = all(r == 1 for _, r in rc)
            t0 = time.time()
            first, bad, vals, pred = T.run_case(fac, i, N)
            w = T.weights(fac, i)
            check('%s, colour %d: (I) summand = I-function (d <= 3); (R) no pole; (V) value = h_(w_i), w_i = {%s}, n <= %d'
                  % (name, i, ', '.join(str(x) for x in w), N), okI and first is None and not bad,
                  '' if (okI and first is None and not bad) else 'I %s, first pole z^%s, mismatches %s'
                  % ([str(r) for _, r in rc], first, bad[:4]))
            if mods:
                wq = w[:nb - 1] + [-x for x in w[nb - 1:]]
                ctl = 'all twist weights with the opposite sign in h'
            else:
                wq = [w[0] + Fr(1, 11)] + w[1:]
                ctl = 'one tangent weight + 1/11 in h'
            f2, bad2, _, _ = T.run_case(fac, i, NC, pred_weights=wq)
            check('   control (%s): value fails, first at z^%s' % (ctl, bad2[0] if bad2 else None), bool(bad2))
            out('   (%.1f s) c_1 = %s, c_3 = %s' % (time.time() - t0, vals[1], vals[3] if N >= 3 else '-'))
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main(); out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
