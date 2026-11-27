"""g5_t6.py -- Stage 5 (general): the hypotheses of the chart argument, one at a time (exact over Q).

The chart argument predicts, for a singular factor Gamma(e* + c*(sigma - s*)/z) and regular Gamma factors
Gamma(e_f + (w_f + c_f (sigma - s*))/z):  removable singularity at s*, and value = the normalised Stirling factors of
the regular factors at the FIRST POLE sigma = s* - e* z/c* of the singular factor, i.e. with shifts
h_f = e_f - c_f e*/c*  (the Euler-integral insertion y_0^(e*) when the restriction G = 0 is solved for the singular
variable).  Tests:
 A  q = oo chart of the formal quintic: singular Gamma(1+5s) (c* = 5, e* = 1) at s* = 0; prediction Ginf (shift 1/5);
    the UNshifted regular tail must fail.
 B  the same for the O(5)-twisted P^4 with general weights: prediction exp(Sum_n beta_n(1/5) z^n Sum_j lambda_j^-n).
 C  local P^2 at q = oo (orbifold point, singular Gamma(3s - mu/z), c* = 3, e* = 0): prediction without shift.
 D  shifted singular factor Gamma(1/3 + (lambda_i - sigma)/z) in the twisted P^4: shifts e_f + c_f/3.
 E  orbifold-type singular factor Gamma((lambda_1 - 2 sigma)/z) (weighted P(1,2)), c* = -2; both colours.
 F  regular part NOT of Gamma type: random WKB data (incl. even z-powers); prediction exp(Sum_n z^n A_n(s*)).
 G  coincident weights (two or three singular factors, M = 2, 3): the chart point is degenerate; prediction: poles;
    for M = 2 the eps^-n coefficient of c_n is the Hankel coefficient ((2n-1)!!)^2/(n! 16^n) of K_0.
Usage: py -I -B g5_t6.py N      Log: logs/g5_t6_<N>.log
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D, g5_theory as T

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_t6_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def expser(lg, N):
    E = [Fr(1)] + [Fr(0)] * N
    for m in range(1, N + 1):
        E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
    return E


def run(groups, sing, N, pred):
    ct = C.engine(D.build_loc(groups, sing, N, N), N, N)
    first, vals, poles = C.laurent_report(ct, N, N)
    bad = [n for n in range(N + 1) if vals[n] != pred[n]]
    return first, bad, vals, ct


def sing_of(c, e, M=1):
    return dict(M=M, cs=Fr(c), es=Fr(e), bstar=[0] + [D.beta(e, n) for n in range(1, N + 1)])


def main():
    B = 2 * N + 2
    out('g5_t6: hypotheses of the chart argument, z^0..z^%d, exact over Q' % N)
    # A: quintic at q = oo
    gr = D.quintic_groups(B, 'qinf'); sg = sing_of(5, 1)
    ginf = expser([Fr(0)] + [(-1) ** (n // 5 + 1) * D.bpoly(n + 1, Fr(1, 5)) / ((n // 5) * (n + 1)) if n % 5 == 0 else
                             Fr(0) for n in range(1, N + 1)], N)
    first, bad, vals, _ = run(gr, sg, N, ginf)
    check('A quintic, L -> 0 (q = oo): no pole and c_n(0) = [z^n] Ginf for n <= %d (Ginf from B_(5l+1)(1/5))' % N,
          first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
    check('A  consistency: Ginf = the shifted-Stirling prediction of g5_data.predicted', ginf == D.predicted(gr, sg, N))
    nos = D.predicted(gr, sg, N, shift=False)
    check('A  control: the UNshifted regular tail is NOT the value (first mismatch at z^%s)'
          % next((n for n in range(N + 1) if nos[n] != vals[n]), None), any(nos[n] != vals[n] for n in range(N + 1)))
    # B: general weights at q = oo
    for name, lam in (('W1', [1, 3, -2, Fr(5, 2), Fr(-7, 3)]), ('W2', [2, -1, 3, -5, 1]), ('W3', [1, Fr(1, 2), Fr(1, 3), Fr(1, 4), Fr(1, 5)])):
        fac = T.base(lam) + [T.convex(5)]
        gr, s0 = T.groups_for(fac, 5, B); sg = sing_of(5, 1)
        pred = expser([Fr(0)] + [D.beta(Fr(1, 5), n) * sum(Fr(l) ** (-n) for l in lam) for n in range(1, N + 1)], N)
        first, bad, vals, _ = run(gr, sg, N, pred)
        check('B %s at q = oo (s* = %s): no pole, value exp(Sum beta_n(1/5) z^n p_-n), n <= %d' % (name, s0, N),
              first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
    # C: local P^2 at q = oo (orbifold point)
    for mu in (Fr(0), Fr(1, 2)):
        fac = T.base([1, 3, -2]) + [T.concave(3, mu)]
        gr, s0 = T.groups_for(fac, 3, B); sg = sing_of(3, 0)
        pred = D.predicted(gr, sg, N)
        first, bad, vals, _ = run(gr, sg, N, pred)
        check('C local P^2 (lift -3H + %s) at q = oo, s* = %s: no pole, unshifted value (e* = 0), n <= %d' % (mu, s0, N),
              first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
    # D: shifted singular factor Gamma(1/3 + (lambda_i - sigma)/z)
    lam = [1, 3, -2, Fr(5, 2), Fr(-7, 3)]
    fac = T.base(lam) + [T.convex(5)]
    for i in (0, 2):
        gr, s0 = T.groups_for(fac, i, B); sg = sing_of(-1, Fr(1, 3))
        pred = D.predicted(gr, sg, N)
        first, bad, vals, _ = run(gr, sg, N, pred)
        nos = D.predicted(gr, sg, N, shift=False)
        check('D singular Gamma(1/3 + (lambda_%d - sigma)/z): no pole, value with shifts e_f + c_f/3, n <= %d; the '
              'unshifted value fails at z^%s' % (i, N, next((n for n in range(N + 1) if nos[n] != vals[n]), None)),
              first is None and not bad and nos != vals, 'first pole %s, mismatches %s' % (first, bad[:4]))
    # E: weighted P(1,2)
    fac = [dict(c=Fr(-1), d=Fr(3), e=Fr(0), kind='base'), dict(c=Fr(-2), d=Fr(-1), e=Fr(0), kind='base')]
    for i in (0, 1):
        gr, s0 = T.groups_for(fac, i, B); sg = sing_of(fac[i]['c'], 0)
        pred = D.predicted(gr, sg, N)
        first, bad, vals, _ = run(gr, sg, N, pred)
        check('E P(1,2), Gamma(3/z - s) Gamma(-1/z - 2s), colour s* = %s (c* = %s): no pole, value = regular tail, n <= %d'
              % (s0, fac[i]['c'], N), first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED so far: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main()
    import g5_t6b
    g5_t6b.main(N, out, check, FAILS)
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)
    out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
