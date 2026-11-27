"""g5_t1.py -- Stage 5 (general), test T1 (+ the quintic controls of T4).

Formal quintic, colour 0.  c(z; L) = normalised formal SP series at sigma = L of the Mellin-Barnes summand with all
six Gamma factors Stirling-expanded.  New method (not used in Stage 4): expansion at L = 1 + eps directly, as a
Laurent series in eps at every order z^n (g5_core: scaled variable T = (sigma - L)/eps, nu = z/eps, Lagrange-
Buermann).  Checks, for every n <= N:  (R) c_n(1 + eps) has no negative power of eps;  (V) [eps^0] c_n = h_n,
h = h_(w_0), w_0 = {1 - zeta^j}_(j=1..4) U {-5}.
Controls (n <= NC): wrong twist sign in h; Stirling series of Gamma(1+5s) dropped; Stirling series of the singular
factor dropped; b_1 -> b_1 + 1/7 everywhere; b_3 -> b_3 + 1/7 in the singular factor only / in the regular factors
only; and the positive control Gamma(5s) in place of Gamma(1+5s) (the chart argument predicts the same q = 0 value).
Usage: py -I -B g5_t1.py N [K] [Q|p] [NC]     (K = eps-truncation, default N; p = mod 2^61-1 and 2^89-1)
Log: logs/g5_t1_<N>_<K>_<mode>.log (relative to this script).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
K = int(sys.argv[2]) if len(sys.argv) > 2 else N
MODE = sys.argv[3] if len(sys.argv) > 3 else 'Q'
NC = int(sys.argv[4]) if len(sys.argv) > 4 else min(N, 12)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_t1_%d_%d_%s.log' % (N, K, MODE)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def short(x):
    s = str(x)
    return s if len(s) <= 44 else s[:20] + '..(%d chars)' % len(s)


def run(N, K, gopt=None, bfun=None, pred_groups=None):
    gr = D.quintic_groups(N + K + 1, 'q0', **(gopt or {}))
    sg = D.sing_base(N, bfun)
    ct = C.engine(D.build_loc(gr, sg, N, K), N, K)
    first, vals, poles = C.laurent_report(ct, N, K)
    h = D.predicted(pred_groups if pred_groups else D.quintic_groups(N + K + 1, 'q0'), D.sing_base(N), N)
    hv = [C.cv(x) for x in h]
    badv = [n for n in range(N + 1) if n <= K and not C.iszero(vals[n] - hv[n])]
    return first, badv, vals, h, poles


def main():
    primes = [None] if MODE == 'Q' else [2 ** 61 - 1, 2 ** 89 - 1]
    out('g5_t1: formal quintic colour 0, c(z; 1 + eps) as Laurent series in eps; z^0..z^%d, eps-truncation %d; field %s'
        % (N, K, 'Q (exact)' if MODE == 'Q' else 'F_p, p = 2^61-1 and 2^89-1'))
    for p in primes:
        C.setfield(p); t0 = time.time()
        first, badv, vals, h, poles = run(N, K)
        tag = 'Q' if p is None else 'p = 2^%d-1' % (p.bit_length())
        if p is None:
            for n in range(min(N, 10) + 1):
                out('   n = %2d: eps-orders -%d..-1 all zero: %s;  [eps^0] c_n = %s;  h_n = %s'
                    % (n, n, n not in poles, short(vals[n]), short(h[n])))
        check('(R) [%s] no negative power of eps in c_n(1+eps), all n <= %d' % (tag, N), first is None,
              '' if first is None else 'first pole at z^%d: %s' % (first, poles[first][:5]))
        check('(V) [%s] [eps^0] c_n = h_n for all n <= %d (incl. (N) at n = 5, 10, ..)' % (tag, min(N, K)), not badv,
              badv[:5])
        out('   (%.1f s)' % (time.time() - t0))
    C.setfield(None)
    out('Controls (exact over Q, z^0..z^%d):' % NC)
    b = lambda n: D.beta(0, n)
    wrong = [D.group_rational(5, 1, [-5], NC + NC + 1), D.quintic_groups(NC + NC + 1, 'q0')[1]]
    ctrls = [('wrong twist sign in h (w contains +5)', dict(pred_groups=wrong), 'V'),
             ('Stirling series of Gamma(1+5s) dropped', dict(gopt={'tw': {'drop': True}}), 'V'),
             ('Stirling series of the singular factor dropped', dict(bfun=lambda n: Fr(0)), 'R'),
             ('b_1 + 1/7 in every factor', dict(gopt={'tw': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 1 else 0)},
                                                   'base': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 1 else 0)}},
                                             bfun=lambda n: b(n) + (Fr(1, 7) if n == 1 else 0)), 'R'),
             ('b_3 + 1/7 in the singular factor only', dict(bfun=lambda n: b(n) + (Fr(1, 7) if n == 5 else 0)), 'R'),
             ('b_3 + 1/7 in the regular factors only', dict(gopt={'tw': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 5 else 0)},
                                                             'base': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 5 else 0)}}), 'V')]
    for name, kw, expect in ctrls:
        first, badv, vals, h, poles = run(NC, NC, **kw)
        fv = badv[0] if badv else None
        out('   control (%s): first pole at z^%s; first value mismatch at z^%s' % (name, first, fv))
        if expect == 'R':
            check('control fails through (R), as predicted: %s' % name, first is not None)
        else:
            check('control keeps (R) and fails through (V), as predicted: %s' % name, first is None and fv is not None)
    # b_3 + 1/7 in the regular factors only: the local Mellin-inversion prediction is the PERTURBED regular tail
    pb = {'tw': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 5 else 0)}, 'base': {'beta': lambda n: b(n) + (Fr(1, 7) if n == 5 else 0)}}
    first, badv, vals, h, poles = run(NC, NC, gopt=pb)
    lg = [Fr(0)] * (NC + 1)
    for g in D.quintic_groups(2 * NC + 1, 'q0'):
        for n in range(1, NC + 1):
            lg[n] += (D.beta(0, n) + (Fr(1, 7) if n == 5 else 0)) * g['ps'][n]
    E = [Fr(1)] + [Fr(0)] * NC
    for m in range(1, NC + 1):
        E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
    check('regular-only b_3 perturbation: value = exp(Sum z^n A_n(1)) of the PERTURBED regular factors, n <= %d' % NC,
          first is None and all(vals[n] == E[n] for n in range(NC + 1)))
    first, badv, vals, h, poles = run(NC, NC, gopt={'tw': {'e': Fr(0)}})
    check('positive control Gamma(5s) for Gamma(1+5s): (R) and (V) still hold, n <= %d (chart: only the elementary '
          'factor changes)' % NC, first is None and not badv)
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main(); out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
