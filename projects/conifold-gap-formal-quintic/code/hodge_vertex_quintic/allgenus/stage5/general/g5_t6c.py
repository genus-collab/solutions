"""g5_t6c.py -- further hypothesis tests (exact over Q).

F3 no regular phase at all (Phi = 0), random A_0, A_n: the proof sketch's Legendre step needs Phi''(s*) != 0, the
   statement should not (the coefficients are polynomial in the Taylor data of Phi).  Prediction: no pole, value
   exp(Sum_n z^n A_n(s*)).
H1 TWO singular factors with different slopes at one point: O(5)-twisted P^4 with lambda_i = 0 at q = oo, where
   Gamma(-sigma/z) (c = -1) and Gamma(1 + 5 sigma/z) (c = 5) both degenerate at sigma = 0 (four critical points
   tend to 0 as q -> oo, sigma^4 ~ 1/q).  The chart is degenerate; prediction: poles.
H2 the elementary part and the Bernoulli part of the singular factor from DIFFERENT Gammas: amplitude of
   Gamma_St(1/2 + Y) with the Bernoulli series of Gamma_St(Y) (= Gamma_St(Y) Y^(1/2), i.e. a factor (sigma - s*)^(1/2)
   in B).  Prediction: poles.  Positive control: amplitude of Gamma_St(1 + Y) with the series of Gamma_St(Y)
   (= Y Gamma_St(Y) = Gamma_St(1 + Y) exactly): no pole.
L0 the eps^0 layer: singular factor alone, ct_n(0) = 0 (and 1/12, 25/288, .. without its Stirling series).
Usage: py -I -B g5_t6c.py N      Log: logs/g5_t6c_<N>.log
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D, g5_theory as T
from g5_t6b import build_raw

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_t6c_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def main():
    K = N
    out('g5_t6c: further hypothesis tests, z^0..z^%d, exact over Q' % N)
    rnd = random.Random(1114)
    rr = lambda: Fr(rnd.randint(-9, 9), rnd.randint(1, 7))
    for trial in range(2):
        A0p = [rr() for _ in range(K + 1)]; An = [None] + [[rr() for _ in range(K + 1)] for _ in range(N)]
        ct = C.engine(build_raw([Fr(0)] * (K + 1), A0p, An, D.sing_base(N), N, K), N, K)
        first, vals, _ = C.laurent_report(ct, N, K)
        lg = [Fr(0)] + [An[n][0] for n in range(1, N + 1)]
        E = [Fr(1)] + [Fr(0)] * N
        for m in range(1, N + 1):
            E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
        bad = [n for n in range(N + 1) if vals[n] != E[n]]
        check('F3 Phi = 0, random A_0, A_n (trial %d): no pole, value exp(Sum z^n A_n(s*)), n <= %d' % (trial, N),
              first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
    lam = [Fr(0), Fr(3), Fr(-2), Fr(5, 2), Fr(-7, 3)]
    regs = [D.group_rational(-1, 0, lam[1:], N + K + 1)]       # base j != 0 at sigma* = 0: w_j = lambda_j
    sg = D.sing_base(N)                                       # Gamma((0 - sigma)/z)
    sg['extra'] = [dict(M=1, cs=5, es=1, bstar=[0] + [D.beta(1, n) for n in range(1, N + 1)])]
    ct = C.engine(D.build_loc(regs, sg, N, K), N, K)
    first, vals, poles = C.laurent_report(ct, N, K)
    orders = [max([n - r for r in poles.get(n, [])] or [0]) for n in range(N + 1)]
    check('H1 lambda_0 = 0 in the O(5)-twisted P^4 at q = oo (two singular factors, c = -1 and 5): poles, as the '
          'degenerate chart predicts', first is not None, 'first pole at z^%s; pole orders %s' % (first, orders[:10]))
    for es, expect_pole, name in ((Fr(1, 2), True, 'H2 amplitude of Gamma_St(1/2 + Y), Bernoulli series of Gamma_St(Y)'),
                                  (Fr(1), False, 'H2 positive control: amplitude of Gamma_St(1 + Y), series of Gamma_St(Y)')):
        gr = D.quintic_groups(N + K + 1, 'q0')
        sg = dict(M=1, cs=-1, es=es, bstar=[0] + [D.beta(0, n) for n in range(1, N + 1)])
        ct = C.engine(D.build_loc(gr, sg, N, K), N, K)
        first, vals, poles = C.laurent_report(ct, N, K)
        check('%s (quintic colour-0 regular part): %s' % (name, 'pole' if expect_pole else 'no pole'),
              (first is not None) == expect_pole, 'first pole at z^%s' % first)
    # L0: the leading layer eps^0 of ct_n is the singular factor alone (no regular part)
    ct = C.engine(D.build_loc([], D.sing_base(N), N, K), N, K)
    check('L0 singular factor alone: ct_n(0) = 0 for 1 <= n <= %d (one-Gamma identity, leading poles)' % N,
          all(ct[n][0] == 0 for n in range(1, N + 1)))
    ct = C.engine(D.build_loc([], D.sing_base(N, lambda n: Fr(0)), N, K), N, K)
    out('   L0 without the singular Stirling series: ct_1(0), ct_2(0), ct_3(0) = %s' % [str(ct[n][0]) for n in (1, 2, 3)])
    check('L0   so the z^1 residue is 1/12 - b_1 = 0', ct[1][0] == Fr(1, 12))
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main(); out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
