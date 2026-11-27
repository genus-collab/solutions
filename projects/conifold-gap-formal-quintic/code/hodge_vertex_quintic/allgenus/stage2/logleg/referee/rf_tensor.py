"""Referee: the variational formula (h,n) -> (h,n+1) as a TENSOR identity in the n old slots, for n <= 2.
Pairs both sides with p^m in every old slot for ALL m <= M, M >= 10(kmax+1) + 10, which determines a differential
with poles of order <= 2 kmax + 2 at the five branch points (expansion at infinity, s.3 of REFEREE.md).
Identity tested (v = 1, homogeneity for d_v|_p):  <m, -log p>_h = e <m>_h + s * Sum_l m_l Sum_q <.., m_l-5-5q, ..>_h,
s = +1 (claim), controls s = -1, 0 must fail.
usage: py -I -B rf_tensor.py PRIME_EXP KM h,n [h,n ...]"""
import sys, time
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
from rf_rec import Curve
from rf_build import build, closure
from rf_pair import leg_pow, leg_log, contract
from rf_tr import inv

def main():
    e2, KM = int(sys.argv[1]), int(sys.argv[2])
    cases = [tuple(int(x) for x in s.split(',')) for s in sys.argv[3:]]
    P = 2 ** e2 - 1; i5 = inv(5, P)
    t0 = time.time()
    cv = Curve(P, KM)
    W = build(cv, closure([(h, n + 1) for h, n in cases]), log=lambda s: None)
    LL = leg_log(cv)
    for (h, n) in cases:
        kmax = 3 * h - 2 + n
        M = 10 * (kmax + 1) + 10
        Lp = {m: leg_pow(cv, m) for m in range(1, M + 1)}
        T0, T1 = W[(h, n)], contract(W[(h, n + 1)], LL, P)
        def pr(T, ms):
            for m in ms: T = contract(T, Lp[m], P)
            return T.get((), 0)
        if n == 1:
            tuples = [(m,) for m in range(1, M + 1)]
        else:
            tuples = [(a, b) for a in range(1, M + 1) for b in range(a, M + 1)]
        # cache single contractions
        C0 = {m: contract(T0, Lp[m], P) for m in range(1, M + 1)} if n == 2 else None
        C1 = {m: contract(T1, Lp[m], P) for m in range(1, M + 1)} if n == 2 else None
        def val(which, ms):
            if n == 1: return pr(T0 if which == 0 else T1, ms)
            C = C0 if which == 0 else C1
            return contract(C[ms[0]], Lp[ms[1]], P).get((), 0)
        fails = {1: 0, -1: 0, 0: 0}; nz = 0
        for ms in tuples:
            lhs = val(1, ms)
            nz += lhs != 0
            e = ((2 - 2 * h - n) + sum(ms) * i5) % P
            xi = 0
            for i, m in enumerate(ms):
                q = 0
                while m - 5 - 5 * q >= 1:
                    ms2 = tuple(sorted(ms[:i] + (m - 5 - 5 * q,) + ms[i + 1:]))
                    xi += m * val(0, ms2); q += 1
            for s in fails:
                if (e * val(0, ms) + s * xi - lhs) % P: fails[s] += 1
        print(f"  [{'PASS' if fails[1] == 0 and fails[-1] and fails[0] else 'FAIL'}] ({h},{n})->({h},{n+1}) tensor identity, "
              f"all m_i <= {M}: {len(tuples)} tuples ({nz} non-zero LHS), failures: + {fails[1]}, controls - {fails[-1]}, no-xi {fails[0]}")
    print(f"  total {time.time()-t0:.1f} s (P = 2^{e2}-1)")

if __name__ == '__main__':
    main()
