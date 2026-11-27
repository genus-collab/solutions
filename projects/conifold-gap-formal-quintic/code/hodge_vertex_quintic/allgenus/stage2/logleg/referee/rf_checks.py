"""Referee: (V1) step, dilaton form (D), controls, and the worker's quoted numbers, on top of rf_run.
usage: py -I -B rf_checks.py PRIME_EXP KM g,n [g,n ...]"""
import sys, time
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
from fractions import Fraction
import rf_run
from rf_pair import leg_pow, leg_log, contract, modq, ratrec
from rf_tr import inv

def pair(cv, W, h, legs):
    """full pairing of omega_{h,n} with the list of leg functionals"""
    T = W[(h, len(legs))]
    for L in legs:
        T = contract(T, L, cv.P)
    return T.get((), 0)

def main():
    cv, W, vals, order = rf_run.main()
    P = cv.P; i5 = inv(5, P)
    # ---- (V1): <m, (-log p)^(k+1)> = e <m,(-log p)^k> + sgn * Sum_i Sum_q m_i <.., m_i-5-5q, ..>
    res = {+1: [0, 0], -1: [0, 0], 0: [0, 0]}   # sign -> [tested, failed]
    for (h, a, k1, ms), lhs in vals.items():
        k = k1 - 1
        if k < 0: continue
        n = a + k
        if (h, n) == (0, 2) or 2 * h - 2 + n <= 0: continue
        if a >= 1 and sum(ms) > 5 * (2 * h - 2 + n) + 5: continue
        base = vals.get((h, a, k, ms))
        if base is None: continue
        e = ((2 - 2 * h - n) + sum(ms) * i5) % P
        xi = 0
        okall = True
        for i, m in enumerate(ms):
            q = 0
            while m - 5 - 5 * q >= 1:
                ms2 = tuple(sorted(ms[:i] + (m - 5 - 5 * q,) + ms[i + 1:]))
                v2 = vals.get((h, a, k, ms2))
                if v2 is None: okall = False
                else: xi += m * v2
                q += 1
        if not okall: continue
        for sg in (1, -1, 0):
            res[sg][0] += 1
            if (e * base + sg * xi - lhs) % P: res[sg][1] += 1
    print(f"  [{'PASS' if res[1][1] == 0 else 'FAIL'}] (V1) with + sign: {res[1][0]} tuples, {res[1][1]} failures")
    print(f"  [{'PASS' if res[-1][1] > 0 else 'FAIL'}] control (V1) with - sign must fail: {res[-1][1]} failures of {res[-1][0]}")
    print(f"  [{'PASS' if res[0][1] > 0 else 'FAIL'}] control (V1) without xi-term must fail: {res[0][1]} failures of {res[0][0]}")
    # ---- dilaton (D): <J,-log p>_h = (2-2h-n) <J>_h - (1/5) <J, p^5>_h,  J = (ms, (-log p)^k)
    nt = nf = 0
    for (h, a, k1, ms), lhs in vals.items():
        if k1 < 1: continue
        k = k1 - 1; n = a + k
        if (h, n) == (0, 2) or 2 * h - 2 + n <= 0: continue
        b0 = vals.get((h, a, k, ms)); m5 = tuple(sorted(ms + (5,)))
        b5 = vals.get((h, a + 1, k, m5))
        if b0 is None or b5 is None: continue
        nt += 1
        if ((2 - 2 * h - n) * b0 - i5 * b5 - lhs) % P: nf += 1
    print(f"  [{'PASS' if nf == 0 else 'FAIL'}] dilaton (D) (EO Thm 4.7, no (Q1)): {nt} tuples, {nf} failures")
    # ---- controls: off-L legs and a log leg based at p = 2
    LL = leg_log(cv); LL2 = leg_log(cv, 2)
    for (h, n) in order:
        chi = 2 * h - 2 + n
        found = []
        for j in (1, 2, 5, 10, 15):
            Lj = leg_pow(cv, -j)
            if n == 1:
                val = pair(cv, W, h, [Lj])
                if val: found.append(((-j,), ratrec(val, P)))
                continue
            for m in range(1, 5 * chi + j):
                if (m - j) % 5: continue
                if m - j >= 5 * chi: break
                legs = [Lj, leg_pow(cv, m)] + [LL] * (n - 2)
                val = pair(cv, W, h, legs)
                if val: found.append(((-j, m), ratrec(val, P)))
        nz2 = 0
        if n >= 2:
            for m in range(1, 5 * chi):
                legs = [LL2, leg_pow(cv, m)] + [LL] * (n - 2)
                if pair(cv, W, h, legs): nz2 += 1
        print(f"  [{'PASS' if found else 'FAIL'}] ({h},{n}) control off-L legs p^-j: {len(found)} non-zero singular, e.g. {found[:3]}"
              + (f"; log leg based at p=2: {nz2} non-zero singular" if n >= 2 else ""))
    # ---- the worker's quoted numbers
    LP = lambda m: leg_pow(cv, m)
    quoted = [((4, [LP(40), LL]), Fraction(-868220928)), ((3, [LP(35), LL, LL]), Fraction(-144609822225, 1024)),
              ((2, [LP(15)]), Fraction(2079, 128)), ((3, [LP(25)]), Fraction(3132675, 1024)),
              ((4, [LP(35)]), Fraction(-45541532421, 32768)), ((4, [LP(-5), LP(25)]), Fraction(-518283675, 16)),
              ((4, [LP(-2), LP(37)]), Fraction(-2084566248765, 1024)), ((1, [LP(-5)]), Fraction(15, 4)),
              ((0, [LP(-1), LP(1), LL]), Fraction(-1)), ((3, [LP(-5)]), Fraction(0)), ((4, [LP(-5)]), Fraction(0)),
              ((0, [LP(1), LP(4), LL]), Fraction(4)), ((0, [LP(5), LL, LL]), Fraction(-5))]
    for (h, legs), want in quoted:
        if (h, len(legs)) not in W: continue
        val = pair(cv, W, h, legs)
        print(f"  [{'PASS' if val == modq(want, P) else 'FAIL'}] worker's number at ({h},{len(legs)}): {want}; referee {ratrec(val, P)}")
    # ---- F_h from its definition (2-2h) F_h = Sum Res Phi omega_{h,1}, Phi = p^5/5 - v log p
    for h in range(2, 5):
        if (h, 1) not in W: continue
        F = ((i5 * pair(cv, W, h, [LP(5)]) + pair(cv, W, h, [LL])) * inv(2 - 2 * h, P)) % P
        print(f"  F_{h}(v=1) from the TR definition = {ratrec(F, P)}; <p^5>_{h} = {ratrec(pair(cv, W, h, [LP(5)]), P)}")

if __name__ == '__main__':
    t0 = time.time(); main(); print(f"  total {time.time()-t0:.1f} s")
