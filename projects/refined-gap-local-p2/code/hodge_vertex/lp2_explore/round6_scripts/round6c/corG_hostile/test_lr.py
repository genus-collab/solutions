# Large-radius side: same code, F = theta log A, against GV (table) and refined BPS (table).
# usage: py test_lr.py g Dmax
import sys, time
from gsum import *

GV = {
 0: [3, -6, 27, -192, 1695, -17064, 188454, -2228160, 27748899, -360012150],
 1: [0, 0, -10, 231, -4452, 80948, -1438086, 25301295, -443384578, 7760515332],
 2: [0, 0, 0, -102, 5430, -194022, 5784837, -155322234, 3894455457, -93050366010],
 3: [0, 0, 0, 15, -3672, 290853, -15363990, 649358826, -23769907110, 786400843911],
 4: [0, 0, 0, 0, 1386, -290400, 29056614, -2003386626, 109496290149, -5094944994204],
}
# refined BPS N^d_{jL,jR}, d<=4 (two j's times 2 as integers)
REF = {
 1: {(0, 2): 1},
 2: {(0, 5): 1},
 3: {(0, 6): 1, (1, 9): 1},
 4: {(0, 5): 1, (0, 9): 1, (0, 13): 1, (1, 8): 1, (1, 10): 1, (1, 12): 1, (2, 11): 1, (3, 14): 1},
}

def ser_mul(a, b, n):
    r = [Fr(0)] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]):
                r[i + j] += x * y
    return r

def sinc_pow(p, n):
    # (sin(x/2)/(x/2))^p, even powers kept, in x, up to x^(n-1)
    s = [Fr(0)] * n
    for m in range(0, (n + 1) // 2):
        if 2 * m < n:
            s[2 * m] = Fr((-1) ** m, factorial(2 * m + 1) * 4 ** m)
    if p == 0:
        return [Fr(1)] + [Fr(0)] * (n - 1)
    base = s
    if p < 0:
        inv = [Fr(0)] * n
        inv[0] = Fr(1)
        for k in range(1, n):
            inv[k] = -sum(s[i] * inv[k - i] for i in range(1, k + 1))
        base = inv
        p = -p
    r = [Fr(1)] + [Fr(0)] * (n - 1)
    for _ in range(p):
        r = ser_mul(r, base, n)
    return r

def N_gd(g, d):
    tot = Fr(0)
    for k in range(1, d + 1):
        if d % k:
            continue
        dd = d // k
        for h in range(0, g + 1):
            if h not in GV:
                raise ValueError('need n^%d' % h)
            n = GV[h][dd - 1]
            if n == 0:
                continue
            c = sinc_pow(2 * h - 2, 2 * g - 2 * h + 1)[2 * g - 2 * h]
            tot += n * Fr(1, k) * Fr(k) ** (2 * g - 2) * c
    return tot

def exp_series(c, n):
    # exp(c x) up to x^(n-1)
    return [c ** i / factorial(i) for i in range(n)]

def refined_GW(g, d, e1, e2):
    """lambda^(2g-2) coefficient, Q^d, of the refined BPS sum, multi-covers included."""
    tot = Fr(0)
    # work with lam series of length 2g+2 after pulling out lam^-2
    n = 2 * g + 1
    for k in range(1, d + 1):
        if d % k:
            continue
        dd = d // k
        if dd not in REF:
            raise ValueError
        # 1/(2 sinh(k e1 lam/2) 2 sinh(k e2 lam/2)) = lam^-2/(k^2 e1 e2) * sinhc(k e1 lam)^-1 sinhc(k e2 lam)^-1
        def shc_inv(a):
            s = [Fr(0)] * n
            for m in range(0, (n + 1) // 2):
                if 2 * m < n:
                    s[2 * m] = Fr(a) ** (2 * m) / (factorial(2 * m + 1) * 4 ** m)
            inv = [Fr(0)] * n
            inv[0] = Fr(1)
            for j in range(1, n):
                inv[j] = -sum(s[i] * inv[j - i] for i in range(1, j + 1))
            return inv
        den = ser_mul(shc_inv(k * e1), shc_inv(k * e2), n)
        for (twojL, twojR), N in REF[dd].items():
            sign = (-1) ** (twojL + twojR)
            # chi_j(u^k), u = exp((e1-e2)lam/2), v=exp((e1+e2)lam/2): chi_j(x)=sum_{m=-j..j} x^{2m}
            def chi(twoj, a):
                r = [Fr(0)] * n
                for mm in range(-twoj, twoj + 1, 2):
                    # x^{2m}, 2m=mm, x=exp(k a lam /2): exp(mm k a lam/2)
                    ex = exp_series(Fr(mm * k) * a / 2, n)
                    r = [p + q for p, q in zip(r, ex)]
                return r
            chiL = chi(twojL, e1 - e2)
            chiR = chi(twojR, e1 + e2)
            s = ser_mul(ser_mul(chiL, chiR, n), den, n)
            tot += sign * N * Fr(1, k) * Fr(1, k * k) / (e1 * e2) * s[2 * g]
    return tot

def gw_from_G(F, tot, su, Dmax):
    """Given tot (series in t), check it is a power series in y=t^3 and invert the mirror map."""
    # negative powers must vanish
    neg = [n for n in range(tot.val, 0) if n < tot.prec and not F.is0(tot.coeff(n))]
    nonpow = [n for n in range(0, min(tot.prec, 3 * Dmax + 3)) if n % 3 and not F.is0(tot.coeff(n))]
    ys = [tot.coeff(3 * d) for d in range(Dmax + 1)]
    # Q(y)=y Qy(y): subtract successively
    Qy = [F.fr(x) for x in su['Qy'][:Dmax + 1]]
    res = []
    work = list(ys)
    # powers of Qy
    pw = [F.one] + [F.zero] * Dmax
    for d in range(Dmax + 1):
        c = work[d]
        res.append(c)
        # subtract c * Q^d = c y^d Qy^d
        for i in range(Dmax + 1 - d):
            work[d + i] = F.red(work[d + i] - c * pw[i])
        # update pw*=Qy
        new = [F.zero] * (Dmax + 1)
        for i in range(Dmax + 1):
            s = F.zero
            for j in range(i + 1):
                s = s + pw[j] * Qy[i - j]
            new[i] = F.red(s)
        pw = new
    return res, neg, nonpow

def solve_mod(A, b, p):
    n = len(A)
    M = [row[:] + [x] for row, x in zip(A, b)]
    for i in range(n):
        piv = next(r for r in range(i, n) if M[r][i] % p)
        M[i], M[piv] = M[piv], M[i]
        inv = pow(M[i][i], -1, p)
        M[i] = [x * inv % p for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] % p:
                f = M[r][i]
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]

if __name__ == '__main__':
    g = int(sys.argv[1])
    Dmax = int(sys.argv[2])
    F = FpField()
    K = needed_K(g)
    samples = [(Fr(3), Fr(7)), (Fr(1), Fr(1)), (Fr(5, 2), Fr(-3, 7)), (Fr(2), Fr(-5, 3)), (Fr(1, 3), Fr(1)), (Fr(4, 5), Fr(1))][:max(g + 1, 3)]
    Nt = 3 * Dmax + 25 + 12 * g
    allgw = []
    for (e1, e2) in samples:
        t0 = time.time()
        D, su = lr_data(F, e1, e2, Dmax, Nt, K, 'hodge')
        tot = D.total(g)
        res, neg, nonpow = gw_from_G(F, tot, su, Dmax)
        allgw.append(res)
        print('eps', e1, e2, 'time %.1f' % (time.time() - t0), 'val', tot.val, 'prec', tot.prec, 'neg t-powers nonzero:', neg[:5], ' non-y powers nonzero:', nonpow[:5])
        print(' GW_{g,d} d=0..4 :', [str(rr(x)) for x in res[:5]])
        okr = []
        for d in range(1, min(Dmax, 4) + 1):
            okr.append(F.is0(F.red(res[d] - F.fr(refined_GW(g, d, e1, e2)))))
        print(' refined BPS match d=1..4:', okr)
    # interpolate
    n = g + 1
    A = []
    for (e1, e2) in samples[:n]:
        E = e1 + e2
        A.append([F.fr(E ** (2 * g - 2 * j) * (e1 * e2) ** j) for j in range(n)])
    for d in range(0, Dmax + 1):
        b = [F.red(F.fr(e1 * e2) * gw[d]) for (e1, e2), gw in zip(samples, allgw)]
        f = solve_mod(A, b, F.p)
        if d == 0:
            print('d=0 f_j (j=0..g), f_g=F_{0,g} const:', [str(rr(x)) for x in f])
        else:
            exp = N_gd(g, d)
            print('d=%d F_{0,%d} = %s  GV N_{%d,%d} = %s  match %s ; others f_j: %s' % (d, g, rr(f[g]), g, d, exp, F.is0(F.red(f[g] - F.fr(exp))), [str(rr(x)) for x in f[:g]]))
    # consistency sample if available
    if len(samples) > n:
        for (e1, e2), gw in zip(samples[n:], allgw[n:]):
            pass
