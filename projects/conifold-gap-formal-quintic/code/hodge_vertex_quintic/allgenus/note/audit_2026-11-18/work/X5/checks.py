from fractions import Fraction as Fr
from math import factorial, comb

def bernoulli(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B

B = bernoulli(30)
N = 12
# a(x) = sum_{k>=1} x^k/(2k+1); compute [x^h] (1+a)^{1-2h} via power series
def series_pow(c, e, N):
    # c: list with c[0]=1; returns (c)^e for rational e (integer here) up to x^N via log/exp-free recursion
    # use J.C.P. Miller recurrence: b[0]=1, b[n] = 1/(n c0) sum_{k=1}^n ((e+1)k - n) c[k] b[n-k]
    b = [Fr(0)] * (N + 1)
    b[0] = Fr(1)
    for n in range(1, N + 1):
        s = Fr(0)
        for k in range(1, n + 1):
            s += ((e + 1) * k - n) * c[k] * b[n - k]
        b[n] = s / n
    return b
c = [Fr(1)] + [Fr(1, 2 * k + 1) for k in range(1, N + 1)]
for h in range(1, 7):
    p = series_pow(c, 1 - 2 * h, N)
    sig = factorial(2 * h - 2) * p[h]
    rhs = -Fr(4) ** h * B[2 * h] / (2 * h)
    # also direct sum over j of (-1)^j (2h-2+j)!/j! [x^h] a^j
    # a^j coefficient
    a = [Fr(0)] + [Fr(1, 2 * k + 1) for k in range(1, N + 1)]
    def mul(u, v):
        w = [Fr(0)] * (N + 1)
        for i in range(N + 1):
            if u[i] == 0: continue
            for j in range(N + 1 - i):
                w[i + j] += u[i] * v[j]
        return w
    aj = [Fr(1)] + [Fr(0)] * N
    direct = Fr(0)
    for j in range(1, h + 1):
        aj = mul(aj, a)
        direct += (-1) ** j * Fr(factorial(2 * h - 2 + j), factorial(j)) * aj[h]
    print('h', h, 'Sigma', sig, 'direct', direct, 'claimed', rhs, sig == rhs == direct)

# c_g values
for g in range(2, 6):
    cg = (-1) ** (g - 1) * B[2 * g] / (2 * g * (2 * g - 2) * 5 ** (g - 1))
    print('c_g', g, cg)
# F^TR_h(v=1) = B_2h/(2h(2h-2))
for h in range(2, 6):
    print('FTR', h, B[2 * h] / (2 * h * (2 * h - 2)))
