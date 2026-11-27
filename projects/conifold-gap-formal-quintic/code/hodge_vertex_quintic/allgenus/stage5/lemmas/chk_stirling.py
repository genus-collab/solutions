# -*- coding: utf-8 -*-
"""chk_stirling.py -- Lemma L3 (Stirling's series = formal SP of Euler's integral), exact over Q.

Phase F(x) = a log x - x, critical point x = a, F''(a) = -1/a (covariance C = a); x = a + t.
x^(a/z+m) e^(-x/z) dx/x = e^(F/z) x^(m-1) dx  and  x^(a/z+m) e^(-x/z) dx = e^(F/z) x^m dx: in both cases the
normalised amplitude is (1 + t/a)^e (e = m-1, resp. e = m) and L3 predicts SP^norm = St_(e+1)(z/a), with
St_m(w) := exp((1/w + m - 1/2) log(1 + m w) - m + Sum_k b_k w^(2k-1) (1 + m w)^(1-2k)),  b_k = B_2k/(2k(2k-1)),
i.e. the series part of z^(a/z+m) Gamma_St(a/z+m).  Computed by Wick's theorem in one variable (the definition).
Also: the recursion St_(m+1)(w) = (1 + m w) St_m(w) (formal functional equation), the a-shift identity
(1+w)^(1/w+1/2) e^-1 St_0(w/(1+w)) = St_0(w) of the formal proof (note s.9), homogeneity in a.
Controls (must fail): b_2 -> -b_2; measure dx/x versus dx at m = 1 (St_1 versus St_2); St_0 used for m = 2.
Usage: py -I -B chk_stirling.py [N] [NA]   (a = 1 through z^N, default 40; a = 3/7, -2 through z^NA, default 12)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from chk_util import Log, first_diff

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NA = int(sys.argv[2]) if len(sys.argv) > 2 else 12
log = Log(HERE, 'chk_stirling_%d_%d' % (N, NA))
log.out('chk_stirling: Lemma L3 through z^%d at a = 1, z^%d at a = 3/7 and -2; exact over Q; log %s' % (N, NA, log.rel))


def bernoulli(n):
    B = [Fr(1)]
    for k in range(1, n + 1):
        B.append(-sum(comb(k + 1, j) * B[j] for j in range(k)) / (k + 1))
    return B


def smul(A, B, n):
    R = [Fr(0)] * (n + 1)
    for i, x in enumerate(A[:n + 1]):
        if x:
            for j, y in enumerate(B[:n + 1 - i]):
                R[i + j] += x * y
    return R


def sexp(A, n):                              # exp of a series with A[0] = 0 (A' E = E')
    E = [Fr(1)] + [Fr(0)] * n; D = [k * A[k] for k in range(n + 1)]
    for k in range(1, n + 1):
        E[k] = sum(D[j] * E[k - j] for j in range(1, k + 1)) / k
    return E


def St(m, n, bk):
    """St_m(w) through w^n, from the Bernoulli numbers (b_k list, bk[k] for k >= 1)"""
    L = [Fr(0)] * (n + 1)
    lg = [Fr(0)] + [Fr((-1) ** (j + 1) * m ** j, j) for j in range(1, n + 2)]          # log(1 + m w)
    for k in range(n + 1):
        L[k] += lg[k + 1] + (m - Fr(1, 2)) * lg[k]                                     # (1/w + m - 1/2) log(1+mw)
    L[0] -= m
    for k in range(1, (n + 1) // 2 + 1):
        p = 2 * k - 1
        for i in range(n + 1 - p):                                                     # w^p (1 + m w)^(1-2k)
            L[p + i] += bk[k] * gbin(-p, i) * m ** i
    assert L[0] == 0
    return sexp(L, n)


def gbin(r, i):                               # generalised binomial coefficient C(r, i), r integer
    out = Fr(1)
    for j in range(i):
        out = out * (r - j) / (j + 1)
    return out


def sp_euler(a, e_list, n):
    """normalised formal SP at x = a of e^((a log x - x)/z) (1 + t/a)^e, z^0..z^n, by Wick: <t^2k> = (2k-1)!! a^k"""
    D = 6 * n + 2
    f3 = [Fr(0)] * 3 + [Fr((-1) ** (k + 1), k) / a ** (k - 1) for k in range(3, D + 1)]
    mom = [Fr(0)] * (D + 1); mom[0] = Fr(1)
    for d in range(2, D + 1, 2):
        mom[d] = mom[d - 2] * (d - 1) * a
    P = [Fr(1)] + [Fr(0)] * D; Ps = []
    for mp in range(2 * n + 1):
        Ps.append(P[:2 * n + 2 * mp + 1])
        P = [x / (mp + 1) for x in smul(P, f3, 2 * n + 2 * (mp + 1))]
    out = {}
    for e in e_list:
        amp = [gbin(e, i) / a ** i for i in range(D + 1)]
        res = [Fr(0)] * (n + 1)
        for mp, Pm in enumerate(Ps):
            for k in range(n + 1):
                d = 2 * (k + mp)                     # f_3^mp has t-degree >= 3 mp
                res[k] += sum(Pm[i] * amp[d - i] for i in range(3 * mp, min(d, len(Pm) - 1) + 1)) * mom[d]
        out[e] = res
    return out


B = bernoulli(2 * N + 4); bk = [None] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, N + 2)]
t = time.time(); E_LIST = [-3, -2, -1, 0, 1, 2, 3, 4]
sp = sp_euler(Fr(1), E_LIST, N)
log.out('     Wick computation at a = 1 through z^%d: %.1f s' % (N, time.time() - t))
Sts = {m: St(m, N, bk) for m in range(-2, 6)}
for e in E_LIST:
    m_dxx, m_dx = e + 1, e                       # dx/x with x^m: e = m - 1;  dx with x^m: e = m
    fd = first_diff(dict(enumerate(sp[e])), dict(enumerate(Sts[e + 1])), N)
    log.check('L3 a = 1: SP of x^(a/z%+d) e^(-x/z) dx/x  and of  x^(a/z%+d) e^(-x/z) dx  = St_%d(z) through z^%d'
              % (m_dxx, m_dx, e + 1, N), fd is None, 'z^1, z^2: %s, %s' % (sp[e][1], sp[e][2]))
for m in range(-2, 5):
    rec = smul(Sts[m], [Fr(1), Fr(m)], N)
    log.check('functional equation St_%d(w) = (1 + %d w) St_%d(w) through w^%d' % (m + 1, m, m, N), rec == Sts[m + 1])
psi = Sts[0]; L = [Fr(0)] + [Fr((-1) ** (j + 1), j) for j in range(1, N + 2)]
ex = [L[k + 1] + Fr(1, 2) * L[k] for k in range(N + 1)]; ex[0] -= 1             # (1/w + 1/2) log(1+w) - 1
lhs = smul(sexp(ex, N), [sum(psi[k] * gbin(-k, i - k) for k in range(i + 1)) for i in range(N + 1)], N)
log.check('a-shift identity (1+w)^(1/w+1/2) e^-1 St_0(w/(1+w)) = St_0(w) through w^%d' % N, lhs == psi)
for a in (Fr(3, 7), Fr(-2)):
    spa = sp_euler(a, [-1, 0, 2], NA)
    ok = all(spa[e][k] == sp[e][k] / a ** k for e in (-1, 0, 2) for k in range(NA + 1))
    log.check('L3 homogeneity a = %s: SP^norm = St_(e+1)(z/a) through z^%d (e = -1, 0, 2)' % (a, NA), ok)
bad = list(bk); bad[2] = -bad[2]
fd = first_diff(dict(enumerate(sp[0])), dict(enumerate(St(1, N, bad))), N)
log.control('b_2 -> -b_2', fd is not None, 'first differs at z^%s' % fd)
fd = first_diff(dict(enumerate(sp[1])), dict(enumerate(Sts[1])), N)
log.control('measure dx/x versus dx at m = 1 (St_1 instead of St_2)', fd is not None, 'first differs at z^%s' % fd)
fd = first_diff(dict(enumerate(sp[1])), dict(enumerate(Sts[0])), N)
log.control('St_0 used for m = 2 (dx/x)', fd is not None, 'first differs at z^%s' % fd)
log.out('     St_0 = exp(Sum b_k w^(2k-1)): w^1..w^4 = %s' % [str(x) for x in Sts[0][1:5]])
log.finish('CHK_STIRLING %d %d' % (N, NA))
