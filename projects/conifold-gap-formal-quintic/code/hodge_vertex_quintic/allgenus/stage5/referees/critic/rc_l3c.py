# -*- coding: utf-8 -*-
"""rc_l3c.py NW NL -- critic of the Stage-5 referees: check of the two NEW lemmas of the LMI referees (not needed by the chain).

Exact over Q, own code (stdlib only; nothing imported from lf_*, lm_*, ra_*, ch_* or the project).
  D  Formal DLMF 5.11.8 (lmi_formal REFEREE s.4.2; lmi_adversary Lemma F5.11.8), as an identity of polynomials in h, w^-1..w^-NW:
     (w + h - 1/2) log(1 + h/w) - h + Sum_k b_k (w + h)^(1-2k) = Sum_n beta_n(h) w^-n,  beta_n(h) = (-1)^(n+1) B_(n+1)(h)/(n(n+1)).
     Also the two inputs of the interpolation proof: the difference equation of the right side, and the value at h = 0.
  L  Lemma L3-C (lmi_adversary s.2.3): N[log t - t, t^(M-1) dt; t = 1](z) = exp(Sum_n beta_n(M) z^n) as polynomials in M through
     z^NL, by a direct Wick sum with SYMBOLIC M; and, for integer M = -3..4, the Pochhammer form Poch_M(z) exp(Sum_k b_k z^(2k-1))
     of (L3) (CP2 s.2.8, FSL4 s.9).  Controls: amplitude t^M (wrong measure); beta_n at h = M + 1/5.
Usage: py -I -B rc_l3c.py 60 10   (log: logs/rc_l3c_NW_NL.log next to this script)."""
import os, sys, time
from fractions import Fraction as Fr
from math import comb, factorial
HERE = os.path.dirname(os.path.abspath(__file__))
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 60
NL = int(sys.argv[2]) if len(sys.argv) > 2 else 10
LOGREL = 'logs/rc_l3c_%d_%d.log' % (NW, NL)
OUT, NFAIL, T0 = [], [0], time.time()


def say(s):
    print(s)
    OUT.append(s)


def check(name, ok, control=False):
    good = (not ok) if control else ok
    say(('PASS ' if good else 'FAIL ') + name + ('  [control: must fail; it fails]' if control and good else
                                                 '  [control: must fail; it HOLDS]' if control else ''))
    NFAIL[0] += (not good)


def pnorm(p):
    while p and p[-1] == 0:
        p = p[:-1]
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return pnorm([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pmul(p, q):
    if not p or not q:
        return []
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return pnorm(out)


def psc(p, c):
    return pnorm([a * c for a in p])


def peval(p, x):
    return sum(a * x ** i for i, a in enumerate(p))


B = [Fr(1)]
for m in range(1, 2 * max(NW, NL) + 5):
    B.append(-sum(comb(m + 1, i) * B[i] for i in range(m)) / (m + 1))      # B_1 = -1/2
bk = [None] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, max(NW, NL) + 2)]


def bernpoly(m):                 # B_m(h) = sum_i C(m, i) B_i h^(m-i)
    p = [Fr(0)] * (m + 1)
    for i in range(m + 1):
        p[m - i] += comb(m, i) * B[i]
    return pnorm(p)


def beta(n):
    return psc(bernpoly(n + 1), Fr((-1) ** (n + 1), n * (n + 1)))


def gbinom(e, j):                # generalised binomial C(e, j) for an integer e
    num = Fr(1)
    for i in range(j):
        num *= (e - i)
    return num / factorial(j)


say('rc_l3c NW = %d, NL = %d (critic of the Stage-5 referees; exact over Q)' % (NW, NL))
# ---- D: formal 5.11.8
okD = True
for n in range(1, NW + 1):
    lhs = [Fr(0)] * (n + 2)
    lhs[n + 1] += Fr((-1) ** n, n + 1)                        # w * log term, m = n + 1
    lhs[n + 1] += Fr((-1) ** (n + 1), n)                      # (h - 1/2) * log term, m = n: h part
    lhs[n] += Fr((-1) ** (n + 1), n) * Fr(-1, 2)              #                              -1/2 part
    for k in range(1, (n + 1) // 2 + 1):
        j = n + 1 - 2 * k
        lhs[j] += bk[k] * gbinom(1 - 2 * k, j)
    okD &= pnorm(lhs) == beta(n)
check('D1 formal 5.11.8: coefficients of w^-1..w^-%d agree as polynomials in h' % NW, okD)
okD2 = all(
    pnorm([sum(c * comb(i, r) for i, c in enumerate(bernpoly(m)) if i >= r) for r in range(m + 1)]) ==
    padd(bernpoly(m), pnorm([Fr(0)] * (m - 1) + [Fr(m)])) for m in range(2, NW + 2))
check('D2 input of the interpolation proof: B_m(h + 1) = B_m(h) + m h^(m-1), m = 2..%d' % (NW + 1), okD2)
check('D3 input: beta_n(0) = b_k for n = 2k-1 and 0 for even n (n <= %d)' % NW,
      all(peval(beta(n), Fr(0)) == (bk[(n + 1) // 2] if n % 2 else 0) for n in range(1, NW + 1)))
# ---- L: Lemma L3-C by a Wick sum with symbolic M
EMAX = 2 * NL
ck = {k: Fr((-1) ** (k + 1), k) for k in range(3, EMAX + 3)}
E = [{0: Fr(1)}]                                              # E_a as {eta power: rational}
for a in range(1, EMAX + 1):
    acc = {}
    for k in range(1, a + 1):
        for p, v in E[a - k].items():
            acc[p + k + 2] = acc.get(p + k + 2, 0) + k * ck[k + 2] * v
    E.append({p: v / a for p, v in acc.items() if v})


def dfact(m):                    # (m-1)!! for even m, Gaussian moment <eta^m>
    if m % 2:
        return 0
    r = 1
    for i in range(1, m, 2):
        r *= i
    return r


def amp_coeffs(shift):           # C(M - 1 + shift, j) as polynomials in M, j = 0..EMAX
    out, p = [], [Fr(1)]
    for j in range(EMAX + 1):
        out.append(psc(p, Fr(1, factorial(j))))
        p = pmul(p, [Fr(-1 + shift - j), Fr(1)])
    return out


def wick(shift=0):
    am = amp_coeffs(shift)
    Nn = []
    for n in range(NL + 1):
        tot = []
        for a in range(2 * n + 1):
            j = 2 * n - a
            mom = sum(v * dfact(p + j) for p, v in E[a].items())
            if mom:
                tot = padd(tot, psc(am[j], Fr(mom)))
        Nn.append(tot)
    return Nn


def expser(bs):                  # exp(sum_n bs[n] z^n), polynomial coefficients
    G = [[Fr(1)]]
    for n in range(1, NL + 1):
        acc = []
        for k in range(1, n + 1):
            acc = padd(acc, psc(pmul(bs[k], G[n - k]), Fr(k)))
        G.append(psc(acc, Fr(1, n)))
    return G


Nw = wick()
Gb = expser([None] + [beta(n) for n in range(1, NL + 1)])
check('L1 Lemma L3-C: N[log t - t, t^(M-1); 1] = exp(Sum beta_n(M) z^n) as polynomials in M through z^%d' % NL, Nw == Gb)
check('L1 degree: [z^n] is a polynomial of degree 2n in M (n <= %d)' % NL, all(len(Nw[n]) == 2 * n + 1 for n in range(NL + 1)))
stir = [Fr(0)] * (NL + 1)
for k in range(1, NL + 1):
    if 2 * k - 1 <= NL:
        stir[2 * k - 1] = bk[k]
G0 = expser([None] + [[stir[n]] if stir[n] else [] for n in range(1, NL + 1)])
okP = True
for m in range(-3, 5):
    poch = [Fr(1)]
    for i in (range(m) if m > 0 else []):
        poch = pmul(poch, [Fr(1), Fr(i)])
    pinv = [Fr(1)] + [Fr(0)] * NL
    for i in (range(1, -m + 1) if m < 0 else []):
        pinv = [sum(pinv[r] * Fr(i) ** (n - r) for r in range(n + 1)) for n in range(NL + 1)]   # times 1/(1 - i z)
    lhs = [peval(Nw[n], Fr(m)) for n in range(NL + 1)]
    pp = pmul(poch, pinv) if m < 0 else poch
    rhs = [sum((pp[r] if r < len(pp) else 0) * (G0[n - r][0] if G0[n - r] else 0) for r in range(n + 1)) for n in range(NL + 1)]
    okP &= lhs == rhs
check('L2 integer M = -3..4: the polynomial form equals Poch_M(z) exp(Sum_k b_k z^(2k-1)) (the (L3) of CP2 s.2.8)', okP)
check('L3 control: amplitude t^M dt in place of t^(M-1) dt (i.e. M -> M + 1) equals exp(Sum beta_n(M) z^n)', wick(1) == Gb, control=True)
Gs = expser([None] + [pnorm([peval(beta(n), Fr(1, 5))]) for n in range(1, NL + 1)])
check('L3 control: at M = 0 the series equals exp(Sum beta_n(1/5) z^n)', [pnorm([peval(Nw[n], Fr(0))]) for n in range(NL + 1)]
      == Gs, control=True)
say('rc_l3c %d %d: %s   (%.1f s)' % (NW, NL, 'ALL CHECKS PASSED' if NFAIL[0] == 0 else '%d FAILURES' % NFAIL[0], time.time() - T0))
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
with open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8') as fh:
    fh.write('# %s (py -I -B rc_l3c.py %d %d)\n' % (LOGREL, NW, NL))
    fh.write('\n'.join(OUT[:-1]) + '\n' + OUT[-1].split('   (')[0] + '\n')
print('log written: ' + LOGREL)
