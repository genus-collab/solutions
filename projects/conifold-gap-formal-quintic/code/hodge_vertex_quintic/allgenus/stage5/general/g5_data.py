"""g5_data.py -- local data for g5_core.engine built from Gamma factors, and the predicted values.

A regular factor Gamma_St(e + (c sigma + d)/z) has value w = c s* + d at the colour point.  Factors with the same (c, e)
form a GROUP, stored through its power sums ps[b] = Sum_(f in group) w_f^-b (b = 0..BMAX).  With x = sigma - s*:
  Phi''(s*+x) = Sum_g c^2 Sum_m (-c)^m ps[m+1] x^m                       (regular phase, (c s+d)log(c s+d) - (c s+d))
  A_0'(s*+x)  = Sum_g (e - 1/2) c Sum_m (-c)^m ps[m+1] x^m                 (from (Y + e - 1/2) log Y)
  A_n(s*+x)   = Sum_g beta_n(e) Sum_m C(-n, m) c^m ps[n+m] x^m,  beta_n(h) = (-1)^(n+1) B_(n+1)(h)/(n(n+1)).
Prediction (chart argument / first pole of the singular factor): the regular factors evaluated at the first pole
sigma = s* - e* z/c* of Gamma(e* + c*(sigma - s*)/z), i.e. shifts h_g = e_g - c_g e*/c*:
  log c(z; s*) = Sum_n z^n Sum_g beta_n(h_g) ps[n].       (For e* = 0 this is Sum_n z^n A_n(s*) = log h_w.)
Group options: 'beta' = function n -> beta_n used in the integrand (default beta_n(e)); 'drop' = True drops its Stirling
corrections from the integrand (controls).  The prediction always uses the unperturbed beta_n(h_g).
"""
from fractions import Fraction as Fr
from math import comb
import g5_core as C


def bernoulli(n):
    B = [Fr(0)] * (n + 1); B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BER = bernoulli(260)


def bpoly(n, h):
    h = Fr(h)
    return sum(comb(n, k) * BER[k] * h ** (n - k) for k in range(n + 1))


def beta(h, n):
    return Fr((-1) ** (n + 1)) * bpoly(n + 1, h) / (n * (n + 1))


def group_rational(c, e, ws, BMAX, **opt):
    """group of factors with the same (c, e) and rational values ws"""
    ps = [sum(Fr(1) / Fr(w) ** b for w in ws) for b in range(BMAX + 1)]
    g = dict(c=Fr(c), e=Fr(e), ps=ps); g.update(opt)
    return g


def newton_psums(a, BMAX):
    """power sums p_b, b = 0..BMAX, of the roots of the monic x^d + a[0] x^(d-1) + ... + a[d-1]"""
    d = len(a); p = [Fr(d)]
    for b in range(1, BMAX + 1):
        s = -sum(a[i - 1] * p[b - i] for i in range(1, min(b - 1, d) + 1))
        if b <= d:
            s -= b * a[b - 1]
        p.append(s)
    return p


def build_loc(groups, sing, N, K):
    """sing = dict(M, cs, es, bstar) with bstar[n] = beta*_n (n = 0..N)."""
    rho = [0] * (K + 1); a0p = [0] * (K + 1)
    for g in groups:
        c, e, ps = g['c'], g['e'], g['ps']
        for m in range(K + 1):
            rho[m] += c * c * (-c) ** m * ps[m + 1]
            a0p[m] += (e - Fr(1, 2)) * c * (-c) ** m * ps[m + 1]
    Phi = [Fr(0), Fr(0)] + [rho[m] / ((m + 1) * (m + 2)) for m in range(K)]
    H0 = [C.cv(0)] + [C.cv(a0p[m] / (m + 1)) for m in range(K)]
    An = [None]
    for n in range(1, N + 1):
        row = [Fr(0)] * (K + 1)
        for g in groups:
            if g.get('drop'):
                continue
            bn = g['beta'](n) if 'beta' in g else beta(g['e'], n)
            if bn == 0:
                continue
            c, ps = g['c'], g['ps']
            for m in range(K + 1):
                row[m] += bn * gb(-n, m) * c ** m * ps[n + m]
        An.append([C.cv(x) for x in row])
    Bk = [C.sexp(H0, K)]
    for n in range(1, N + 1):
        acc = [0] * (K + 1)
        for i in range(1, n + 1):
            pr = C.smul(An[i], Bk[n - i], K)
            for t in range(K + 1):
                acc[t] += i * pr[t]
        Bk.append([C.red(x * C.finv(C.cv(n))) for x in acc])
    loc = dict(sing); loc['Phi'] = Phi; loc['Bk'] = Bk
    return loc


def gb(n, m):
    """C(n, m) for integer n (possibly negative)"""
    r = Fr(1)
    for i in range(m):
        r = r * (n - i) / (i + 1)
    return r


def predicted(groups, sing, N, shift=True):
    """[z^n] of the predicted value at eps = 0, exact over Q (list n = 0..N)"""
    cs, es = Fr(sing['cs']), Fr(sing['es'])
    lg = [Fr(0)] * (N + 1)
    for g in groups:
        h = g['e'] - g['c'] * es / cs if shift else g['e']
        for n in range(1, N + 1):
            lg[n] += beta(h, n) * g['ps'][n]
    E = [Fr(0)] * (N + 1); E[0] = Fr(1)
    for m in range(1, N + 1):
        E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
    return E


def sing_base(N, bfun=None):
    """singular factor Gamma_St((lambda_i - sigma)/z): M = 1, c* = -1, e* = 0"""
    bf = bfun if bfun else (lambda n: beta(0, n))
    return dict(M=1, cs=-1, es=0, bstar=[0] + [bf(n) for n in range(1, N + 1)])


def quintic_groups(BMAX, colour_point='q0', **opt):
    """formal quintic. 'q0': colour 0 at s* = 1 (singular factor j = 0); regular: Gamma(1+5s) (w = 5) and
    j = 1..4 (w = zeta^j - 1, power sums by Newton from v^4 + 2v^3 + 2v^2 + v + 1/5, v = 1/(zeta^j - 1)).
    'qinf': s* = 0, singular factor Gamma(1+5s); regular: j = 0..4 (w = zeta^j, ps[b] = 5 [5 | b])."""
    if colour_point == 'q0':
        tw = group_rational(5, 1, [5], BMAX); tw.update(opt.get('tw', {}))
        pb = newton_psums([Fr(2), Fr(2), Fr(1), Fr(1, 5)], BMAX)
        base = dict(c=Fr(-1), e=Fr(0), ps=pb); base.update(opt.get('base', {}))
        return [tw, base]
    base = dict(c=Fr(-1), e=Fr(0), ps=[Fr(5) if b % 5 == 0 else Fr(0) for b in range(BMAX + 1)])
    base.update(opt.get('base', {}))
    return [base]
