"""g5_theory.py -- Gamma-factor descriptions of twisted equivariant projective spaces, the Mellin-Barnes summand, the
check of the summand against the twisted I-function, and one generic test run (removable singularity + value).

A theory = list of factors {c, d, e, kind}: the summand is  prod_f Gamma(e_f + c_f s + d_f/z) * Q^s,  sigma = z s.
  base P^(N-1), weight lambda_j:          Gamma(lambda_j/z - s)          c = -1, d = lambda_j, e = 0
  convex O(l), c_1 = l H + mu (Euler):    Gamma(1 + l s + mu/z)          c = l,  d = mu,       e = 1
  concave O(-l), c_1 = -l H + mu (e^-1):  Gamma(l s - mu/z)              c = l,  d = -mu,      e = 0
Q = (-1)^(N + l_conc) z^(Sum_f c_f) q, l_conc = Sum of the concave degrees l  (derived in GENERAL.md s.2 and checked by residue_check below).
Colour i: s* = lambda_i, singular factor = base factor i.  Prediction: h_(w_i), w_i = {-(c_f lambda_i + d_f)}_(f != i).
"""
from fractions import Fraction as Fr
import g5_core as C, g5_data as D


def base(lams):
    return [dict(c=Fr(-1), d=Fr(l), e=Fr(0), kind='base') for l in lams]


def convex(l, mu=0):
    return dict(c=Fr(l), d=Fr(mu), e=Fr(1), kind='convex O(%s)' % l)


def concave(l, mu=0):
    return dict(c=Fr(l), d=Fr(-mu), e=Fr(0), kind='concave O(-%s)' % l)


def groups_for(factors, star, BMAX, opt=None):
    """groups of the regular factors at the colour point of the singular factor 'star' (index)"""
    fs = factors[star]; sstar = -fs['d'] / fs['c']
    bykey = {}
    for k, f in enumerate(factors):
        if k == star:
            continue
        w = f['c'] * sstar + f['d']
        if w == 0:
            raise ValueError('degenerate: a regular factor vanishes at the colour point (factor %d)' % k)
        bykey.setdefault((f['c'], f['e']), []).append(w)
    gr = []
    for (c, e), ws in sorted(bykey.items()):
        g = D.group_rational(c, e, ws, BMAX)
        if opt and (c, e) in opt:
            g.update(opt[(c, e)])
        gr.append(g)
    return gr, sstar


def weights(factors, star):
    """the predicted QRR weight set w_i = {-(c_f s* + d_f)}"""
    fs = factors[star]; sstar = -fs['d'] / fs['c']
    return [-(f['c'] * sstar + f['d']) for k, f in enumerate(factors) if k != star]


def run_case(factors, star, N, K=None, opt=None, bfun=None, pred_weights=None):
    """returns (first pole or None, list of n with value mismatch, values, prediction)"""
    K = N if K is None else K
    gr, sstar = groups_for(factors, star, N + K + 1, opt)
    fs = factors[star]
    sg = dict(M=1, cs=fs['c'], es=fs['e'],
              bstar=[0] + [(bfun(n) if bfun else D.beta(fs['e'], n)) for n in range(1, N + 1)])
    ct = C.engine(D.build_loc(gr, sg, N, K), N, K)
    first, vals, poles = C.laurent_report(ct, N, K)
    if pred_weights is None:
        pg, _ = groups_for(factors, star, N + K + 1)
        pred = D.predicted(pg, sg, N)
    else:
        lg = [Fr(0)] * (N + 1)
        for n in range(1, N + 1):
            lg[n] = -D.beta(0, n) * sum(Fr(a) ** (-n) for a in pred_weights) if n % 2 else Fr(0)
        pred = [Fr(1)] + [Fr(0)] * N
        for m in range(1, N + 1):
            pred[m] = sum(i * lg[i] * pred[m - i] for i in range(1, m + 1)) / m
    bad = [n for n in range(min(N, K) + 1) if not C.iszero(vals[n] - C.cv(pred[n]))]
    return first, bad, vals, pred


def poch(x, n):
    r = Fr(1)
    for k in range(n):
        r *= (x + k)
    return r


def residue_check(factors, i, z, dmax, imod):
    """compare residues of the summand at s = lambda_i/z + d (ratio to d = 0) with the I-function coefficient
    imod(d) at H = lambda_i (imod = product of the twist modification factors), exact, rational z.
    Returns the list of (d, ratio summand/I) -- all equal to 1 iff the summand and Q are right."""
    lam = [f['d'] for f in factors if f['kind'] == 'base']
    Nb = len(lam); csum = sum(f['c'] for f in factors)
    lconc = sum(int(f['c']) for f in factors if f['kind'].startswith('concave'))
    out = []
    for d in range(1, dmax + 1):
        res = Fr((-1) ** d) / poch(Fr(1), d)                  # Gamma(lambda_i/z - s) at its d-th pole, over d = 0
        for k, f in enumerate(factors):
            if f['kind'] == 'base' and f['d'] == lam[i] and k == i:
                continue
            x = f['e'] + f['c'] * lam[i] / z + f['d'] / z        # argument at d = 0
            n = int(f['c'] * d)
            res *= poch(x, n) if n >= 0 else 1 / poch(x + n, -n)
        res *= Fr((-1) ** (Nb + lconc)) ** d * z ** (csum * d)   # Q^d / q^d, Q = (-1)^(N + l_conc) z^(Sum c) q
        ival = imod(d)
        for j in range(Nb):
            for k in range(1, d + 1):
                ival /= (lam[i] - lam[j] + k * z)
        out.append((d, res / ival))
    return out
