# -*- coding: utf-8 -*-
"""pb_res.py NU -- Stage-5 'pb' audit: the residue identity V = C0 (Y^C_Q)_0 (POLY s.6 step 6 with kappa = C0), exact over Q.

V := Res_(x=0) e^(-W/z) at dx = sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)] (e^(-G/z) at),   Lam = -log(1-u),
G(x) = sum_m x^(5m-1)/(m(5m-1))   (from S(1/x) = -W, re-derived in PB_AUDIT s.2),
log at = -(1/2) log(1-x^5) + sum_k b_k (xz)^(2k-1) [5^(1-2k) - 5 sum_P C(2k-2+5P,5P) x^(5P)].
Y^C_Q: GI Lemma Q, y^(0) = -e_1, y^(n+1)_a = int_0^u b_a y^(n)_(a+1); (Y^C_Q)_0 = sum_n z^-n y^(n)_0; C0 = varpi0 of the frame.
Checks every coefficient u^1..u^NU, every power of z, frames (0,0) and (1/3,-2/7).  Controls: no Stirling terms (b_k = 0);
a wrong vanishing period (gamma -> 1/2).  Also: [z^-1] V = -g (the vanishing period, frame free).
Log: logs/pb_res_<NU>.log (relative paths only).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C
from pb_core import LS
from fractions import Fraction as Fr
from math import comb

LOG = []


def log(s=''):
    print(s); LOG.append(s)


def lp_add(a, b, k=1):
    for e, c in b.items():
        a[e] = a.get(e, 0) + k * c
    return a


def lp_mul(a, b):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            r[e1 + e2] = r.get(e1 + e2, 0) + c1 * c2
    return r


def xexp(H, nx):
    """exp of an x-series H (H[0] = 0) with Laurent-polynomial coefficients in z: E_j, j < nx."""
    E = [{0: Fr(1)}]
    for j in range(1, nx):
        acc = {}
        for i in range(1, j + 1):
            if H[i]:
                lp_add(acc, lp_mul(H[i], E[j - i]), i)
        E.append({e: c / j for e, c in acc.items() if c})
    return E


def residue_V(NU, bern, stirling=True):
    nx = NU
    Gx = [Fr(0)] * nx
    for m in range(1, nx):
        if 5 * m - 1 < nx:
            Gx[5 * m - 1] = Fr(1, m * (5 * m - 1))
    H1 = [({-1: -Gx[j]} if Gx[j] else {}) for j in range(nx)]
    E1 = xexp(H1, nx)
    H2 = [{} for _ in range(nx)]
    for m in range(1, nx):
        if 5 * m < nx:
            lp_add(H2[5 * m], {0: Fr(1, 2 * m)})
    if stirling:
        for k in range(1, nx):
            r = 2 * k - 1
            if r >= nx:
                break
            bk = bern[2 * k] / (2 * k * (2 * k - 1))
            for P in range(0, nx):
                j = r + 5 * P
                if j >= nx:
                    break
                c = -5 * comb(2 * k - 2 + 5 * P, 5 * P) + (Fr(1, 5 ** r) if P == 0 else 0)
                lp_add(H2[j], {r: bk * c})
    E2 = xexp(H2, nx)
    prod = [{} for _ in range(nx)]
    for i in range(nx):
        for j in range(nx - i):
            lp_add(prod[i + j], lp_mul(E1[i], E2[j]))
    Lam = LS.from_list([Fr(0)] + [Fr(1, n) for n in range(1, NU + 1)], 0, NU + 1)
    V = {}
    Lk = LS.const(1, NU + 1); fact = 1
    for k in range(1, NU + 1):
        Lk = Lk.mul(Lam); fact *= k
        A = prod[k - 1]
        for n in range(NU + 1):
            cu = Lk.coef(n)
            if cu:
                for e, c in A.items():
                    key = (n, e - k)
                    V[key] = V.get(key, 0) + (-1) ** k * cu * c / fact
    return {k: v for k, v in V.items() if v}


def vanishing_row0(K, NU):
    """C0 (Y^C_Q)_0 as {(u-power, z-power): coeff}, u-powers <= NU."""
    N = NU + 8
    one_m = LS.from_list([1, -1], 0, N + 2)
    w0 = K.varpi0.trunc(N + 2); tcp = K.tc.trunc(N + 3).deriv()
    b1 = w0.mul(w0).mul(tcp).mul(tcp).mul(one_m).mul(one_m).mul(one_m).shift(1).inv().neg()
    b3 = w0.mul(one_m.inv()).neg()
    b = [tcp, b1, tcp, b3, b3]
    y = [LS.const(0, N), LS.const(-1, N), LS.const(0, N), LS.const(0, N), LS.const(0, N)]
    out = {}
    nmax = (5 * NU - 1) // 4 + 1
    C0 = K.varpi0.trunc(NU + 1)
    for n in range(0, nmax + 1):
        row0 = C0.mul(y[0])
        for e in range(NU + 1):
            if e < row0.prec:
                c = row0.coef(e)
                if c:
                    out[(e, -n)] = out.get((e, -n), 0) + c
            else:
                raise ValueError('precision')
        y = [b[a].mul(y[(a + 1) % 5]).integ() for a in range(5)]
    return out


def first_diff(A, B, NU):
    keys = sorted(set(A) | set(B))
    bad = [k for k in keys if 1 <= k[0] <= NU and A.get(k, 0) != B.get(k, 0)]
    return (min(bad) if bad else None), len([k for k in keys if 1 <= k[0] <= NU])


def main():
    NU = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    C.set_field(0)
    bern = C.bernoulli(2 * NU + 4)
    t0 = time.time()
    log('pb_res.py %d   (exact over Q; every coefficient u^1..u^%d, all powers of z)' % (NU, NU))
    V = residue_V(NU, bern)
    log('  V built: %d non-zero coefficients, %.1f s; positive or zero z-powers present: %s'
        % (len(V), time.time() - t0, any(e >= 0 for (_, e) in V)))
    log('  V at u^1: %s;  [u^2 z^-1] V = %s (phase -1/2 plus Stirling b_1 = 1/12 gives -1/5: -7/10 = -gamma)'
        % ({e: str(c) for (n, e), c in V.items() if n == 1}, V.get((2, -1))))
    K = C.Conifold(NU + 16)
    ok = True
    for fr in ((0, 0), (Fr(1, 3), Fr(-2, 7))):
        K.set_frame(*fr)
        W = vanishing_row0(K, NU)
        d, nk = first_diff(V, W, NU)
        log('  frame (%s,%s): V = C0 (Y^C_Q)_0 on all %d keys (u^1..u^%d): %s'
            % (fr[0], fr[1], nk, NU, 'PASS' if d is None else 'FAIL first at (u,z) = %s' % (d,)))
        ok = ok and d is None
        g = K.g
        zm1 = all(V.get((n, -1), 0) == -g.coef(n) for n in range(1, NU + 1))
        log('    [z^-1] V = -g (vanishing period u + (7/10)u^2 + ...) for u^1..u^%d: %s' % (NU, 'PASS' if zm1 else 'FAIL'))
        ok = ok and zm1
    # controls
    K.set_frame(0, 0)
    W = vanishing_row0(K, NU)
    V0 = residue_V(NU, bern, stirling=False)
    d, _ = first_diff(V0, W, NU)
    log('  control [no Stirling terms in at]: first mismatch %s (predicted u^2, z^-1): %s'
        % (d, 'PASS (fails as it must)' if d == (2, -1) else 'CONTROL DID NOT BEHAVE'))
    ok = ok and d == (2, -1)
    gwrong = K.f[1].add(K.f[2].scale(Fr(1, 2)))
    K.g = gwrong; K.set_frame(0, 0)
    W2 = vanishing_row0(K, NU)
    d, _ = first_diff(V, W2, NU)
    log('  control [vanishing period with gamma = 1/2]: first mismatch %s (predicted u^2, z^-1): %s'
        % (d, 'PASS (fails as it must)' if d == (2, -1) else 'CONTROL DID NOT BEHAVE'))
    ok = ok and d == (2, -1)
    log('total %.1f s' % (time.time() - t0))
    log('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    name = 'pb_res_%d.log' % NU
    with open(os.path.join(HERE, 'logs', name), 'w') as f:
        f.write('\n'.join(LOG) + '\n')
    print('log written to logs/' + name)


if __name__ == '__main__':
    main()
