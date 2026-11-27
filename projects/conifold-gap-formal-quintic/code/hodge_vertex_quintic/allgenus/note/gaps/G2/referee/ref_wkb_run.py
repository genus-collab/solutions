# -*- coding: utf-8 -*-
"""ref_wkb_run.py -- driver of ref_wkb.py (see its docstring).  Usage: py -I -B ref_wkb_run.py [nz] [nq]."""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_wkb import (check, FAIL, zc, zpow, zadd, zneg, zmul, zinv, zscal, sadd, smul, sD, L_series, pf5)

BK = {1: Fr(1, 12), 2: Fr(-1, 360), 3: Fr(1, 1260)}          # B_2k/(2k(2k-1))


def h_series(w, nz):
    """h_w(z) = exp(-Sum_k bk_k z^(2k-1) Sum_(a in w) a^(1-2k)) through z^nz, over Q(zeta)"""
    lg = [zc(0)] * (nz + 1)
    for k, bk in BK.items():
        if 2 * k - 1 <= nz:
            s = zc(0)
            for a in w:
                ia = zinv(a); p = ia
                for _ in range(2 * k - 2):
                    p = zmul(p, ia)
                s = zadd(s, p)
            lg[2 * k - 1] = zscal(-bk, s)
    out = [zc(1)] + [zc(0)] * nz                              # exp via n e_n = Sum_k k l_k e_(n-k)
    for n in range(1, nz + 1):
        acc = zc(0)
        for k in range(1, n + 1):
            acc = zadd(acc, zscal(k, zmul(lg[k], out[n - k])))
        out[n] = zscal(Fr(1, n), acc)
    return out


def weights(i):
    w = [zadd(zpow(i), zneg(zpow(j))) for j in range(5) if j != i]
    return w + [zscal(-5, zpow(i))]


def wkb(i, consts, nz, nq, pert=False):
    """amplitude psi_i = Sum_n z^n L g_n with g_n(0) = consts[n]; returns list of q-series and the no-log record"""
    L = L_series(nq)
    lam = [zmul(zpow(i), a) for a in L]
    piv = zscal(5, zpow(4 * i))
    psi = []; nolog = True
    for n in range(nz + 1):
        F = psi + [[zc(0)] * nq] * (2)                         # z^0..z^(n+1), psi_n = 0 for the residual
        F = F[:n + 2]
        r = pf5(F, lam, pert)[n + 1]
        rhs = [zneg(a) for a in r]
        nolog = nolog and not any(rhs[0])
        g = [consts[n]] + [zmul(rhs[k], zinv(zscal(k, piv))) for k in range(1, nq)]
        psi.append(smul(L, g))
    return psi, nolog


def transport_ok(i, nq):
    """[z^1] of PF5^(u_i)[L g] equals 5 zeta^(4i) D g, tested on g = 1 and g = q"""
    L = L_series(nq); lam = [zmul(zpow(i), a) for a in L]
    ok = True
    for g in ([zc(1)] + [zc(0)] * (nq - 1), [zc(0), zc(1)] + [zc(0)] * (nq - 2)):
        r = pf5([smul(L, g), [zc(0)] * nq], lam)
        want = [zscal(5, zmul(zpow(4 * i), a)) for a in sD(g)]
        ok = ok and not any(any(a) for a in r[0]) and all(x == y for x, y in zip(r[1][:nq - 1], want[:nq - 1]))
    return ok


def covariant(psi_i, psi_0, i):
    return all(psi_i[n][k] == zmul(zpow(-i * n), psi_0[n][k]) for n in range(len(psi_0)) for k in range(len(psi_0[n]) - 1))


def main():
    nz = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nq = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    print("ref_wkb_run.py: WKB amplitudes of PF5 at large radius, z^0..z^%d, q^0..q^%d, exact over Q(zeta)" % (nz, nq - 2))
    h = [h_series(weights(i), nz) for i in range(5)]
    check("h_(w_0) = 1 - (3/20) z + (9/800) z^2 + ... is rational", all(all(c == 0 for c in h[0][n][1:]) for n in range(nz + 1)),
          "[%s]" % ", ".join(str(h[0][n][0]) for n in range(nz + 1)))
    check("lem:set-hrot: h_(w_i)(z) = h_(w_0)(zeta^-i z), i = 1..4 (from the weights)",
          all(h[i][n] == zmul(zpow(-i * n), h[0][n]) for i in range(1, 5) for n in range(nz + 1)))
    check("transport operator T(L g) = 5 zeta^(4i) D g and symbol zero, i = 0..4", all(transport_ok(i, nq) for i in range(5)))
    psi0, nl0 = wkb(0, h[0], nz, nq)
    L = L_series(nq)
    check("colour 0: no logarithm at any order; psi_(0,0) = L", nl0 and psi0[0][:nq - 1] == L[:nq - 1])
    lam0 = L
    res = pf5(psi0 + [[zc(0)] * nq], lam0)
    check("colour 0: PF5^(u_0)[psi_0] = 0 at z^0..z^%d (direct substitution)" % nz,
          all(not any(res[m][k]) for m in range(nz + 1) for k in range(nq - 1)))
    for i in range(1, 5):
        psii, nli = wkb(i, h[i], nz, nq)
        check("colour %d: no logarithm; psi_%d(q, z) = psi_0(q, zeta^-%d z) (all z^n, q^k in range)" % (i, i, i),
              nli and covariant(psii, psi0, i))
        bad, _ = wkb(i, h[0], nz, nq)
        check("colour %d: control: normalised by h_(w_0) instead of h_(w_%d), the covariance fails" % (i, i), not covariant(bad, psi0, i))
    p0, _ = wkb(0, h[0], nz, nq, pert=True)
    p1, _ = wkb(1, h[1], nz, nq, pert=True)
    check("control: PF5 + q z^2 (not invariant under z -> zeta z): the covariance fails", not covariant(p1, p0, 1))
    print("  psi_(0,1) = %s + %s q + %s q^2 + ..." % tuple(str(psi0[1][k][0]) for k in range(3)))
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
