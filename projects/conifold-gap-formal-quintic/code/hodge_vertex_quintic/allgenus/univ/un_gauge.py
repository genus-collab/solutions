# -*- coding: utf-8 -*-
"""un_gauge.py -- Theorem-U gauges for FICTITIOUS self-dual cyclic systems with one simple pole (E2), and the pure-1/z solutions.

System:  th Y_a = b_a(u) Y_{a+1} (a in Z/5, th = z d/du),   b = (f, B/(f^2 g^2), f, g, g),   B = prod b_a = -1/(u (1-u)^5),
         f = t_c'(u) F(u),   g = -(varpi_0(u)/(1-u)) Gm(u)      (F = Gm = 1: the formal quintic in the conifold frame (alpha, beta)).
The shape (f, ., f, g, g) is self-duality (the pairing Sum_a Y_a(z) eta^{a,3-a} Y_{3-a}(-z) is u-constant iff b_{a-1} = b_{3-a});
B is kept, so the model variable v(u) = u tau^(-5/4) and tau = 1/2F1(1, 4/5; 9/5; u) are those of the quintic.
The gauge is computed by natural_form/theory/ur_gauge.gauge (imported; its module-level function `coefficients` is replaced IN MEMORY
by the perturbed one for the duration of the call -- no file of the project is touched).

API: system(NU, F, Gm, ab) -> (bq, bt, tau, v, C, f);  pert_gauge(K, NU, F, Gm, ab) -> res (ur_gauge format, + 'f', 'T');
     rows_from_gauge(res, Rh, R, KZ, K) -> Qh (plain u-lists);  yc_model(v, nmax, NU), check_G(res, kmax, NU) (hypothesis (G)).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from fq_fp import FpRing, QRing, P2
from fq_series import Ser
from fq_conifold import Conifold
import ur_gauge
from ur_gauge import CEXP, mu
Q = QRing()


def system(NU, F=(1,), Gm=(1,), ab=(0, 0)):
    C = Conifold(FpRing(P2), Nu=NU, alpha=Fr(ab[0]), beta=Fr(ab[1])); N = NU
    var0 = C.varpi0.truncate(N + 2); tc = C.tc.truncate(N + 2)
    one_m = Ser.from_list(Q, [1, -1], 0, N + 2)
    Fs = Ser.from_list(Q, [Fr(x) for x in F], 0, N + 2); Gs = Ser.from_list(Q, [Fr(x) for x in Gm], 0, N + 2)
    f = tc.deriv().mul(Fs); g = var0.mul(one_m.inv()).neg().mul(Gs)
    B = Ser.mono(Q, 1, 1, N + 2).mul(one_m.pow(5)).inv().neg()
    b1 = B.mul(f.mul(f).mul(g).mul(g).inv())
    bq = [f, b1, f, g, g]
    tau = Ser.from_list(Q, [Fr(4, 5 * n + 4) for n in range(N + 2)], 0, N + 2).inv()
    v = tau.log().scale(Fr(-5, 4)).exp().shift(1); vp = v.deriv()
    bt = [vp, vp.mul(v.inv()).neg(), vp, vp.neg(), vp.neg()]
    return bq, bt, tau, v, C, f


def pert_gauge(K, NU, F=(1,), Gm=(1,), ab=(0, 0), log=None):
    sysd = system(NU, F, Gm, ab); orig = ur_gauge.coefficients
    ur_gauge.coefficients = lambda NU_: sysd[:5]
    try: res = ur_gauge.gauge(K, NU, log)
    finally: ur_gauge.coefficients = orig
    f = sysd[5]; res['f'] = f
    res['T'] = Ser(Q, 1, [f.coef(n) / (n + 1) for n in range(f.prec)])        # T = Int_0^u f
    return res


def rows_from_gauge(res, Rh, R, KZ, K):
    """Qh'[a][i][n] = Sum_{k <= n, b = a - k} Ghat_(k,b)(u) tau(u)^(n-k) Rh[b][i][n-k]   (plain lists of K residues)"""
    p = R.p; hat = ur_gauge.hatted(res); tau = res['tau']
    def red(s):
        assert s.prec >= K, "gauge precision O(u^%d) < %d" % (s.prec, K)
        assert (not s.c) or s.v >= 0, "pole in a hatted gauge entry"
        return [R.of(s.coef(n)) for n in range(K)]
    def mul(a, b):
        out = [0] * K
        for i in range(K):
            if a[i]:
                for j in range(K - i): out[i + j] += a[i] * b[j]
        return [x % p for x in out]
    tp = [[1] + [0] * (K - 1)]; t1 = red(tau)
    for m in range(1, KZ + 1): tp.append(mul(tp[-1], t1))
    H = {key: red(s) for key, s in hat.items() if key[0] <= KZ}
    Qh = [[[[0] * K for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (k, b), h in H.items():
        if not any(h): continue
        a = (b + k) % 5
        for n in range(k, KZ + 1):
            ht = mul(h, tp[n - k])
            for i in range(5):
                c = Rh[b][i][n - k]
                if c:
                    q = Qh[a][i][n]
                    for j in range(K): q[j] = (q[j] + ht[j] * c) % p
    return Qh


def yc_model(v, nmax, NU):
    """the pure-1/z solution of the model, Y^C = (I_0, -I_-1, -I_3, I_2, I_1)(v(u)), I_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!:
    {n: (component, Ser in u)} for the coefficient of z^-n, n = 0..nmax  (one component per n: a = (1 - n) mod 5)"""
    import math
    comp = {0: (0, 1), 1: (-1, -1), 2: (3, -1), 3: (2, 1), 4: (1, 1)}       # a -> (m, sign)
    out = {}
    vpow = [Ser.const(Q, 1, NU)]
    for e in range(1, nmax + 2): vpow.append(vpow[-1].mul(v).truncate(NU))
    for a, (m, sg) in comp.items():
        i = 0
        while True:
            e = m + 4 * i + 1; n = e + i                      # (-v)^e z^-(e+i) (-1/4)^i / (i! e!)
            if n > nmax: break
            if e >= 0:
                c = Fr((-1) ** e) * Fr(-1, 4) ** i / (math.factorial(i) * math.factorial(e)) * sg
                assert n not in out and (1 - n) % 5 == a
                out[n] = (a, vpow[e].scale(c))
            i += 1
    return out


def yc_system(bq, nmax, NU):
    """the pure-1/z solution of th Y_a = b_a Y_{a+1} with y_0 = -e_1 and zero integration constants: {n: (component, Ser)}"""
    out = {0: (1, Ser.const(Q, -1, NU))}
    for n in range(1, nmax + 1):
        a = (1 - n) % 5; integrand = bq[a].mul(out[n - 1][1])
        assert (not integrand.c) or integrand.v >= 0, "logarithm in y_%d" % n
        out[n] = (a, Ser(Q, 1, [integrand.coef(k) / (k + 1) for k in range(min(integrand.prec, NU))]))
    return out


def check_G(res, kmax, J):
    """hypothesis (G): P Y^C_model(v(u)) = c(z) Y^C_system with c(z) = P_11(0, z), compared mod u^J on every coefficient z^m with
    -nJ <= m <= kmax - nJ, nJ = max{n : ord_u y_n = n - floor(n/5) < J} (the range in which the truncation k <= kmax is exact).
    Returns (ok, c, bad, number of z-powers compared, (m_min, m_max))."""
    nJ = max(n for n in range(0, 5 * J + 5) if n - n // 5 < J); NU = J + 2
    p = res['p']; ym = yc_model(res['v'], nJ, NU); ys = yc_system(res['bq'], nJ, NU)
    c = {k: p[k][1].coef(0) for k in range(0, kmax + 1) if k % 5 == 0}      # P_11(0, z) = Sum_k p_{k,1}(0) z^k (graded: k = 0 mod 5)
    bad = []; nchk = 0
    for m in range(-nJ, kmax - nJ + 1):
        acc = Ser.zero(Q, NU)
        for n in range(max(0, -m), nJ + 1):                    # coefficient of z^m of P Y^C: Sum_{k - n = m} p_{k, a(n)} y_n
            a, y = ym[n]; acc = acc.add(p[m + n][a].mul(y))
        ref = Ser.zero(Q, NU)
        for k5, ck in c.items():                               # c(z) Y^C_system at z^m: Sum_{k5 - n = m} c_{k5} ys_n
            n = k5 - m
            if 0 <= n <= nJ and ck != 0: ref = ref.add(ys[n][1].scale(ck))
        d = acc.sub(ref); nchk += 1
        assert d.prec >= J, "precision O(u^%d) < %d in check_G" % (d.prec, J)
        if any(d.coef(j) != 0 for j in range(min(d.v, 0), J)): bad.append((m, d.fmt(d.v + 2)))
    return (not bad), c, bad, nchk, (-nJ, kmax - nJ)
