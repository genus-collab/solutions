# T1: singular row at generic CY weights against the closed form; constants; leaf structure.
# Evidence only. Run: py -I t1_rows.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cylib import *
from fractions import Fraction as Fr

SAMPLES = [  # (a0,a1,a2,e1,e2); first ones CY
    (Fr(1), Fr(2), Fr(4), Fr(3), Fr(4)),
    (Fr(3, 5), Fr(-7, 3), Fr(11, 2), Fr(5, 7), None),
    (Fr(-2), Fr(9), Fr(1, 3), Fr(13, 4), None),
    (Fr(0), Fr(5, 3), Fr(2, 7), Fr(5, 3), Fr(2, 7)),   # fibrewise torus, control of the code
    (Fr(0), Fr(5, 3), Fr(2, 7), Fr(1, 3), None),       # a0 = 0, a_i != e_i
]


def mk(sample):
    a0, a1, a2, e1, e2 = sample
    if e2 is None:
        e2 = a0 + a1 + a2 - e1
    return weights(a0, a1, a2, e1, e2)


def zmul(a, b, N):
    r = [0] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x:
            for j, y in enumerate(b[:N + 1 - i]):
                r[i + j] += x * y
    return [v % P for v in r]


def check_identities(w, R):
    """Delta_1 lead, rho, L_1 pole, L sum/product, discriminant"""
    L = R.L
    E = w.Ea
    ok = []
    # Delta_1 t^2 at t=0 = - e1 e2 E^3
    ok.append(('Delta1 lead = -e1 e2 E^3', lcoef(R.D1, -2) == (-w.e1 * w.e2 % P * pow(E, 3, P)) % P))
    ok.append(('t L_1(0) = E', lcoef(R.L1, -1) == E))
    # L_1+L_2 = E/tau - a0, L_1 L_2 = a1 a2/tau  (compare Laurent)
    itau = (-1, R.itau_reg)
    S = ladd(R.L1, R.L2, L)
    Sx = ladd(lscal(E, itau), (0, const(P - w.a0, L)), L)
    ok.append(('L1+L2 = E/tau - a0', all(lcoef(S, k) == lcoef(Sx, k) for k in range(-1, L - 3))))
    Pr = lmul(R.L1, R.L2, L)
    Px = lscal(w.a1 * w.a2 % P, itau)
    ok.append(('L1 L2 = a1 a2/tau', all(lcoef(Pr, k) == lcoef(Px, k) for k in range(-1, L - 3))))
    # (13.2): Delta_i = e1 e2 (a0 a1 a2 + 2 a1 a2 L_i - sig L_i^2)
    for nm, Li, Di in (('1', R.L1, R.D1), ('2', R.L2, R.D2)):
        t = ladd(lscal(2 * w.a1 * w.a2 % P, Li), (0, const(w.a0 * w.a1 % P * w.a2 % P, L)), L)
        t = ladd(t, lscal((-E) % P, lmul(Li, Li, L)), L)
        t = lscal(w.e1 * w.e2 % P, t)
        ok.append(('Delta_%s two ways' % nm, all(lcoef(t, k) == lcoef(Di, k) for k in range(-2, L - 4))))
    ok.append(('L_2(0) = a1 a2/E', lcoef(R.L2, 0) == w.a1 * w.a2 % P * inv(E) % P))
    ok.append(('Delta_2(0) = e1 e2 a1 a2 (a0+a1)(a0+a2)/E',
               lcoef(R.D2, 0) == w.e1 * w.e2 % P * w.a1 % P * w.a2 % P * w.B % P * w.A % P * inv(E) % P))
    return ok


def closed_form_G(w, R, N, K):
    """(E/(E+z)) (tau/t) F Ephi as [t^k] -> z-series 0..N; computed in z-Laurent, returns G[k][m]."""
    L = R.L
    E = w.Ea
    A = w.A
    a2 = w.a2
    tau = R.tau
    # dphi_t = (1-tau) a2 s*/(1-tau s*)   (envelope)
    onemt = sub(const(1, L), tau)
    dphidt = mul(mul(onemt, scal(a2, R.s), L), sinv(R.U1, L), L)
    dphi = sint(dphidt, L)
    # expand in tau rather than t: F_n rational in z with z^-n; use Laurent in z of length Nz with offset K
    off = K + 2
    Nz = N + off + 3 * K + 8

    def laur_one():
        r = [0] * Nz
        r[off] = 1
        return r

    def laur_mul(a, b):
        # a,b lists with offset 'off' each; result offset 'off'
        r = [0] * Nz
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    k = i + j - off
                    if y and 0 <= k < Nz:
                        r[k] += x * y
        return [v % P for v in r]
    # F as a t-series with coefficient Laurent in z.
    # F_n = prod_{j=1}^n (a2+jz)(A+jz) / (n! z^n prod_{j=2}^{n+1}(E+jz))
    Fn = []
    cur = laur_one()
    for n in range(0, K + 1):
        if n > 0:
            # multiply by (a2+nz)(A+nz)/(n z (E+(n+1)z))
            lin1 = [0] * Nz
            lin1[off] = a2
            lin1[off + 1] = n % P
            lin2 = [0] * Nz
            lin2[off] = A
            lin2[off + 1] = n % P
            cur = laur_mul(cur, lin1)
            cur = laur_mul(cur, lin2)
            # divide by z
            cur = cur[1:] + [0]
            # divide by n
            cur = scal(inv(n), cur)
            # geometric 1/(E+(n+1)z)
            ge = [0] * Nz
            iE = inv(E)
            for i in range(0, Nz - off):
                ge[off + i] = iE * pow((-(n + 1) * iE) % P, i, P) % P
            cur = laur_mul(cur, ge)
        Fn.append(cur)
    # tau^n as t series; multiply by Ephi: exp(-dphi/z): as t-series with coefficients Laurent in z
    G = [[0] * (N + 1) for _ in range(K + 1)]
    # build (tau/t) tau^n in t
    tt = tau[1:] + [0]
    # Ephi coefficient in t: E_k(z) = sum_l (-1)^l/l! z^-l [t^k] dphi^l
    Eph = [[0] * Nz for _ in range(K + 1)]
    dp = const(1, L)
    lf = 1
    for l in range(0, K + 1):
        if l > 0:
            dp = mul(dp, dphi, L)
            lf = lf * l % P
        co = (-1) ** l * inv(lf) % P
        for k in range(K + 1):
            if off - l >= 0:
                Eph[k][off - l] = (Eph[k][off - l] + co * dp[k]) % P
    # F as t-series
    Ft = [[0] * Nz for _ in range(K + 1)]
    tn = tt[:]  # (tau/t) tau^0
    for n in range(K + 1):
        for k in range(K + 1):
            if tn[k]:
                for i in range(Nz):
                    Ft[k][i] = (Ft[k][i] + tn[k] * Fn[n][i]) % P
        tn = mul(tn, tau, L)
    pref = [0] * Nz
    iE = inv(E)
    for i in range(0, Nz - off):
        pref[off + i] = E * iE * pow((-iE) % P, i, P) % P
    for k in range(K + 1):
        s = [0] * Nz
        for a in range(k + 1):
            t = laur_mul(Ft[a], Eph[k - a])
            for i in range(Nz):
                s[i] = (s[i] + t[i]) % P
        s = laur_mul(s, pref)
        # no negative powers expected
        neg = any(s[i] for i in range(off))
        if neg:
            raise RuntimeError('pole at z=0 at t^%d' % k)
        for m in range(N + 1):
            G[k][m] = s[off + m]
    return G


def main():
    N = 12
    K = 6
    L = K + 8
    for sm in SAMPLES:
        w = mk(sm)
        print('=== weights', [str(x) for x in w.raw], 'CY' if w.cy else 'nonCY')
        R = build_rows(w, N, L)
        for nm, b in check_identities(w, R):
            print('  ', nm, b)
        E = w.Ea
        # constants at t=0
        pt0 = [R.r0pt1[m][0] for m in range(N + 1)]
        geo = [E * inv(E) * pow((-inv(E)) % P, m, P) % P for m in range(N + 1)]
        hb = h_series([E, (-w.a1) % P, (-w.A) % P], N)
        print('   r^pt_{0,1}(z;0) = h_(E,-a1,-A) E/(E+z):', pt0 == zmul(hb, geo, N))
        hee = h_series([E, (-w.e1) % P, (-w.e2) % P], N)
        hnw = h_series([w.a1, w.A, (-w.e1) % P, (-w.e2) % P], N)
        true0 = zmul(hnw, pt0, N)
        print('   h_nw r^pt_{0,1}(z;0) = h_(E,-e1,-e2) E/(E+z):', true0 == zmul(hee, geo, N))
        print('   control: h_nw r^pt_{0,1}(z;0) = h_(E,-a1,-A) E/(E+z):', true0 == zmul(hb, geo, N))
        # regular colour at t=0
        reg0 = [R.r0pt2[m][0] for m in range(N + 1)]
        hr = h_series([w.a2, (-E) % P, w.B], N)
        print('   r^pt_{0,2}(z;0) = h_(a2,-E,a0+a1):', reg0 == hr)
        # frame B rows
        rows, fac = framed_rows(R, [E, (-w.e1) % P, (-w.e2) % P])
        r0 = rows[1][0]
        print('   frame-B r_{0,1}(z;0) = E/(E+z):', [lcoef(r0[m], 0) for m in range(N + 1)] == geo)
        # rho = (Delta_1 t^2/(-e1 e2 E^3))^(-1/2), root 1 at 0
        D1 = R.D1[1]
        lead = (-w.e1 * w.e2 % P * pow(E, 3, P)) % P
        rho_m2 = scal(inv(lead), D1)
        rm1 = ssqrt(rho_m2, 1, L)
        rho = sinv(rm1, L)
        # G from the Wick rows
        Gw = [[0] * (N + 1) for _ in range(K + 1)]
        for m in range(N + 1):
            ser = mul(rho, r0[m][1], L)
            for k in range(K + 1):
                Gw[k][m] = ser[k]
        Gc = closed_form_G(w, R, N, K)
        print('   rho r_{0,1} = (E/(E+z))(tau/t) F Ephi through t^%d z^%d:' % (K, N), Gw == Gc)
        # leaf test: [t^k](z G) * prod_{j<=k+1}(E+jz) polynomial of degree <= k (k>=1), no constant term
        res = []
        for k in range(1, K + 1):
            ser = [0] + Gw[k][:]  # zG, z^0..z^(N+1)
            NN = len(ser) - 1
            for j in range(1, k + 2):
                lin = [E, j % P]
                ser = [(E * ser[i] + (j * ser[i - 1] if i > 0 else 0)) % P for i in range(NN + 1)]
            deg = max([i for i, x in enumerate(ser) if x] + [-1])
            res.append((k, deg, ser[0]))
        print('   leaf test k -> degree of [t^k](zG) prod(E+jz) (expect <= k), constant term 0:', res)
        # control: wrong scale E' = E+1
        res = []
        Ep = (E + 1) % P
        for k in range(1, 3):
            ser = [0] + Gw[k][:]
            NN = len(ser) - 1
            for j in range(1, k + 2):
                ser = [(Ep * ser[i] + (j * ser[i - 1] if i > 0 else 0)) % P for i in range(NN + 1)]
            res.append((k, max([i for i, x in enumerate(ser) if x] + [-1])))
        print('   control wrong scale E+1 (expect full degree):', res)
        # closed form of the leaf coefficients: A_j = (-1)^j C_j t^-1 zeta^j 2F1(j B/E, j a1/E; j; tau)
        # with zeta = tau exp(dphi/E), C_j = P_{j-1}/(j! (j-1)! E^(2j-4)), P_k = prod_{i=1}^k ((k+1) a2 - i E)((k+1) A - i E)
        tau = R.tau
        onemt = sub(const(1, L), tau)
        dphidt = mul(mul(onemt, scal(w.a2, R.s), L), sinv(R.U1, L), L)
        dphi = sint(dphidt, L)
        okA = True
        Ahat = {}
        for j in range(1, K + 2):
            Pk = 1
            for i in range(1, j):
                Pk = Pk * (j * w.a2 - i * E) % P * (j * w.A - i * E) % P
            fj = 1
            for q in range(2, j + 1):
                fj = fj * q % P
            fj1 = 1
            for q in range(2, j):
                fj1 = fj1 * q % P
            Cj = Pk * inv(fj * fj1 % P) % P * pow(inv(E), 2 * j - 4, P) % P if j >= 2 else E * E % P
            # 2F1
            aa = j * w.B % P * inv(E) % P
            bb = j * w.a1 % P * inv(E) % P
            F2 = [0] * L
            term = const(1, L)
            tp = const(1, L)
            co = 1
            for n in range(L):
                F2 = add(F2, scal(co, tp))
                co = co * (aa + n) % P * (bb + n) % P * inv((j + n) * (n + 1) % P) % P
                tp = mul(tp, tau, L)
            zeta = mul(tau, sexp(scal(inv(E), dphi), L), L)
            sA = mul(spow(zeta, j, L), F2, L)
            sA = scal((-1) ** j * Cj % P, sA)
            Ahat[j] = sA  # t A_j
        # compare partial fractions: [t^k] zG = Sum_j [t^(k+1)]Ahat_j/(E+jz), k>=1; at k=0: E - E^2/(E+z)
        for k in range(0, K + 1):
            for m in range(0, N):
                pred = 0
                for j in range(1, K + 2):
                    co = Ahat[j][k + 1]
                    pred += co * inv(E) * pow((-j * inv(E)) % P, m + 1, P)
                pred %= P
                if pred != Gw[k][m]:
                    okA = False
        print('   closed form A_j = (-1)^j C_j t^-1 zeta^j 2F1(j(a0+a1)/E, j a1/E; j; tau) reproduces [t^k] zG:', okA)
        print('   P_1 = -u? ', (2 * w.a2 - E) * (2 * w.A - E) % P == (-(E - 2 * w.a1) * (E - 2 * w.a2)) % P)


if __name__ == '__main__':
    main()
