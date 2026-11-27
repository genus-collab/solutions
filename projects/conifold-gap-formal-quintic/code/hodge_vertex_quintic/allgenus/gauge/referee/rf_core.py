# -*- coding: utf-8 -*-
"""rf_core.py -- REFEREE of task `gauge`: own exact tools (stdlib only; imports NO gi_* module and NO project module).

Power series in u: lists of Fractions [c_0, .., c_{N-1}] = Sum c_k u^k + O(u^N).
  periods(N, al, be)  own Frobenius solution of the quintic Picard-Fuchs operator at u = 1 - 5^5 q = 0 (recursion on monomials:
                      P = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u), D = q d/dq = -(1-u) d/du);
                      the vanishing period g from the vanishing of the u^-3..u^-1 obstruction of g log u + h;
                      varpi_0 = f_0 + al g + be f_2, t_c = g/varpi_0.
  systems(N, al, be)  b_a (quintic) and bt_a (model pulled back by v = u tau^(-5/4)); index 1 is stored as u*b_1 (a unit).
  model_vc / quintic_vc   the two vanishing solutions {n: [5 series]} (coefficient of z^-n), as in the brief.
"""
from fractions import Fraction as Fr
CEXP = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)


def mul(a, b, N):
    out = [Fr(0)] * N
    for i, x in enumerate(a[:N]):
        if x:
            for j, y in enumerate(b[:N - i]):
                if y: out[i + j] += x * y
    return out


def add(a, b, N, c=1):
    return [(a[k] if k < len(a) else 0) + c * (b[k] if k < len(b) else 0) for k in range(N)]


def scale(a, c): return [c * x for x in a]


def inv(a, N):
    assert a[0] != 0
    out = [Fr(1) / a[0]] + [Fr(0)] * (N - 1)
    for n in range(1, N):
        s = sum(a[j] * out[n - j] for j in range(1, min(n, len(a) - 1) + 1))
        out[n] = -s / a[0]
    return out


def deriv(a): return [k * a[k] for k in range(1, len(a))]


def integ(a): return [Fr(0)] + [a[k] / (k + 1) for k in range(len(a))]


def powf(a, al, N):
    """a^al for a[0] = 1 (recursion n f_n = Sum_{j>=1} (al j - n + j) a_j f_{n-j})"""
    assert a[0] == 1
    f = [Fr(1)] + [Fr(0)] * (N - 1)
    for n in range(1, N):
        f[n] = sum((al * j - n + j) * a[j] * f[n - j] for j in range(1, min(n, len(a) - 1) + 1)) / n
    return f


def ldict_mul(A, B):
    out = {}
    for i, x in A.items():
        for j, y in B.items(): out[i + j] = out.get(i + j, 0) + x * y
    return {k: v for k, v in out.items() if v != 0}


def ldict_D(A):
    """D = -(1-u) d/du on a Laurent polynomial {power: coefficient}"""
    out = {}
    for m, c in A.items():
        if m:
            out[m - 1] = out.get(m - 1, 0) - m * c; out[m] = out.get(m, 0) + m * c
    return {k: v for k, v in out.items() if v != 0}


PCOEF = {4: {1: Fr(1)}, 3: {0: Fr(-2), 1: Fr(2)}, 2: {0: Fr(-7, 5), 1: Fr(7, 5)}, 1: {0: Fr(-2, 5), 1: Fr(2, 5)},
         0: {0: Fr(-24, 625), 1: Fr(24, 625)}}


def L_on(A):
    """the Picard-Fuchs operator on a Laurent polynomial"""
    jets = [A]
    for k in range(4): jets.append(ldict_D(jets[-1]))
    out = {}
    for k in range(5):
        for m, c in ldict_mul(PCOEF[k], jets[k]).items(): out[m] = out.get(m, 0) + c
    return {k: v for k, v in out.items() if v != 0}


def frob(N, b012):
    """power-series solution with (b_0, b_1, b_2) = b012, O(u^N)"""
    Lm = [L_on({n: Fr(1)}) for n in range(N + 4)]
    for n in range(N + 4):
        assert min(Lm[n]) >= n - 3 and Lm[n].get(n - 3, 0) == n * (n - 1) ** 2 * (n - 2)
    b = [Fr(x) for x in b012] + [Fr(0)] * (N - 3)
    for m in range(0, N - 3):
        s = sum(b[n] * Lm[n].get(m, 0) for n in range(max(m - 1, 0), m + 3))
        b[m + 3] = -s / Lm[m + 3][m]
    return b


def log_obstruction(g):
    """u^-3, u^-2, u^-1 coefficients of K(g) = P(g log u) - log u P(g) (Leibniz for the derivation D; w = D log u = 1 - 1/u)"""
    G = {k: x for k, x in enumerate(g[:10]) if x}
    w = [{0: Fr(1), -1: Fr(-1)}]
    for j in range(4): w.append(ldict_D(w[-1]))
    gj = [G]
    for j in range(4): gj.append(ldict_D(gj[-1]))
    binom = lambda k, j: [1, k, k * (k - 1) // 2, k * (k - 1) * (k - 2) // 6, 1][j] if j < 4 else 1
    out = {}
    for k in range(1, 5):
        for j in range(1, k + 1):
            t = ldict_mul(PCOEF[k], ldict_mul(gj[k - j], w[j - 1]))
            for m, c in t.items(): out[m] = out.get(m, 0) + binom(k, j) * c
    return [out.get(m, Fr(0)) for m in (-3, -2, -1)]


def periods(N, al=0, be=0):
    f0 = frob(N + 2, (1, 0, 0)); f1 = frob(N + 2, (0, 1, 0)); f2 = frob(N + 2, (0, 0, 1))
    o0, o1, o2 = log_obstruction(f0), log_obstruction(f1), log_obstruction(f2)
    assert o1[0] == o1[1] == 0 and o2[0] == o2[1] == 0 and o2[2] != 0
    gam = -o1[2] / o2[2]
    g = add(f1, f2, N + 2, gam)
    assert log_obstruction(g) == [0, 0, 0]
    v0 = add(add(f0, g, N + 2, Fr(al)), f2, N + 2, Fr(be))
    tc = mul(g, inv(v0, N + 2), N + 2)
    return dict(f0=f0, f1=f1, f2=f2, g=g, gamma=gam, obst_f0=o0, varpi0=v0[:N + 1], tc=tc[:N + 1])


def systems(N, al=0, be=0):
    """b = [b0, u*b1, b2, b3, b4], bt likewise (power series O(u^N)); also tau, v, varpi0, tc"""
    per = periods(N + 3, al, be); v0 = per['varpi0']; tc = per['tc']; M = N + 2
    tcp = deriv(tc)[:M]; one_m = [Fr(1), Fr(-1)] + [Fr(0)] * (M - 2)
    om3 = mul(mul(one_m, one_m, M), one_m, M)
    ub1 = scale(inv(mul(mul(mul(v0, v0, M), om3, M), mul(tcp, tcp, M), M), M), -1)
    b3 = scale(mul(v0, inv(one_m, M), M), -1)
    b = [tcp[:N], ub1[:N], tcp[:N], b3[:N], b3[:N]]
    tau = inv([Fr(4, 5 * n + 4) for n in range(M)], M)
    vs = powf(tau, Fr(-5, 4), M); v = [Fr(0)] + vs[:M - 1]; vp = deriv(v)
    ubt1 = scale(mul(vp, inv(vs, M), M), -1)
    bt = [vp[:N], ubt1[:N], vp[:N], scale(vp[:N], -1), scale(vp[:N], -1)]
    return dict(b=b, bt=bt, tau=tau[:N], v=v[:N], varpi0=v0[:N], tc=tc[:N], per=per)


def fact(n):
    out = 1
    for k in range(2, n + 1): out *= k
    return out


def ic_terms(m, nmax):
    """I^C_m = Sum_i (-1/(4z))^i/i! (-v/z)^j/j!, j = m+4i+1 >= 0: list of (n = i + j, j, coefficient of v^j z^-n)"""
    out = []
    for i in range(0, nmax + 2):
        j = m + 4 * i + 1
        if j < 0: continue
        if i + j > nmax: break
        out.append((i + j, j, Fr((-1) ** i, 4 ** i * fact(i)) * Fr((-1) ** j, fact(j))))
    return out


def model_vc(v, N, nmax):
    """Y^C(v(u), z) = (I_0, -I_-1, -I_3, I_2, I_1): {n: [5 series O(u^N)]}, coefficient of z^-n"""
    vp = [[Fr(1)] + [Fr(0)] * (N - 1)]
    for j in range(1, 2 * nmax + 2): vp.append(mul(vp[-1], v, N))
    out = {}
    for b in range(5):
        for (n, j, c) in ic_terms(CEXP[b], nmax):
            vec = out.setdefault(n, [[Fr(0)] * N for _ in range(5)])
            vec[b] = add(vec[b], vp[j], N, c * EPS[b])
    return out


def quintic_vc(b, N, nmax):
    """y^(0) = -e_1, y^(n)_a = Int_0^u b_a y^(n-1)_{a+1} (b_1 = (u b_1)/u: the integrand is divided by u, asserting no 1/u term)"""
    y = [[Fr(0)] * N for _ in range(5)]; y[1][0] = Fr(-1); out = {0: y}
    for n in range(1, nmax + 1):
        new = []
        for a in range(5):
            t = mul(b[a], y[(a + 1) % 5], N)
            if a == 1:
                assert t[0] == 0, "log term in the iterated integral at n = %d" % n
                t = t[1:] + [Fr(0)]
            new.append(integ(t)[:N])
        y = new; out[n] = y
    return out
