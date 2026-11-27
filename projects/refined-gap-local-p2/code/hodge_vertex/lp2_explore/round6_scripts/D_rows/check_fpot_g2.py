"""s.1.5 / (H0) at genus 2 with the actual Hodge class Omega_w, exact; evidence only.
Fresh Hodge integrals: Witten-Kontsevich by DVV, lambda_1 by Mumford, lambda_2 and lambda_1 lambda_2 by the lambda_g, lambda_g lambda_{g-1} formulas.
Tests: known values (controls); Xi_2((2))=0; Xi_2((3)) nonzero; Xi_2(1^p)=0; (S_P2); (L0_P2); Z[2,0] vs -beta_2/(e1 e2 E^3);
(H0) and the overlap of the Taylor and dilaton-shifted definitions, with the true q_t expanded in t_c."""
import sys, itertools, time
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial
sys.setrecursionlimit(100000)
from rows_core import *
from check_conifold import Ser, powser_1mx, setup

def dfact(n):
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out

@lru_cache(maxsize=None)
def wk(g, d):
    """<tau_d1...tau_dn>_g, d sorted tuple"""
    n = len(d)
    if g < 0: return Fr(0)
    if 2 * g - 2 + n <= 0: return Fr(0)
    if sum(d) != 3 * g - 3 + n: return Fr(0)
    if g == 0 and n == 3: return Fr(1)
    if g == 1 and n == 1: return Fr(1, 24)
    d = list(d)
    if d[0] == 0:
        rest = d[1:]
        tot = Fr(0)
        for j in range(len(rest)):
            if rest[j] > 0:
                r = rest[:]; r[j] -= 1
                tot += wk(g, tuple(sorted(r)))
        return tot
    if 1 in d:
        i = d.index(1)
        r = d[:i] + d[i + 1:]
        return (2 * g - 2 + len(r)) * wk(g, tuple(r))
    # DVV on the largest
    k = d[-1] - 1
    rest = d[:-1]
    tot = Fr(0)
    for j in range(len(rest)):
        r = rest[:]; r[j] += k
        tot += Fr(dfact(2 * k + 2 * rest[j] + 1), dfact(2 * rest[j] - 1)) * wk(g, tuple(sorted(r)))
    for r_ in range(0, k):
        s_ = k - 1 - r_
        cf = dfact(2 * r_ + 1) * dfact(2 * s_ + 1)
        tot += Fr(cf, 2) * wk(g - 1, tuple(sorted(rest + [r_, s_])))
        idx = list(range(len(rest)))
        for mask in range(1 << len(rest)):
            I = [rest[i] for i in idx if mask >> i & 1]
            J = [rest[i] for i in idx if not mask >> i & 1]
            # dimension fixes g1
            need = sum(I) + r_ - len(I) - 1 + 3
            if need % 3: continue
            g1 = need // 3
            g2 = g - g1
            if g1 < 0 or g2 < 0: continue
            tot += Fr(cf, 2) * wk(g1, tuple(sorted(I + [r_]))) * wk(g2, tuple(sorted(J + [s_])))
    return tot / dfact(2 * k + 3)

@lru_cache(maxsize=None)
def lam1(g, a):
    """Int_{Mbar_{g,n}} lambda_1 psi^a, Mumford: lambda_1=(kappa_1 - Sum psi + delta)/12"""
    n = len(a)
    if sum(a) != 3 * g - 3 + n - 1: return Fr(0)
    a = list(a)
    tot = wk(g, tuple(sorted(a + [2])))
    for i in range(n):
        b = a[:]; b[i] += 1
        tot -= wk(g, tuple(sorted(b)))
    tot += Fr(1, 2) * wk(g - 1, tuple(sorted(a + [0, 0])))
    idx = range(n)
    for mask in range(1 << n):
        I = [a[i] for i in idx if mask >> i & 1]
        J = [a[i] for i in idx if not mask >> i & 1]
        for g1 in range(0, g + 1):
            g2 = g - g1
            if 2 * g1 - 2 + len(I) + 1 <= 0 or 2 * g2 - 2 + len(J) + 1 <= 0: continue
            tot += Fr(1, 2) * wk(g1, tuple(sorted(I + [0]))) * wk(g2, tuple(sorted(J + [0])))
    return tot / 12

def multinom(n, a):
    out = factorial(n)
    for x in a:
        out //= factorial(x)
    return out

def lam2(a):
    """g=2: Int lambda_2 psi^a = multinomial(1+n; a) * 7/5760, Sum a = 1 + n"""
    n = len(a)
    if sum(a) != 1 + n: return Fr(0)
    return multinom(1 + n, a) * Fr(7, 5760)

@lru_cache(maxsize=None)
def lam12(a):
    """g=2: Int lambda_1 lambda_2 psi^a, Sum a = n"""
    n = len(a)
    if sum(a) != n: return Fr(0)
    if 0 in a:
        a = list(a); a.remove(0)
        tot = Fr(0)
        for j in range(len(a)):
            if a[j] > 0:
                b = a[:]; b[j] -= 1
                tot += lam12(tuple(sorted(b)))
        return tot
    den = 1
    for x in a:
        den *= dfact(2 * x - 1)
    return Fr(factorial(1 + n)) * Fr(1, 30) / (8 * 24 * den)

class G2:
    def __init__(self, e1, e2):
        self.e1, self.e2 = e1, e2
        self.E = e1 + e2
        u = [1 / self.E, -1 / e1, -1 / e2]
        self.U1 = sum(u)
        self.U2 = u[0] * u[1] + u[0] * u[2] + u[1] * u[2]
        self.U3 = u[0] * u[1] * u[2]
    def integ(self, a):
        a = tuple(sorted(a))
        n = len(a); s = sum(a)
        dim = 3 + n
        if s == dim: return wk(2, a)
        if s == dim - 1: return -self.U1 * lam1(2, a)
        if s == dim - 2: return self.U1 ** 2 * lam2(a)
        if s == dim - 3: return -(self.U1 * self.U2 - self.U3) * lam12(a)
        return Fr(0)

def lamser(E, j, D):
    """lambda_j(psi) coefficients through psi^D"""
    return [Fr(j) ** n / E ** (n + 1) for n in range(D + 1)]

def Tser(E, D):
    # T = -psi^2/(E-psi)
    return [Fr(0), Fr(0)] + [-(1 / E) ** (n - 1) for n in range(2, D + 1)]

def integrate(G, series_list):
    """Int Omega prod_i g_i(psi_i), g_i coefficient lists.  Exponents cut by dim."""
    n = len(series_list)
    dim = 3 + n
    # enumerate
    tot = Fr(0)
    def rec(i, cur, used, coef):
        nonlocal tot
        if i == n:
            if used >= dim - 3:
                tot += coef * G.integ(cur)
            return
        ser = series_list[i]
        for e, c in enumerate(ser):
            if used + e > dim: break
            if c == 0: continue
            rec(i + 1, cur + [e], used + e, coef * c)
    rec(0, [], 0, Fr(1))
    return tot

def frakd(G, gs, T):
    """frak d^s (g_1,...,g_s) ; T = series"""
    s = len(gs)
    tot = Fr(0)
    for m in range(0, 3 + s + 1):
        tot += integrate(G, gs + [T] * m) / factorial(m)
    return tot

def dilshift(G, fs, q):
    """dilaton-shifted d^n Fpot(q)[f]: q=-psi+..., q+psi=O(psi^2)"""
    assert q[0] == 0 and q[1] == -1
    Tq = [q[0], Fr(0)] + q[2:]
    n = len(fs)
    tot = Fr(0)
    for m in range(0, 3 + n + 1):
        tot += integrate(G, fs + [Tq] * m) / factorial(m)
    return tot

def main():
    # controls
    print("controls: <tau_4>_2 =", wk(2, (4,)), "(1/1152); <tau_2 tau_3>_2 =", wk(2, (2, 3)), "(29/5760); <tau_2^3>_2 =", wk(2, (2, 2, 2)), "(7/240)")
    print("  <tau_3 lam1>_2 =", lam1(2, (3,)), "(1/480); <tau_2^2 lam1>_2 =", lam1(2, (2, 2)), "(5/576); <lam1^3>: 2*lam1lam2 =", 2 * lam12(()), "(1/2880); <tau_1 lam1 lam2> =", lam12((1,)), "(dilaton 2*1/5760)")
    print("  lam1 on Mbar_2,0 check: lam1(2,(1,1)) via lam1 =", lam1(2, (1, 1)) , " lam1(2,(0,2)) =", lam1(2, (0, 2)))
    for (e1, e2) in [(Fr(-2, 7), Fr(9, 7)), (Fr(2, 3), Fr(5, 7)), (Fr(-11, 5), Fr(3, 2))]:
        E = e1 + e2
        G = G2(e1, e2)
        D = 12
        T = Tser(E, D)
        Z0 = frakd(G, [], T)
        beta2 = (7 * e1 ** 4 + 10 * e1 ** 2 * e2 ** 2 + 7 * e2 ** 4) / (5760 * e1 * e2)
        print(f"e=({e1},{e2}): Z[2,0] = {Z0}; -beta_2/(e1 e2 E^3) = {-beta2/(e1*e2*E**3)}; equal {Z0 == -beta2/(e1*e2*E**3)}")
        X = lambda *R: frakd(G, [lamser(E, j, D) for j in R], T)
        print("  Xi_2((2)) =", X(2), " Xi_2((3)) =", X(3), " Xi_2((2,2)) =", X(2, 2), " Xi_2(1) =", X(1), " Xi_2(1,1) =", X(1, 1))
        # string: Xi(1,R)=(J_R/E^2)Xi(R)
        print("  (S_P2) R=(3):", X(1, 3) == Fr(3) / E ** 2 * X(3), " R=(2,3):", X(1, 2, 3) == Fr(5) / E ** 2 * X(2, 3))
        # L0: Xi(0,R)=((D_R-2h+2)/E^2)Xi(R)
        print("  (L0_P2) R=(3):", X(0, 3) == Fr(2 - 4 + 2) / E ** 2 * X(3), " R=empty:", X(0) == Fr(-2) / E ** 2 * Z0, " R=(3,4):", X(0, 3, 4) == Fr(5 - 2) / E ** 2 * X(3, 4))
    # (H0) with the true q_t.  Work mod t_c^(K+1).
    e1, e2 = Fr(2, 3), Fr(5, 7); E = e1 + e2
    G = G2(e1, e2)
    K, D = 3, 14
    Ser.M = K + 8
    x, wc, tc, theta, y = setup()
    xt = Ser([Fr(0), Fr(1)]); Tv = Ser([Fr(0), Fr(1)])
    for _ in range(Ser.M + 1):
        xt = Tv - (tc - x).compose(xt)
    v = ((1 - powser_1mx(Fr(-1, 3))) * (Fr(3) / E)).compose(xt)
    rho = (v * (-E)).shift(-1)
    S = series_S(E, D, None)
    # q_t(psi) = -psi S(v,-psi): list over psi^d of Ser in t_c
    vp = [Ser.const(1)]
    for k in range(1, D + 2):
        vp.append(vp[-1] * v)
    qt = [Ser.const(0)]
    for d in range(1, D + 1):
        n = d - 1
        val = Ser.const(0)
        for e_, c in S[n].items():
            val = val + vp[e_] * c
        qt.append(val * Fr((-1) ** d))
    # coefficient table: qt[d].c[k]
    def tc_coeff(serlist, k):
        return [s.c[k] for s in serlist]
    q0 = [Fr(0)] + [-(1 / E) ** (d - 1) for d in range(1, D + 1)]
    print("q_t at t_c^0 = q_0:", tc_coeff(qt, 0) == q0, "; Dil(1) = -1 at every t_c order:", qt[1].c[0] == -1 and all(c == 0 for c in qt[1].c[1:]))
    f = lamser(E, 2, D)
    T = Tser(E, D)
    # polys in t_c: represent as lists length K+1
    def tmul(a, b):
        o = [Fr(0)] * (K + 1)
        for i, x_ in enumerate(a):
            for j, y_ in enumerate(b):
                if i + j <= K: o[i + j] += x_ * y_
        return o
    # Taylor definition: sum_r (1/r!) d^(n+r)(f, delta^r), delta = Sum_k>=1 t^k delta_k.
    def taylor(fs, deltas):
        """deltas: dict k -> series (k>=1). returns list t_c^0..K"""
        out = [Fr(0)] * (K + 1)
        keys = sorted(deltas)
        for r in range(0, K + 1):
            for combo in itertools.product(keys, repeat=r):
                w = sum(combo)
                if w > K: continue
                val = frakd(G, fs + [deltas[k] for k in combo], T) / factorial(r)
                out[w] += val
        return out
    # delta for q_t
    deltas = {k: tc_coeff(qt, k) for k in range(1, K + 1)}
    t0 = time.time()
    LHS_taylor = taylor([f], deltas)
    print("  Taylor form done", round(time.time() - t0, 1))
    # dilaton-shifted form, coefficientwise in t_c: expand multilinearly
    def dil_form(fs):
        # q+psi = T + delta; identical expansion would be the same as Taylor, so compute directly by evaluating at
        # t_c-polynomial coefficients: Sum_m 1/m! Int prod (T+delta)
        out = [Fr(0)] * (K + 1)
        n = len(fs)
        Tq = {0: T}
        for k, ser in deltas.items():
            Tq[k] = ser
        keys = sorted(Tq)
        for m in range(0, 3 + n + 1):
            for combo in itertools.product(keys, repeat=m):
                w = sum(combo)
                if w > K: continue
                out[w] += integrate(G, fs + [Tq[k] for k in combo]) / factorial(m)
        return out
    LHS_dil = dil_form([f])
    print("overlap of the two definitions, d^1 Fpot(q_t)[lam_2], t_c^0..K:", LHS_taylor == LHS_dil, [str(x_)[:18] for x_ in LHS_dil])
    # rho q_t
    rq = [s * rho for s in qt]
    assert tc_coeff(rq, 0) == q0
    deltas_r = {k: tc_coeff(rq, k) for k in range(1, K + 1)}
    print("rho q_t has psi-linear delta:", [deltas_r[k][1] for k in deltas_r], " (=-[t^k] rho)", [-rho.c[k] for k in deltas_r])
    # RHS: rho^(2h-2) * d^1 Fpot(rho q_t)[rho f] = rho^3 * Taylor
    fs_t = taylor([f], deltas_r)
    rho3 = [(rho * rho * rho).c[k] for k in range(K + 1)]
    RHS = tmul(rho3, fs_t)
    print("(H0)/(7.1) n=1, lambda=rho:", RHS == LHS_dil)
    # n=0
    L0 = taylor([], deltas)
    R0 = tmul([(rho * rho).c[k] for k in range(K + 1)], taylor([], deltas_r))
    print("(H0) n=0: Fpot(q_t) = rho^2 Fpot(rho q_t):", L0 == R0, [str(x_)[:18] for x_ in L0])
    # control: wrong exponent
    R0bad = tmul([(rho).c[k] for k in range(K + 1)], taylor([], deltas_r))
    print("  control wrong exponent:", L0 == R0bad)
    # the polar structure: Fpot(rho q_t) t_c^k coefficient, D0 check: k=1 should be Taylor at Lambda_2 leaves with no constant
    print("(L1) [t_c^k](rho q_t) leaves, no constant: k=1:", [str(c)[:10] for c in deltas_r[1][:4]])
    # expected for genus 2 (Theorem R at single vertex, M_irr(2)): t_c^-1 coefficient of Delta_0 Fpot(q_t) vanishes => [t^1]Fpot(rho q_t)=0
    print("gap at genus 2: [t_c^1] Fpot_2(rho q_t) =", taylor([], deltas_r)[1], " [t_c^0] =", taylor([], deltas_r)[0])

if __name__ == "__main__":
    main()
