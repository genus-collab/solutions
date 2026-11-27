"""Expansions at the conifold, all exact in Q(E) with rational E; evidence only.
Fresh code: Lemma 1.2, Cor 2.9 (F2, L1, Z-infinity), Lemma D0, s.0.2 F1-F3, L1, L3, K=1. Controls: off-shell c1."""
import sys, time
from fractions import Fraction as Fr
from rows_core import *

# ---- truncated power series in one variable, list of coefficients, order < M
class Ser:
    M = 30
    __slots__ = ("c",)
    def __init__(self, c=None):
        c = list(c) if c else []
        self.c = (c + [Fr(0)] * Ser.M)[:Ser.M]
    @staticmethod
    def const(a):
        return Ser([Fr(a)])
    def __add__(s, o):
        if not isinstance(o, Ser): o = Ser.const(o)
        return Ser([a + b for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return Ser([-a for a in s.c])
    def __sub__(s, o):
        if not isinstance(o, Ser): o = Ser.const(o)
        return Ser([a - b for a, b in zip(s.c, o.c)])
    def __rsub__(s, o): return Ser.const(o) - s
    def __mul__(s, o):
        if not isinstance(o, Ser):
            return Ser([a * o for a in s.c])
        out = [Fr(0)] * Ser.M
        for i, a in enumerate(s.c):
            if a == 0: continue
            for j, b in enumerate(o.c):
                if i + j >= Ser.M: break
                if b != 0:
                    out[i + j] += a * b
        return Ser(out)
    __rmul__ = __mul__
    def inv(s):
        assert s.c[0] != 0
        out = [Fr(0)] * Ser.M
        out[0] = 1 / s.c[0]
        for n in range(1, Ser.M):
            t = Fr(0)
            for k in range(1, n + 1):
                t += s.c[k] * out[n - k]
            out[n] = -t * out[0]
        return Ser(out)
    def __truediv__(s, o):
        if not isinstance(o, Ser):
            return Ser([a / o for a in s.c])
        return s * o.inv()
    def d(s):
        return Ser([(i + 1) * s.c[i + 1] for i in range(Ser.M - 1)])
    def integ(s):
        return Ser([Fr(0)] + [s.c[i] / (i + 1) for i in range(Ser.M - 1)])
    def shift(s, n):
        """multiply by x^n, n may be negative if low terms vanish"""
        if n >= 0:
            return Ser([Fr(0)] * n + s.c)
        assert all(a == 0 for a in s.c[:-n]), "shift not divisible"
        return Ser(s.c[-n:])
    def compose(s, g):
        """s(g(x)) with g(0)=0"""
        out = Ser.const(0)
        p = Ser.const(1)
        for a in s.c:
            if a != 0:
                out = out + p * a
            p = p * g
        return out

def powser_1mx(alpha):
    """(1-x)^alpha"""
    c = [Fr(1)]
    for n in range(1, Ser.M):
        c.append(c[-1] * (Fr(n - 1) - alpha) / n)
    return Ser(c)

def wc_series():
    # solve x(x-1)f''+(2x-1)f'+(2/9)f=0 from L2b directly, with L2b written fresh in x
    a = [Fr(1)]
    for n in range(0, Ser.M - 1):
        # coefficient of x^n : n(n-1)a_n -(n+1)n a_{n+1} + 2n a_n -(n+1)a_{n+1} + (2/9)a_n
        a.append((Fr(n * (n - 1) + 2 * n) + Fr(2, 9)) * a[n] / ((n + 1) ** 2))
    return Ser(a)

def setup():
    x = Ser([Fr(0), Fr(1)])
    wc = wc_series()
    # check against definition of L2b using theta=(x-1)d/dx, 27y = x-1
    theta = lambda f: (x - 1) * f.d()
    y = (x - 1) / 27
    res = x * theta(theta(wc)) + 27 * y * theta(wc) + 6 * y * wc
    assert all(a == 0 for a in res.c[:Ser.M - 3]), "wc does not solve L2b"
    tc = (wc * (1 - x).inv()).integ()   # dt_c/dx = -wc/(x-1)
    Ct = theta(tc)
    assert all(a + b == 0 for a, b in list(zip(Ct.c, wc.c))[:Ser.M - 2]), "theta t_c != -wc"
    return x, wc, tc, theta, y

def z_inf(Mx=61):
    Ser.M = Mx
    x, wc, tc, theta, y = setup()
    r = rec_r(Mx + 2)
    Ev = (1 - powser_1mx(Fr(-1, 3))) * 3    # Ev does not depend on E
    # Ev has zero constant term
    Sinf = Ser.const(0)
    p = Ser.const(1)
    fact = Fr(1)
    for l in range(0, Mx):
        p = p * Ev
        fact = fact * (l + 1)
        Sinf = Sinf + p * (r[l] / fact)
    ok = all(a + b == 0 for a, b in list(zip(Sinf.c, tc.c))[:Mx - 1])
    # L2b(theta Sinf)=0 check
    th = theta(Sinf)
    res = x * theta(theta(th)) + 27 * y * theta(th) + 6 * y * th
    print(f"Z-infinity: S_inf + t_c = 0 through x^{Mx-2}: {ok}")
    # control: wrong recurrence value
    r2 = rec_r(Mx + 2, (4, Fr(1, 10 ** 6)))
    Sbad = Ser.const(0); p = Ser.const(1); fact = Fr(1)
    for l in range(0, Mx):
        p = p * Ev; fact *= (l + 1)
        Sbad = Sbad + p * (r2[l] / fact)
    print("  control (r_4 perturbed): first order of failure", next(i for i, (a, b) in enumerate(zip(Sbad.c, tc.c)) if a + b != 0))

def bivar_tc(E, K, Z, offshell=None):
    """Return dict with sigma_a^(k)[n] = coefficient of t_c^k z^n in rho P^^(0)(z), k=-2..K, n=0..Z.
    offshell: None or series function multiplying c1 (control)."""
    Ser.M = K + 12
    x, wc, tc, theta, y = setup()
    # reversion: express x in t_c
    # solve by Newton-free iteration
    xt = Ser([Fr(0), Fr(1)])
    T = Ser([Fr(0), Fr(1)])
    for _ in range(Ser.M + 1):
        # x = T - (tc(x)-x)
        xt = T - (tc - x).compose(xt)
    assert all(a == b for a, b in zip(tc.compose(xt).c[:Ser.M - 1], T.c[:Ser.M - 1]))
    v_x = (1 - powser_1mx(Fr(-1, 3))) * (Fr(3) / E)
    v = v_x.compose(xt)
    c1_x = -wc
    if offshell is not None:
        c1_x = c1_x * offshell(x)
    XX_x = theta(c1_x) * c1_x.inv()
    C1i = c1_x.inv().compose(xt)
    XX = XX_x.compose(xt)
    assert all(a == 0 for a in (v.c[0],))
    rho = (v * (-E)).shift(-1)
    assert rho.c[0] == 1
    # d2 = E(C1i^2-1)/(2x); as series in t_c. xt = t*(unit)
    d2_x = ((c1_x.inv() * c1_x.inv() - 1) * (E / 2)).shift(-1)
    d2 = d2_x.compose(xt)
    a, b, Thb = abc_vs(E, Z, None)
    # precompute powers of v as Laurent series: v = t*u, write v^k = t^k u^k ; store list of (shift, Ser)
    u = v.shift(-1)
    upow = {}
    uinv = u.inv()
    def vpow(k):
        if k not in upow:
            if k == 0: upow[k] = Ser.const(1)
            elif k > 0: upow[k] = vpow(k - 1) * u
            else: upow[k] = vpow(k + 1) * uinv
        return upow[k]
    LOW = 3
    TOP = K
    def lp_to_laurent(lp):
        """Laurent polynomial in v -> dict k->Fraction t_c^k truncated to k<=TOP+2"""
        out = {}
        for e, cf in lp.items():
            if e > TOP + 3:
                continue
            s = vpow(e)
            for i, co in enumerate(s.c):
                if co != 0 and e + i <= TOP + 3:
                    out[e + i] = out.get(e + i, Fr(0)) + cf * co
        return out
    def ldict_mul_ser(d, s, sh=0):
        out = {}
        for e, cf in d.items():
            for i, co in enumerate(s.c):
                if co != 0 and e + i + sh <= TOP + 3:
                    out[e + i + sh] = out.get(e + i + sh, Fr(0)) + cf * co
        return out
    def ldict_add(p, q, c=Fr(1)):
        out = dict(p)
        for k, val in q.items():
            out[k] = out.get(k, Fr(0)) + c * val
        return out
    # rows: P^_0=a/v, P^_1=C1i b/v, P^_2=(x/C1i)[Thb - z XX b - d2 b]/v
    xs = xt.shift(-1)  # x = t * xs
    sig = [[None] * (Z + 1) for _ in range(3)]
    for n in range(Z + 1):
        a_v = lp_shift(a[n], -1)
        b_v = lp_shift(b[n], -1)
        Th_v = lp_shift(Thb[n], -1)
        A0 = lp_to_laurent(a_v)
        A1 = ldict_mul_ser(lp_to_laurent(b_v), C1i)
        inner = lp_to_laurent(Th_v)
        inner = ldict_add(inner, ldict_mul_ser(lp_to_laurent(b_v), d2), Fr(-1))
        if n >= 1:
            inner = ldict_add(inner, ldict_mul_ser(lp_to_laurent(lp_shift(b[n - 1], -1)), XX), Fr(-1))
        A2 = ldict_mul_ser(ldict_mul_ser(inner, xs, 1), C1i.inv() )
        sig[0][n] = ldict_mul_ser(A0, rho)
        sig[1][n] = ldict_mul_ser(A1, rho)
        sig[2][n] = ldict_mul_ser(A2, rho)
    return sig, rho

def in_z_at_k(sig, a, k, Z, sign=-1):
    """coefficient of t_c^k, as series in psi (z=-psi): list psi^n"""
    return [sig[a][n].get(k, Fr(0)) * (sign ** n) for n in range(Z + 1)]

def leaf_fit(ser, E, m, allow_const=False):
    """Fit ser(psi) = Sum_{j=1}^m beta_j/(E-j psi) (+ c). Return betas or None"""
    Z = len(ser) - 1
    # unknowns beta_1..beta_m (and c); equations psi^n: Sum beta_j j^n/E^(n+1) (+c delta_n0)
    unknown = m + (1 if allow_const else 0)
    rows = []
    for n in range(Z + 1):
        row = [Fr(j) ** n / E ** (n + 1) for j in range(1, m + 1)]
        if allow_const:
            row.append(Fr(1) if n == 0 else Fr(0))
        rows.append(row + [ser[n]])
    # gaussian
    piv = []
    r = 0
    for c in range(unknown):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                p = i; break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        f = rows[r][c]
        rows[r] = [x / f for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                g = rows[i][c]
                rows[i] = [x - g * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    for i in range(r, len(rows)):
        if rows[i][-1] != 0:
            return None
    sol = [Fr(0)] * unknown
    for i, c in enumerate(piv):
        sol[c] = rows[i][-1]
    return sol

def main():
    z_inf(61)
    for (e1, e2) in [(Fr(2, 3), Fr(5, 7)), (Fr(-11, 5), Fr(3, 2))]:
        E = e1 + e2
        K, Z = 7, 26
        t0 = time.time()
        sig, rho = bivar_tc(E, K, Z)
        print("E =", E, "time", round(time.time() - t0, 1))
        # poles
        for a in range(3):
            lows = sorted(set(k for n in range(Z + 1) for k, val in sig[a][n].items() if val != 0))[:3]
            print(f"  component {a}: lowest t_c powers present {lows}")
        print("  sigma^(-1) z^0:", [sig[a][0].get(-1, 0) for a in range(3)], " expected (0,E,-E^2)=", (0, E, -E * E))
        print("  sigma^(-1) z^n, n>=1 vanish:", all(sig[a][n].get(-1, Fr(0)) == 0 for a in range(3) for n in range(1, Z + 1)))
        s00 = leaf_fit(in_z_at_k(sig, 0, 0, Z), E, 2)
        s10 = leaf_fit(in_z_at_k(sig, 1, 0, Z), E, 2)
        s20 = in_z_at_k(sig, 2, 0, Z)
        print("  sigma^(0) comps: beta_0", s00, " expected [E,0]=", E, " beta_1", s10, " expected [E^2/9,0]=", E * E / 9, " comp2 zero:", all(x == 0 for x in s20))
        # L1 for sigma_0, Z-infinity
        for k in range(0, K + 1):
            ser = in_z_at_k(sig, 0, k, Z)
            betas = leaf_fit(ser, E, k + 1)
            smaller = leaf_fit(ser, E, k) if k >= 1 else None
            s = sum(b / (j + 1) for j, b in enumerate(betas)) if betas else None
            print(f"  k={k}: sigma_0^(k) in Lambda_{k+1}: {betas is not None}; in Lambda_{k}: {smaller is not None}; Sum beta_j/j = {s} (expected {E if k==0 else 0})")
        # H2(b) form for sigma_1, sigma_2
        for k in range(0, K + 1):
            for a in (1, 2):
                ser = in_z_at_k(sig, a, k, Z)
                bb = leaf_fit(ser, E, k + 2)
                bb2 = leaf_fit(ser, E, k + 1)
                print(f"  k={k}: sigma_{a}^(k) in Lambda_{k+2}: {bb is not None}, in Lambda_{k+1}: {bb2 is not None}")
    # control off shell : c1 -> c1*(1+x)
    E = Fr(29, 21)
    sig, rho = bivar_tc(E, 4, 16, offshell=lambda x: 1 + x)
    print("off-shell control c1->c1(1+x):")
    for k in range(0, 4):
        for a in range(3):
            ser = in_z_at_k(sig, a, k, 16)
            print(f"  k={k} a={a} in Lambda_{k+2}: {leaf_fit(ser, E, k + 2) is not None}; with constant: {leaf_fit(ser, E, k + 2, True) is not None}")
    return

if __name__ == "__main__":
    main()
