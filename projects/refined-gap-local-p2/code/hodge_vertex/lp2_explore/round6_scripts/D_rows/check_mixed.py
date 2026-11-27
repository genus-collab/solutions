"""(F3), (L2) for mixed edges in Q(omega), low t_c-orders; regular-colour factors analytic; K=1 mixed.  Evidence only."""
from fractions import Fraction as Fr
from rows_core import *
from check_conifold import Ser, powser_1mx, setup, bivar_tc, leaf_fit
from check_edges import edge_from_num

def regular_rows(E, j, K, Z):
    Ser.M = K + 12
    x, wc, tc, theta, y = setup()
    xt = Ser([Fr(0), Fr(1)]); T = Ser([Fr(0), Fr(1)])
    for _ in range(Ser.M + 1):
        xt = T - (tc - x).compose(xt)
    om_inv = OMEGA ** (-j)
    p = powser_1mx(Fr(-1, 3))
    vjx = Ser([Fr(3) / E * (1 if n == 0 else 0) - om_inv * (3 * p.c[n] / E) for n in range(Ser.M)])
    v = vjx.compose(xt)
    assert not Qw._c(v.c[0]).iszero()
    C1i = (-wc).inv().compose(xt)
    XX = ((theta(-wc)) * (-wc).inv()).compose(xt)
    d2 = (((-wc).inv() * (-wc).inv() - 1) * (E / 2)).shift(-1).compose(xt)
    a, b, Thb = abc_vs(E, Z)
    # v powers
    vinv = v.inv()
    cache = {0: Ser.const(1)}
    def vp(e):
        if e not in cache:
            cache[e] = vp(e - 1) * v if e > 0 else vp(e + 1) * vinv
        return cache[e]
    def ev(lp):
        out = Ser.const(0)
        for e, c in lp.items():
            out = out + vp(e) * c
        return out
    rows = [[None] * (Z + 1) for _ in range(3)]
    xfac = xt  # Delta_c as series in t_c
    C1inv = C1i.inv()
    for n in range(Z + 1):
        rows[0][n] = ev(lp_shift(a[n], -1))
        rows[1][n] = ev(lp_shift(b[n], -1)) * C1i
        inner = ev(lp_shift(Thb[n], -1)) - ev(lp_shift(b[n], -1)) * d2
        if n >= 1:
            inner = inner - ev(lp_shift(b[n - 1], -1)) * XX
        rows[2][n] = inner * xfac * C1inv
    return rows

def main():
    e1, e2 = Fr(2, 3), Fr(5, 7); E = e1 + e2
    K, D = 5, 14
    n00 = e1 * e2 * E ** 3 / 9; n11 = -3 * e1 * e2 * E; n12 = -3 * e1 * e2
    for j in (1, 2):
        sig, rho = bivar_tc(E, K + 2, D + 1)
        rows = regular_rows(E, j, K + 2, D + 1)
        def s(a, k, m):
            return sig[a][m].get(k, Fr(0)) * ((-1) ** m)
        def r(a, k, n):
            if k < 0 or k >= Ser.M: return Fr(0)
            return rows[a][n].c[k] * ((-1) ** n)
        for tot in range(-1, K + 1):
            N = {}
            for m in range(D + 1):
                for n in range(D + 1 - m):
                    val = 0
                    for k1 in range(-1, tot + 1):
                        k2 = tot - k1
                        val = val + n00 * s(0, k1, m) * r(0, k2, n) + n11 * s(1, k1, m) * r(1, k2, n) + n12 * (s(1, k1, m) * r(2, k2, n) + s(2, k1, m) * r(1, k2, n))
                    N[(m, n)] = -val
            anti = [sum(((-1) ** m) * N[(m, n - m)] for m in range(n + 1)) for n in range(D + 1)]
            antizero = all(Qw._c(x).iszero() for x in anti)
            ed = edge_from_num(N, D)
            if tot < 0:
                print(f"j={j} t_c^-1: rho Ehat_0j all zero {all(Qw._c(v).iszero() for v in ed.values())}; N(-w,w)=0: {antizero}")
                continue
            res = []
            for b in range(0, 4):
                ser = [ed[(m, b)] for m in range(0, D - b)]
                fit = leaf_fit(ser, E, tot + 1)
                small = leaf_fit(ser, E, tot) if tot >= 1 else None
                res.append((fit is not None, small is not None))
            print(f"j={j} t_c^{tot}: N(-w,w)=0 {antizero}; [w^b] in Lambda_{tot+1}, smaller span, b=0..3: {res}")

if __name__ == "__main__":
    main()
