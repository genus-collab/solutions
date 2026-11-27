# T0: hand identities at generic weights, exact rationals at random points. Evidence only.
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cylib import weights, build_rows, framed_rows, lcoef, h_series, P, inv, mul, const
random.seed(20261102)


def rnd():
    return Fr(random.randint(-30, 30) or 7, random.randint(1, 17))


def amplitude_identity():
    """generic (3.16): tau(1-tau) Dam^2 am0 + (gamma-(alpha+beta+1)tau) Dam am0 - alpha beta am0 = -alpha (kE' + f' kE/z)"""
    ok = True
    for _ in range(50):
        a0, a1, a2, z = rnd(), rnd(), rnd(), rnd()
        s, t = rnd(), rnd()
        A = a0 + a2
        B = a0 + a1
        E = a0 + a1 + a2
        al = 1 + a2 / z
        be = 1 + A / z
        ga = 2 + E / z
        try:
            am0 = 1 / (1 - t * s)
            D1 = al * s * am0 ** 2
            D2 = al * (1 + al) * s ** 2 * am0 ** 3
            lhs = t * (1 - t) * D2 + (ga - (al + be + 1) * t) * D1 - al * be * am0
            kE = s * (1 - s) * am0 ** 2
            kEp = (1 - 2 * s) * am0 ** 2 + 2 * t * s * (1 - s) * am0 ** 3
            fp_ = A / s - a1 / (1 - s) + a2 * t / (1 - t * s)
            rhs = -al * (kEp + fp_ * kE / z)
            # critical polynomial
            cr = B * t * s ** 2 - (E + a0 * t) * s + A
            ok2 = fp_ * s * (1 - s) * (1 - t * s) == cr
        except ZeroDivisionError:
            continue
        ok = ok and lhs == rhs and ok2
    return ok


def operator_identity():
    """generic (4.4) on tau^n: (Dd+omega)^2 - (E/tau - a0)(Dd+omega) + a1 a2/tau = ((1-tau) z^2/tau) Gauss"""
    ok = True
    for _ in range(50):
        a0, a1, a2, z = rnd(), rnd(), rnd(), rnd()
        t = rnd()
        if t in (0, 1):
            continue
        A = a0 + a2
        E = a0 + a1 + a2
        if a2 + z == 0:
            continue
        for n in (-2, -1, 0, 1, 2, 5):
            f = t ** n
            ft = n * t ** (n - 1)
            ftt = n * (n - 1) * t ** (n - 2)
            om = (E + z) / t - (A + z)
            omt = -(E + z) / t ** 2
            # Dd f = z (1-t) f_t
            Df = z * (1 - t) * ft
            DDf = z * (1 - t) * (-z * ft + z * (1 - t) * ftt)
            lhs = DDf + 2 * om * Df + z * (1 - t) * omt * f + om * om * f - (E / t - a0) * (Df + om * f) + a1 * a2 / t * f
            al, be, ga = 1 + a2 / z, 1 + A / z, 2 + E / z
            gauss = t * (1 - t) * ftt + (ga - (al + be + 1) * t) * ft - al * be * f
            rhs = (1 - t) * z ** 2 / t * gauss
            ok = ok and lhs == rhs
    return ok


def frobenius():
    ok = True
    for _ in range(50):
        a0, a1, a2, e1 = rnd(), rnd(), rnd(), rnd()
        E = a0 + a1 + a2
        e2 = E - e1
        if 0 in (e1, e2, a1, a2, a0 + a1, a0 + a2, E):
            continue
        eta = [[E / (a1 * (a0 + a1) * a2 * (a0 + a2) * e1 * e2), 1 / ((a0 + a1) * (a0 + a2) * e1 * e2)],
               [1 / ((a0 + a1) * (a0 + a2) * e1 * e2), -a0 / ((a0 + a1) * (a0 + a2) * e1 * e2)]]
        ei = [[e1 * e2 * a0 * a1 * a2, e1 * e2 * a1 * a2], [e1 * e2 * a1 * a2, -e1 * e2 * E]]
        pr = [[sum(eta[i][k] * ei[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        ok = ok and pr == [[1, 0], [0, 1]]
        u = (E - 2 * a1) * (E - 2 * a2)
        ok = ok and 2 * a0 * E + 4 * a1 * a2 == u + E ** 2
        ok = ok and (a1 - a2) ** 2 == a0 ** 2 - u
        # discriminant of Disc in Q
        b = 2 * a0 * (a1 + a2) + 4 * a1 * a2
        c = (a1 - a2) ** 2
        ok = ok and b ** 2 - 4 * a0 ** 2 * c == 16 * a1 * a2 * (a0 + a1) * (a0 + a2)
    return ok


def scaling_a0zero():
    """a0 = 0: frame-B rows = fibrewise rows at x = a1/a2, z -> z/a2"""
    N = 8
    L = 8
    a1, a2 = Fr(5, 3), Fr(2, 7)
    w = weights(0, a1, a2, Fr(1, 3), a1 + a2 - Fr(1, 3))
    R = build_rows(w, N, L)
    rows, _ = framed_rows(R, [w.Ea, (-w.e1) % P, (-w.e2) % P])
    x = a1 / a2
    wf = weights(0, x, 1, x, 1)
    Rf = build_rows(wf, N, L)
    rowsf, _ = framed_rows(Rf, [wf.Ea, (-wf.e1) % P, (-wf.e2) % P])
    ia = inv(w.a2)
    ok = True
    for c in (1, 2):
        for m in range(N + 1):
            for k in range(0, 4):
                ok = ok and lcoef(rows[c][0][m], k) == lcoef(rowsf[c][0][m], k) * pow(ia, m, P) % P
    return ok


def gamma_connection():
    """tails: h_{wt_+} h_{(-a0,E,-a1,-a2)} = h_{(E,-e1,-e2)}; h_{wt_-} h_{(a0,E,-a0-a1,-a0-a2)} = same;
       regular solution: h_{wt_+} h_{(-a0,-E,a0+a1,a0+a2)} = h_nw h_{(a2,-E,a0+a1)}"""
    N = 12
    w = weights(Fr(3, 5), Fr(-7, 3), Fr(11, 2), Fr(5, 7), Fr(3, 5) + Fr(-7, 3) + Fr(11, 2) - Fr(5, 7))
    E = w.Ea
    a0, a1, a2, e1, e2, A, B = w.a0, w.a1, w.a2, w.e1, w.e2, w.A, w.B
    n = lambda x: (-x) % P
    t = h_series([E, n(e1), n(e2)], N)
    c1 = h_series([a0, a1, a2, n(e1), n(e2), n(a0), E, n(a1), n(a2)], N) == t
    c2 = h_series([n(a0), B, A, n(e1), n(e2), a0, E, n(B), n(A)], N) == t
    c3 = h_series([a0, a1, a2, n(e1), n(e2), n(a0), n(E), B, A], N) == h_series([a1, A, n(e1), n(e2), a2, n(E), B], N)
    return c1, c2, c3


print('amplitude identity (generic 3.16):', amplitude_identity())
print('operator identity (generic 4.4):', operator_identity())
print('eta^-1, u, Disc discriminant:', frobenius())
print('a0 = 0 frame-B rows are the fibrewise rows at x = a1/a2, z/a2:', scaling_a0zero())
print('Gamma-connection tails:', gamma_connection())
