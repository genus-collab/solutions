# T4: low-order hand values at generic CY weights. Evidence only.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cylib import *
from t1_rows import SAMPLES, mk


def main():
    N = 10
    L = 12
    for sm in SAMPLES:
        w = mk(sm)
        R = build_rows(w, N, L)
        E = w.Ea
        u = (E - 2 * w.a1) * (E - 2 * w.a2) % P
        rows, fac = framed_rows(R, [E, (-w.e1) % P, (-w.e2) % P])
        lead = (-w.e1 * w.e2 % P * pow(E, 3, P)) % P
        rm1 = ssqrt(scal(inv(lead), R.D1[1]), 1, L)
        rho = sinv(rm1, L)
        print('weights', [str(x) for x in w.raw])
        print('  Delta_1 = -e1 e2 E^3 t^-2 (1 - u t/E^2 + ...):', lcoef(R.D1, -1) == lead * (-u * inv(E * E)) % P % P)
        print('  rho_1 = u/(2 E^2):', rho[1] == u * inv(2 * E * E) % P)
        print("  Delta_2'(0) = 0 (fibrewise accident):", lcoef(R.D2, 1) == 0)
        # q^(1) = (u/2)(lam_1 - lam_2): [t^1] z G at z^(m+1): (u/2)(1/(E+z) - 1/(E+2z))
        G1 = [mul(rho, rows[1][0][m][1], L)[1] for m in range(N + 1)]
        ok = True
        iE = inv(E)
        for m in range(N):
            # zG coefficient z^(m+1) = G1[m]
            pred = (u * inv(2)) % P * iE % P * (pow((-iE) % P, m + 1, P) - pow((-2 * iE) % P, m + 1, P)) % P
            ok = ok and pred == G1[m]
        print('  [t^1](rho q_t) = (u/2)(lam_1 - lam_2):', ok)
        # sigma^(0) = (E lam_1, a1 a2 lam_1): [t^0] rho r_{1,1}(z) = a1 a2/(E+z)
        r1 = rows[1][1]
        ok = True
        for m in range(N):
            ser = lmul((0, rho), r1[m], L)
            pred = w.a1 * w.a2 % P * iE % P * pow((-iE) % P, m, P) % P
            ok = ok and lcoef(ser, 0) == pred
            if m > 0:
                ok = ok and lcoef(ser, -1) == 0
        ok = ok and lcoef(lmul((0, rho), r1[0], L), -1) == E
        print('  rho r_{1,1} = E/t + a1 a2/(E+z) + O(t):', ok)
        # E_11 at t=0 = -(e1 e2 a1 a2 A B/2) lam_1 lam_1'
        E11 = edge_table(R, rows, 1, 1, 5)
        co = (-w.e1 * w.e2 % P * w.a1 % P * w.a2 % P * w.A % P * w.B % P * inv(2)) % P
        ok = all(lcoef(E11[(a, b)], 0) == co * pow(iE, a + b + 2, P) % P for (a, b) in E11)
        print('  E_11|t=0 = -(e1 e2 a1 a2 (a0+a1)(a0+a2)/2) lam_1 lam_1\':', ok)
        # E_12 at t=0 = lam_1(psi) chi(w)
        E12 = edge_table(R, rows, 1, 2, 5)
        ok = all(lcoef(E12[(a, b)], 0) == lcoef(E12[(0, b)], 0) * pow(iE, a, P) % P for (a, b) in E12 if a + b <= 5)
        print('  E_12|t=0 rank one with lam_1 at the singular end:', ok)
        h2 = h_series([(-w.a1) % P, (-w.A) % P, (-w.a2) % P, (-w.B) % P, E, E], 5)
        print()


if __name__ == '__main__':
    main()
