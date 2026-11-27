"""F3, L2, L3 at low t_c-orders, exact (evidence only), fresh code.  Also Delta_j(0), Dil^_j analytic.
Controls: off-shell c1 gives the pole quoted in s.3.3; L3 with a smaller span fails."""
import time
from fractions import Fraction as Fr
from rows_core import *
from check_conifold import Ser, powser_1mx, setup, bivar_tc, leaf_fit

def edge_from_num(N, D):
    """N: dict (m,n)->value ; returns Ed(m,n)=Sum_k (-1)^k N(m+k+1,n-k), m+n<=D-1"""
    out = {}
    for m in range(D):
        for n in range(D - m):
            s = 0
            for k in range(n + 1):
                s = s + ((-1) ** k) * N.get((m + k + 1, n - k), 0)
            out[(m, n)] = s
    return out

def solve_span(data, basis):
    """data: dict key->value; basis: list of dict key->value ; least: exact solve; return sol or None"""
    keys = list(data.keys())
    rows = [[b.get(key, Fr(0)) for b in basis] + [data[key]] for key in keys]
    r = 0
    piv = []
    nb = len(basis)
    for c in range(nb):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                p = i
                break
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
    sol = [Fr(0)] * nb
    for i, c in enumerate(piv):
        sol[c] = rows[i][-1]
    return sol

def test_L3(e1, e2, K, D, offshell=None):
    E = e1 + e2
    sig, rho = bivar_tc(E, K + 2, D + 1, offshell)
    n00 = e1 * e2 * E ** 3 / 9
    n11 = -3 * e1 * e2 * E
    n12 = -3 * e1 * e2
    # sigma at z=-psi: coefficient of psi^m is (-1)^m [z^m]
    def s(a, k, m):
        return sig[a][m].get(k, Fr(0)) * ((-1) ** m)
    res = {}
    for tot in range(-2, K + 1):
        N = {}
        for m in range(D + 1):
            for n in range(D + 1 - m):
                val = Fr(0)
                for k1 in range(-1, tot + 2):
                    k2 = tot - k1
                    val += n00 * s(0, k1, m) * s(0, k2, n)
                    val += n11 * s(1, k1, m) * s(1, k2, n)
                    val += n12 * (s(1, k1, m) * s(2, k2, n) + s(2, k1, m) * s(1, k2, n))
                N[(m, n)] = -val
        # P(-w,w) coefficient check: Sum_{m} (-1)^m N(m, n-m) for n>0 must vanish
        anti = [sum(((-1) ** m) * N[(m, n - m)] for m in range(n + 1)) for n in range(D + 1)]
        res[tot] = (edge_from_num(N, D), anti)
    return res, E

def lam_basis(E, D, i, l):
    b = {}
    for m in range(D):
        for n in range(D - m):
            b[(m, n)] = Fr(i) ** m / E ** (m + 1) * Fr(l) ** n / E ** (n + 1)
    return b

def main():
    t0 = time.time()
    for (e1, e2) in [(Fr(2, 3), Fr(5, 7)), (Fr(-11, 5), Fr(3, 2))]:
        K, D = 6, 16
        res, E = test_L3(e1, e2, K, D)
        c0 = 3 * e1 * e2 * E ** 3
        print("weights", e1, e2)
        for tot in sorted(res):
            ed, anti = res[tot]
            antiok = {tot: anti[0]} if tot == -2 else None
            nz = [i for i, x in enumerate(anti) if x != 0]
            if tot < 0:
                print(f"  t_c^{tot}: rho^2 Ehat_00 all zero: {all(v == 0 for v in ed.values())}; N(-w,w) nonzero at w^{nz}, value {anti[0]} (K=1: -c0 = {-c0} at t_c^-2)")
                continue
            basis = []
            idx = []
            for i in range(1, tot + 2):
                for l in range(1, tot + 2):
                    if (i - 1) + (l - 1) <= tot:
                        basis.append(lam_basis(E, D, i, l))
                        idx.append((i, l))
            sol = solve_span(ed, basis)
            # control: smaller span
            basis2 = [b for b, (i, l) in zip(basis, idx) if (i - 1) + (l - 1) < tot]
            sol2 = solve_span(ed, basis2) if basis2 else None
            sym = all(ed[(m, n)] == ed[(n, m)] for (m, n) in ed if (n, m) in ed)
            print(f"  t_c^{tot}: N(-w,w) zero: {nz==[]}; L3 span OK: {sol is not None}; smaller span OK: {sol2 is not None}; symmetric: {sym}")
            if tot == 0 and sol is not None:
                print("     leading coefficient:", sol[0], " expected rho^2 Ehat(0) = -e1e2E^4/27 =", -e1 * e2 * E ** 4 / 27)
        print("  time", round(time.time() - t0, 1))
    # off-shell control
    e1, e2 = Fr(2, 3), Fr(5, 7)
    res, E = test_L3(e1, e2, 3, 12, offshell=lambda x: 1 + x)
    ed, anti = res[-1]
    print("off-shell c1->c1(1+x): t_c^-1 of rho^2 Ehat_00 at psi^0 psi'^0 =", ed[(0, 0)], " s.3.3 predicts 3 e1 e2 E^2 (XX(0)+2/9) with XX(0)+2/9=-1:", -3 * e1 * e2 * E ** 2)
    print("   constant in psi,psi' ?", all(v == 0 for k, v in ed.items() if k != (0, 0)))
    ed0, anti0 = res[-2]
    print("   t_c^-2 vanishes:", all(v == 0 for v in ed0.values()))
    ed1, _ = res[0]
    print("   L3 span at t_c^0 off shell:", solve_span(ed1, [lam_basis(E, 12, 1, 1)]) is not None)

if __name__ == "__main__":
    main()
