# T2: frame-B edges at generic CY weights: F3, symplecticity, flow, L2, L3. Evidence only.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cylib import *
from t1_rows import SAMPLES, mk


def minval(A):
    v, s = A
    for i, x in enumerate(s):
        if x:
            return v + i
    return None


def reg(A, L):
    assert all(lcoef(A, k) == 0 for k in range(A[0], 0))
    return [lcoef(A, k) for k in range(L)]


def main():
    N = 13
    K = 5
    L = K + 8
    for sm in SAMPLES:
        w = mk(sm)
        print('=== weights', [str(x) for x in w.raw])
        R = build_rows(w, N, L)
        E = w.Ea
        for name, cls in (('frame B (E,-e1,-e2)', [E, (-w.e1) % P, (-w.e2) % P]),
                          ('control class (E,-a1,-A)', [E, (-w.a1) % P, (-w.A) % P]),
                          ('control point frame', None)):
            rows, fac = framed_rows(R, cls)
            M = 8
            E11 = edge_table(R, rows, 1, 1, M)
            E12 = edge_table(R, rows, 1, 2, M)
            E21 = edge_table(R, rows, 2, 1, M)
            E22 = edge_table(R, rows, 2, 2, M)
            mv = {}
            for nm, T in (('E11', E11), ('E12', E12), ('E21', E21), ('E22', E22)):
                vals = [minval(x) for x in T.values() if minval(x) is not None]
                mv[nm] = min(vals)
            dil_min = min(minval(rows[c][0][m]) for c in (1, 2) for m in range(N) if minval(rows[c][0][m]) is not None)
            print('  ', name, 'min t-order of edges:', mv, 'dil:', dil_min)
            if name.startswith('frame B'):
                # symmetry and transposition
                sym = all(all(lcoef(E11[(a, b)], k) == lcoef(E11[(b, a)], k) for k in range(0, 5)) for (a, b) in E11)
                tr = all(all(lcoef(E21[(b, a)], k) == lcoef(E12[(a, b)], k) for k in range(0, 5)) for (a, b) in E12)
                print('   E11 symmetric', sym, ' E21(b,a)=E12(a,b)', tr)
                # K_i: <v_i(z),v_i(-z)> = Delta_i, <v_1(z),v_2(-z)>=0, orders z^0..z^N
                okK = True
                for (i, j) in ((1, 1), (2, 2), (1, 2)):
                    for m in range(N + 1):
                        tot = (0, [0] * L)
                        for p in range(m + 1):
                            q = m - p
                            s = pair(w, (rows[i][0][p], rows[i][1][p]), (rows[j][0][q], rows[j][1][q]), L)
                            tot = ladd(tot, lscal((-1) ** q % P, s), L)
                        if i == j and m == 0:
                            tot = ladd(tot, lneg(R.D1 if i == 1 else R.D2), L)
                        for k in range(-2, 5):
                            if lcoef(tot, k):
                                okK = False
                print('   <v_i(z),v_j(-z)> = delta_ij Delta_i through z^%d: %s' % (N, okK))
                # Row 2: (L_i + z D_i) r_1 = -(a1 a2/tau) r_0 + (E/tau - a0) r_1
                okR = True
                itau = (-1, R.itau_reg)
                for i, Li, Di in ((1, R.L1, R.D1), (2, R.L2, R.D2)):
                    vD, sD = Di
                    dl = lscal(inv(2), lmul(lder(Di), (-vD, sinv(sD, L)), L))
                    r0, r1 = rows[i]
                    for m in range(N + 1):
                        lhs = lmul(Li, r1[m], L)
                        if m >= 1:
                            lhs = ladd(lhs, ladd(lder(r1[m - 1]), lneg(lmul(dl, r1[m - 1], L)), L), L)
                        rhs = lscal((-w.a1 * w.a2) % P, lmul(itau, r0[m], L))
                        rhs = ladd(rhs, lmul(ladd(lscal(E, itau), (0, const(P - w.a0, L)), L), r1[m], L), L)
                        d = ladd(lhs, lneg(rhs), L)
                        for k in range(-3, 4):
                            if lcoef(d, k):
                                okR = False
                print('   Row 2 of the QDE, both colours, generic Euler operator:', okR)
                # flow: [t^k](rho^2 E_11) against (7.5) with Theta = Delta_2/(L_2-L_1)
                lead = (-w.e1 * w.e2 % P * pow(E, 3, P)) % P
                rm1 = ssqrt(scal(inv(lead), R.D1[1]), 1, L)
                rho = sinv(rm1, L)
                rho2 = mul(rho, rho, L)
                # g_k[a]: [psi^a][t^k] rho r_{0,1}(-psi) = (-1)^a [t^k] rho r_{0,1}[a]
                g = [[(-1) ** a * mul(rho, rows[1][0][a][1], L)[k] % P for a in range(N + 1)] for k in range(L)]
                # Theta
                dL = ladd(R.L2, lneg(R.L1), L)
                vd, sd = dL
                Th = lmul(R.D2, (-vd, sinv(sd, L)), L)  # series with val 1
                Theta = [lcoef(Th, k) for k in range(L - 2)]
                okF = Theta[0] == 0
                for a in range(0, 7):
                    for b in range(0, 7 - a):
                        ser = mul(rho2, reg(E11[(a, b)], L), L)
                        for k in range(0, K):
                            pred = 0
                            for k0 in range(1, k + 2):
                                for k1 in range(0, k + 2 - k0):
                                    k2 = k + 1 - k0 - k1
                                    pred += Theta[k0] * (k1 + 1) * (k2 + 1) * g[k1][a] * g[k2][b]
                            pred = pred * inv(k + 2) % P
                            if pred != ser[k]:
                                okF = False
                print('   flow (7.5) with Theta = Delta_2/(L_2-L_1), Theta_0 = 0:', okF)
                print('   Theta_1 = -Delta_2(0)/E :', Theta[1] == (-lcoef(R.D2, 0) * inv(E)) % P)
                # L3 support: [t^k](rho^2 E11) * prod_{i<=k+1}(E - i psi) polynomial in psi of degree <= k at each psi'-order
                sup = True
                M2 = 8
                for k in range(0, 3):
                    for b in range(0, 2):
                        col = [mul(rho2, reg(E11[(a, b)], L), L)[k] for a in range(M2 + 1 - b)]
                        ser = col[:]
                        for j in range(1, k + 2):
                            ser = [(E * ser[i] - (j * ser[i - 1] if i > 0 else 0)) % P for i in range(len(ser))]
                        deg = max([i for i, x in enumerate(ser) if x] + [-1])
                        if deg > k:
                            sup = False
                print('   L3 pole support through t^2 (m+n<=8):', sup)
                # L2 for j=2: [t^k](rho f_{2,b}): needs higher psi order, use edge table with bigger M
                M3 = 11
                E12big = edge_table(R, rows, 1, 2, M3)
                res = []
                for k in range(0, 4):
                    for b in range(0, 3):
                        col = [mul(rho, reg(E12big[(a, b)], L), L)[k] for a in range(M3 + 1 - b)]
                        if None in col:
                            res.append('val')
                            continue
                        ser = col[:]
                        for j in range(1, k + 2):
                            ser = [(E * ser[i] - (j * ser[i - 1] if i > 0 else 0)) % P for i in range(len(ser))]
                        res.append(max([i for i, x in enumerate(ser) if x] + [-1]))
                print('   L2 j=2: degree of [t^k](rho f_{2,b}) prod(E-j psi) (expect <= k), k=0..3,b=0..2:', res)
                # v_2(z;0) proportional to (1, a1 a2/E)
                hreg = h_series([w.a1, w.A, w.a2, w.B, (-E) % P, (-E) % P], N)
                ok2 = all(lcoef(rows[2][0][m], 0) == hreg[m] and lcoef(rows[2][1][m], 0) == hreg[m] * w.a1 % P * w.a2 % P * inv(E) % P for m in range(N - 1))
                print('   v_2(z;0) = h_(a1,A,a2,a0+a1,-E,-E)(z) (1, a1 a2/E):', ok2)
                # sigma^(0)
                l0 = lcoef(R.L1, 0)
                print('   L_1 = E/t + l_0: l_0 =', l0 == ((E * inv(2) - w.a0 * inv(2) - (2 * w.a0 * E + 4 * w.a1 * w.a2) * inv(4 * E) + 0) % P), '(E/2 - a0/2 - (a0 E/2 + a1 a2)/E)')


if __name__ == '__main__':
    main()
