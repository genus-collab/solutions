# ref5_item3_main.py -- driver of ref5_item3.py (cor:mai-natfact on actual data); see the header there.  Evidence only.
# Usage: py -I -B ref5_item3_main.py [U] [K]
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref5_item3 import p, inv, fr, ZETA, EPS, umul, uinv, uadd, rows_colour0, model_rows, frob, logcheck
from ref5_item1 import ginf

def zp(k): return pow(ZETA, k % 5, p)
def zero(U): return [0] * U
def one(U): return [1] + [0] * (U - 1)

def mmul(X, Y, U, K):                                   # X[i][j][k] = u-series of the zhat^k coefficient
    Z = [[[zero(U) for _ in range(K)] for _ in range(5)] for _ in range(5)]
    for i in range(5):
        for j in range(5):
            for l in range(5):
                for k1 in range(K):
                    if any(X[i][l][k1]):
                        for k2 in range(K - k1):
                            if any(Y[l][j][k2]):
                                Z[i][j][k1 + k2] = uadd(Z[i][j][k1 + k2], umul(X[i][l][k1], Y[l][j][k2], U))
    return Z

def trm(X):                                             # transpose and zhat -> -zhat
    K = len(X[0][0])
    return [[[[(-1) ** k * c % p for c in X[j][i][k]] for k in range(K)] for j in range(5)] for i in range(5)]

def first_fail_one(X, U, K):
    for k in range(K):
        for n in range(U):
            for i in range(5):
                for j in range(5):
                    want = 1 if (i == j and k == 0 and n == 0) else 0
                    if X[i][j][k][n] != want:
                        return (n, k)
    return None

def minv(X, U, K):                                      # inverse of 1 + O(zhat)
    N = [[[one(U) if i == j else zero(U) for j in range(5)] for i in range(5)]]
    for k in range(1, K):
        acc = [[zero(U) for _ in range(5)] for _ in range(5)]
        for j in range(1, k + 1):
            for a in range(5):
                for c in range(5):
                    if any(X[a][c][j]):
                        for b in range(5):
                            acc[a][b] = uadd(acc[a][b], umul(X[a][c][j], N[k - j][c][b], U), -1)
        N.append(acc)
    return [[[N[k][i][j] for k in range(K)] for j in range(5)] for i in range(5)]

def build(U, K, free, alpha, beta, gamma=Fr(7, 10), use_tau=True):
    Q0, C0, C1, closing = rows_colour0(U, K, free, alpha, beta, gamma)
    Qh = [[[[zp(i * (a + 1) - i * k) * c % p for c in Q0[a][k]] for k in range(K)] for i in range(5)] for a in range(5)]
    Psi = [[Qh[b][j][0] for j in range(5)] for b in range(5)]
    w = [uinv(C0, U), uinv(umul(C0, C1, U), U), umul(C0, C1, U), C0, one(U)]
    lem3A = all(Psi[a][j] == [zp(j * (a + 1)) * c % p for c in w[a]] for a in range(5) for j in range(5))
    i5 = inv(5)
    def pair(X, Y):                                     # (1/5) sum_a X_(3-a) Y_a  (eta^-1 pairing of two flat vectors of u-series)
        r = zero(U)
        for a in range(5):
            r = uadd(r, umul(X[(3 - a) % 5], Y[a], U))
        return [x * i5 % p for x in r]
    ortho = all(pair([Psi[b][j] for b in range(5)], [Psi[a][i] for a in range(5)]) == (one(U) if i == j else zero(U))
                for i in range(5) for j in range(5))
    Rt = [[[pair([Psi[b][j] for b in range(5)], [Qh[a][i][k] for a in range(5)]) for k in range(K)] for i in range(5)] for j in range(5)]
    rows, mclose = model_rows(K)
    taus = uinv([fr(Fr(4, 4 + 5 * n)) for n in range(U)], U) if use_tau else one(U)
    Rm = [[[zp(i * (b + 1) - i * k) * fr(rows[b][k]) % p for k in range(K)] for i in range(5)] for b in range(5)]
    sym_b = all(Qh[b][i][k][0] == Rm[b][i][k] for b in range(5) for i in range(5) for k in range(K))      # (b) at u = 0
    M0 = [[[sum(zp(j * (4 - a)) * EPS[(3 - a) % 5] * Rm[a][i][k] for a in range(5)) * i5 % p for k in range(K)]
           for i in range(5)] for j in range(5)]
    tk = [one(U)]
    for k in range(1, K):
        tk.append(umul(tk[-1], taus, U))
    M = [[[[M0[j][i][k] * c % p for c in tk[k]] for k in range(K)] for i in range(5)] for j in range(5)]
    return dict(Psi=Psi, Rt=Rt, M=M, lem3A=lem3A, ortho=ortho, sym_b=sym_b, closing=closing, mclose=mclose)

def poly_report(G, U, K):                               # violations of deg_zhat [u^n] G <= floor((5n+4)/4)
    bad = []
    for n in range(U):
        for k in range((5 * n + 4) // 4 + 1, K):
            if any(G[i][j][k][n] for i in range(5) for j in range(5)):
                bad.append((n, k)); break
    return bad

def natural(d, U, K, left=False):
    Minv = minv(d['M'], U, K)
    return mmul(Minv, d['Rt'], U, K) if left else mmul(d['Rt'], Minv, U, K)

def main():
    U = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    Gi = ginf(K); free = {5 * m: Gi[5 * m] for m in range(1, (K - 1) // 5 + 1)}
    print('ref5_item3: u^0..u^%d, zhat^0..zhat^%d, mod 2^61-1; free coefficients [z^5m] Ginf (D = 1)' % (U - 1, K - 1))
    f0, f1, f2 = frob(U + 2)
    print('  log solution exists (no negative powers) for gamma = 7/10: %s; for gamma = 7/10 + 1/3: %s (must be False)'
          % (not logcheck(f1, f2, Fr(7, 10), U), not logcheck(f1, f2, Fr(7, 10) + Fr(1, 3), U)))
    for (al, be) in ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7))):
        d = build(U, K, free, al, be); G = natural(d, U, K)
        print('[frame (%s, %s)] rows close (4 -> 0): %s; model T-system closes: %s; Lemma 3.A: %s; Psi^T eta^-1 Psi = 1: %s'
              % (al, be, d['closing'], d['mclose'], d['lem3A'], d['ortho']))
        g0 = all(G[i][j][0] == (one(U) if i == j else zero(U)) for i in range(5) for j in range(5))
        gu = all(G[i][j][k][0] == (1 if (i == j and k == 0) else 0) for i in range(5) for j in range(5) for k in range(K))
        print('   (b) Qhat(0, zhat) = Rhat^mod(zhat): %s; G(u,0) = 1: %s; G(0,zhat) = 1: %s' % (d['sym_b'], g0, gu))
        print('   (c) degree violations (n, k) =', poly_report(G, U, K), '(must be empty); first failure of Rt(-zh)^T Rt(zh) = 1:',
              first_fail_one(mmul(trm(d['Rt']), d['Rt'], U, K), U, K), '; of G(-zh)^T G(zh) = 1:', first_fail_one(mmul(trm(G), G, U, K), U, K))
        if al != 0:
            continue
        print('   control M_0 on the left: violations =', poly_report(natural(d, U, K, left=True), U, K)[:3], '(must be non-empty)')
        d2 = build(U, K, free, al, be, use_tau=False)
        print('   control M_0(zhat), no tau: violations =', poly_report(natural(d2, U, K), U, K)[:3], '(must be non-empty)')
        d3 = build(U, K, free, al, be, gamma=Fr(7, 10) + Fr(1, 3))
        print('   control gamma = 7/10 + 1/3: rows close: %s; violations =' % d3['closing'], poly_report(natural(d3, U, K), U, K)[:3])
        Ps = d['Psi']; i5 = inv(5)
        ou3 = True
        for a in range(5):
            for c in range(5):
                r = zero(U)
                for b in range(5):
                    r = uadd(r, umul(Ps[a][(3 - b) % 5], Ps[c][b], U))
                ou3 &= [x * i5 % p for x in r] == ([i5] + [0] * (U - 1) if (a + c) % 5 == 3 else zero(U))
        print('   control: OU3 relation Psi eta^-1 Psi^T = eta^-1 holds: %s (must be False)' % ou3)
        fx = {k: (v + Fr(1, 7) if k == 5 else v) for k, v in free.items()}
        d5 = build(U, K, fx, al, be); G5 = natural(d5, U, K)
        gu5 = all(G5[i][j][k][0] == (1 if (i == j and k == 0) else 0) for i in range(5) for j in range(5) for k in range(K))
        print('   control f_1 + 1/7 (D = 1 + d z^5): violations =', poly_report(G5, U, K), '; G(0,zhat) = 1:', gu5,
              '; first failure of G(-zh)^T G(zh) = 1 at (u^n, zh^k) =', first_fail_one(mmul(trm(G5), G5, U, K), U, K))

if __name__ == '__main__':
    main()
