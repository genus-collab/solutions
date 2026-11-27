# Independent check of eq:J=I and the QDE at generic alpha on T_max, by exact fractions.
# J side: direct graph sum of virtual localisation on Mbar_{0,n}(Z,d), Z = K_P2 x C^2.
# I side: BS eq:Ifun at H = -alpha_i.  Finite-order evidence only.
from fractions import Fraction as Fr
from math import factorial as fac
import sys
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4


class T:
    def __init__(s, al, nu):
        s.al = [Fr(a) for a in al]
        s.nu = [Fr(n) for n in nu]
        a = s.al
        s.W = [(a[(i+1) % 3]-a[i])*(a[(i+2) % 3]-a[i])*3*a[i]*s.nu[0]*s.nu[1] for i in range(3)]
        s.memo = {}
        s.ins = None

    def edge(s, p, c, e):
        a = s.al
        om = (a[c]-a[p])/e
        num = Fr(1)
        for m in range(1, 3*e):
            num *= 3*a[p] + m*om
        den = s.nu[0]*s.nu[1]
        for k in range(3):
            for m in range(e+1):
                if (k == p and m == 0) or (k == c and m == e):
                    continue
                den *= a[k]-a[p]-m*om
        return num/den/e

    def tuples(s, c, D, m):
        """ordered tuples of branches out of colour c, total degree D, m extra markings"""
        if D == 0:
            if m == 0:
                yield []
            return
        for c2 in range(3):
            if c2 == c:
                continue
            for e in range(1, D+1):
                for D1 in range(e, D+1):
                    for m1 in range(m+1):
                        v = s.branch(c, e, c2, D1, m1)
                        w = (s.al[c2]-s.al[c])/e
                        for rest in s.tuples(c, D-D1, m-m1):
                            yield [(w, v)] + rest

    def branch(s, p, e, c, D, m):
        key = (p, e, c, D, m, s.ins)
        if key in s.memo:
            return s.memo[key]
        om = (s.al[p]-s.al[c])/e
        tot = Fr(0)
        for mh in range(m+1):
            for tup in s.tuples(c, D-e, m-mh):
                k = len(tup)
                n = 1+k+mh
                val = s.W[c]**k / fac(k) * (1/om)
                sm = 1/om
                for w, v in tup:
                    val *= v/w
                    sm += 1/w
                val *= sm**(n-3)
                if mh:
                    val *= (-s.al[c])**s.ins
                tot += val
        res = s.edge(p, c, e)*tot
        s.memo[key] = res
        return res

    def root(s, i, d, z=None, a=None):
        """<[H^a,] phi_i/(z-psi)>_{0,n,d}; phi_i restricts to W_i at P_i and 0 elsewhere.
        z=None means the psi-free insertion phi_i."""
        s.ins = a
        M = 0 if a is None else 1
        tot = Fr(0)
        for mh in range(M+1):
            for tup in s.tuples(i, d, M-mh):
                k = len(tup)
                n = 1+k+mh
                val = s.W[i]**k/fac(k)
                sm = Fr(0)
                for w, v in tup:
                    val *= v/w
                    sm += 1/w
                if z is not None:
                    val *= 1/z
                    sm += 1/z
                val *= sm**(n-3)
                if mh:
                    val *= (-s.al[i])**a
                tot += val
        return tot

    def eta(s, f, g):
        return sum(f(i)*g(i)/s.W[i] for i in range(3))

    def corr1(s, b, d):
        return sum((-s.al[i])**b/s.W[i]*s.root(i, d) for i in range(3))

    def corr2(s, a, b, d):
        return sum((-s.al[i])**b/s.W[i]*s.root(i, d, a=a) for i in range(3))


def smul(f, g, N=DMAX):
    h = [Fr(0)]*(N+1)
    for i, x in enumerate(f):
        if x == 0:
            continue
        for j, y in enumerate(g):
            if i+j <= N:
                h[i+j] += x*y
    return h


def sexp(f, N=DMAX):
    res = [Fr(1)]+[Fr(0)]*N
    term = res[:]
    for k in range(1, N+1):
        term = [t/k for t in smul(term, f, N)]
        res = [r+t for r, t in zip(res, term)]
    return res


def sexp_c(f, c, N=DMAX):
    return sexp([c*x for x in f], N)


g = [Fr(0)]+[3*Fr(fac(3*d-1), fac(d)**3)*(-1)**d for d in range(1, DMAX+1)]


def I_coeffs(t, i, z):
    a = t.al
    H = -a[i]
    out = []
    for d in range(DMAX+1):
        num = Fr(1)
        for m in range(3*d):
            num *= -3*H-m*z
        den = Fr(1)
        for j in range(3):
            for m in range(1, d+1):
                den *= H+a[j]+m*z
        out.append(z*num/den)
    return out


def check_JI(t, z, g_use=g, label=''):
    ok = True
    for i in range(3):
        H = -t.al[i]
        I = I_coeffs(t, i, z)
        J = [Fr(0)]*(DMAX+1)
        J[0] = z
        for d in range(1, DMAX+1):
            Jd = t.root(i, d, z=z)
            s = sexp_c(g_use, Fr(d))
            for n in range(DMAX+1-d):
                J[d+n] += Jd*s[n]
        J = smul(sexp_c(g_use, H/z), J)
        if any(x != y for x, y in zip(I, J)):
            ok = False
        if label and i == 0:
            print('    ', label, 'I and J at P_0, y^1:', I[1], J[1])
    return ok


def Ian_check(t):
    """1/z coefficient: Sum_d Qcal^d <phi_i>_{0,1,d} against eq:Ian"""
    e1 = sum(t.al)
    ok = True
    I11 = [Fr(0)]*(DMAX+1)
    I21 = [Fr(0)]*(DMAX+1)
    for d in range(1, DMAX+1):
        Hd = sum(Fr(1, k) for k in range(1, d+1))
        Hdiff = sum(Fr(1, k) for k in range(d+1, 3*d))
        c = Fr(fac(3*d-1), fac(d)**3)*(-1)**d
        I11[d] = -3*e1*c*Hd
        I21[d] = 9*c*Hdiff
    g2 = smul(g, g)
    for i in range(3):
        H = -t.al[i]
        lhs = [Fr(0)]*(DMAX+1)
        for d in range(1, DMAX+1):
            s = sexp_c(g, Fr(d))
            v = t.root(i, d)
            for n in range(DMAX+1-d):
                lhs[d+n] += v*s[n]
        rhs = [I11[k]*H+(I21[k]-g2[k]/2)*H*H for k in range(DMAX+1)]
        if lhs != rhs:
            ok = False
    return ok


def solve3(A, b):
    n = 3
    M = [list(A[r])+[b[r]] for r in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        M[c] = [x/M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x-f*y for x, y in zip(M[r], M[c])]
    return [M[r][n] for r in range(n)]


def qde_check(t, z):
    """z theta_Q S(H^a) = Sum_c C^c_{1a} S(H^c), a = 0,1,2, at each P_i, with
    S(H^a)|_i = Q^{H/z}(H^a + Sum_d Q^d <H^a, phi_i/(z-psi)>_{0,2,d}) from two-point
    localisation, and C from three-point invariants through the divisor equation."""
    al = t.al
    ok = True
    G = [[t.eta(lambda i: (-al[i])**a, lambda i: (-al[i])**b) for b in range(3)] for a in range(3)]
    C = {}
    for a in range(3):
        for d in range(DMAX+1):
            rhsv = []
            for b in range(3):
                if d == 0:
                    rhsv.append(t.eta(lambda i: (-al[i])**(a+1), lambda i: (-al[i])**b))
                else:
                    rhsv.append(d*t.corr2(a, b, d))
            C[(a, d)] = solve3(G, rhsv)
    for i in range(3):
        H = -al[i]
        S = {}
        for a in range(3):
            S[a] = [H**a]+[t.root(i, d, z=z, a=a) for d in range(1, DMAX+1)]
        for a in range(3):
            lhs = [(H+z*d)*S[a][d] for d in range(DMAX+1)]
            rhs = [Fr(0)]*(DMAX+1)
            for d in range(DMAX+1):
                for n in range(DMAX+1-d):
                    for c in range(3):
                        rhs[d+n] += C[(a, d)][c]*S[c][n]
            if lhs != rhs:
                ok = False
    return ok


def SJ_check(t, z):
    """S(1)=J/z and S(H)=theta_Q J: two-point localisation against the one-point sum"""
    ok = True
    for i in range(3):
        H = -t.al[i]
        for d in range(1, DMAX+1):
            Jd = t.root(i, d, z=z)
            if t.root(i, d, z=z, a=0) != Jd/z:
                ok = False
            if t.root(i, d, z=z, a=1) != (H/z+d)*Jd:
                ok = False
    return ok


tests = [(T(['2/7', '5/11', '-3/13'], ['7/5', '-9/4']), Fr(3, 17)),
         (T(['1/3', '-4/5', '11/2'], ['-2/9', '5/6']), Fr(-7, 19)),
         # CY subtorus: e1(alpha) = eps1 + eps2, nu = -eps
         (T(['1/7', '3/5', '11/13'], ['-1/3', str(Fr(-1,7)+Fr(-3,5)+Fr(-11,13)+Fr(1,3))]), Fr(5, 23))]
for t, z in tests:
    print('alpha', [str(a) for a in t.al], 'nu=-eps', [str(n) for n in t.nu], 'z', z)
    a = t.al
    e1 = sum(a)
    e3 = a[0]*a[1]*a[2]
    f = 1/(3*t.nu[0]*t.nu[1])
    B = [[f/e3, 0, 0], [0, 0, -f], [0, -f, f*e1]]
    print('  eta_ab against BS eq:etaflat:',
          all(t.eta(lambda i: (-a[i])**p, lambda i: (-a[i])**q) == B[p][q] for p in range(3) for q in range(3)))
    Nd = [t.corr1(1, d)/d*t.nu[0]*t.nu[1] for d in range(1, DMAX+1)]
    print('  GW_{0,d} eps1 eps2 = <H>_{0,1,d} eps1 eps2 / d:', [str(x) for x in Nd])
    print('  <H^2>_{0,1,d} eps1 eps2 / e1(alpha):',
          [str(t.corr1(2, d)*t.nu[0]*t.nu[1]/e1) for d in range(1, DMAX+1)])
    print('  J = I through y^%d at P_0,P_1,P_2:' % DMAX, check_JI(t, z, label='d=1'))
    print('  eq:Ian 1/z coefficients against localisation:', Ian_check(t))
    print('  control, mirror map dropped (must be False):', check_JI(t, z, g_use=[Fr(0)]*(DMAX+1)))
    print('  S(1) = J/z, S(H) = theta_Q J from two-point sums:', SJ_check(t, z))
    print('  QDE z theta_Q S(H^a) = Sum C^c_{1a} S(H^c), a=0,1,2 through Q^%d:' % DMAX, qde_check(t, z))
