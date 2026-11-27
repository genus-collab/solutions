# Independent refutation test (evidence only) of Lemmas 2.3, 2.4, 2.5(i), the y=0 value and pi(Phi)=Theta(-z/E) S.
# Arithmetic mod a large prime. Phi is computed directly from the jets of Lemma 2.3, with coefficients
# truncated series in delta (Y = Y0 + delta), so that theta acts through theta Y = -Y prod q / W.
import sys, math
from fractions import Fraction as Fr
p = (1 << 89) - 1
def inv(a): return pow(a % p, p - 2, p)
def fr(a, b=1): return a % p * inv(b) % p
ORD = 5
class TS:
    __slots__ = ('c',)
    def __init__(s, c=None):
        s.c = c if c is not None else [0] * ORD
    @staticmethod
    def const(a):
        c = [0] * ORD; c[0] = a % p; return TS(c)
    def __add__(s, o):
        o = o if isinstance(o, TS) else TS.const(o)
        return TS([(a + b) % p for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return TS([-a % p for a in s.c])
    def __sub__(s, o):
        o = o if isinstance(o, TS) else TS.const(o)
        return TS([(a - b) % p for a, b in zip(s.c, o.c)])
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if not isinstance(o, TS):
            o = o % p; return TS([a * o % p for a in s.c])
        r = [0] * ORD
        for i, a in enumerate(s.c):
            if a:
                for j in range(ORD - i):
                    r[i + j] += a * o.c[j]
        return TS([x % p for x in r])
    __rmul__ = __mul__
    def inv(s):
        a0 = inv(s.c[0]); r = [0] * ORD; r[0] = a0
        for n in range(1, ORD):
            t = 0
            for i in range(1, n + 1): t += s.c[i] * r[n - i]
            r[n] = (-t * a0) % p
        return TS(r)
    def deriv(s):
        r = [0] * ORD
        for i in range(1, ORD): r[i - 1] = i * s.c[i] % p
        return TS(r)
def var(Y0):
    c = [0] * ORD; c[0] = Y0 % p; c[1] = 1; return TS(c)
def monos(deg):
    for a in range(deg + 1):
        for b in range(deg + 1 - a):
            yield (a, b, deg - a - b)
def pmul(f, g, maxd):
    r = {}
    for m1, c1 in f.items():
        d1 = sum(m1)
        for m2, c2 in g.items():
            if d1 + sum(m2) > maxd: continue
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            x = c1 * c2
            r[m] = r[m] + x if m in r else x
    return r
def lap(f, Hinv):
    r = {}
    for m, c in f.items():
        for a in range(3):
            for b in range(3):
                mm = list(m); fac = 1
                if mm[a] == 0: continue
                fac *= mm[a]; mm[a] -= 1
                if mm[b] == 0: continue
                fac *= mm[b]; mm[b] -= 1
                mt = tuple(mm)
                x = c * (Hinv[a][b] * fac)
                r[mt] = r[mt] + x if mt in r else x
    return r
def phi_coeffs(alpha, Y, N):
    q = [Y + a for a in alpha]
    W = (sum(alpha) % p) * (Y * Y) + (2 * (alpha[0] * alpha[1] + alpha[0] * alpha[2] + alpha[1] * alpha[2]) % p) * Y + (3 * alpha[0] * alpha[1] * alpha[2] % p)
    third = fr(1, 3)
    HW = [[(q[a] if a == b else TS.const(0)) - Y * third for b in range(3)] for a in range(3)]
    def cof(a, b):
        i = [x for x in range(3) if x != a]; j = [x for x in range(3) if x != b]
        d = HW[i[0]][j[0]] * HW[i[1]][j[1]] - HW[i[0]][j[1]] * HW[i[1]][j[0]]
        return d if (a + b) % 2 == 0 else -d
    det = W * third
    dinv = det.inv()
    Hinv_W = [[cof(b, a) * dinv for b in range(3)] for a in range(3)]
    Hinv = [[-Hinv_W[a][b] for b in range(3)] for a in range(3)]  # inverse Hessian of phi = -W
    maxd = 2 * N + 2
    Wng = {}
    for d in range(3, maxd + 1):
        for m in monos(d):
            mf = math.factorial(m[0]) * math.factorial(m[1]) * math.factorial(m[2])
            coeff = (-Y) * fr(1, pow(3, d - 1) * mf)
            for j in range(3):
                if m[j] == d:
                    coeff = coeff + q[j] * fr(1, math.factorial(d))
            Wng[m] = coeff
    powers = [{(0, 0, 0): TS.const(1)}]
    Kmax = 2 * N
    for k in range(1, Kmax + 1):
        powers.append(pmul(powers[-1], Wng, 2 * (N + k)))
    out = [TS.const(0) for _ in range(N + 1)]
    for k in range(0, Kmax + 1):
        cur = powers[k]
        mmax = N + k
        for m in range(0, mmax + 1):
            n = m - k
            if 0 <= n <= N and 3 * k <= 2 * m:
                c = cur.get((0, 0, 0))
                if c is not None:
                    out[n] = out[n] + c * fr(1, pow(2, m) * math.factorial(m) * math.factorial(k))
            if m < mmax:
                cur = lap({mm: cc for mm, cc in cur.items() if sum(mm) <= 2 * (mmax - m)}, Hinv)
    return out, W
def pf_test(alpha, Y0, N):
    Y = var(Y0)
    Phi, W = phi_coeffs(alpha, Y, N)
    q = [Y + a for a in alpha]
    thY = -(Y * q[0] * q[1] * q[2]) * W.inv()
    th = lambda f: f.deriv() * thY
    y = -(q[0] * q[1] * q[2]) * (Y * Y * Y * 27).inv()
    half = fr(1, 2)
    thlogW = th(W) * W.inv()
    def T(f):
        r = [TS.const(0) for _ in range(len(f) + 1)]
        for i, c in enumerate(f):
            r[i] = r[i] + Y * c
            r[i + 1] = r[i + 1] + th(c) - half * thlogW * c
        return r
    A = list(Phi)
    for j in range(3):
        tA = T(A)
        A = [tA[i] + (A[i] * alpha[j] if i < len(A) else TS.const(0)) for i in range(len(tA))]
    B = list(Phi)
    for j in range(3):
        tb = T(B)
        C = [x * 3 for x in tb]
        for i, c in enumerate(B):
            C[i + 1] = C[i + 1] + c * j
        B = C
    tot = [A[i] + y * B[i] for i in range(len(A))]
    return tot, Phi, W
def bern(n):
    B = [Fr(0)] * (n + 1); B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(Fr(math.comb(m + 1, k)) * B[k] for k in range(m)) / (m + 1)
    return B
def Theta_series(x, N):
    B = bern(N + 3)
    l = [0] * (N + 1)
    for m in range(1, N + 2):
        k = 2 * m - 1
        if k <= N:
            f = B[2 * m] / (2 * m * (2 * m - 1))
            l[k] = fr(f.numerator, f.denominator) * pow(x, k, p) % p
    e = [0] * (N + 1); e[0] = 1
    for n in range(1, N + 1):
        s = 0
        for k in range(1, n + 1):
            s += k * l[k] * e[n - k]
        e[n] = s % p * inv(n) % p
    return e
def smul(a, b, N):
    r = [0] * (N + 1)
    for i, x in enumerate(a):
        for j, y_ in enumerate(b):
            if i + j <= N: r[i + j] += x * y_
    return [x % p for x in r]
def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    alpha = [fr(5, 7), fr(11, 3), fr(-2, 13)]
    tot, Phi, W = pf_test(alpha, fr(17, 19), N)
    print("generic alpha: Phi_0 == 1:", Phi[0].c[0] == 1)
    # usable delta orders: 4 - number of thetas; test the delta^0 and delta^1 coefficients where valid
    for i, c in enumerate(tot[:N + 2]):
        print(" z^%d coefficient of W^{1/2} P W^{-1/2} Phi vanishes, delta^0,1:" % i, [x == 0 for x in c.c[:2]])
    for i in range(3):
        P, W0 = phi_coeffs(alpha, TS.const(-alpha[i]), N)
        j, k = [x for x in range(3) if x != i]
        a = [3 * alpha[i] % p, (alpha[j] - alpha[i]) % p, (alpha[k] - alpha[i]) % p]
        pr = [1] + [0] * N
        for ai in a:
            pr = smul(pr, Theta_series((-inv(ai)) % p, N), N)
        print(" colour", i, "Phi(alpha,-alpha_i) = prod Theta(-z/a):", [P[n].c[0] == pr[n] for n in range(N + 1)])
    E = fr(23, 5); a = E * fr(1, 3) % p
    alpha_s = [a, a, a]
    Y0 = fr(9, 11)
    P, W1 = phi_coeffs(alpha_s, TS.const(Y0), N)
    r = [1]
    for m in range(0, N + 3):
        rp = (3 * (3 * m * m + 3 * m + 1) * r[m] - m ** 3 * (r[m - 1] if m >= 1 else 0)) % p * inv(27 * (m + 1)) % p
        r.append(rp)
    v = inv((Y0 + a) % p)
    S = [0] * (N + 1)
    for l in range(N + 1):
        den = [1] + [0] * N
        for jj in range(1, l + 2):
            g = [pow((-jj * inv(E)) % p, n, p) for n in range(N + 1)]
            den = smul(den, g, N)
        for n in range(N + 1 - l):
            S[n + l] = (S[n + l] + r[l] * pow(v, l, p) * den[n]) % p
    TH = smul(Theta_series((-inv(E)) % p, N), S, N)
    print(" T_F Phi = Theta(-z/E) S:", [P[n].c[0] == TH[n] for n in range(N + 1)])
    tot, _, _ = pf_test(alpha_s, Y0, N)
    print(" T_F PF:", [c.c[0] == 0 for c in tot[:N + 2]])
main()
