"""Referee: leg functionals on the dxi-basis and contraction of symmetric tensors (independent of worker code)."""
from fractions import Fraction
from rf_tr import inv, smul, sinv, sderiv

def leg_pow(cv, m):
    """<p^m, dxi_{b,k}> = Res_{c_b} p^m dxi_{b,k} = [zeta^(2k+1)] (c_b + t_b(zeta))^m  (m may be negative)."""
    P, out = cv.P, {}
    for b in range(5):
        pz = cv.tser[b][:]; pz[0] = (pz[0] + cv.c[b]) % P
        base = pz if m >= 0 else sinv(pz, cv.N, P)
        r = [1] + [0] * (cv.N - 1)
        for _ in range(abs(m)):
            r = smul(r, base, cv.N, P)
        for k in range(cv.KM + 1):
            out[cv.idx(b, k)] = r[2 * k + 1]
    return out

def leg_log(cv, a=0):
    """<-log(p-a), dxi_{b,k}> = Int_a^oo dxi_{b,k} = -(1/(2k+1)) [zeta^(2k)] t'/(c_b - a + t).  a = 0: the log leg."""
    P, out = cv.P, {}
    for b in range(5):
        t = cv.tser[b]; td = sderiv(t, P)
        d = t[:]; d[0] = (d[0] + cv.c[b] - a) % P
        f = smul(td, sinv(d, cv.N, P), cv.N, P)
        for k in range(cv.KM + 1):
            out[cv.idx(b, k)] = (-f[2 * k] * inv(2 * k + 1, P)) % P
    return out

def leg_log_ebasis(cv):
    """Second route: dxi_{b,k} = Sum_m T[b,k][m] dp/(p-c_b)^(m+1), Int_0^oo dp/(p-c)^(m+1) = (-1)^m/(m c^m)."""
    P, out = cv.P, {}
    for b in range(5):
        c = cv.c[b]
        for k in range(cv.KM + 1):
            r = 0
            for m, cf in cv.T[b, k].items():
                r += cf * (-1) ** m * inv(m * pow(c, m, P), P)
            out[cv.idx(b, k)] = r % P
    return out

def contract(T, vec, P):
    """Symmetric tensor T (sorted tuples) contracted in one slot with the functional vec."""
    out = {}
    for S, w in T.items():
        last = None
        for q, x in enumerate(S):
            if x == last: continue
            last = x
            f = vec.get(x, 0)
            if f:
                R = S[:q] + S[q + 1:]
                out[R] = (out.get(R, 0) + w * f) % P
    return {k: v for k, v in out.items() if v}

def bern(n):
    B = [Fraction(1)]
    from math import comb
    for m in range(1, n + 1):
        B.append(-sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1))
    return B[n]

def Fderiv(h, k):
    """k-th v-derivative at v = 1 of the TR-normalised F_h: B_2h/(2h(2h-2)) v^(2-2h) (h>=2),
    -(1/12) log v (h=1, k>=1), (1/2) v^2 log v (h=0, k>=3)."""
    if h >= 2:
        c = bern(2 * h) / (2 * h * (2 * h - 2)); e = 2 - 2 * h
        for _ in range(k): c *= e; e -= 1
        return c
    if h == 1:
        c = Fraction(-1, 12); e = -1          # first derivative -1/(12 v)
        for _ in range(k - 1): c *= e; e -= 1
        return c
    c = Fraction(1); e = -1                   # third derivative 1/v
    for _ in range(k - 3): c *= e; e -= 1
    return c

def modq(fr, P):
    return fr.numerator % P * inv(fr.denominator, P) % P

def ratrec(a, P):
    """Rational reconstruction of a mod P (small heights)."""
    a %= P
    r0, r1, s0, s1 = P, a, 0, 1
    import math
    bound = math.isqrt(P // 2)
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)
