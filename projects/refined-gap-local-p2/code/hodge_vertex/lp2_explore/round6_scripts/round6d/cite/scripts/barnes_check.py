# Independent high-precision check of Barnes's second lemma (form of FRAMEB C Step 6).
# Stdlib only (decimal). Evidence, not a proof step.
# Trapezoid rule on a vertical line: exponentially convergent for an analytic integrand
# decaying like exp(-c|t|); step h gives error about exp(-2 pi d / h), d the distance to the
# nearest pole. Run with two steps to measure the truncation.
from decimal import Decimal as Dd, getcontext
from fractions import Fraction as Fr
import sys, time

getcontext().prec = 60
PI = Dd("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")

def cadd(a, b): return (a[0]+b[0], a[1]+b[1])
def csub(a, b): return (a[0]-b[0], a[1]-b[1])
def cmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cinv(a):
    n = a[0]*a[0]+a[1]*a[1]
    return (a[0]/n, -a[1]/n)

def dsin_cos(x):
    # reduce mod 2 pi, then Taylor
    twopi = 2*PI
    k = (x/twopi).to_integral_value()
    x = x-k*twopi
    s, c, term = Dd(0), Dd(0), Dd(1)
    n = 0
    eps = Dd(10)**(-getcontext().prec-5)
    while abs(term) > eps or n < 4:
        if n % 4 == 0: c += term
        elif n % 4 == 1: s += term
        elif n % 4 == 2: c -= term
        else: s -= term
        n += 1
        term = term*x/n
    return s, c

def datan(x):
    # atan by argument halving, then Taylor
    neg = x < 0
    x = abs(x)
    inv = x > 1
    if inv: x = 1/x
    k = 0
    while x > Dd("0.05"):
        x = x/(1+(1+x*x).sqrt())
        k += 1
    s, p, n = Dd(0), x, 1
    x2 = x*x
    eps = Dd(10)**(-getcontext().prec-5)
    while abs(p) > eps:
        s += p/n if (n % 4 == 1) else -p/n
        p *= x2
        n += 2
    s *= 2**k
    if inv: s = PI/2-s
    return -s if neg else s

def clog(a):
    re = ((a[0]*a[0]+a[1]*a[1]).ln())/2
    x, y = a
    if x > 0: th = datan(y/x)
    elif y >= 0: th = PI/2-datan(x/y) if y != 0 else PI
    else: th = -PI/2-datan(x/y)
    return (re, th)

def cexp(a):
    r = a[0].exp()
    s, c = dsin_cos(a[1])
    return (r*c, r*s)

def bernoulli(nmax):
    B = [Fr(0)]*(nmax+1)
    A = [Fr(0)]*(nmax+1)
    for m in range(nmax+1):
        A[m] = Fr(1, m+1)
        for j in range(m, 0, -1):
            A[j-1] = j*(A[j-1]-A[j])
        B[m] = A[0]
    return B

BN = bernoulli(70)
STIR = [(k, Dd(BN[2*k].numerator)/Dd(BN[2*k].denominator)/(2*k*(2*k-1))) for k in range(1, 35)]
HALF_LOG_2PI = (2*PI).ln()/2

def clgamma(z):
    # log Gamma by upward shift, then Stirling
    shift = 40
    acc = (Dd(0), Dd(0))
    w = z
    # accumulate log of product z(z+1)...(z+shift-1) by summing logs (branch: continuous sum)
    for _ in range(shift):
        acc = cadd(acc, clog(w))
        w = (w[0]+1, w[1])
    lw = clog(w)
    r = cmul((w[0]-Dd("0.5"), w[1]), lw)
    r = csub(r, w)
    r = (r[0]+HALF_LOG_2PI, r[1])
    inv = cinv(w)
    inv2 = cmul(inv, inv)
    p = inv
    for k, c in STIR:
        r = cadd(r, (p[0]*c, p[1]*c))
        p = cmul(p, inv2)
    return csub(r, acc)

def cgamma(z): return cexp(clgamma(z))
def rgamma(x): return cgamma((Dd(x), Dd(0)))[0]

def gamma_cross():
    # control: Gamma(1/3)Gamma(2/3) = 2 pi / sqrt 3
    t = Dd(1)/Dd(3)
    return rgamma(t)*rgamma(2*t), 2*PI/Dd(3).sqrt()

def lhs_sigma(a, b, c, d, e, sig, h, T):
    """(1/2 pi i) Int Gamma(a+s)Gamma(b+s)Gamma(c+s)Gamma(1-d-s)Gamma(-s)/Gamma(e+s) ds on Re s = sig"""
    tot = (Dd(0), Dd(0))
    n = int(T/h)
    for k in range(-n, n+1):
        t = h*k
        s = (sig, t)
        lg = (Dd(0), Dd(0))
        for p in (a, b, c):
            lg = cadd(lg, clgamma((p+s[0], s[1])))
        lg = cadd(lg, clgamma((1-d-s[0], -s[1])))
        lg = cadd(lg, clgamma((-s[0], -s[1])))
        lg = csub(lg, clgamma((e+s[0], s[1])))
        tot = cadd(tot, cexp(lg))
    # ds = i dt, so (1/2 pi i) Int f ds = (1/2 pi) Int f dt
    return tot[0]*h/(2*PI), tot[1]*h/(2*PI)

def rhs(a, b, c, d, e):
    g = rgamma
    return (g(a)*g(b)*g(c)*g(1-d+a)*g(1-d+b)*g(1-d+c))/(g(e-a)*g(e-b)*g(e-c))

def run(a, b, c, d, sig, h, T, label):
    e = a+b+c-d+1
    t0 = time.time()
    v = lhs_sigma(a, b, c, d, e, sig, h, T)
    r = rhs(a, b, c, d, e)
    print(label, "sigma line", str(sig)[:12], "h", h, "T", T, "e'", str(e)[:12])
    print("  LHS re", str(v[0])[:45])
    print("  LHS im", str(v[1])[:25])
    print("  RHS   ", str(r)[:45])
    print("  LHS-RHS", str(v[0]-r)[:25], " time", round(time.time()-t0, 1))
    return v[0], r

if __name__ == "__main__":
    a, b = gamma_cross()
    print("control Gamma(1/3)Gamma(2/3) - 2pi/sqrt3 =", str(a-b)[:25])
    print("control Gamma(5)=", str(rgamma(5))[:30], " Gamma(1/2)^2-pi =", str(rgamma(Dd("0.5"))**2-PI)[:20])
    th = Dd(1)/Dd(3)
    # C Step 6: (a',b',c',d',e') = (2/3,2/3,2/3,4/3,5/3), contour -2/3 < Re sigma < -1/3
    h = Dd(1)/Dd(20)
    run(2*th, 2*th, 2*th, 4*th, -Dd(1)/Dd(2), h, Dd(16), "CStep6 h=1/20")
    run(2*th, 2*th, 2*th, 4*th, -Dd(1)/Dd(2), Dd(1)/Dd(40), Dd(16), "CStep6 h=1/40")
    print("  target (2 pi/sqrt 3)^3 =", str((2*PI/Dd(3).sqrt())**3)[:45])
    # Generic control of the printed lemma with unrelated parameters; contour max(-a,-b,-c,1-d) < ... < 0
    run(Dd("0.3"), Dd("0.45"), Dd("0.7"), Dd("0.2"), -Dd("0.15"), Dd(1)/Dd(40), Dd(16), "generic h=1/40")
    # Negative control: e' != a'+b'+c'-d'+1 must fail
    a0, b0, c0, d0 = 2*th, 2*th, 2*th, 4*th
    v = lhs_sigma(a0, b0, c0, d0, Dd(2), -Dd(1)/Dd(2), Dd(1)/Dd(20), Dd(16))
    print("negative control e'=2 (wrong): LHS", str(v[0])[:20], " formula", str(rhs(a0, b0, c0, d0, Dd(2)))[:20])
    # T(1) value in the s-variable, K_B times the integral, against 2 pi / sqrt 3
    KB = (Dd(3).sqrt()/(2*PI))**2
    tot = (Dd(0), Dd(0))
    h = Dd(1)/Dd(40); n = int(Dd(16)/h); sig = -Dd(1)/Dd(6)
    for k in range(-n, n+1):
        t = h*k
        lg = (Dd(0), Dd(0))
        for _ in range(3):
            lg = cadd(lg, clgamma((-sig, -t)))
        lg = cadd(lg, clgamma((th+sig, t)))
        lg = cadd(lg, clgamma((2*th+sig, t)))
        lg = csub(lg, clgamma((1-sig, -t)))
        tot = cadd(tot, cexp(lg))
    val = tot[0]*h/(2*PI)*KB
    print("T(1) direct on Re s=-1/6:", str(val)[:45], " 2pi/sqrt3:", str(2*PI/Dd(3).sqrt())[:45], " diff", str(val-2*PI/Dd(3).sqrt())[:25])
