# check4.py -- evidence only. Controls on the finite formula (STAR); closed unstable table; empty profile through M_0 = -A_h sigma^(1-h).
from fractions import Fraction as Fr
from math import factorial
from engine import *
from schur import *
from check3 import lag, sig, SchurXi, series_from_M

def star(h, R, tau, dal=0, autfac=True, sgn=1):
    r = len(R); J = sum(R); N = 3*h-3+r
    P = 2*h-2+r+N+J+N
    C = ConnSeries(tau, P)
    s = sig(tau)
    tot = Fr(0)
    for a in range(N+1):
        mu = tuple(sorted(R+(1,)*a, reverse=True))
        l = len(mu)
        tot += Fr(sgn)**a * Fr(-1)**a * s**(-a) * lag(N-a, 2*h-2+r+a+dal, -J-a) * C.coeff(mu, 2*h-2+l)
    pref = (aut(R) if autfac else 1) / (s**(h+r-1) * PiMu(R, tau))
    return pref*tot

tau = Fr(2,7)
print("(STAR) in the a!-free form, true formula, in range (expected 0), and first point outside:")
for h,R in ((2,(2,)),(3,(3,)),(3,(2,2)),(3,(4,)),(4,(4,)),(3,(5,)),(2,(3,))):
    print("  ", h, R, star(h,R,tau), " direct:", Xi_direct(h,R,ts_of(tau)) if h<4 else "-")
print("Controls at in-range points (must be nonzero):")
for h,R in ((2,(2,)),(3,(3,)),(3,(2,2))):
    print("   Laguerre index off by one:", h, R, star(h,R,tau,dal=-1) != 0, "; no (-1)^a:", star(h,R,tau,sgn=-1) != 0)
print("Control, |Aut R| dropped at an outside point (2,(2,2)): differs from direct:", star(2,(2,2),tau,autfac=False) != Xi_direct(2,(2,2),ts_of(tau)))
# empty profile: Xi_h(empty) with M_0 = -A_h sigma^(1-h), rest from the Schur side
def Ag(g):
    return Fr(1, 2*factorial(2*g-2))*abs(bern(2*g-2))/(2*g-2)*abs(bern(2*g))/(2*g)
for h in (2,3):
    for tau in (Fr(2,7),Fr(-7,3)):
        N = 3*h-3
        C = ConnSeries(tau, 2*h-2+N+N)
        s = sig(tau)
        tot = lag(N, 2*h-2, Fr(0))*(-Ag(h)*s**(1-h))
        for a in range(1,N+1):
            tot += Fr(-1)**a * lag(N-a, 2*h-2+a, Fr(-a)) * C.coeff((1,)*a, 2*h-2+a) * s**(-(h+a-1))
        print(f"  Xi_{h}(empty) at tau={tau}: from MV + M_0 = {tot}, direct = {Xi_direct(h,(),ts_of(tau))}", "match" if tot==Xi_direct(h,(),ts_of(tau)) else "MISMATCH")
