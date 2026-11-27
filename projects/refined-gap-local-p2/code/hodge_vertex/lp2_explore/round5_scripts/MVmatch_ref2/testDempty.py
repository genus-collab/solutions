# (D) empty-profile form, literal; + tau=0 all-genus closed form check of (A) at g up to 12; + one-part closed form control.
from fractions import Fraction as F
from testD import *
from hodge import Xi_def,BB
from math import factorial

def A(g):
    return F(1,2*factorial(2*g-2))*abs(BB[2*g-2])/(2*g-2)*abs(BB[2*g])/(2*g)

for h in [2,3,4,5]:
    for tau in [F(5,3),F(-1,2),F(-9,7)]:
        sig=tau*(1+tau);N=3*h-3
        P=2*h-2+N+N+1
        L=Line((),tau,N,P)
        tot=lag(N,2*h-2,0)*(-A(h)*sig**(1-h))
        ctl=F(0)
        for a in range(1,N+1):
            co=L.coeff(a,2*h-2+a)
            term=(-1)**a*sig**(-a)*lag(N-a,2*h-2+a,-a)*co/sig**(h-1)
            tot+=term;ctl+=term
        d=Xi_def(h,tau,()) if h<=4 else None
        print("h",h,"tau",tau,"Xi(D-form)",tot,"def",d,"OK" if d==tot else ("no def" if d is None else "DIFF"),"| drop a=0 equal:",ctl==d,flush=True)
