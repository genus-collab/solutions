# High precision check of Lemma (t_cf). Evidence only.
from decimal import Decimal as D, getcontext
from fractions import Fraction as Fr
import math
getcontext().prec=60
PI=D("3.14159265358979323846264338327950288419716939937510582097494459")
s3=D(3).sqrt()
def toD(fr): return D(fr.numerator)/D(fr.denominator)
def H(n): return sum(Fr(1,k) for k in range(1,n+1))
def series_y(y,N=400):
    # g=3 Sum (3d-1)!/(d!)^3 (-y)^d, g2 = 9 Sum (3d-1)!/(d!)^3 (psi(3d)-psi(d+1)) (-y)^d, A, h
    Y=-y   # positive
    g=D(0); g2=D(0); A=D(1); h=D(0)
    # a_d=(3d)!/(d!)^3 in Decimal recursively
    a=D(1); pw=D(1); h3=Fr(0); h1=Fr(0)
    for d in range(1,N+1):
        a=a*D((3*d-2)*(3*d-1)*3*d)/D(d**3)
        pw=pw*Y
        term=a*pw
        h1+=Fr(1,d)
        h3+=(Fr(1,3*d-3) if d>1 else Fr(0))+Fr(1,3*d-2)+Fr(1,3*d-1)   # H_{3d-1}
        # psi(3d)-psi(d+1)=H_{3d-1}-H_d
        g+=term*3/D(3*d)*3/3 if False else term/D(d)
        g2+=term/D(3*d)*9*toD(h3-h1)
        A+=term
        h+=term*3*toD(h3+Fr(1,3*d)-h1)
        if abs(term)<D(10)**(-55): break
    return g,g2,A,h
def wc_tc(xi,N=800):
    # wc=2F1(1/3,2/3;1;xi), t_c=Int_0^xi wc/(1-x)
    c=D(1); pw=D(1); s=D(0); tc=D(0); w=D(0)
    for n in range(N):
        w+=c*pw
        s+=c
        tc+=s*pw*xi/(n+1)
        c=c*(D(1)/3+n)*(D(2)/3+n)/D((n+1)**2)
        pw*=xi
    return w,tc
for y in (Fr(-1,54),Fr(-1,81),Fr(-1,270),Fr(-1,40)):
    yd=toD(y); xi=1+27*yd
    l=(-yd).ln()
    g,g2,A,h=series_y(yd)
    tcf=(s3/(2*PI))*(l*l/2+g*l+g2-PI*PI/2)
    # other branch: log y = l - i pi : Re and Im
    # I12 = (L^2)/2+gL+g2, I01 = L+g; expression I12 - i pi I01 - pi^2
    re2=l*l/2-PI*PI/2+g*l+g2 - PI*PI - PI*PI
    im2=-PI*l-PI*g-PI*(l+g)
    w,tc=wc_tc(xi)
    Acf=(s3/(2*PI))*(A*l+h)
    print(f"y={y}: xi={xi:.6f} t_c={tc:.30f}\n   tcf-t_c={tcf-tc:.3e}; A_cf(=sqrt3/2pi (A l+h)) + wc = {Acf+w:.3e}; other branch re={float((s3/(2*PI))*re2):.5f} im={float((s3/(2*PI))*im2):.5f}")
    # Barnes constant check: 2 pi/sqrt3 and (log G)''(0)
# Wronskian-free check of theta t_cf=-w: dt/dtheta done above via A_cf + wc=0
