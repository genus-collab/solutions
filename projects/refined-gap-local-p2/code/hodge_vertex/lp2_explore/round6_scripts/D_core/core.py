# Rows of Definition 1.3 from the closed form only, expanded in t_c.  No use of Lemma 4.1, H2, Prop 2.10 in the construction.
from fractions import Fraction as Fr
from alg import *

def rec_r(N,perturb=None):
    r=[Fr(1)]
    for m in range(N+1):
        rm=r[m]; rm1=r[m-1] if m>=1 else Fr(0)
        r.append((3*(3*m*m+3*m+1)*rm-m**3*rm1)/(27*(m+1)))
    if perturb:
        for k,v in perturb.items(): r[k]=r[k]+v
    return r

def compose(fs,D,prec):
    out=LS({},prec)
    for c in reversed(fs):
        out=(out*D)+LS({0:c},prec)
    return out

class Base:
    """z-free series in t_c at the conifold, on the real branch."""
    def __init__(s,E,prec,fsec=None):
        s.E=E; s.prec=prec
        n=prec+2
        wc=[Fr(1)]
        for k in range(n): wc.append(wc[k]*(k*k+k+Fr(2,9))/((k+1)**2))
        # dt_c/dDelta = wc/(1-Delta)
        dt=[]; acc=Fr(0)
        for k in range(n+1):
            acc+=wc[k]; dt.append(acc)
        tcofD=[Fr(0)]+[dt[k]/(k+1) for k in range(n+1)]
        tc=LS({1:Fr(1)},n)
        D=tc
        for it in range(n):
            D=tc-compose([Fr(0),Fr(0)]+tcofD[2:],D,n)
        s.Delta=D
        s.tc=tc
        s.wc=compose(wc,D,n)
        s.C1cf=-s.wc
        s.g=s.C1cf
        # (1-Delta)^(-1/3)
        u=[Fr(1)]
        for k in range(n): u.append(u[k]*(Fr(1,3)+k)/(k+1))
        s.u=compose(u,D,n)
        s.v0=(s.u*(-3/E))+LS({0:3/E},n)
        s.vinv0=s.v0.inv()
        s.Y0=s.vinv0-LS({0:E/3},n)
        s.rho=(s.v0*(-E)).shift(-1)
        s.rhoinv=s.rho.inv()
        f=fsec if fsec is not None else [Fr(1)]
        s.f=compose(f,D,n)
        s.c1=s.C1cf*s.f
        s.C1i=s.c1.inv()
        c12=s.c1*s.c1; C12=s.C1i*s.C1i
        s.c12=c12
        # d2=E(C1i^2-1)/(2 Delta_c)
        num=(C12-LS({0:Fr(1)},n))*(E/2)
        s.d2=num*s.Delta.inv()
        s.y=(s.Delta-LS({0:Fr(1)},n))*Fr(1,27)
    def theta(s,F):
        return s.g*F.deriv()

def S_colour0(B,r,sgn=-1,L=None):
    E=B.E; L=L if L is not None else B.prec+1
    out=LS({},B.prec)
    vp=LS({0:Fr(1)},B.prec)
    den=R([1])
    zp=R([1])
    z=R([0,Fr(sgn)])
    for l in range(L+1):
        # z=sgn*psi; term r_l (v z)^l / prod_{j<=l+1}(1+j z/E)
        den=den*R([1,Fr(sgn*(l+1))/E])
        term=zp/den*r[l]
        out=out+vp*term
        vp=vp*B.v0
        zp=zp*z
    return out

def rows0(B,r,sgn=-1,xxshift=None):
    """a=(a0,a1,a2) = sqrt(c0) P^B at colour 0, sigma = rho a/v = (-E/t_c) a, z = sgn*psi"""
    E=B.E
    S=S_colour0(B,r,sgn)
    z=R([0,Fr(sgn)])
    a0=B.v0*S
    a1=B.C1i*((B.Y0*a0)+(B.theta(a0)*z))
    ins=(B.Y0*a1)+(B.theta(a1)*z)-(B.d2*a1)
    if xxshift is not None:
        # XX -> XX + delta without changing c1: theta(C1i Q) gains -delta C1i Q, i.e. -z delta P1
        ins=ins-(a1*(z*xxshift))
    a2=(B.Delta*B.c12)*ins
    fac=LS({-1:-E},10**6)
    sig=[fac*a0,fac*a1,fac*a2]
    return S,(a0,a1,a2),sig

def stripped_PF(B,a0,z):
    """Delta Th^3 + (E+27 y z) Th^2 + (E^2/3+6 y z^2) Th + E^3/27 applied to a0"""
    E=B.E
    Th=lambda F:(B.Y0*F)+(B.theta(F)*z)
    t1=Th(a0); t2=Th(t1); t3=Th(t2)
    A=B.Delta*t3
    C=((B.y*27)*t2)*z
    D=(t2*E)
    G=((B.y*6)*t1)*(z*z)
    H=t1*(E*E/3)
    I=a0*(E**3/27)
    return A+C+D+G+H+I

def eta_of(e1,e2):
    E=e1+e2
    return {(0,0):e1*e2*E**3/9,(1,1):-3*e1*e2*E,(1,2):-3*e1*e2,(2,1):-3*e1*e2}
