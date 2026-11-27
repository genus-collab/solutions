#!/usr/bin/env python
"""Mechanical check of the grading: expand (*) in monomials x^a w^b from D=(w-1)w^2 d_w, Taylor shifts,
with Y = DF built from G_0, G_1 exact and G_h = w^{2h-2}(Z_h + c w + ...), Z_h from beta formula (dual numbers in tau).
Grade(x^a w^b)=b-a. Claim: grades<=0 vanish identically; grade 1 vanishes for Z_h=beta and does not see c;
perturbing one Z_h breaks grade 1; grade 2 does see c. Evidence only."""
from fractions import Fraction as Fr
import sys, random, math
from vi_check import Dual, beta_dual
amax=int(sys.argv[1]) if len(sys.argv)>1 else 12
def add(d,k,v):
    d[k]=d.get(k,Dual(0))+v
def D(Y):
    out={}
    for (a,b),c in Y.items():
        if b!=0:
            add(out,(a,b+2),c*b); add(out,(a,b+1),c*(-b))
    return out
def dtau(Y): return {k:Dual(c.b,0) for k,c in Y.items()}   # first derivative only: see note
def lin(*terms):
    out={}
    for f,Y in terms:
        for k,c in Y.items(): add(out,k,c*f)
    return out
def mulx(Y,n,cut=None):
    return {(a+n,b):c for (a,b),c in Y.items() if a+n<=amax+3}
def trunc(Y): return {k:c for k,c in Y.items() if k[0]<=amax+3}
def run(tau0,pert=None,cterm=False):
    t=Fr(tau0)
    # coefficients are Dual(val, dval/dtau) -- for T1 we need d/dtau only once; use val part elsewhere
    tau=Dual(t,1); e=1+tau; sig=tau*e
    Y={}
    # G_0: DG_0=(1-w^-2)/2
    add(Y,(-2,0),Dual(1,0)/sig*Fr(1,2)); add(Y,(-2,-2),Dual(1,0)/sig*Fr(-1,2))
    # G_1: DG_1=(w-1)(w-1-1/sig)/24 = (w^2 -(2+1/sig) w + (1+1/sig))/24
    add(Y,(0,2),Dual(Fr(1,24))); add(Y,(0,1),(2+1/sig)*Fr(-1,24)); add(Y,(0,0),(1+1/sig)*Fr(1,24))
    h=2
    while 2*h-2<=amax+3:
        z=beta_dual(h,-tau,e)    # sigma^{h-1} Z_h
        if pert==h: z=z+Fr(1)
        n=2*h-2
        # D w^n = n w^{n+2}-n w^{n+1}
        add(Y,(n,n+2),z*n); add(Y,(n,n+1),z*(-n))
        if cterm:
            cc=Fr(random.randint(1,50),7)*(tau+1)  # grade 1 term c w^{n+1}, tau dependent
            add(Y,(n,n+3),cc*(n+1)); add(Y,(n,n+2),cc*(-(n+1)))
        h+=1
    val={k:Dual(c.a,0) for k,c in Y.items()}
    dY={k:Dual(c.b,0) for k,c in Y.items()}
    sg=sig.a; sp=sig.b
    # Taylor: E = sum beta^k D^k/k!
    powers=[val]; 
    Kmax=amax+6
    for k in range(1,Kmax+1): powers.append(trunc(D(powers[-1])))
    pd=[dY]
    for k in range(1,Kmax+1): pd.append(trunc(D(pd[-1])))
    bfac=lambda k:(t**k)/math.factorial(k)
    # S := e^{beta D}Y - Y
    S={}
    for k in range(1,Kmax+1):
        for (a,b),c in powers[k].items():
            if a+k<=amax+3: add(S,(a+k,b),c*bfac(k))
    # T1 = -(2/x) sum_{k>=1} beta^k/k! D^{k-1}(dtau Y + (sp/sg) D Y)
    T={}
    for k in range(1,Kmax+1):
        for (a,b),c in pd[k-1].items():
            if a+k-1<=amax+3: add(T,(a+k-1,b),c*bfac(k)*(-2))
        for (a,b),c in powers[k].items():
            if a+k-1<=amax+3: add(T,(a+k-1,b),c*bfac(k)*(-2)*(sp/sg))
    # T2=(D-1)(2Y+S)
    two=lin((Fr(2),val),(Fr(1),S))
    for k,c in trunc(D(two)).items(): add(T,k,c)
    for k,c in two.items(): add(T,k,c*(-1))
    # T3=S^2
    items=list(S.items())
    for (a,b),c in items:
        for (a2,b2),c2 in items:
            if a+a2<=amax+1: add(T,(a+a2,b+b2),c*c2)
    res={}
    for (a,b),c in T.items():
        if a<=amax and c.a!=0:
            g=b-a
            res.setdefault(g,[]).append((a,b))
    return res
for tau0 in (Fr(2,7),Fr(-5,3)):
    r=run(tau0); print("tau",tau0,"Z=beta: nonzero by grade:",{g:len(v) for g,v in sorted(r.items())}, "lowest grade", min(r) if r else None)
    r=run(tau0,cterm=True); print("  + random grade-1 terms: nonzero by grade:",{g:len(v) for g,v in sorted(r.items())})
    r=run(tau0,pert=3); print("  Z_3 perturbed: nonzero by grade:",{g:len(v) for g,v in sorted(r.items())}, "lowest", min(r) if r else None)
