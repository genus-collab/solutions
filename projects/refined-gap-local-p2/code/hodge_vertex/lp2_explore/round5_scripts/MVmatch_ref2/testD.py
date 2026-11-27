# (D): Schur-only Xi_h(R) to high genus, all tau; compare with literal definition at h<=5 where engine can reach;
# test Lemma-A range zeros and the Laurent-in-sigma claim (E); controls that must fail.
from fractions import Fraction as F
from hodge import Xi_def,Qk,ch,omega_ts,omega_monos
from schur import *
from math import comb,factorial
from collections import Counter
import time,sys

def lag(n,al,x):
    # L_n^{(al)}(x)=sum_i (-1)^i C(n+al,n-i) x^i/i!
    if n<0: return F(0)
    s=F(0)
    for i in range(n+1):
        # binomial with integer n+al>=0 assumed
        s+=(-1)**i*comb(n+al,n-i)*F(x)**i/factorial(i)
    return s

def chi_mixed(nu,R,a):
    return chi(nu,tuple(sorted(R,reverse=True))+(1,)*a)

class Line:
    """Compute [p_S p_1^a] C for S sublists of R via series in t=p_1 and nilpotents for R (labelled)."""
    def __init__(s,R,tau,amax,Pext):
        s.R=tuple(R);s.tau=F(tau);s.amax=amax
        s.J=sum(R)
        # relative precision: coefficient of t^a in masks
        s.P=Pext
        l=len(R)
        s.l=l
        # Z[mask][a]: series; (labelled |Aut R| conv) : coefficient of t^a/ a!-less
        s.Z={}
        s.cum={}
        for m in range(1<<l):
            sub=tuple(R[i] for i in range(l) if m>>i&1)
            pr=1
            for q in sub: pr*=q
            rows=[]
            for a in range(amax+1):
                n=sum(sub)+a
                r=[F(0)]*s.P
                for nu in partitions(n):
                    c=chi_mixed(nu,sub,a) if n>0 else 1
                    if c==0: continue
                    y=Y(nu,s.tau,s.P)
                    for i in range(s.P): r[i]+=c*y[i]
                f=F(1,pr*factorial(a))
                rows.append([v*f for v in r])
            s.Z[m]=rows
        # log: c_m(t) with Z_m=sum_{B containing low} c_B Z_{m^B}; c_0 = log Z_0 handled: for m=0 use t d/dt recursion
        # Valuation convention: coefficient t^a p_S has x^{-(|S|+a)}*series
        c0=[None]*(amax+1)
        Z0=s.Z[0]
        for a in range(1,amax+1):
            # a Z_a = sum_{b=1}^a b c_b Z_{a-b}
            r=[a*v for v in Z0[a]]
            for b in range(1,a):
                pp=mul(c0[b],Z0[a-b],s.P)
                for i in range(s.P): r[i]-=b*pp[i]
            c0[a]=[v/a for v in r]
        c0[0]=[F(0)]*s.P
        s.cum[0]=c0
        for m in range(1,1<<l):
            low=m&-m
            rest=m^low
            rows=[list(v) for v in s.Z[m]]
            for sp in submasks(rest):
                B=low|sp
                if B==m: continue
                other=m^B
                for a in range(amax+1):
                    for b in range(a+1):
                        pp=mul(s.cum[B][b],s.Z[other][a-b],s.P)
                        for i in range(s.P): rows[a][i]-=pp[i]
            # divide by Z_0(t): c_m Z_0 = rows
            out=[None]*(amax+1)
            for a in range(amax+1):
                r=list(rows[a])
                for b in range(a):
                    pp=mul(out[b],Z0[a-b],s.P)
                    for i in range(s.P): r[i]-=pp[i]
                out[a]=r
            s.cum[m]=out
    def coeff(s,a,e):
        """ |Aut R| [p_R p_1^a] C at x^e """
        m=(1<<s.l)-1
        i=e+s.J+a
        if i<0: return F(0)
        return s.cum[m][a][i]

def Xi_schur(h,R,tau,variant=None):
    r=len(R);J=sum(R);N=3*h-3+r
    tau=F(tau);sig=tau*(1+tau)
    Pext=2*h-2+r+N+J+N+1
    L=Line(R,tau,N,Pext)
    piR=F(1)
    for j in R: piR*=Pi(j,tau)
    tot=F(0)
    for a in range(N+1):
        e=2*h-2+r+a
        co=L.coeff(a,e)
        if variant=="lagidx":
            lg=lag(N-a,2*h-3+r+a,-J-a)
        elif variant=="nosign":
            lg=lag(N-a,2*h-2+r+a,-J-a)*(-1)**a
        elif variant=="J":
            lg=lag(N-a,2*h-2+r+a,-J)
        else:
            lg=lag(N-a,2*h-2+r+a,-J-a)
        tot+=(-1)**a*sig**(-a)*lg*co
    return tot/(sig**(h+r-1)*piR)

def D(R): return sum(R)-len(R)

if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="cmp":
        for h,Rs in [(2,[(2,),(3,),(4,),(2,2),(3,2),(5,)]),(3,[(3,),(5,),(6,),(2,2),(4,2)]),(4,[(2,),(7,),(8,),(3,3)])]:
            for tau in [F(5,3),F(-9,7)]:
                for R in Rs:
                    t0=time.time()
                    a=Xi_schur(h,R,tau)
                    b=Xi_def(h,tau,R)
                    print(h,R,tau,"D",D(R),"schur",a,"def",b,"OK" if a==b else "DIFF",round(time.time()-t0,1),flush=True)
                    if a!=0:
                        for v in ["lagidx","J","nosign"]:
                            print("   ctl",v,Xi_schur(h,R,tau,v)==b)
    if mode=="hi":
        # Schur-only high genus
        for h,Rs in [(5,[(2,),(3,),(4,),(5,),(8,),(9,),(2,2),(3,2),(4,2),(3,3),(2,2,2)]),(6,[(2,),(3,),(5,),(11,),(2,2)]),(7,[(2,),(4,),(13,)])]:
            for R in Rs:
                t0=time.time()
                tau=F(3,2)
                v=Xi_schur(h,R,tau)
                print("h",h,"R",R,"D",D(R),"inrange",D(R)<=2*h-3,"Xi",v,round(time.time()-t0,1),flush=True)
