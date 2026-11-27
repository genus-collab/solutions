# Test: for ARBITRARY B_s and arbitrary higher-grade terms, the grade-1 part of (*) in monomials x^a w^b
# equals L[B_s] - RHS(E1), and grade<=0 vanish.  No Hodge data, no beta. Controls: wrong a/log/quadratic coefficient.
from fractions import Fraction as Fr
import random, sys
from math import factorial
from collections import defaultdict
def run(tau,X,seed,dgen1=Fr(0),dlog=Fr(0),dg0=Fr(0),verbose=True):
    random.seed(seed)
    e=1+tau; sig=tau*e; sp=1+2*tau
    K=X+6
    rnd=lambda: Fr(random.randint(-9,9),random.randint(1,7))
    # Phi = DF, Psi = F_tau|_y mod constants, as dict (a,b)->coef, x^a w^b
    Phi=defaultdict(Fr); Psi=defaultdict(Fr)
    def add(d,a,b,c):
        if a<=X+4: d[(a,b)]+=c
    def D(d):
        o=defaultdict(Fr)
        for (a,b),c in d.items():
            if b: o[(a,b+2)]+=b*c; o[(a,b+1)]-=b*c
        return o
    # genus 0: F=x^-2/sig*(-v^2/4-v^3/6) , v=w^-1
    g0=defaultdict(Fr); 
    g0[(-2,-2)]=(Fr(-1,4)+dg0)/sig; g0[(-2,-3)]=Fr(-1,6)/sig
    # genus 1: (1/24)(-log v+(1+1/sig) v): Dlog v: D(-log v)= -(1-1/v)/v = (w-1)... compute directly
    g1=defaultdict(Fr); g1[(0,-1)]=(1+1/sig)/24+dgen1
    for k,v in D(g0).items(): Phi[k]+=v
    for k,v in D(g1).items(): Phi[k]+=v
    # D(-(1/24+dlog) log v): log v=-log w ; D log w=(w-1)w => D(-log v)=D(log w)=w^2-w
    Phi[(0,2)]+=Fr(1,24)+dlog; Phi[(0,1)]-=Fr(1,24)+dlog
    dinv=-sp/sig**2 # d/dtau 1/sig
    Psi[(-2,-2)]=(Fr(-1,4)+dg0)*dinv; Psi[(-2,-3)]=Fr(-1,6)*dinv
    # derivative of g1 coefficient: (1/24)*dinv (control shift treated tau independent)
    Psi[(0,-1)]=dinv/24
    # arbitrary B_s
    bk={k:rnd() for k in range(1,K)}; bkp={k:rnd() for k in range(1,K)}
    Fh=defaultdict(Fr)
    for k in range(1,K):
        Fh[(k,k)]=bk[k]; Psi[(k,k)]+=bkp[k]
    # arbitrary higher grade terms, a>=-1, grade 1..3
    for a in range(-1,K):
        for g in (1,2,3):
            if a+g==0: continue
            Fh[(a,a+g)]+=rnd(); Psi[(a,a+g)]+=rnd()
    for k,v in D(Fh).items(): Phi[k]+=v
    def trunc(d,mx):
        return defaultdict(Fr,{k:v for k,v in d.items() if v!=0 and k[0]<=mx})
    Phi=trunc(Phi,X+4); Psi=trunc(Psi,X+4)
    def shift(d,s,co=Fr(1)):
        o=defaultdict(Fr)
        for (a,b),c in d.items(): o[(a+s,b)]+=c*co
        return o
    def addd(*ds):
        o=defaultdict(Fr)
        for d in ds:
            for k,v in d.items(): o[k]+=v
        return o
    def exp_minus1(d):
        # (e^{tau x D}-1) d
        o=defaultdict(Fr); cur=d
        for n in range(1,X+10):
            cur=trunc(D(cur),X+4)
            amin=min((k[0] for k in cur),default=None)
            if amin is None: break
            # term (tau x)^n/n! D^n
            co=tau**n/factorial(n)
            for (a,b),c in cur.items():
                if a+n<=X+4: o[(a+n,b)]+=c*co
        return o
    Ftau=addd(Psi,shift(Phi,0,sp/sig))
    T1=shift(exp_minus1(Ftau),-1,Fr(-2))
    em=exp_minus1(Phi)
    Sum=addd(Phi,Phi,em)  # (1+e^{bD})Phi = 2Phi+em
    T2=addd(D(Sum),shift(Sum,0,Fr(-1)))
    # T3 square of -em
    T3=defaultdict(Fr)
    items=list(em.items())
    for (a,b),c in items:
        for (a2,b2),c2 in items:
            if a+a2<=X: T3[(a+a2,b+b2)]+=c*c2
    LHS=addd(T1,T2,T3)
    # residual by grade
    gr=defaultdict(dict)
    for (a,b),c in LHS.items():
        if c!=0 and a<=X: gr[b-a][(a,b)]=c
    # expected grade-1: N^-b coefficient
    # series in 1/N: dict p->coef of N^-p (p may be negative)
    P=X+8
    def ser_E(p):
        # (N+tau)^-p expanded in 1/N, p>=1
        o={}
        for j in range(0,P+2-p):
            co=Fr(1)
            for i in range(j): co*=Fr(-p-i,i+1)
            o[p+j]=co*tau**j
        return o
    L=defaultdict(Fr)
    for k in range(1,K):
        # -2 b'_k Delta N^-k
        for q,c in ser_E(k).items(): L[q]+=-2*bkp[k]*c
        L[k]+=2*bkp[k]
        # B' = -k b_k N^-k-1
        cp=-k*bk[k]; q0=k+1
        # -(1+E)B'
        for q,c in ser_E(q0).items(): L[q]+=-cp*c
        L[q0]+=-cp
        # -((M+N)/e) Delta B'
        dl=defaultdict(Fr)
        for q,c in ser_E(q0).items(): dl[q]+=cp*c
        dl[q0]-=cp
        for q,c in dl.items():
            # (2N+tau) N^-q
            L[q-1]+=-2*c/e; L[q]+=-tau*c/e
    # RHS (E1)
    R=defaultdict(Fr)
    for j in range(1,P+4):
        c=Fr((-1)**(j-1),j)*tau**j  # l coefficient N^-j
        # -l*(MN - sp/12)= -l*(N^2+tau N - sp/12)
        R[j-2]+=-c; R[j-1]+=-c*tau; R[j]+=c*sp/12
    # tau(M+N)/2= tau(2N+tau)/2
    R[-1]+=tau; R[0]+=tau*tau/2
    for q,c in ser_E(1).items(): R[q]+=-(tau*(1+2*sig)/24)*c
    R[1]+=-(tau*(1+2*sig)/24)
    bad=0; n=0
    # compare grade 1: x^a w^(a+1) <-> N^-(a+1)
    exp=defaultdict(Fr)
    for q,c in L.items(): exp[q]+=c
    for q,c in R.items(): exp[q]-=c/sig**2
    for b in range(-6,X+2):
        a=b-1
        if a>X: continue
        n+=1
        if gr[1].get((a,b),Fr(0))!=exp.get(b,Fr(0)): bad+=1
    neg=sum(len(v) for g,v in gr.items() if g<=0)
    if verbose: print(f"tau={tau} X={X} seed={seed} dgen1={dgen1} dlog={dlog} dg0={dg0}: nonzero monomials grade<=0: {neg}; grade 1 compared {n}, mismatches {bad}; grade-2 nonzero (not constrained) {len(gr[2])}")
    return neg,bad
for tau in [Fr(3,11),Fr(-7,4)]:
    for s in (1,2):
        run(tau,10,s)
run(Fr(3,11),10,1,dgen1=Fr(1,7))
run(Fr(3,11),10,1,dlog=Fr(1,7))
run(Fr(3,11),10,1,dg0=Fr(1,7))
