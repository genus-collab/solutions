# Referee own Schur side: Murnaghan-Nakayama on partitions (rim hooks by beta numbers), hook products by
# exp(log series). No Hodge input.
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from hodge import BB

@lru_cache(maxsize=None)
def partitions(n,m=None):
    if m is None: m=n
    if n==0: return (( ),)
    out=[]
    for p in range(min(n,m),0,-1):
        for r in partitions(n-p,p):
            out.append((p,)+r)
    return tuple(out)

def hooks(nu):
    conj=[sum(1 for a in nu if a>j) for j in range(nu[0])] if nu else []
    return [nu[i]-j+conj[j]-i-1 for i in range(len(nu)) for j in range(nu[i])]

def kappa(nu):
    return sum(a*(a-2*(i+1)+1) for i,a in enumerate(nu))

@lru_cache(maxsize=None)
def chi(nu,mu):
    if not mu: return 1 if not nu else 0
    k=mu[0]; rest=mu[1:]
    l=len(nu)
    beta=[nu[i]+l-1-i for i in range(l)]
    bs=set(beta)
    tot=0
    for i,b in enumerate(beta):
        c=b-k
        if c>=0 and c not in bs:
            ht=sum(1 for z in beta if c<z<b)
            nb=sorted([z for z in beta if z!=b]+[c],reverse=True)
            nn=tuple(v for v in (nb[j]-(l-1-j) for j in range(l)) if v>0)
            tot+=(-1)**ht*chi(nn,rest)
    return tot

class Sers:
    pass

def mul(a,b,P):
    r=[F(0)]*P
    for i,x in enumerate(a):
        if x==0 or i>=P: continue
        for j,y in enumerate(b):
            if i+j>=P: break
            if y!=0: r[i+j]+=x*y
    return r

def exp_series(L,P):
    # L[0]==0
    e=[F(0)]*P; e[0]=F(1)
    for n in range(1,P):
        s=F(0)
        for k in range(1,n+1):
            if k<len(L) and L[k]!=0:
                s+=k*L[k]*e[n-k]
        e[n]=s/n
    return e

_ycache={}
def Y(nu,tau,P):
    """x^{|nu|} * e^{(tau+1/2)kappa x/2} prod 1/(2 sinh(h x/2)) as power series to O(x^P)"""
    key=(nu,tau,P)
    if key in _ycache: return _ycache[key]
    hs=hooks(nu) if nu else []
    L=[F(0)]*P
    if P>1:
        L[1]=(F(tau)+F(1,2))*kappa(nu)/2
    k=1
    while 2*k<P:
        L[2*k]=-BB[2*k]/(2*k*factorial(2*k))*sum(h**(2*k) for h in hs)
        k+=1
    e=exp_series(L,P)
    den=1
    for h in hs: den*=h
    r=[c/den for c in e]
    _ycache[key]=r
    return r

def Z_label(mu,tau,P):
    """ |Aut mu| [p_mu] Z = sum_nu chi_nu(mu)/prod mu_i Y_nu, relative series of x^{|mu|}*... """
    n=sum(mu)
    mu=tuple(sorted(mu,reverse=True))
    pr=1
    for m in mu: pr*=m
    r=[F(0)]*P
    for nu in partitions(n):
        c=chi(nu,mu)
        if c==0: continue
        y=Y(nu,tau,P)
        for i in range(P): r[i]+=c*y[i]
    return [v/pr for v in r]

def hatc(mu,tau,P):
    """hat c_mu := |Aut mu| [p_mu] C, returned as list a with x^{-|mu|+i} coefficient a[i], i<P.
    cumulants over labelled parts."""
    l=len(mu)
    Zs={}
    cs={}
    def Zm(m):
        if m not in Zs:
            sub=tuple(mu[i] for i in range(l) if m>>i&1)
            Zs[m]=Z_label(sub,tau,P) if m else [F(1)]+[F(0)]*(P-1)
        return Zs[m]
    for m in range(1,1<<l):
        low=m&-m
        rest=m^low
        r=list(Zm(m))
        sub=rest
        while True:
            # proper B = low|sub' with sub' subset of rest, complement non-empty
            sp=rest^sub  # sp subset of rest ; B=low|sp
            B=low|sp
            if B!=m:
                pass
            if sub==0: break
            sub=(sub-1)&rest
        # enumerate properly
        for sp in submasks(rest):
            B=low|sp
            if B==m: continue
            pp=mul(cs[B],Zm(m^B),P)
            for i in range(P): r[i]-=pp[i]
        cs[m]=r
    return cs[(1<<l)-1]

def submasks(m):
    s=m
    while True:
        yield s
        if s==0: break
        s=(s-1)&m

def Pi(m,tau):
    tau=F(tau)
    r=F(1)
    for a in range(1,m):
        r*=(m*tau+a)
    return r/factorial(m-1)
