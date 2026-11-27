# Own MV-side engine: F=log Sum_nu e^{tau x f2} s_nu(q^rho) s_nu(p), restricted to monomials p_{R'} p_1^a.
# Mode: modulus=None -> exact Fractions; modulus=p -> arithmetic mod p.
from common import *
from itertools import product
import time

class Fld:
    def __init__(self,p=None):
        self.p=p
    def c(self,x):
        x=Fr(x)
        if self.p is None: return x
        return x.numerator%self.p*pow(x.denominator,-1,self.p)%self.p
    def zero(self):
        return Fr(0) if self.p is None else 0
    def inv(self,x):
        if self.p is None: return 1/x
        return pow(x,-1,self.p)
    def smul(self,a,b,P):
        """truncated product of series a,b (lists length P)."""
        p=self.p
        if p is None:
            r=[Fr(0)]*P
            for i,x in enumerate(a):
                if x:
                    for j in range(P-i):
                        y=b[j]
                        if y: r[i+j]+=x*y
            return r
        r=[0]*P
        for i,x in enumerate(a):
            if x:
                for j in range(P-i):
                    r[i+j]+=x*b[j]
        return [v%p for v in r]

def W_series(nu,tau,P,K):
    """coefficients of x^{-N+i}, i<P, of e^{tau x f2(nu)} s_nu(q^rho) via hook form (Prop 2.1)"""
    N=sum(nu)
    f2=contents_sum(nu)
    hs=hooks(nu)
    h0=1
    for h in hs: h0*=h
    # log coefficients
    l=[K.zero()]*P
    if P>1: l[1]=K.c((tau+Fr(1,2))*f2)
    for lidx in range(1,P//2+1):
        ii=2*lidx
        if ii>=P: break
        ps=sum(h**ii for h in hs)
        l[ii]=K.c(-BERN[ii]/(ii*factorial(ii))*ps)
    e=[K.zero()]*P
    e[0]=K.c(Fr(1,h0))
    nzi=[i for i in range(1,P) if l[i]]
    invs=INV(P,K)
    if K.p is None:
        for n in range(1,P):
            s=Fr(0)
            for i in nzi:
                if i>n: break
                s+=i*l[i]*e[n-i]
            e[n]=s/n
    else:
        p=K.p
        for n in range(1,P):
            s=0
            for i in nzi:
                if i>n: break
                s+=i*l[i]%p*e[n-i]
            e[n]=s%p*invs[n]%p
    return e

_inv={}
def INV(P,K):
    if K.p is None: return None
    key=K.p
    t=_inv.get(key,[0,1])
    while len(t)<P+1:
        t.append(pow(len(t),-1,K.p))
    _inv[key]=t
    return t

def compute_F(R,amax,P,tau,K,verbose=False):
    """Returns dict (Rsub exponents key tuple sorted parts, a) -> list F coefficient, index i = power x^{-W+i}, W=J'+a.
    R: tuple of parts >=2."""
    R=tuple(sorted(R,reverse=True))
    items=sorted(Counter(R).items(),reverse=True)
    # exponent vectors
    ranges=[range(m+1) for _,m in items]
    subs=[]
    for ev in product(*ranges):
        sub=()
        for (v,_),e in zip(items,ev):
            sub+=(v,)*e
        subs.append((ev,sub))
    W={}
    Wneeded=sorted(set(sum(sub)+a for ev,sub in subs for a in range(amax+1)))
    t0=time.time()
    for N in Wneeded:
        for nu in partitions(N):
            W[nu]=W_series(nu,tau,P,K)
    if verbose: print("W done",len(W),time.time()-t0,flush=True)
    Z={}
    for ev,sub in subs:
        for a in range(amax+1):
            mu=sub+(1,)*a
            N=sum(mu)
            zz=zee(mu)
            acc=[K.zero()]*P
            if N==0:
                acc[0]=K.c(1)
                Z[(ev,a)]=acc
                continue
            for nu in partitions(N):
                ch=chi(nu,mu)
                if ch:
                    w=W[nu]
                    cc=K.c(Fr(ch,zz))
                    if K.p is None:
                        for i in range(P):
                            acc[i]+=cc*w[i]
                    else:
                        for i in range(P):
                            acc[i]+=cc*w[i]
            if K.p is not None:
                acc=[v%K.p for v in acc]
            Z[(ev,a)]=acc
    if verbose: print("Z done",time.time()-t0,flush=True)
    F={}
    # order by weight
    keys=sorted(Z.keys(),key=lambda t:(sum(t[0])+t[1]))
    for ev,a in keys:
        wt=sum(ev)+a
        if wt==0: continue
        res=[K.c(wt)*x for x in Z[(ev,a)]]
        for ev1 in product(*[range(e+1) for e in ev]):
            for a1 in range(a+1):
                if (sum(ev1)+a1)==0: continue
                ev2=tuple(x-y for x,y in zip(ev,ev1)); a2=a-a1
                if sum(ev2)+a2==0: continue
                w1=sum(ev1)+a1
                pr=K.smul(F[(ev1,a1)],Z[(ev2,a2)],P)
                cw=K.c(w1)
                for i in range(P):
                    res[i]-=cw*pr[i]
        iw=K.inv(K.c(wt))
        if K.p is None:
            F[(ev,a)]=[x*iw for x in res]
        else:
            F[(ev,a)]=[x%K.p*iw%K.p for x in res]
    if verbose: print("F done",time.time()-t0,flush=True)
    return F,items

def key_for(items,sub):
    c=Counter(sub)
    return tuple(c.get(v,0) for v,_ in items)

def coeff(F,items,sub,a,power):
    ev=key_for(items,sub)
    W=sum(sub)+a
    i=power+W
    L=F[(ev,a)]
    if i<0: return 0
    return L[i]

def laguerre(n,al,x):
    """L_n^{(al)}(x)=Sum_i (-1)^i binom(n+al,n-i) x^i/i!"""
    s=Fr(0)
    for i in range(n+1):
        s+=(-1)**i*comb(n+al,n-i)*Fr(x)**i/factorial(i)
    return s

def Xi_from_F(F,items,R,h,tau,K):
    """s.8.2: Xi_h(R)"""
    R=tuple(R); r=len(R); J=sum(R); N=3*h-3+r
    sigma=tau*(1+tau)
    PiR=Fr(1)
    for j in R: PiR*=Pi(j,tau)
    pref=K.c(Fr(aut(R))/(sigma**(h+r-1)*PiR))
    tot=K.zero()
    sg=K.c(1/sigma)
    for a in range(N+1):
        L=K.c(laguerre(N-a,2*h-2+r+a,-J-a))
        c=coeff(F,items,R,a,2*h-2+r+a)
        term=K.c((-1)**a*(1/sigma)**a)*L*c
        tot+=term
    r2=tot*pref
    return r2 if K.p is None else r2%K.p
