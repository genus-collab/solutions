# Fresh Hodge integrals: Witten-Kontsevich by DVV, ch classes by Mumford (Faber's algorithm). stdlib only.
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial
import sys
sys.setrecursionlimit(100000)

def dfact(n):
    r=1
    while n>1:
        r*=n;n-=2
    return r

@lru_cache(maxsize=None)
def wk(g,a):
    a=tuple(sorted(a)); n=len(a)
    if g<0: return Fr(0)
    if 2*g-2+n<=0: return Fr(0)
    if sum(a)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and n==1: return Fr(1,24)
    if a[0]==0:
        rest=a[1:]; s=Fr(0)
        for i,x in enumerate(rest):
            if x>0:
                s+=wk(g,tuple(sorted(rest[:i]+(x-1,)+rest[i+1:])))
        return s
    if 1 in a:
        i=a.index(1); rest=a[:i]+a[i+1:]
        return (2*g-2+len(rest))*wk(g,rest)
    k=a[-1]; rest=a[:-1]; n1=len(rest)
    s=Fr(0)
    for i,x in enumerate(rest):
        oth=rest[:i]+rest[i+1:]
        s+=Fr(dfact(2*k+2*x-1),dfact(2*x-1))*wk(g,tuple(sorted(oth+(k+x-1,))))
    for p in range(0,k-1):
        q=k-2-p
        s+=Fr(1,2)*dfact(2*p+1)*dfact(2*q+1)*wk(g-1,tuple(sorted(rest+(p,q))))
    idx=list(range(n1))
    for p in range(0,k-1):
        q=k-2-p
        for mask in range(1<<n1):
            I=tuple(rest[i] for i in idx if mask>>i&1)
            J=tuple(rest[i] for i in idx if not mask>>i&1)
            dI=sum(I)+p
            x=dI+2-len(I)
            if x%3: continue
            g1=x//3; g2=g-g1
            if g1<0 or g2<0: continue
            v1=wk(g1,tuple(sorted(I+(p,))))
            if v1==0: continue
            v2=wk(g2,tuple(sorted(J+(q,))))
            s+=Fr(1,2)*dfact(2*p+1)*dfact(2*q+1)*v1*v2
    return s/dfact(2*k+1)

def bern(n):
    A=[Fr(0)]*(n+1)
    for m in range(n+1):
        A[m]=Fr(1,m+1)
        for j in range(m,0,-1):
            A[j-1]=j*(A[j-1]-A[j])
    return A[0]

@lru_cache(maxsize=None)
def H(g,a,K):
    """Int_{Mbar_{g,n}} prod psi^a prod_{k in K} ch_k. a, K sorted tuples, K odd integers."""
    a=tuple(sorted(a)); K=tuple(sorted(K)); n=len(a)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(a)+sum(K)!=3*g-3+n: return Fr(0)
    if not K: return wk(g,a)
    if g==0: return Fr(0)
    if K[-1]>=2*g: return Fr(0)
    k=K[-1]; Kr=K[:-1]; l=(k+1)//2
    pref=bern(2*l)/factorial(2*l)
    s=H(g,a+(k+1,),Kr)
    for i,x in enumerate(a):
        s-=H(g,tuple(sorted(a[:i]+(x+k,)+a[i+1:])),Kr)
    bd=Fr(0)
    for p in range(0,k):
        q=k-1-p; sg=(-1)**p
        bd+=sg*H(g-1,tuple(sorted(a+(p,q))),Kr)
        n_=len(a); nk=len(Kr)
        for mask in range(1<<n_):
            I=tuple(a[i] for i in range(n_) if mask>>i&1)
            J=tuple(a[i] for i in range(n_) if not mask>>i&1)
            for kmask in range(1<<nk):
                K1=tuple(Kr[i] for i in range(nk) if kmask>>i&1)
                K2=tuple(Kr[i] for i in range(nk) if not kmask>>i&1)
                x=sum(I)+p+sum(K1)+2-len(I)
                if x%3: continue
                g1=x//3; g2=g-g1
                if g1<0 or g2<0: continue
                v1=H(g1,tuple(sorted(I+(p,))),tuple(sorted(K1)))
                if v1==0: continue
                bd+=sg*v1*H(g2,tuple(sorted(J+(q,))),tuple(sorted(K2)))
    s+=Fr(1,2)*bd
    return pref*s

def odd_multisets(total,maxk):
    def rec(rem,mx):
        if rem==0:
            yield ()
            return
        k=mx if mx%2==1 else mx-1
        while k>=1:
            if k<=rem:
                for r in rec(rem-k,k):
                    yield (k,)+r
            k-=2
    yield from rec(total,maxk)
