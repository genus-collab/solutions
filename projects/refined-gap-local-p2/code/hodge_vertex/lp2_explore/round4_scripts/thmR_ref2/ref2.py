from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import itertools
E=F(7,3)
# trivial CohFT correlators genus 0,1 via string/dilaton
@lru_cache(None)
def corr(g,ds):
    ds=tuple(sorted(ds)); n=len(ds)
    if sum(ds)!=3*g-3+n or any(d<0 for d in ds): return F(0)
    if g==0 and n==3: return F(1)
    if g==1 and ds==(1,): return F(1,24)
    if 2*g-2+n<=0: return F(0)
    if ds[0]==0:
        rest=list(ds[1:]); s=F(0)
        for i in range(len(rest)):
            if rest[i]>0:
                r=rest[:]; r[i]-=1; s+=corr(g,tuple(r))
        return s
    if ds[0]==1:
        return (2*g-2+n-1)*corr(g,ds[1:])
    raise Exception("unreduced",g,ds)
def lamc(j,a): return F(j**a)/E**(a+1)
def Tc(a): return F(0) if a<2 else -F(1)/E**(a-1)
def Xi(h,R,extra=None):
    # extra: list of coefficient functions for extra insertions
    ins=[ (lambda a,j=j: lamc(j,a)) for j in R]+(extra or [])
    r=len(ins); tot=F(0)
    for m in range(0,40):
        n=r+m; dim=3*h-3+n
        if dim<0: continue
        if 2*m>dim and m>0: break
        s=F(0)
        def rec(i,rem,acc,ds):
            nonlocal s
            if i==n:
                if rem==0: s+=acc*corr(h,tuple(ds))
                return
            for a in range(rem+1):
                c= ins[i](a) if i<r else Tc(a)
                if c: rec(i+1,rem-a,acc*c,ds+[a])
        rec(0,dim,F(1),[])
        tot+=s/factorial(m)
    return tot
ok=True
for h in (0,1):
    for R in [(2,3,5),(2,2,4),(3,3),(2,5,1),(4,2,2),(2,),(3,)]:
        if h==0 and len(R)<3: continue
        X=Xi(h,R); J=sum(R); D=sum(j-1 for j in R)
        X1=Xi(h,(1,)+R); X0=Xi(h,(0,)+R)
        s1 = X1==F(J)/E**2*X; s0 = X0==F(D-2*h+2)/E**2*X
        c1 = X1==F(J)/E*X;  c0 = X0==F(D)/E**2*X if h==0 else X0==F(D+2)/E**2*X
        print(h,R,"X=",X,"string",s1,"lam0",s0,"controls(should be False unless X=0)",c1,c0)
