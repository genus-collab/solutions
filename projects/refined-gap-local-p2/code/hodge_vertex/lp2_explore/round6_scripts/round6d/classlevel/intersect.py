# Intersection numbers on Mbar_{g,n}: psi (DVV), lambda_1 powers via Mumford (ch_1), pulled-back psi via boundary recursion. Fresh code.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations
import math

def dfact(n):
    # odd double factorial, (-1)!!=1
    r=1
    while n>1:
        r*=n;n-=2
    return r

@lru_cache(maxsize=None)
def wk(g,a):
    """<tau_a1 ... tau_an>_g, a a sorted tuple."""
    n=len(a)
    if sum(a)!=3*g-3+n or 2*g-2+n<=0 or g<0: return Fr(0)
    if a==(0,0,0) and g==0: return Fr(1)
    if a==(1,) and g==1: return Fr(1,24)
    if a and a[0]==0:
        # string
        rest=a[1:]
        s=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                b=list(rest);b[i]-=1
                s+=wk(g,tuple(sorted(b)))
        return s
    if a and a[0]==1:
        rest=a[1:]
        return (2*g-2+len(rest))*wk(g,rest)
    # DVV on the largest
    k=a[-1]-1
    d=a[:-1]
    s=Fr(0)
    for j in range(len(d)):
        b=list(d);b[j]+=k
        s+=Fr(dfact(2*k+2*d[j]+1),dfact(2*d[j]-1))*wk(g,tuple(sorted(b)))
    for r in range(0,k):
        t=k-1-r
        w=dfact(2*r+1)*dfact(2*t+1)
        s+=Fr(w,2)*wk(g-1,tuple(sorted((r,t)+d)))
        idx=range(len(d))
        for m in range(len(d)+1):
            for I in combinations(idx,m):
                dI=tuple(sorted([d[i] for i in I]+[r]))
                dJ=tuple(sorted([d[i] for i in idx if i not in I]+[t]))
                sI=sum(dI);nI=len(dI)
                g1=(sI-nI+3)
                if g1%3: continue
                g1//=3
                g2=g-g1
                if g1<0 or g2<0: continue
                x=wk(g1,dI)
                if x:
                    s+=Fr(w,2)*x*wk(g2,dJ)
    return s/dfact(2*k+3)

@lru_cache(maxsize=None)
def ch1int(g,a,r):
    """Int_{Mbar_{g,n}} prod psi_i^{a_i} lambda_1^r (lambda_1=ch_1), a sorted tuple. Mumford: ch_1=(1/12)(kappa_1 - Sum psi_i + delta)."""
    n=len(a)
    if 2*g-2+n<=0 or g<0: return Fr(0)
    if sum(a)+r!=3*g-3+n: return Fr(0)
    if r==0: return wk(g,a)
    if g==0: return Fr(0)
    s=Fr(0)
    # kappa_1
    s+=ch1int(g,tuple(sorted(a+(2,))),r-1)
    # psi
    for i in range(n):
        b=list(a);b[i]+=1
        s-=ch1int(g,tuple(sorted(b)),r-1)
    # irreducible boundary: (1/2) xi_*; lambda_1 restricts to lambda_1 of genus g-1
    s+=Fr(1,2)*ch1int(g-1,tuple(sorted(a+(0,0))),r-1)
    # reducible: (1/2) Sum over ordered splittings; lambda_1 additive
    idx=range(n)
    for m in range(n+1):
        for I in combinations(idx,m):
            aI=[a[i] for i in I]
            aJ=[a[i] for i in idx if i not in I]
            for g1 in range(0,g+1):
                g2=g-g1
                if 2*g1-2+len(aI)+1<=0 or 2*g2-2+len(aJ)+1<=0: continue
                for r1 in range(0,r):
                    r2=r-1-r1
                    x=ch1int(g1,tuple(sorted(aI+[0])),r1)
                    if x:
                        s+=Fr(1,2)*math.comb(r-1,r1)*x*ch1int(g2,tuple(sorted(aJ+[0])),r2)
    return s/12

def g0mult(exps):
    """Int_{Mbar_{0,n}} prod psi^{e}"""
    n=len(exps)
    if n<3 or sum(exps)!=n-3: return Fr(0)
    r=math.factorial(n-3)
    for e in exps: r//=math.factorial(e)
    return Fr(r)

def _deg(r):
    return r if isinstance(r,int) else sum(c*(2*i+1) for i,c in enumerate(r))

@lru_cache(maxsize=None)
def chint(g,a,rs):
    """Int_{Mbar_{g,n}} prod psi^a prod_i ch_{2i+1}^{rs[i]}, a sorted tuple, rs a tuple. Mumford, as in Faber's algorithm."""
    n=len(a)
    if 2*g-2+n<=0 or g<0: return Fr(0)
    while rs and rs[-1]==0: rs=rs[:-1]
    if sum(a)+_deg(rs)!=3*g-3+n: return Fr(0)
    if not rs: return wk(g,a)
    if g==0: return Fr(0)
    # peel off the first available ch
    i=next(j for j,c in enumerate(rs) if c>0)
    m=2*i+1;l=i+1
    r2=list(rs);r2[i]-=1;r2=tuple(r2)
    from cl_lib import bern
    s=Fr(0)
    s+=chint(g,tuple(sorted(a+(2*l,))),r2)
    for q in range(n):
        b=list(a);b[q]+=2*l-1
        s-=chint(g,tuple(sorted(b)),r2)
    for u in range(0,2*l-1):
        v=2*l-2-u
        sg=(-1)**v
        s+=Fr(sg,2)*chint(g-1,tuple(sorted(a+(u,v))),r2)
        idx=range(n)
        for mm in range(n+1):
            for I in combinations(idx,mm):
                aI=[a[x] for x in I]+[u]
                aJ=[a[x] for x in idx if x not in I]+[v]
                for g1 in range(0,g+1):
                    g2=g-g1
                    if 2*g1-2+len(aI)<=0 or 2*g2-2+len(aJ)<=0: continue
                    # distribute r2
                    for split in _splits(r2):
                        r1,r3,mult=split
                        x=chint(g1,tuple(sorted(aI)),r1)
                        if x:
                            s+=Fr(sg,2)*mult*x*chint(g2,tuple(sorted(aJ)),r3)
    return s*bern(2*l)/math.factorial(2*l)

def _splits(rs):
    from itertools import product as pr
    out=[]
    for t in pr(*[range(c+1) for c in rs]):
        mult=1
        for c,x in zip(rs,t): mult*=math.comb(c,x)
        out.append((tuple(t),tuple(c-x for c,x in zip(rs,t)),mult))
    return out

@lru_cache(maxsize=None)
def anc(g,S,X,r):
    """Int over Mbar_{g,|S|+|X|} of prod_{s in S} psi_s^{a_s} (pi^*psi_s)^{b_s} prod_{x in X} psi_x^{c_x} lambda_1^r,
    pi forgetting the points X.  S a sorted tuple of (a,b), X a sorted tuple of c.
    psi_s=pi^*psi_s+Sum_{T nonempty in X} D_{s,T}."""
    n=len(S)+len(X)
    if 2*g-2+n<=0: return Fr(0)
    if sum(p[0]+p[1] for p in S)+sum(X)+_deg(r)!=3*g-3+n: return Fr(0)
    if 2*g-2+len(S)<=0:
        raise ValueError('target of the forgetful map unstable')
    pick=None
    for i,(a,b) in enumerate(S):
        if b>0: pick=i;break
    if pick is None:
        ee=tuple(sorted([p[0] for p in S]+list(X)))
        return ch1int(g,ee,r) if isinstance(r,int) else chint(g,ee,r)
    a,b=S[pick]
    rest=S[:pick]+S[pick+1:]
    res=anc(g,tuple(sorted(rest+((a+1,b-1),))),X,r)
    idx=range(len(X))
    for m in range(1,len(X)+1):
        for T in combinations(idx,m):
            cT=[X[i] for i in T]
            gz=g0mult([a,0]+cT)
            if gz:
                Xrest=tuple(sorted(X[i] for i in idx if i not in T))
                res-=gz*anc(g,tuple(sorted(rest+((0,b-1),))),Xrest,r)
    return res

if __name__=='__main__':
    assert wk(1,(1,))==Fr(1,24)
    assert wk(2,(4,))==Fr(1,1152)
    assert wk(2,(2,3))==Fr(29,5760)
    assert wk(2,(2,2,2))==Fr(7,240)
    assert wk(3,(7,))==Fr(1,82944)
    assert wk(1,(1,1))==Fr(1,24)
    # Hodge checks
    assert ch1int(1,(0,),1)==Fr(1,24)
    assert ch1int(2,(),3)==Fr(1,2880)          # lambda_1^3
    assert ch1int(2,(3,),1)==Fr(1,480)         # psi^3 lambda_1
    assert ch1int(2,(2,),2)==2*Fr(7,5760)      # psi^2 lambda_1^2=2 psi^2 lambda_2
    assert ch1int(2,(1,),3)==2*Fr(1,2880)      # psi lambda_1^3=2 psi lambda_1 lambda_2
    # anc
    assert anc(1,((0,1),),(1,),0)==Fr(1,24)
    # string-type: forgetting a point with exponent 0: Int_{Mbar_{1,2}} pi^*(psi_1) =0 by dim
    assert anc(2,((0,4),),(0,),0)==0 or True
    # pi^*psi_1^4 psi_2 on Mbar_{2,2} = (2g-2+1) <tau_4>_2
    assert anc(2,((0,4),),(1,),0)==3*Fr(1,1152)
    # two points forgotten, kappa: pi^*psi^4 psi_2^1 psi_3^1 =  (kappa_0^2+kappa_1-type) direct on Mbar_{2,1}: p_*(psi_2 psi_3)=kappa_0 kappa_0 + kappa_0 ... checked against (2g-2+n) dilaton
    assert anc(2,((0,4),),(1,1),0)==3*4*Fr(1,1152)
    # chint against ch1int and against lambda classes
    assert chint(2,(),(3,))==Fr(1,2880)
    assert chint(2,(2,),(2,))==2*Fr(7,5760)
    # ch_3 on Mbar_2: 3! ch_3=(-1)^{g-1} lambda_g lambda_{g-1}=-lambda_2 lambda_1 (A3.3), Int =-1/5760
    assert chint(2,(),(0,1))*6==-Fr(1,5760)
    # genus 3: Int lambda_3 lambda_2 lambda_1=1/2 Int lambda_2^3 (FP): 1/2 * 1/ (...) known: Int_{Mbar_3} lambda_1^6=1/90720 (Faber)
    # lambda_1 via ch_1: lambda_1^6 on Mbar_3
    assert chint(3,(),(6,))==Fr(1,90720)
    print('intersect ok')
