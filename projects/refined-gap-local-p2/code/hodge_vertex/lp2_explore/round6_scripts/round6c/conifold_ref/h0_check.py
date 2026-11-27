# Evidence only: (H0) of s.13.0 at h = 2, exact Fractions, truncated in t; also Teleman-free hand items.
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import product
from math import factorial
import sys
sys.setrecursionlimit(10000)
def dfact(n):
    r=1
    while n>1: r*=n; n-=2
    return r
@lru_cache(None)
def wk(g,ds):
    # <tau_d1...tau_dn>_g by DVV; ds sorted tuple
    n=len(ds)
    if 2*g-2+n<=0 or sum(ds)!=3*g-3+n or (n and min(ds)<0): return Fr(0)
    if n==0: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    ds=list(ds)
    k=ds.pop()  # largest
    res=Fr(0)
    for j in range(len(ds)):
        d=ds[j]; rest=ds[:j]+ds[j+1:]
        res+=Fr(dfact(2*k+2*d-1),dfact(2*d-1))*wk(g,tuple(sorted(rest+[d+k-1])))
    if g>=1:
        for a in range(k-1):
            b=k-2-a
            res+=Fr(1,2)*dfact(2*a+1)*dfact(2*b+1)*wk(g-1,tuple(sorted(ds+[a,b])))
    n1=len(ds)
    for mask in range(1<<n1):
        I=[ds[i] for i in range(n1) if mask>>i&1]; J=[ds[i] for i in range(n1) if not mask>>i&1]
        for a in range(k-1):
            b=k-2-a
            for g1 in range(g+1):
                v1=wk(g1,tuple(sorted(I+[a]))); 
                if v1==0: continue
                v2=wk(g-g1,tuple(sorted(J+[b])))
                res+=Fr(1,2)*dfact(2*a+1)*dfact(2*b+1)*v1*v2
    return res/dfact(2*k+1)
assert wk(2,(4,))==Fr(1,1152) and wk(2,(2,3))==Fr(29,5760) and wk(2,(2,2,2))==Fr(7,240)
h=2; c=Fr(3,2); Nt=3
# t-series: list of length Nt+1
def tmul(a,b):
    r=[Fr(0)]*(Nt+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if i+j<=Nt: r[i+j]+=x*y
    return r
# psi-series with t-series coefficients: dict d -> tseries
PM=14
def lam(j): return {d:[Fr(j)**d/c**(d+1)]+[Fr(0)]*Nt for d in range(PM+1)}
T={d:[-c**(1-d)]+[Fr(0)]*Nt for d in range(2,PM+1)}
q0={d:list(T.get(d,[Fr(0)]*(Nt+1))) for d in range(PM+1)}
q0[1]=[Fr(-1)]+[Fr(0)]*Nt
def add(a,b):
    r={d:list(v) for d,v in a.items()}
    for d,v in b.items():
        r[d]=[x+y for x,y in zip(r.get(d,[Fr(0)]*(Nt+1)),v)]
    return r
def sc(ts,a): return {d:tmul(ts,v) for d,v in a.items()}
# class Omega = 1 (a pull-back; (H0) uses only the dilaton equation)
def dd(slots):
    s=len(slots); tot=[Fr(0)]*(Nt+1)
    dim0=3*h-3+s
    for m in range(0,dim0+1):
        dim=dim0+m
        # enumerate degrees of T insertions as sorted partitions with weight
        def recT(k,left,minimum,cur):
            if k==0: yield list(cur); return
            for d in range(minimum,left+1):
                cur.append(d); yield from recT(k-1,left-d,d,cur); cur.pop()
        for Tdeg in recT(m,dim,2,[]):
            sT=sum(Tdeg)
            if sT>dim: continue
            # weight: 1/m! * number of orderings * prod T coeff
            from collections import Counter
            cnt=Counter(Tdeg); mult=factorial(m)
            for v in cnt.values(): mult//=factorial(v)
            w=Fr(mult,factorial(m))
            for d in Tdeg: w*=T[d][0]
            # slots degrees
            def recS(i,left,val,degs):
                nonlocal tot
                if i==s:
                    if left==0:
                        I=wk(h,tuple(sorted(degs+Tdeg)))
                        if I:
                            for j in range(Nt+1): tot[j]+=w*I*val[j]
                    return
                for d,v in slots[i].items():
                    if d<=left and any(v):
                        recS(i+1,left-d,tmul(val,v),degs+[d])
            recS(0,dim-sT,[Fr(1)]+[Fr(0)]*Nt,[])
    return tot
def dF(q,fs):
    delta=add(q,sc([Fr(-1)]+[Fr(0)]*Nt,q0))
    assert all(v[0]==0 for v in delta.values())
    tot=[Fr(0)]*(Nt+1)
    for r in range(0,Nt+1):
        v=dd(list(fs)+[delta]*r)
        for j in range(Nt+1): tot[j]+=v[j]/factorial(r)
    return tot
def tpow(lamb,e):
    # lamb=1+O(t), integer e possibly negative
    inv=[Fr(0)]*(Nt+1); inv[0]=Fr(1)
    for k in range(1,Nt+1):
        inv[k]=-sum(lamb[i]*inv[k-i] for i in range(1,k+1))
    b=lamb if e>=0 else inv
    r=[Fr(1)]+[Fr(0)]*Nt
    for _ in range(abs(e)): r=tmul(r,b)
    return r
lamb=[Fr(1),Fr(2,3),Fr(-1,5),Fr(7,11)]
# q with a psi-linear t-term and a leaf at t^1, t^2
q=dict((d,list(v)) for d,v in q0.items())
q=add(q,sc([Fr(0),Fr(1,2),Fr(0),Fr(0)],lam(2)))
q=add(q,sc([Fr(0),Fr(0),Fr(-3,7),Fr(0)],lam(3)))
q[1][1]+=Fr(5,9)
for n,fs in [(0,[]),(1,[lam(3)])]:
    L=dF(sc(lamb,q),fs); R=tmul(tpow(lamb,2-2*h-n),dF(q,fs))
    print("n=",n,"H0 residual",[a-b for a,b in zip(L,R)])
    print("  values",L[:2])
