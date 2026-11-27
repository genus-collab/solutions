# Evidence only: fresh check of (4.2),(4.3) with Z_0,Z_1 built from Murnaghan-Nakayama characters.
from fractions import Fraction as Fr
from itertools import combinations
import sys
sys.setrecursionlimit(10000)
def partitions(n,m=None):
    m=n if m is None else m
    if n==0: yield (); return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k): yield (k,)+r
def beta(nu,L): return [ (nu[i] if i<len(nu) else 0)+L-1-i for i in range(L)]
def chi(nu,mu):
    # M-N by removing rim hooks with beta numbers
    if not mu: return 1 if sum(nu)==0 else 0
    L=sum(nu)+1
    b=beta(nu,L); s=set(b); k=mu[0]; t=0
    for x in b:
        y=x-k
        if y>=0 and y not in s:
            h=sum(1 for z in b if y<z<x)
            nb=sorted([z for z in b if z!=x]+[y],reverse=True)
            nn=tuple(v for v in (nb[i]-(L-1-i) for i in range(L)) if v>0)
            t+=(-1)**h*chi(nn,mu[1:])
    return t
def f2(nu): return sum(j-i for i,r in enumerate(nu) for j in range(r))
def zmu(mu):
    z=1
    for k in set(mu):
        m=mu.count(k)
        for i in range(1,m+1): z*=i
        z*=k**m
    return z
# polys in t_1..t_W : dict exponent tuple -> Fr. Use monomials in p_mu -> t_mu
W=9
def mono(mu):
    e=[0]*W
    for k in mu: e[k-1]+=1
    return tuple(e)
def schur(nu):
    d={}
    for mu in partitions(sum(nu)):
        c=Fr(chi(nu,mu),zmu(mu))
        for k in mu: c*=k   # p_k=k t_k
        if c: d[mono(mu)]=c
    return d
def add(d,e,c=1):
    for k,v in e.items():
        d[k]=d.get(k,0)+c*v
        if d[k]==0: del d[k]
def mul(a,b,maxw):
    r={}
    for k1,v1 in a.items():
        w1=sum((i+1)*e for i,e in enumerate(k1))
        for k2,v2 in b.items():
            w2=sum((i+1)*e for i,e in enumerate(k2))
            if w1+w2>maxw: continue
            k=tuple(x+y for x,y in zip(k1,k2))
            r[k]=r.get(k,0)+v1*v2
    return r
def der(a,i):
    r={}
    for k,v in a.items():
        if k[i-1]>0:
            e=list(k); e[i-1]-=1
            r[tuple(e)]=v*k[i-1]
    return r
def build(A,pp,n,wrong=False):
    Z={}
    def sp(nu):
        s=Fr(0)
        for mu in partitions(sum(nu)):
            c=Fr(chi(nu,mu),zmu(mu))
            for k in mu: c*=pp[k]
            s+=c
        return s
    for d in range(0,W+1):
        for nu in partitions(d):
            c=(A**(f2(nu)+n*d))*sp(nu)
            add(Z,schur(nu),c)
    return Z
A=Fr(3,2); pp={k:Fr(2*k+1,k+3) for k in range(1,W+1)}
Z0=build(A,pp,0); Z1=build(A,pp,1)
Zbad=build(A,pp,-1)
def hirota(alpha,X,Y,maxw):
    r={}
    n=len(alpha)
    for m in range(n+1):
        for S in combinations(range(n),m):
            a=X; b=Y
            for i in range(n):
                if i in S: a=der(a,alpha[i])
                else: b=der(b,alpha[i])
            sign=(-1)**(n-m)
            add(r,mul(a,b,maxw),sign)
    return r
def c_alpha(al,k):
    c=Fr(1)
    for i in set(al):
        m=al.count(i)
        c*=Fr(-1,i)**m
        for j in range(1,m+1): c/=j
    if al==(k,): c+=1
    return c
for k in range(1,8):
    maxw=W-k
    for name,Y in (('Z1',Z1),('Z_{-1}',Zbad)):
        tot={}
        for al in partitions(k):
            add(tot,hirota(al,Z0,Y,maxw),c_alpha(al,k))
        print('P_%d'%k,name,'nonzero coeffs through weight',maxw,':',len(tot))
