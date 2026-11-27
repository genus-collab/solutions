# Independent Hodge-side engine: Witten-Kontsevich via DVV, ch-classes via Mumford (Faber algorithm).
# Nothing here uses Marino-Vafa, the Toda hierarchy, or the project code.
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, factorial
from itertools import product
sys.setrecursionlimit(100000)

def dfact(n):
    r=1
    while n>1:
        r*=n; n-=2
    return r

def bern(n,_c={}):
    if n in _c: return _c[n]
    if n==0: r=Fr(1)
    else:
        r=Fr(0)
        for k in range(n):
            r-=comb(n+1,k)*bern(k)
        r/= (n+1)
    _c[n]=r
    return r

def submultisets(t):
    """t sorted tuple; yield (sub, complement, multiplicity) over labelled subsets grouped by sub-multiset"""
    vals=sorted(set(t)); cnt=[t.count(v) for v in vals]
    for choice in product(*[range(c+1) for c in cnt]):
        mult=1; a=[];b=[]
        for v,c,k in zip(vals,cnt,choice):
            mult*=comb(c,k); a+= [v]*k; b+=[v]*(c-k)
        yield tuple(a),tuple(b),mult

@lru_cache(maxsize=None)
def wk(g,ds):
    n=len(ds)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds[0]==0: # string
        rest=ds[1:]; s=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0 and (i==0 or rest[i]!=rest[i-1]):
                c=rest.count(rest[i])
                l=list(rest); l[i]-=1
                s+=c*wk(g,tuple(sorted(l)))
        return s
    if 1 in ds: # dilaton
        l=list(ds); l.remove(1)
        return (2*g-2+n-1)*wk(g,tuple(l))
    # DVV on largest
    d=ds[-1]; k=d-1; rest=ds[:-1]
    s=Fr(0)
    seen=set()
    for i,di in enumerate(rest):
        if di in seen: continue
        seen.add(di); c=rest.count(di)
        l=list(rest); l[i]=di+k
        s+=c*Fr(dfact(2*k+2*di+1),dfact(2*di-1))*wk(g,tuple(sorted(l)))
    for a in range(k):
        b=k-1-a
        co=Fr(dfact(2*a+1)*dfact(2*b+1),2)
        s+=co*wk(g-1,tuple(sorted(rest+(a,b))))
        for I,J,m in submultisets(rest):
            for g1 in range(g+1):
                g2=g-g1
                if 2*g1-2+len(I)+1<=0 or 2*g2-2+len(J)+1<=0: continue
                if sum(I)+a!=3*g1-3+len(I)+1: continue
                v1=wk(g1,tuple(sorted(I+(a,))))
                if v1==0: continue
                s+=co*m*v1*wk(g2,tuple(sorted(J+(b,))))
    return s/dfact(2*k+3)

@lru_cache(maxsize=None)
def corr(g,ds,chs):
    """<prod tau_d prod ch_k>_g"""
    n=len(ds)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(ds)+sum(chs)!=3*g-3+n: return Fr(0)
    if not chs: return wk(g,ds)
    if g==0: return Fr(0)
    if chs[-1]>=2*g or any(c%2==0 for c in chs): return Fr(0)
    k=chs[-1]; rest=chs[:-1]; l=(k+1)//2
    s=corr(g,tuple(sorted(ds+(2*l,))),rest)
    seen=set()
    for i,di in enumerate(ds):
        if di in seen: continue
        seen.add(di); c=ds.count(di)
        L=list(ds); L[i]=di+2*l-1
        s-=c*corr(g,tuple(sorted(L)),rest)
    for a in range(2*l-1):
        b=2*l-2-a; sg=(-1)**a
        t=corr(g-1,tuple(sorted(ds+(a,b))),rest)
        for I,J,m in submultisets(ds):
            for C,Cc,mc in submultisets(rest):
                for g1 in range(g+1):
                    g2=g-g1
                    if 2*g1-2+len(I)+1<=0 or 2*g2-2+len(J)+1<=0: continue
                    if sum(I)+a+sum(C)!=3*g1-3+len(I)+1: continue
                    v1=corr(g1,tuple(sorted(I+(a,))),C)
                    if v1==0: continue
                    t+=m*mc*v1*corr(g2,tuple(sorted(J+(b,))),Cc)
        s+=Fr(sg,2)*t
    return s*bern(2*l)/factorial(2*l)

def partitions(n,maxp=None):
    if maxp is None: maxp=n
    if n==0:
        yield (); return
    for p in range(min(n,maxp),0,-1):
        for q in partitions(n-p,p):
            yield (p,)+q

def odd_parts(n):
    for p in partitions(n):
        if all(x%2==1 for x in p): yield tuple(sorted(p))

def a_coef(k,tau):
    # Omega_g = exp(Sum_{k odd} a_k ch_k), a_k = -(k-1)! (1 + tau^-k + (-1-tau)^-k)
    return -factorial(k-1)*(1+tau**(-k)+(-1-tau)**(-k))

def ch_weight(mon,tau):
    w=Fr(1)
    for k in set(mon):
        m=mon.count(k)
        w*=a_coef(k,tau)**m/factorial(m)
    return w

@lru_cache(maxsize=None)
def table(g,dtuple):
    """sorted dtuple of probe psi-powers: dict chmon -> Sum_m 1/m! <tau_d ch prod_leaves (-psi^2/(1-psi))>"""
    r=len(dtuple); N=3*g-3+r-sum(dtuple)
    res={}
    if N<0: return res
    for c in range(0,N+1):
        for mon in odd_parts(c):
            if mon and (g==0 or max(mon)>=2*g): continue
            val=Fr(0)
            for lp in partitions(N-c):  # lp parts = l-1 >=1
                ls=tuple(sorted(x+1 for x in lp))
                co=Fr((-1)**len(ls))
                for v in set(ls):
                    co/=factorial(ls.count(v))
                val+=co*corr(g,tuple(sorted(dtuple+ls)),mon)
            if val!=0: res[mon]=val
    return res

@lru_cache(maxsize=None)
def table_noleaf(g,dtuple):
    r=len(dtuple); N=3*g-3+r-sum(dtuple); res={}
    if N<0: return res
    for mon in odd_parts(N):
        if mon and (g==0 or max(mon)>=2*g): continue
        v=corr(g,dtuple,mon)
        if v!=0: res[mon]=v
    return res

def comps(r,N):
    if r==0:
        if N>=0: yield ()
        return
    for d in range(N+1):
        for q in comps(r-1,N-d):
            yield (d,)+q

def probe_sum(g,U,tau,tab):
    """Sum_{d} prod u_i^{d_i} * tab(g,sorted d) with Omega at tau"""
    r=len(U); N=3*g-3+r
    s=Fr(0)
    for d in comps(r,N):
        co=Fr(1)
        for u,di in zip(U,d): co*=u**di
        t=tab(g,tuple(sorted(d)))
        for mon,v in t.items():
            s+=co*ch_weight(mon,tau)*v
    return s

_X={}
def X(g,U,tau):
    """Xi_g(U) with leaves (Hodge side). Direct, string NOT used."""
    U=tuple(sorted(U)); key=(g,U,tau)
    if key in _X: return _X[key]
    v=probe_sum(g,U,tau,table)
    _X[key]=v
    return v

_M={}
def M(g,U,tau):
    U=tuple(sorted(U)); key=(g,U,tau)
    if key in _M: return _M[key]
    if g==0:
        v=Fr(sum(U))**(len(U)-3)
    else:
        v=probe_sum(g,U,tau,table_noleaf)
    _M[key]=v
    return v
