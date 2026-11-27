# Own Hodge engine: DVV + Mumford (Faber's algorithm). No MV anywhere.
import sys
from common import *
from functools import lru_cache
from collections import Counter
sys.setrecursionlimit(100000)

def dfact(n):
    r=1
    while n>1:
        r*=n; n-=2
    return r

def submultisets(ms):
    """ms: sorted tuple. yields (sub, complement, multiplicity)"""
    items=sorted(Counter(ms).items())
    def rec(i):
        if i==len(items):
            yield ((),(),1)
            return
        v,m=items[i]
        for (a,b,c) in rec(i+1):
            for t in range(m+1):
                yield ((v,)*t+a,(v,)*(m-t)+b,c*comb(m,t))
    return rec(0)

@lru_cache(maxsize=None)
def wk(g,ds):
    """<tau_ds>_g; ds sorted tuple"""
    n=len(ds)
    if 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds and ds[0]==0:
        rest=ds[1:]
        tot=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                l=list(rest); l[i]-=1
                tot+=wk(g,tuple(sorted(l)))
        return tot
    if ds and 1 in ds:
        l=list(ds); l.remove(1)
        return (2*g-2+n-1)*wk(g,tuple(l))
    # DVV on the last (largest) with k=d-1>=1
    k=ds[-1]-1
    rest=ds[:-1]
    tot=Fr(0)
    for j in range(len(rest)):
        l=list(rest); dj=l[j]; l[j]=dj+k
        tot+=Fr(dfact(2*k+2*dj+1),dfact(2*dj-1))*wk(g,tuple(sorted(l)))
    for a in range(k):
        b=k-1-a
        cf=dfact(2*a+1)*dfact(2*b+1)
        s=Fr(0)
        if g>=1:
            s+=wk(g-1,tuple(sorted(rest+(a,b))))
        for I,J,mult in submultisets(rest):
            sa=sum(I)+a
            # 3 g1 -3 + |I|+1 = sa
            g1n=sa+2-len(I)
            if g1n%3: continue
            g1=g1n//3
            g2=g-g1
            if g2<0 or g1<0: continue
            v1=wk(g1,tuple(sorted(I+(a,))))
            if not v1: continue
            s+=mult*v1*wk(g2,tuple(sorted(J+(b,))))
        tot+=Fr(cf,2)*s
    return tot/dfact(2*k+3)

@lru_cache(maxsize=None)
def chint(g,chs,ds):
    """Int_{Mbar_{g,n}} prod ch_l prod psi^d, chs sorted tuple of odd l, ds sorted tuple."""
    n=len(ds)
    if 2*g-2+n<=0: return Fr(0)
    if sum(chs)+sum(ds)!=3*g-3+n: return Fr(0)
    if not chs: return wk(g,ds)
    if g==0: return Fr(0)
    if chs[-1]>2*g-1: return Fr(0)
    if ds and ds[0]==0 and 2*g-2+n-1>0:
        rest=ds[1:]
        tot=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                l=list(rest); l[i]-=1
                tot+=chint(g,chs,tuple(sorted(l)))
        return tot
    if ds and 1 in ds and 2*g-2+n-1>0:
        l=list(ds); l.remove(1)
        return (2*g-2+n-1)*chint(g,chs,tuple(l))
    # Mumford on the largest class
    l=chs[-1]; rest_ch=chs[:-1]
    kk=(l+1)//2
    pref=BERN[2*kk]/factorial(2*kk)
    tot=chint(g,rest_ch,tuple(sorted(ds+(2*kk,))))
    for i in range(n):
        li=list(ds); li[i]+=l
        tot-=chint(g,rest_ch,tuple(sorted(li)))
    bd=Fr(0)
    for a in range(2*kk-1):
        b=2*kk-2-a
        sg=(-1)**a
        s=Fr(0)
        if g>=1:
            s+=chint(g-1,rest_ch,tuple(sorted(ds+(a,b))))
        for I,J,mi in submultisets(ds):
            for A,B,mc in submultisets(rest_ch):
                dimL=sum(I)+a+sum(A)
                g1n=dimL+2-len(I)
                if g1n%3: continue
                g1=g1n//3; g2=g-g1
                if g2<0: continue
                v1=chint(g1,tuple(sorted(A)),tuple(sorted(I+(a,))))
                if not v1: continue
                s+=mi*mc*v1*chint(g2,tuple(sorted(B)),tuple(sorted(J+(b,))))
        bd+=sg*s
    tot+=bd/2
    return pref*tot

def odd_partitions(n,maxp):
    """multisets of odd parts summing to n with parts<=maxp, returned sorted ascending"""
    if n==0:
        yield ()
        return
    p=maxp if maxp%2==1 else maxp-1
    while p>=1:
        if p<=n:
            for r in odd_partitions(n-p,p):
                yield r+(p,)
        p-=2

# polynomials in z as lists
def padd(a,b):
    n=max(len(a),len(b)); r=[Fr(0)]*n
    for i,x in enumerate(a): r[i]+=x
    for i,x in enumerate(b): r[i]+=x
    return r
def pmul(a,b):
    if not a or not b: return []
    r=[Fr(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                r[i+j]+=x*y
    return r
def pscal(c,a): return [c*x for x in a]
def peval(a,z):
    r=Fr(0)
    for c in reversed(a): r=r*z+c
    return r
def ptrim(a):
    a=list(a)
    while a and not a[-1]: a.pop()
    return a

# c_l(z): power sums of (-1, a, b) with a+b=z, ab=z where z=-1/sigma, a=-1/tau, b=1/(1+tau).
_pl=[[Fr(2)],[Fr(0),Fr(1)]]
def plab(l):
    while len(_pl)<=l:
        x=padd(pmul([Fr(0),Fr(1)],_pl[-1]),pmul([Fr(0),Fr(-1)],_pl[-2]))
        _pl.append(x)
    return _pl[l]
def cl(l):
    """power sum t1^l+t2^l+t3^l for odd l, polynomial in z"""
    return padd([Fr(-1)],plab(l))

@lru_cache(maxsize=None)
def omega_coef(chs):
    """coefficient poly in z of ch-monomial chs in Omega=exp(Sum_{l odd}(l-1)! c_l ch_l)"""
    r=[Fr(1)]
    for l,m in Counter(chs).items():
        base=pscal(Fr(factorial(l-1)),cl(l))
        pw=[Fr(1)]
        for i in range(m): pw=pmul(pw,base)
        r=pmul(r,pscal(Fr(1,factorial(m)),pw))
    return tuple(r)

def z_of_tau(tau):
    return -1/(tau*(1+tau))

@lru_cache(maxsize=None)
def int_omega_poly_c(g,ds):
    n=len(ds)
    deg=3*g-3+n-sum(ds)
    if deg<0: return ()
    if g==0: return (wk(0,ds),) if deg==0 else ()
    tot=[]
    for chs in odd_partitions(deg,2*g-1):
        v=chint(g,chs,ds)
        if v:
            tot=padd(tot,pscal(v,list(omega_coef(chs))))
    return tuple(ptrim(tot))

def int_omega_poly(g,ds):
    return list(int_omega_poly_c(g,tuple(sorted(ds))))

def int_omega(g,ds,z):
    return peval(int_omega_poly(g,ds),z)

def partitions_ge(n,maxp):
    if n==0:
        yield ()
        return
    for k in range(min(n,maxp),0,-1):
        for r in partitions_ge(n-k,k):
            yield (k,)+r

def partitions_ge_len(n,maxlen,maxp=None):
    if maxp is None: maxp=n
    if n==0:
        yield ()
        return
    if maxlen==0: return
    for k in range(min(n,maxp),0,-1):
        for r in partitions_ge_len(n-k,maxlen-1,k):
            yield (k,)+r

_cacheE={}
def tp_sum(g,e):
    """for fixed leaf degrees e (sorted): dict m -> poly in z of
    Sum_{d_i>=2 multisets of size m} (1/prod mult!) <Omega tau_e tau_d>.  (This is 1/m! Int Omega psi^e Tp^m.)"""
    key=(g,e)
    if key in _cacheE: return _cacheE[key]
    r=len(e)
    M=3*g-3+r-sum(e)
    res={}
    if M<0:
        _cacheE[key]=res
        return res
    for M1 in range(0,M+1):
        # Sum (d_k-1) = M1, ch degree = M-M1
        for par in partitions_ge(M1,M1 if M1>0 else 1):
            m=len(par)
            ds=tuple(sorted(e+tuple(p+1 for p in par)))
            w=Fr(1)
            for v,c in Counter(par).items():
                w/=factorial(c)
            poly=int_omega_poly(g,ds)
            if poly:
                res[m]=padd(res.get(m,[]),pscal(w,poly))
    _cacheE[key]=res
    return res

def leaf_groups(U,N):
    """Yield (sorted e tuple, weight) with Sum weight = Sum over ordered e, sum e<=N, prod j^e. Equal parts grouped. U may hold Fractions."""
    items=sorted(Counter(U).items(),key=lambda t:Fr(t[0]))
    def rec(i,left):
        if i==len(items):
            yield ((),Fr(1))
            return
        j,m=items[i]
        for tot in range(left+1):
            for par in partitions_ge_len(tot,m):
                pad=par+(0,)*(m-len(par))
                c=factorial(m)
                for v,cnt in Counter(pad).items():
                    c//=factorial(cnt)
                if tot>0 and j==0:
                    continue
                w=Fr(c)*(Fr(j)**tot if tot>0 else Fr(1))
                for e2,w2 in rec(i+1,left-tot):
                    yield (pad+e2,w*w2)
    grouped={}
    for e,w in rec(0,N):
        k=tuple(sorted(e))
        grouped[k]=grouped.get(k,Fr(0))+w
    return grouped

_Qcache={}
def Q_poly(g,U):
    """dict m -> poly in z: coefficient of s^m in Q_g(s;U). U tuple of numbers."""
    U=tuple(sorted(U,key=Fr)); r=len(U)
    key=(g,U)
    if key in _Qcache: return _Qcache[key]
    assert 2*g-2+r>0
    N=3*g-3+r
    res={}
    for e,w in leaf_groups(U,N).items():
        ts=tp_sum(g,e)
        for m,poly in ts.items():
            res[m]=padd(res.get(m,[]),pscal(w,poly))
    _Qcache[key]=res
    return res

def Xi_poly(g,U):
    """Xi_g(U)=Q_g(-1;U) as polynomial in z=-1/sigma."""
    q=Q_poly(g,U)
    tot=[]
    for m,poly in q.items():
        tot=padd(tot,pscal(Fr((-1)**m),poly))
    return ptrim(tot)

_Mcache={}
def M_poly(g,mu):
    """Int Omega_g prod lam_{mu_i}(psi_i) as poly in z, with M_0 conv."""
    mu=tuple(sorted(mu,key=Fr)); l=len(mu)
    key=(g,mu)
    if key in _Mcache: return _Mcache[key]
    if g==0:
        if l>=1:
            res=[Fr(sum(mu))**(l-3)]
        else:
            raise ValueError
    else:
        N=3*g-3+l
        res=[]
        for e,w in leaf_groups(mu,N).items():
            p=int_omega_poly(g,e)
            if p:
                res=padd(res,pscal(w,p))
    _Mcache[key]=res
    return res
