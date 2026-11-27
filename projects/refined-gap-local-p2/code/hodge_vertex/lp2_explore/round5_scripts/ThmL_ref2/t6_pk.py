# P_k(D) tau_0 . tau_1 = 0 for generic content-product tau (random B, random p'), k up to Kmax, own code.
# Also: a literal fermionic check of (F2)-(F4) on wedges in a finite mode window.
import sys,random,time
from fractions import Fraction as Fr
from math import factorial,comb
from functools import lru_cache
from hodge import partitions
from itertools import product
random.seed(20261102)
Kmax=int(sys.argv[1]); dout=int(sys.argv[2])
W=Kmax+dout
B=Fr(7,5); pp={k:Fr(random.randint(-9,9),random.randint(1,7)) or Fr(1) for k in range(1,W+1)}
def content_sum(nu): return sum(j-i for i,r in enumerate(nu) for j in range(r))
@lru_cache(maxsize=None)
def char(nu,mu):
    if not mu: return 1 if not nu else 0
    k=mu[0]; s=0; n=len(nu)
    beta=[nu[i]+n-1-i for i in range(n)]; bs=set(beta)
    for b in beta:
        nb=b-k
        if nb<0 or nb in bs: continue
        ht=sum(1 for x in beta if nb<x<b)
        newb=sorted([x for x in beta if x!=b]+[nb],reverse=True)
        nn=tuple(y for y in (newb[i]-(n-1-i) for i in range(n)) if y>0)
        s+=(-1)**ht*char(nn,mu[1:])
    return s
def zmu(mu):
    z=1
    for v in set(mu):
        m=mu.count(v); z*=v**m*factorial(m)
    return z
t0=time.time()
def build(shift,sgn=1):
    """tau = Sum_nu B^{f2+shift|nu|} s_nu(p') s_nu(p): dict mu->coef"""
    T={():Fr(1)}
    for n in range(1,W+1):
        parts=list(partitions(n))
        sv={}
        for nu in parts:
            s=Fr(0)
            for la in parts:
                c=char(nu,la)
                if c:
                    v=Fr(c,zmu(la))
                    for k in la: v*=pp[k]
                    s+=v
            sv[nu]=s*B**(sgn*content_sum(nu)+shift*n)
        for mu in parts:
            T[mu]=sum(Fr(char(nu,mu),zmu(mu))*sv[nu] for nu in parts)
    return T
print("building",flush=True)
t_0=build(0); t_1=build(1); t_m1=build(-1)
print("built",round(time.time()-t0,1),flush=True)
def deriv(mu,i,j):
    """(i d/dp_i)^j p_mu -> (coef, newmu) with D_i = i d/dp_i since p_i=i t_i"""
    m=mu.count(i)
    if j>m: return 0,None
    c=i**j
    for a in range(j): c*=(m-a)
    l=list(mu)
    for _ in range(j): l.remove(i)
    return c,tuple(l)
def apply_d(T,cnt):
    """T dict; cnt dict i->j"""
    R={}
    for mu,v in T.items():
        c=v; cur=mu
        for i,j in cnt.items():
            c0,cur2=deriv(cur,i,j)
            if c0==0: c=0;break
            c*=c0; cur=cur2
        if c!=0: R[cur]=R.get(cur,0)+c
    return R
def bymsum(T):
    L={}
    for mu,v in T.items():
        L.setdefault(sum(mu),{})[mu]=v
    return L
def Pcoef(alpha,k):
    r=Fr(1)
    for i in set(alpha):
        m=alpha.count(i); r*=Fr(-1,i)**m/factorial(m)
    if alpha==(k,): r+=1
    return r
def hirota(k,Ta,Tb,d,coefun=Pcoef):
    """coefficients at output weight d of P_k(D) Ta.Tb"""
    out={}
    for alpha in partitions(k):
        ca=coefun(tuple(alpha),k)
        vals=sorted(set(alpha)); mult=[alpha.count(v) for v in vals]
        for js in product(*[range(m+1) for m in mult]):
            co=ca
            ca_d={};cb_d={}
            for v,m,j in zip(vals,mult,js):
                co*=comb(m,j)*(-1)**(m-j)
                if j: ca_d[v]=j
                if m-j: cb_d[v]=m-j
            wa=sum(v*j for v,j in ca_d.items()); wb=sum(v*j for v,j in cb_d.items())
            # levels: a+b = d+k
            for la in range(0,d+k+1):
                lb=d+k-la
                if la<wa or lb<wb: continue
                A=apply_d({m:Ta[m] for m in partitions(la)},ca_d)
                Bb=apply_d({m:Tb[m] for m in partitions(lb)},cb_d)
                for m1,v1 in A.items():
                    for m2,v2 in Bb.items():
                        key=tuple(sorted(m1+m2,reverse=True))
                        out[key]=out.get(key,0)+co*v1*v2
    return out
tot=0;bad=0;c1=[0,0];c2=[0,0];c3=[0,0]
for k in range(2,Kmax+1):
    for d in range(0,dout+1):
        o=hirota(k,t_0,t_1,d)
        tot+=len(o); bad+=sum(1 for v in o.values() if v!=0)
        o2=hirota(k,t_1,t_0,d); c1[1]+=len(o2); c1[0]+=sum(1 for v in o2.values() if v!=0)
        o3=hirota(k,t_0,t_m1,d); c2[1]+=len(o3); c2[0]+=sum(1 for v in o3.values() if v!=0)
        def wrong(al,kk):
            r=Pcoef(al,kk)
            if al==(kk,): r+=Fr(1,3)
            return r
        o4=hirota(k,t_0,t_1,d,wrong); c3[1]+=len(o4); c3[0]+=sum(1 for v in o4.values() if v!=0)
    print("k",k,"so far: coefficients",tot,"nonzero",bad,"controls swap",c1,"tau_-1",c2,"wrong c_(k)",c3,"t",round(time.time()-t0,1),flush=True)
