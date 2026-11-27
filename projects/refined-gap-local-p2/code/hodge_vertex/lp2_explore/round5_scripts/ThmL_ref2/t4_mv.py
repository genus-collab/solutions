# Form (MV) of s.1: own computation of log Z (Schur sum, Murnaghan-Nakayama characters) vs own Hodge integrals, genus up to Gmax.
# Report compared g <= 2 directly on the Hodge side; here g up to 5.
import sys,time
from fractions import Fraction as Fr
from math import factorial
from functools import lru_cache
import hodge
from hodge import partitions

tau=Fr(sys.argv[1]); N=int(sys.argv[2]); Gmax=int(sys.argv[3])
L=N+2*Gmax+N  # truncation: series in x of Z coefficient of y_mu, y_k=x^k p_k

def mul(a,b):
    r=[Fr(0)]*(L+1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b):
            if i+j>L: break
            if bj!=0: r[i+j]+=ai*bj
    return r

# series: 1/[h] * h x = (hx/2)/sinh(hx/2) ; exp(c x)
def ser_inv_sinhc(h):
    # S(t)=t/sinh(t) with t=hx/2 -> coefficient in x
    # sinh(t)/t = Sum t^{2n}/(2n+1)!
    a=[Fr(0)]*(L+1)
    for n in range(0,L//2+1):
        a[2*n]=Fr(h,2)**(2*n)/factorial(2*n+1)
    # invert
    b=[Fr(0)]*(L+1); b[0]=Fr(1)
    for i in range(1,L+1):
        s=Fr(0)
        for j in range(1,i+1):
            s+=a[j]*b[i-j]
        b[i]=-s
    return b
_inv={}
def inv(h):
    if h not in _inv: _inv[h]=ser_inv_sinhc(h)
    return _inv[h]

def exp_c(c):
    return [Fr(c)**i/factorial(i) for i in range(L+1)]

def conj(nu):
    return [sum(1 for x in nu if x>i) for i in range(nu[0])] if nu else []

def hooks(nu):
    c=conj(nu); out=[]
    for i,r in enumerate(nu):
        for j in range(r):
            out.append(r-j-1+c[j]-i-1+1)
    return out
def content_sum(nu):
    return sum(j-i for i,r in enumerate(nu) for j in range(r))

@lru_cache(maxsize=None)
def char(nu,mu):
    # MN rule; nu tuple partition; mu tuple cycle type
    if not mu: return 1 if not nu else 0
    k=mu[0]; s=0
    # beta numbers
    n=len(nu)
    beta=[nu[i]+n-1-i for i in range(n)]
    bs=set(beta)
    for idx,b in enumerate(beta):
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
def autmu(mu):
    z=1
    for v in set(mu):
        z*=factorial(mu.count(v))
    return z

# Z_n : dict mu-> series in x for coefficient of y_mu
t0=time.time()
Z={0:{():[Fr(1)]+[Fr(0)]*L}}
for n in range(1,N+1):
    Zn={mu:[Fr(0)]*(L+1) for mu in partitions(n)}
    for nu in partitions(n):
        hs=hooks(nu)
        f2=content_sum(nu)
        ser=exp_c((tau+Fr(1,2))*f2)
        pre=Fr(1)
        for h in hs:
            ser=mul(ser,inv(h)); pre/=h
        ser=[pre*a for a in ser]
        for mu in partitions(n):
            ch=char(tuple(nu),tuple(mu))
            if ch==0: continue
            f=Fr(ch,zmu(mu))
            Zn[mu]=[a+f*b for a,b in zip(Zn[mu],ser)]
    Z[n]=Zn
    print("Z level",n,round(time.time()-t0,1),flush=True)

def pmul(A,B):
    r={}
    for a,va in A.items():
        for b,vb in B.items():
            c=tuple(sorted(a+b,reverse=True))
            s=mul(va,vb)
            if c in r: r[c]=[x+y for x,y in zip(r[c],s)]
            else: r[c]=s
    return r
Log={}
for n in range(1,N+1):
    R={k:[n*x for x in v] for k,v in Z[n].items()}
    for m in range(1,n):
        P=pmul(Log[m],Z[n-m])
        for k,v in P.items():
            R[k]=[x-m*y for x,y in zip(R[k],v)]
    Log[n]={k:[x/n for x in v] for k,v in R.items()}
    print("log level",n,round(time.time()-t0,1),flush=True)

sigma=tau*(1+tau)
def Pi(k):
    r=Fr(1)
    for a in range(1,k): r*=(k*tau+a)
    return r/factorial(k-1)
ok=0;bad=[];zero_low=0;nzl=0
for n in range(1,N+1):
    for mu in partitions(n):
        l=len(mu); ser=Log[n][mu]
        # coefficient of p_mu x^m = ser[m+|mu|]; expect m=2g-2+l
        for idx in range(0,L+1):
            m=idx-n
            pred=None
            if (m-l)%2==0 and m>=l-2:
                g=(m-l+2)//2
                if g>Gmax: continue
                pm=Fr(1)
                for k in mu: pm*=Pi(k)
                pred=sigma**(g+l-1)*pm*hodge.M(g,tuple(mu),tau)/autmu(mu)
            else:
                pred=Fr(0)
                if idx>L-2: continue
                if ser[idx]!=0: nzl+=1
                else: zero_low+=1
                continue
            if ser[idx]==pred: ok+=1
            else: bad.append((mu,g,ser[idx],pred))
print("tau",tau,"N",N,"Gmax",Gmax,"matches",ok,"mismatches",len(bad),bad[:4],"shape zeros",zero_low,"shape violations",nzl)
# control: same test with wrong framing sign must fail
b2=0;c2=0
for n in range(2,N+1):
    for mu in partitions(n):
        l=len(mu); g=1
        pm=Fr(1)
        for k in mu: pm*=Pi(k)
        tt=-1-tau # the symmetric framing: M invariant? LLZ has tau<->-1-tau symmetry up to sign; print ratio only
        ctrl=sigma**(g+l-1)*pm*hodge.M(g,tuple(mu),tau+1)/autmu(mu)
        c2+=1
        if Log[n][mu][n+2*g-2+l]!=ctrl: b2+=1
print("control (M at shifted tau, genus 1) mismatches",b2,"/",c2)
