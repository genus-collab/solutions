# Main-session independent check (evidence only). Exact arithmetic, stdlib only.
from fractions import Fraction as Fr
from math import comb, factorial
from functools import lru_cache
import sys
sys.setrecursionlimit(10000)
tau=Fr(2,7); sigma=tau*(1+tau)
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/(m+1)
    return B[n]
ORD=44
# Laurent series in x: (val, list of coeffs length ORD)
def smul(a,b):
    v=a[0]+b[0]; c=[Fr(0)]*ORD
    for i,x in enumerate(a[1]):
        if x==0: continue
        for j,y in enumerate(b[1]):
            if i+j>=ORD: break
            c[i+j]+=x*y
    return (v,c)
def sadd(a,b):
    if a is None: return b
    if b is None: return a
    v=min(a[0],b[0]); c=[Fr(0)]*ORD
    for i in range(ORD):
        if i-(a[0]-v)>=0 and i-(a[0]-v)<ORD: c[i]+=a[1][i-(a[0]-v)]
        if i-(b[0]-v)>=0 and i-(b[0]-v)<ORD: c[i]+=b[1][i-(b[0]-v)]
    return (v,c)
def sscale(a,k): return (a[0],[k*x for x in a[1]])
def sget(a,n):
    i=n-a[0]
    return a[1][i] if 0<=i<ORD else Fr(0)
def sinv(a):
    i0=next(i for i,x in enumerate(a[1]) if x!=0)
    v=a[0]+i0; c=a[1][i0:]+[Fr(0)]*i0
    r=[Fr(0)]*ORD; r[0]=1/c[0]
    for n in range(1,ORD):
        r[n]=-sum(c[k]*r[n-k] for k in range(1,n+1))/c[0]
    return (-v,r)
@lru_cache(None)
def invsinh(h):  # 1/(2 sinh(h x/2))
    c=[Fr(0)]*ORD
    for k in range(ORD):
        if 2*k+1>=0 and k<ORD:
            pass
    s=[Fr(0)]*ORD
    for k in range(ORD):
        # 2 sinh(hx/2)= sum_{k} 2 (hx/2)^{2k+1}/(2k+1)!
        if 2*k+1<ORD: s[2*k+1]=2*Fr(h,2)**(2*k+1)/factorial(2*k+1)
    return sinv((0,s))
@lru_cache(None)
def expser(k):  # exp(k x)
    return (0,[k**i/Fr(factorial(i)) for i in range(ORD)])
def partitions(n,m=None):
    m=m or n
    if n==0: yield (); return
    for k in range(min(n,m),0,-1):
        for p in partitions(n-k,k): yield (k,)+p
def hooks(nu):
    conj=[sum(1 for x in nu if x>j) for j in range(nu[0])] if nu else []
    return [nu[i]-j+conj[j]-i-1 for i in range(len(nu)) for j in range(nu[i])]
def kappa(nu): return sum(n*(n-2*(i+1)+1) for i,n in enumerate(nu))
def beta_set(nu,L):
    nu=list(nu)+[0]*(L-len(nu)); return [nu[i]+L-1-i for i in range(L)]
@lru_cache(None)
def mn(nu,mu):
    # Murnaghan-Nakayama chi_nu(mu)
    if not mu: return 1 if not nu else 0
    k=mu[0]; L=len(nu)+k
    b=beta_set(nu,L); bs=set(b); tot=0
    for i,x in enumerate(b):
        y=x-k
        if y>=0 and y not in bs:
            h=sum(1 for z in b if y<z<x)
            nb=sorted([z for z in b if z!=x]+[y],reverse=True)
            nn=tuple(v for v in (nb[j]-(L-1-j) for j in range(L)) if v>0)
            tot+=(-1)**h*mn(nn,mu[1:])
    return tot
@lru_cache(None)
def Anu(nu):
    s=expser((tau+Fr(1,2))*kappa(nu)/2)
    for h in hooks(nu): s=smul(s,invsinh(h))
    return s
def zee(mu):
    from collections import Counter
    z=1
    for k,m in Counter(mu).items(): z*=k**m*factorial(m)
    return z
# slice series: list by p1-degree
def Zslice(dm):
    out=[]
    for d in range(dm+1):
        s=None
        mu=(1,)*d
        for nu in partitions(d):
            s=sadd(s,sscale(Anu(nu),Fr(mn(nu,mu),zee(mu))))
        out.append(s)
    return out
def Wslice(R,dm):
    # coefficient of p_R (as monomial) by p1-degree in Z, for R=tuple of parts>=2
    out=[]
    JR=sum(R)
    for a in range(dm+1):
        mu=tuple(sorted(R,reverse=True))+(1,)*a
        s=None
        for nu in partitions(JR+a):
            s=sadd(s,sscale(Anu(nu),Fr(mn(nu,mu),zee(mu))))
        out.append(s)
    return out
def polymul(a,b,dm):
    out=[None]*(dm+1)
    for i,x in enumerate(a):
        if x is None: continue
        for j,y in enumerate(b):
            if i+j>dm or y is None: continue
            out[i+j]=sadd(out[i+j],smul(x,y))
    return out
def polyinv(a,dm):
    # a[0]=1
    out=[None]*(dm+1); out[0]=(0,[Fr(1)]+[Fr(0)]*(ORD-1))
    for n in range(1,dm+1):
        s=None
        for k in range(1,n+1):
            s=sadd(s,smul(a[k],out[n-k]))
        out[n]=sscale(s,Fr(-1))
    return out
def laguerre(n,al,x):
    return sum((-1)**i*comb(n+al,n-i)*x**i/Fr(factorial(i)) for i in range(n+1)) if n>=0 else Fr(0)
def Pi(m):
    p=Fr(1)
    for a in range(1,m): p*=(m*tau+a)
    return p/factorial(m-1)
def Xi_R(h,j):
    # R=(j), r=1
    r=1; N=3*h-3+r; J=j
    Z=Zslice(N); W=Wslice((j,),N)
    # dF/dp_j on slice = W/Z ; F coefficient [p_j p1^a] 
    Q=polymul(W,polyinv(Z,N),N)
    tot=Fr(0)
    for a in range(N+1):
        tot+=(-1)**a*sigma**(-a)*laguerre(N-a,2*h-2+r+a,Fr(-J-a))*sget(Q[a],2*h-2+r+a)
    return tot/(sigma**(h+r-1)*Pi(j))
def logslice(dm):
    Z=Zslice(dm)
    # log via n a_n Z.. : F' = Z'/Z
    F=[None]*(dm+1)
    for n in range(1,dm+1):
        s=sscale(Z[n],Fr(n))
        for k in range(1,n):
            s=sadd(s,sscale(smul(F[k],Z[n-k]),Fr(-k)))
        F[n]=sscale(s,Fr(1,n))
    return F
def Ag(g):
    return abs(bern(2*g-2))*abs(bern(2*g))/(Fr(2*factorial(2*g-2))*(2*g-2)*(2*g))
def Xi_empty(h):
    N=3*h-3; F=logslice(N)
    tot=Fr(comb(N+2*h-2,N))*(-Ag(h)*sigma**(1-h))
    for a in range(1,N+1):
        tot+=(-1)**a*sigma**(-a-(h-1))*laguerre(N-a,2*h-2+a,Fr(-a))*sget(F[a],2*h-2+a)
    return tot
def beta(n,a,b):
    # (n-1)! [t^n] 1/(4 sinh(at/2) sinh(bt/2))
    M=n+4
    S=lambda c:[(bern(2*k)*(2-2**(2*k))*Fr(c)**(2*k)/ (Fr(2)**(2*k)*factorial(2*k)))*1 if False else None for k in range(0)]
    def ss(c):
        # (c t/2)/sinh(c t/2) = sum (1-2^{2k-1}) ... 2 B_{2k} (ct/2)^{2k}/(2k)! *(-1)..
        out=[Fr(0)]*(M)
        for k in range(M//2):
            out[2*k]=-(2**(2*k)-2)*bern(2*k)*(c/2)**(2*k)/factorial(2*k)
        return out
    A=ss(a);B=ss(b);tot=Fr(0)
    for i in range(n+3):
        tot+=A[i]*B[n+2-i]
    return factorial(n-1)*tot/(a*b)
print("sigma=",sigma)
for h,j,exp0 in [(2,2,True),(3,2,True),(3,3,True),(3,4,True),(4,2,True),(2,3,False),(3,5,False)]:
    v=Xi_R(h,j)
    print(f"Xi_{h}(({j})) D={j-1} range={'in' if j-1<=2*h-3 else 'OUT'}: {v}",flush=True)
for h in [2,3,4]:
    a1=-tau;a2=1+tau
    L=sigma**(h-1)*Xi_empty(h); R=beta(2*h-2,a1,a2)
    print(f"VI h={h}: sigma^(h-1) Xi={L}  beta={R}  equal={L==R}",flush=True)
# Prop 2.12
n=24
a=[Fr(factorial(3*d),factorial(d)**3) for d in range(n)]
H=[Fr(0)]*(3*n+1)
for i in range(1,3*n+1): H[i]=H[i-1]+Fr(1,i)
def smulq(x,y):
    r=[Fr(0)]*n
    for i in range(n):
        for j in range(n-i): r[i+j]+=x[i]*y[j]
    return r
def invq(x):
    r=[Fr(0)]*n; r[0]=1/x[0]
    for m in range(1,n): r[m]=-sum(x[k]*r[m-k] for k in range(1,m+1))/x[0]
    return r
sgn=lambda d:(-1)**d
A=[a[d]*sgn(d) for d in range(n)]
G=[-a[d]*H[d]*sgn(d) for d in range(n)]
hh=[3*a[d]*(H[3*d]-H[d])*sgn(d) for d in range(n)]
th=lambda x:[d*x[d] for d in range(n)]
Ainv=invq(A)
# P_I/A = log y + h/A; D=theta(P_I/A)=1+theta(h/A)
hA=smulq(hh,Ainv); D=th(hA); D[0]+=1
E=th(smulq(G,Ainv))
X=[Fr((-27)**d) for d in range(n)]
print("A^2 D = X:",smulq(smulq(A,A),D)==X)
print("2E = D-X:",[2*e for e in E]==[D[i]-X[i] for i in range(n)])
