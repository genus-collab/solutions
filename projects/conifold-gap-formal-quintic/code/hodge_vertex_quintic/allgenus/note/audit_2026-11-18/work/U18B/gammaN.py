# Numeric formal stationary phase for c(z;L) of def:gam-c, to order z^N; fit c_n in Q[L] and evaluate at L=1, compare with h_n.
import cmath, math
from fractions import Fraction as Fr
zeta=cmath.exp(2j*cmath.pi/5)
N=4
TD=2*(N+2*N)+2   # t-degree bound
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(math.comb(m+1,k)*B[k] for k in range(m))/Fr(m+1)
    return B
B=bern(2*N+4)
bk=[None]+[float(B[2*k]/(2*k*(2*k-1))) for k in range(1,N+2)]
def binom(a,r):
    out=1.0
    for i in range(r): out*= (a-i)/(i+1)
    return out
def p(k,L): return sum((zeta**j-L)**(-k) for j in range(5))
# series in (t,z): dict {(dt,dz):coef}, dz may be negative
def mul(A,Bs):
    C={}
    for (a1,b1),x in A.items():
        for (a2,b2),y in Bs.items():
            a=a1+a2; b=b1+b2
            if a>TD: continue
            # weight a/2+b must be <= N ; also prune
            if a/2+b>N+1e-9: continue
            C[(a,b)]=C.get((a,b),0)+x*y
    return C
def expser(A):  # A has only terms of weight >= 1/2 (no constant)
    res={(0,0):1.0}; term={(0,0):1.0}
    for m in range(1,4*N+2):
        term=mul(term,A); term={k:v/m for k,v in term.items()}
        if not term: break
        for k,v in term.items(): res[k]=res.get(k,0)+v
    return res
def dfact(n):  # (n-1)!! for even n
    r=1
    for i in range(n-1,0,-2): r*=i
    return r
def c_series(L):
    C=-L*(1-L**5)/5
    F={}
    for k in range(2,TD):  # t^{k+1}
        F[(k+1,-1)]=(5*(-1)**(k+1)*L**(-k)+p(k,L))/(k*(k+1))
    S1={}
    for k in range(1,TD+1):
        S1[(k,0)]=((-1)**(k+1)*L**(-k)+p(k,L))/(2*k)
    S2={}
    for kk in range(1,N+1):
        for r in range(0,TD+1):
            S2[(r,2*kk-1)]=bk[kk]*binom(1-2*kk,r)*(5**(1-2*kk)*L**(1-2*kk-r)+(-1)**r*p(2*kk-1+r,L))
    # S2 has z^1 constant term (weight 1) ok
    tot={}
    for d in (F,S1,S2):
        for k,v in d.items():
            if k[0]/2+k[1]<=N+1e-9 and k[0]<=TD: tot[k]=tot.get(k,0)+v
    E=expser(tot)
    c=[0]*(N+1)
    for (a,b),v in E.items():
        if a%2: continue
        n=b+a//2
        if 0<=n<=N: c[n]+=v*dfact(a)*C**(a//2)
    return c
# h_n
def hser():
    T=[None]+[5**(1-2*k)+sum((zeta**j-1)**(1-2*k) for j in range(1,5)) for k in range(1,N+2)]
    lg=[0]*(N+1)
    for k in range(1,N+1):
        if 2*k-1<=N: lg[2*k-1]+=bk[k]*T[k]
    # exp of series
    h=[0]*(N+1); h[0]=1
    for n in range(1,N+1):
        h[n]=sum(m*lg[m]*h[n-m] for m in range(1,n+1))/n
    return h
h=hser()
print("h:",[x.real for x in h])
# fit c_n(L)= sum_{k} a_k L^{4n-5k}, k=0..floor(4n/5)
import numpy as np
Ls=[0.5,0.7,1.6,2.1,0.4+0.3j]
vals=[c_series(L) for L in Ls]
for n in range(1,N+1):
    ks=list(range(0,4*n//5+1))
    A=np.array([[L**(4*n-5*k) for k in ks] for L in Ls])
    y=np.array([v[n] for v in vals])
    sol,res,rk,sv=np.linalg.lstsq(A,y,rcond=None)
    print(n,"coeffs",np.round(sol,10),"resid",res,"c_n(1)=",sum(sol),"h_n=",h[n])
