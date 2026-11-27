# Independent sanity checks for unit U06 (lem:gau-poly, lem:gau-wkb, thm:gau-R0).
# Exact arithmetic with fractions. No project script is used.
from fractions import Fraction as Fr
import cmath, math

# ---------- polynomials in u: dict {k: coeff} ----------
def padd(a,b,cb=Fr(1)):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,Fr(0))+cb*v
        if r[k]==0: del r[k]
    return r
def pmul(a,b):
    r={}
    for k1,v1 in a.items():
        for k2,v2 in b.items():
            r[k1+k2]=r.get(k1+k2,Fr(0))+v1*v2
    return {k:v for k,v in r.items() if v!=0}
def pscal(c,a): return {k:c*v for k,v in a.items() if c*v!=0}
def theta(a): return {k:k*v for k,v in a.items() if k*v!=0}
ONE={0:Fr(1)}; U={1:Fr(1)}; OMU={0:Fr(1),1:Fr(-1)}
def upow(n): return {n:Fr(1)}

def hpf(N, free=None):
    """solve (4n-5theta)ph_n = -r_n for n<=N; free[m] = coefficient of u^{4m} in ph_{5m}. returns ph list and resonance values"""
    free = free or {}
    ph=[ONE]
    res={}
    for n in range(1,N+1):
        # Phi^{[j]}_m for m<=n, j=0..5 with ph_n := 0
        P=[[dict() for m in range(n+1)] for j in range(6)]
        for m in range(n): P[0][m]=ph[m]
        P[0][n]={}
        for j in range(1,6):
            for m in range(n+1):
                t=dict(P[j-1][m])
                if m>=1:
                    prev=P[j-1][m-1]
                    op=padd(pscal(Fr(4*(m-1)+j,5),prev), theta(prev), Fr(-1))
                    t=padd(t,pmul(OMU,op))
                P[j][m]=t
        r={}
        for j in range(1,6):
            f=P[j-1][n]
            r=padd(r, padd(pscal(Fr(4*n+j,5),f), theta(f), Fr(-1)))
        r=padd(r,P[4][n],Fr(-3))
        if n>=1: r=padd(r,pmul(U,P[3][n-1]),Fr(-17,5))
        if n>=2: r=padd(r,pmul(upow(2),P[2][n-2]),Fr(-9,5))
        if n>=3: r=padd(r,pmul(upow(3),P[1][n-3]),Fr(-274,625))
        if n>=4: r=padd(r,pmul(upow(4),ph[n-4]),Fr(-24,625))
        new={}
        for k,v in r.items():
            if 5*k==4*n:
                res[n]=v
            else:
                new[k]=-v/(4*n-5*k)
        if n%5==0:
            res.setdefault(n,Fr(0))
            fm=free.get(n//5,Fr(0))
            if fm!=0: new[4*n//5]=fm
        ph.append(new)
    return ph,res

N=11
ph,res=hpf(N)
print("ph_1 =",ph[1]); print("ph_2 =",ph[2]); print("ph_3 =",ph[3])
print("degrees:",[max(p.keys()) if p else None for p in ph])
print("resonances (free=0):",res)
ph1,res1=hpf(N,{1:Fr(7,3),2:Fr(-5,11)})
print("resonances (free f1=7/3,f2=-5/11):",res1)

# lem:gau-poly (iii): Phi*k(z), k=1+k5 z^5, z^5=u^4 zh^5, is the solution with free coefficient shifted
k5=Fr(7,3)
shift=[padd(ph[n], pscal(k5,pmul(upow(4),ph[n-5])) if n>=5 else {}) for n in range(10)]
ph2,_=hpf(9,{1:Fr(7,3)})
print("Phi°*(1+k5 z^5) == solution with f_1=k5 (n<10):", all(shift[n]==ph2[n] for n in range(10)))

# ---------- T-system of colour 0: F_0 ----------
def Tk(k,F,M):
    # T_k t^n = -t^n + (k-4n)/5 t^{n+1}
    r=[Fr(0)]*(M+1)
    for n,c in enumerate(F):
        if c==0: continue
        r[n]-=c
        if n+1<=M: r[n+1]+=Fr(k-4*n,5)*c
    return r
M=12
F0=[Fr(0)]*(M+1); F0[0]=Fr(1)
def closing(F0):
    F1=Tk(-1,F0,M); F2=[-x for x in Tk(-2,F1,M)]; F3=Tk(2,F2,M); F4=[-x for x in Tk(1,F3,M)]
    G=[-x for x in Tk(0,F4,M)]
    return G,(F1,F2,F3,F4)
for n in range(1,M):
    G,_=closing(F0)
    d=G[n+1]-F0[n+1]   # must vanish; depends linearly on F0[n] with slope 4n
    # determine slope numerically
    F0[n]+=1
    G2,_=closing(F0)
    slope=(G2[n+1]-F0[n+1])-d
    F0[n]-=1
    assert slope==4*n, (n,slope)
    F0[n]-=d/slope
print("F_0 =",F0[:5])
print("ph_n(0) == [t^n]F_0 for n<=%d:"%min(N,M-1), all(ph[n].get(0,Fr(0))==F0[n] for n in range(min(N,M-1)+1)))

# ---------- h(z) = h_{w_0}: numeric check ph_n(1) vs h_n ----------
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(Fr(math.comb(m+1,j))*B[j] for j in range(m))/(m+1)
    return B
B=bern(12)
z5=cmath.exp(2j*cmath.pi/5)
w0=[1-z5**j for j in range(1,5)]+[-5]
K=5
ex=[0]*(2*K)
for k in range(1,K+1):
    bk=B[2*k]/(2*k*(2*k-1))
    s=sum(a**(1-2*k) for a in w0)
    ex[2*k-1]=-float(bk)*s
# exponentiate series
h=[1]+[0]*(2*K-1)
for n in range(1,2*K):
    h[n]=sum(m*ex[m]*h[n-m] for m in range(1,n+1))/n
for n in range(0,8):
    v=float(sum(ph[n].values()))
    v1=float(sum(ph1[n].values()))
    print("n=%d  ph°_n(1)=% .12f  (other free coeffs: % .12f)  h_n=% .12f %+.1ej"%(n,v,v1,h[n].real if isinstance(h[n],complex) else h[n], h[n].imag if isinstance(h[n],complex) else 0))
