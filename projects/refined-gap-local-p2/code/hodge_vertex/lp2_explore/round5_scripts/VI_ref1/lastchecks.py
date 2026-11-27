import os
import sys; sys.argv=['x']; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
import faber, importlib
# 1. (1.1)/(1.2) slice normalisation against MV-free Faber integrals: M_h(1^a) vs Schur, exact
def M(h,a,tau):
    def P(k): return 1+tau**(-k)+(-1-tau)**(-k)
    dim=3*h-3+a; tot=Fr(0); odd=set(range(1,dim+1,2))
    for j in range(dim+1):
        for chm in faber.partitions(j,None,odd):
            co=Fr(1)
            for k in chm: co*=-factorial(k-1)*P(k)
            co/=faber.mult(chm); ch=tuple(sorted(chm))
            d=dim-j
            # distribute d into a legs: 1/(1-psi) => all compositions
            def comps(n,m):
                if m==1: yield (n,); return
                for i in range(n+1):
                    for r in comps(n-i,m-1): yield (i,)+r
            seen={}
            for c in comps(d,a):
                s=tuple(sorted(c)); seen[s]=seen.get(s,0)+1
            for s,mu in seen.items():
                tot+=co*mu*faber.hod(h,s,ch)
    return tot
# Schur side exact, p_1^a, x series
def partitions(n,m=None):
    if m is None:m=n
    if n==0: yield ();return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k): yield (k,)+r
def hooks(nu):
    conj=[sum(1 for t in nu if t>j) for j in range(nu[0])] if nu else []
    return [nu[i]-j+conj[j]-i-1 for i in range(len(nu)) for j in range(nu[i])]
def kappa(nu): return sum(t*(t-2*(i+1)+1) for i,t in enumerate(nu))
def mul(a,b,L):
    o=[Fr(0)]*(L+1)
    for i,x in enumerate(a[:L+1]):
        if x:
            for j,y in enumerate(b[:L+1-i]):
                if y:o[i+j]+=x*y
    return o
def ser(tau,nu,L):
    # x^{|nu|} * e^{(tau+1/2)kappa x/2} prod 1/(2 sinh(hx/2)) / prod h   (coefficient of p_1^n)
    o=[Fr(0)]*(L+1);o[0]=Fr(1)
    c=(tau+Fr(1,2))*Fr(kappa(nu),2)
    e=[c**i/factorial(i) for i in range(L+1)]
    o=mul(o,e,L)
    for h in hooks(nu):
        # (hx)/(2 sinh(hx/2)) inverse of sum (hx/2)^{2k}/(2k+1)!
        ii=[Fr(0)]*(L+1)
        for i in range(0,L+1,2): ii[i]=(Fr(h,2))**i/factorial(i+1)
        s=[Fr(0)]*(L+1);s[0]=Fr(1)
        for i in range(1,L+1):
            s[i]=-sum(ii[j]*s[i-j] for j in range(1,i+1))
        o=mul(o,s,L)
        o=[t/(h*h) for t in o]
    return o
def F(tau,nmax,L):
    Z=[None]*(nmax+1)
    for n in range(nmax+1):
        t=[Fr(0)]*(L+1)
        for nu in partitions(n):
            s=ser(tau,nu,L); t=[a+b for a,b in zip(t,s)]
        Z[n]=t
    f=[None]*(nmax+1)
    for n in range(1,nmax+1):
        s=[n*t for t in Z[n]]
        for k in range(1,n):
            m=mul(f[k],Z[n-k],L)
            s=[a-k*b for a,b in zip(s,m)]
        f[n]=[t/n for t in s]
    return f
for tau in [Fr(3,11),Fr(-7,4)]:
    sig=tau*(1+tau)
    f=F(tau,4,16)
    for h,a in [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2)]:
        sch=f[a][2*h-2+2*a]   # [x^{2h-2+a}] F_a, F_a=x^-a f[a]
        mm=M(h,a,tau)
        ok=sch==sig**(h+a-1)*mm/factorial(a)
        print(f"tau={tau} M_{h}(1^{a}): Schur==sigma^(h+a-1)M/a! {ok}; wrong-sign control {sch==-sig**(h+a-1)*mm/factorial(a)}")
# 2. L[B] vs RHS(E1) with true beta via exact dual numbers; kernel control
class D:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o):
        o=o if isinstance(o,D) else D(o); return D(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o):
        o=o if isinstance(o,D) else D(o); return D(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return D(o)-s
    def __neg__(s): return D(-s.a,-s.b)
    def __mul__(s,o):
        o=o if isinstance(o,D) else D(o); return D(s.a*o.a,s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def inv(s): return D(1/s.a,-s.b/s.a**2)
    def __truediv__(s,o):
        o=o if isinstance(o,D) else D(o); return s*o.inv()
    def __rtruediv__(s,o): return D(o)*s.inv()
    def pw(s,n):
        r=D(1)
        for _ in range(n): r=r*s
        return r
def fstar(tau,T):
    # coefficients of t^j, j=-2..T of 1/(4 sinh(e1 t/2) sinh(e2 t/2)), e1=-tau,e2=1+tau, in dual numbers
    e1=-tau;e2=1+tau
    def S(c):
        ii=[(c*Fr(1,2)).pw(i)*Fr(1,factorial(i+1)) if i%2==0 else D(0) for i in range(T+3)]
        s=[D(0)]*(T+3);s[0]=D(1)
        for i in range(1,T+3):
            acc=D(0)
            for j in range(1,i+1): acc=acc+ii[j]*s[i-j]
            s[i]=-acc
        return s
    A=S(e1);B=S(e2);o={}
    for n in range(0,T+3):
        acc=D(0)
        for i in range(n+1): acc=acc+A[i]*B[n-i]
        o[n-2]=acc/(e1*e2)
    return o
def test(tauv,T,ker=None):
    tau=D(tauv,1); t0=tau.a; e=1+t0; sg=t0*e; sp=1+2*t0
    fs=fstar(tau,T)
    co={j:fs[j] for j in fs if j>=2}   # f_reg=sum co[j] t^j
    if ker is not None:
        n,C=ker  # add kernel element f_* (e t)^{n+2} * C
        et=(1+tau)
        for j in list(fs):
            jj=j+n+2
            if jj>=2 and jj<=T:
                co[jj]=co.get(jj,D(0))+fs[j]*et.pw(n+2)*C
    # B=LL[f/t]: co[j] t^{j-1} -> co[j] (j-1)! N^{-j}
    P=T
    def serE(p):
        o={}
        for j in range(0,P+3-p):
            c=Fr(1)
            for i in range(j): c*=Fr(-p-i,i+1)
            o[p+j]=c*t0**j
        return o
    L={}
    def ad(q,c): L[q]=L.get(q,Fr(0))+c
    for j,cc in co.items():
        b=cc.a*factorial(j-1); bp=cc.b*factorial(j-1)
        for q,c in serE(j).items(): ad(q,-2*bp*c)
        ad(j,2*bp)
        cp=-j*b; q0=j+1
        for q,c in serE(q0).items():
            ad(q,-cp*c); ad(q-1,-2*cp*c/e); ad(q,-t0*cp*c/e)
        ad(q0,-cp); ad(q0-1,2*cp/e); ad(q0,t0*cp/e)
    R={}
    def ar(q,c): R[q]=R.get(q,Fr(0))+c
    for j in range(1,P+4):
        c=Fr((-1)**(j-1),j)*t0**j
        ar(j-2,-c); ar(j-1,-c*t0); ar(j,c*sp/12)
    ar(-1,t0); ar(0,t0*t0/2)
    for q,c in serE(1).items(): ar(q,-(t0*(1+2*sg)/24)*c)
    ar(1,-(t0*(1+2*sg)/24))
    bad=[q for q in range(-1,T-1) if L.get(q,Fr(0))!=R.get(q,Fr(0))/sg**2]
    return bad
for tv in [Fr(3,11),Fr(-7,4),Fr(13,5)]:
    print("tau",tv,"L[beta series]-RHS mismatches up to N^-(T-2), T=36:",test(tv,36))
print("kernel control: f_reg + C f_*(et)^(n+2), n=2 (breaks tau-symmetry) still solves (E1):",test(Fr(3,11),30,ker=(2,Fr(5,3))))
print("control: f_reg with coefficient rescaled must fail:")
import copy
