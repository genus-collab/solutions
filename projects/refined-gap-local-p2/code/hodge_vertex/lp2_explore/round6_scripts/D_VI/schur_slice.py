# Fresh code. Slice data from Schur functions at rational tau; no Hodge input, no MV.
# Z = Sum_nu exp(tau x f2) s_nu(q^rho) s_nu(p); slice p=(p1,0,...), s_nu(p)= p1^n / H_nu.
# With Y = p1/x:  Z = Sum_n z_n(x) Y^n, z_n = Sum_nu (1/H^2) exp((tau+1/2) f2 x) Prod_box 1/S(h x), S(z)=sinh(z/2)/(z/2).
# z_n is even in x after pairing nu, nu'. Series in U=x^2.
from fractions import Fraction as Fr
import math, sys, pickle, time
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(math.comb(m+1,k)*B[k] for k in range(m))/(m+1)
    return B
def partitions(n,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    for k in range(min(n,mx),0,-1):
        for r in partitions(n-k,k): yield (k,)+r
def conj(nu):
    return tuple(sum(1 for a in nu if a>i) for i in range(nu[0])) if nu else ()
def hooks(nu):
    c=conj(nu); return [nu[i]-j+c[j]-i-1 for i in range(len(nu)) for j in range(nu[i])]
def f2(nu): return sum(j-i for i in range(len(nu)) for j in range(nu[i]))
def smul(a,b,L):
    r=[Fr(0)]*L
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b):
            if i+j>=L: break
            r[i+j]+=ai*bj
    return r
def sexp(lg,L):
    # exp of series with lg[0]=0
    e=[Fr(0)]*L; e[0]=Fr(1)
    for n in range(1,L):
        e[n]=sum(k*lg[k]*e[n-k] for k in range(1,n+1))/n
    return e
def slice_series(tau,A,L):
    """returns z[n], zt[n]: U-series of length L for Z and for x^{?}: dZ/dtau coefficient.
    Z_tau = Sum_nu x f2 exp(..) ... ; in Y-expansion zt_n(x) is even in x as well."""
    B=bern(2*L+2)
    # log(1/S(z)) = -log(sinh(z/2)/(z/2)) ; d/dz log sinh(z/2) = (1/2)coth(z/2) = Sum B_{2n} z^{2n-1}/(2n)!
    # log(sinh(z/2)/(z/2)) = Sum_{n>=1} B_{2n} z^{2n}/(2n (2n)!)
    l=[Fr(0)]*L
    for j in range(1,L):
        l[j]=-B[2*j]/(2*j*math.factorial(2*j))
    c=tau+Fr(1,2)
    z=[[Fr(0)]*L for _ in range(A+1)]; zt=[[Fr(0)]*L for _ in range(A+1)]
    z[0][0]=Fr(1)
    for n in range(1,A+1):
        for nu in partitions(n):
            nc=conj(nu)
            if nu<nc: continue
            hs=hooks(nu); f=f2(nu)
            H=1
            for h in hs: H*=h
            lg=[Fr(0)]*L
            for j in range(1,L):
                lg[j]=l[j]*sum(h**(2*j) for h in hs)
            e=sexp(lg,L)
            mult=Fr(1,H*H)
            if nu==nc:   # self-conjugate: f2=0
                for i in range(L): z[n][i]+=mult*e[i]
            else:
                # exp(c f x)+exp(-c f x)=2cosh ; x f(exp(c f x)-exp(-c f x)) = 2 x f sinh(c f x) (even)
                cf=c*f
                ch=[2*cf**(2*j)/math.factorial(2*j) for j in range(L)]
                sh=[Fr(0)]*L
                for j in range(1,L):
                    sh[j]=2*f*cf**(2*j-1)/math.factorial(2*j-1)
                a=smul(e,ch,L); b=smul(e,sh,L)
                for i in range(L):
                    z[n][i]+=mult*a[i]; zt[n][i]+=mult*b[i]
    return z,zt
def logseries(z,zt,A,L):
    # F = log Z in Y, F_tau = Z_tau/Z
    F=[[Fr(0)]*L for _ in range(A+1)]
    Ft=[[Fr(0)]*L for _ in range(A+1)]
    for a in range(1,A+1):
        acc=[a*v for v in z[a]]
        for k in range(1,a):
            p=smul(F[k],z[a-k],L)
            for i in range(L): acc[i]-=k*p[i]
        F[a]=[v/a for v in acc]
    for a in range(1,A+1):
        acc=list(zt[a])
        for k in range(1,a):
            p=smul(Ft[k],z[a-k],L)
            for i in range(L): acc[i]-=p[i]
        Ft[a]=acc
    return F,Ft
def build(tau,A,hmax):
    """Return dictionaries FF[a][h]=[x^{2h-2}y^a] F, FT[a][h]=[x^{2h-2}y^a] F_tau, h=0..hmax, checks lower orders vanish."""
    sig=tau*(1+tau)
    L=A+hmax   # U-power needed: F_a starts at U^{a-1}, need up to U^{a-1+hmax}
    t0=time.time()
    z,zt=slice_series(tau,A,L)
    F,Ft=logseries(z,zt,A,L)
    FF={};FT={}
    for a in range(1,A+1):
        # y^a coefficient = F_a /(sig^a x^{2a}) ; x^{2h-2} <-> U^{h-1+a}
        for i in range(a-1):
            assert F[a][i]==0,("low order F",a,i)
            assert Ft[a][i]==0,("low order Ft",a,i)
        FF[a]={h:F[a][h-1+a]/sig**a for h in range(0,hmax+1) if h-1+a<L}
        FT[a]={h:Ft[a][h-1+a]/sig**a for h in range(0,hmax+1) if h-1+a<L}
    return FF,FT
if __name__=="__main__":
    tau=Fr(sys.argv[1]); A=int(sys.argv[2]); hmax=int(sys.argv[3])
    t0=time.time()
    FF,FT=build(tau,A,hmax)
    pickle.dump((tau,A,hmax,FF,FT),open(f"slice_{str(tau).replace('/','_')}_{A}_{hmax}.pkl","wb"))
    print("done",tau,A,hmax,time.time()-t0)
