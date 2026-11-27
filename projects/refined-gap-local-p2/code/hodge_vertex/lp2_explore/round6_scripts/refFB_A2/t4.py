# Exact T_F limit by polynomial interpolation in delta (tests regularity, A1.3) and comparison with refined BPS.
from gsum import *
from engine import bern
def lagr0(xs,ys):
    tot=Fr(0)
    for i,(x,y) in enumerate(zip(xs,ys)):
        t=y
        for j,xx in enumerate(xs):
            if j!=i: t*=Fr(0)-xx; t/=x-xx
        tot+=t
    return tot
def lagrat(xs,ys,x0):
    tot=Fr(0)
    for i,(x,y) in enumerate(zip(xs,ys)):
        t=y
        for j,xx in enumerate(xs):
            if j!=i: t*=x0-xx; t/=x-xx
        tot+=t
    return tot

def tf_value(eps,g,d):
    E=eps[0]+eps[1]
    b=(Fr(1),Fr(2),Fr(-3))
    # d=0: 3 e3 e1e2 GW is polynomial; d>0: e1e2 GW polynomial of degree 2g (claim Stage B(b)).
    deg=2*g if d>0 else 2*g+3
    xs=[Fr(k,37) for k in range(1,deg+3)]
    ys=[]
    for x in xs:
        al=tuple(E/3+x*bb for bb in b)
        v=loc(al,eps,g,d)*eps[0]*eps[1]
        if d==0: v*=3*al[0]*al[1]*al[2]
        ys.append(v)
    # extra point test of polynomiality
    pred=lagrat(xs[:-1],ys[:-1],xs[-1])
    poly=(pred==ys[-1])
    val=lagr0(xs[:-1],ys[:-1])/(eps[0]*eps[1])
    if d==0: val/= (3*(E/3)**3)
    return val,poly

def series_exp_sinh(eps,N,sgn=1):
    pass

def refined(eps,g,d,spins,conv=1):
    """coefficient of degree 2g-2 in eps of sum_k sum_spins (-1)^{2jL+2jR} N chi_jL(u^k) chi_jR(v^k)/(k*4 sinh(k e1/2) sinh(k e2/2)),
       conv=1: real exponent; conv=-1: e->i e. u=exp((e1-e2)/2), v=exp((e1+e2)/2)."""
    e1,e2=eps
    # scale variable s: e->s e; expand in s to order 2g-2 relative to s^-2
    M=2*g+1
    def mul(a,b):
        r=[Fr(0)]*(M+1)
        for i,x in enumerate(a):
            for j,w in enumerate(b):
                if i+j<=M: r[i+j]+=x*w
        return r
    def inv(a):
        r=[Fr(0)]*(M+1); r[0]=1/a[0]
        for n in range(1,M+1):
            r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
        return r
    def sinhc(x):
        # (2 sinh(x s/2))/(x s) in s
        r=[Fr(0)]*(M+1)
        for n in range(0,M+1,2):
            r[n]=x**n/(2**n*factorial(n+1))
        return r
    def cosh(x):
        r=[Fr(0)]*(M+1)
        for n in range(0,M+1,2):
            r[n]=x**n/factorial(n)
        return r
    tot=Fr(0)
    for k in range(1,d+1):
        if d%k: continue
        dk=d//k
        if dk not in spins: continue
        for (jL2,jR2,N) in spins[dk]:  # twice spins
            sign=(-1)**(jL2+jR2)
            def chi(j2,x):
                # sum over m=-j2,-j2+2,...,j2 of exp(m x s)
                r=[Fr(0)]*(M+1)
                for m in range(-j2,j2+1,2):
                    c=cosh(m*x)
                    r=[a+b for a,b in zip(r,c)]
                return r
            f=mul(chi(jL2,Fr(k)*(e1-e2)/2),chi(jR2,Fr(k)*(e1+e2)/2))
            f=mul(f,inv(sinhc(k*e1))); f=mul(f,inv(sinhc(k*e2)))
            # overall 1/(k * k e1 k e2 s^2)
            c=f[2*g]*sign*N/(k*k*e1*k*e2)
            tot+=c*(conv**(g-1))
    return tot

spins={1:[(0,2,1)],2:[(0,5,1)],3:[(0,6,1),(1,9,1)]}
if __name__=="__main__":
    for eps in [(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))]:
        for g in (1,2,3):
            for d in (1,2,3):
                if g==3 and d==3: continue
                if g==1:
                    pass
                v,p=tf_value(eps,g,d)
                print("eps",*eps,"g",g,"d",d,"T_F value",v,"poly?",p,"refined(+)",refined(eps,g,d,spins,1)==v,"refined(i)",refined(eps,g,d,spins,-1)==v,flush=True)
        for g in (2,3):
            v,p=tf_value(eps,g,0)
            print("eps",*eps,"g",g,"d=0 value",v,"poly?",p,"vs E^2/640 (g=2):",(eps[0]+eps[1])**2/640,flush=True)
