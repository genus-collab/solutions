# Independent of MV: Hodge integrals by Mumford/Faber recursion + DVV.  Tests sigma^(h-1) Z_h = beta_{2h-2}(-tau,1+tau).
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial, comb
import sys, itertools
sys.setrecursionlimit(100000)
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/(m+1)
    return B
BB=bern(40)
def dfact(n):
    r=1
    while n>1: r*=n; n-=2
    return r
@lru_cache(None)
def wk(g,t):
    n=len(t)
    if n==0: return Fr(0)
    if sum(t)!=3*g-3+n or 2*g-2+n<=0: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and t==(1,): return Fr(1,24)
    t=tuple(sorted(t))
    if t[0]==0:
        r=t[1:]; s=Fr(0)
        for i in range(len(r)):
            if r[i]>0: s+=wk(g,tuple(sorted(r[:i]+(r[i]-1,)+r[i+1:])))
        return s
    if t[0]==1:
        r=t[1:]; return (2*g-2+len(r))*wk(g,r)
    # DVV on largest
    k=t[-1]; r=list(t[:-1]); s=Fr(0)
    for j in range(len(r)):
        rr=r[:j]+r[j+1:]
        s+=Fr(dfact(2*k+2*r[j]-1),dfact(2*r[j]-1))*wk(g,tuple(sorted(rr+[k+r[j]-1])))
    for a in range(k-1):
        b=k-2-a
        c=dfact(2*a+1)*dfact(2*b+1)
        s+=Fr(c,2)*wk(g-1,tuple(sorted(r+[a,b]))) if g>=1 else 0
        for mask in range(1<<len(r)):
            A=[r[i] for i in range(len(r)) if mask>>i&1]; B=[r[i] for i in range(len(r)) if not mask>>i&1]
            sa=sum(A)+a; # g1 from dimension
            na=len(A)+1
            g1n=sa-na+3
            if g1n%3: continue
            g1=g1n//3; g2=g-g1
            if g1<0 or g2<0: continue
            v=wk(g1,tuple(sorted(A+[a])))
            if v==0: continue
            s+=Fr(c,2)*v*wk(g2,tuple(sorted(B+[b])))
    return s/dfact(2*k+1)
wk_base=wk
@lru_cache(None)
def hod(g,t,c):
    # <prod tau_t ; prod ch_c>_g, c = sorted tuple of odd ints
    n=len(t)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(t)+sum(c)!=3*g-3+n: return Fr(0)
    if not c: return wk(g,t)
    if g==0: return Fr(0)
    k=c[-1]; rest=c[:-1]; l=(k+1)//2
    pref=BB[2*l]/factorial(2*l)
    tl=list(t); s=hod(g,tuple(sorted(tl+[k+1])),rest)
    for i in range(n):
        s-=hod(g,tuple(sorted(tl[:i]+[tl[i]+k]+tl[i+1:])),rest)
    for a in range(k):
        b=k-1-a; sg=(-1)**a
        ss=hod(g-1,tuple(sorted(tl+[a,b])),rest)
        for mt in range(1<<n):
            A=[tl[i] for i in range(n) if mt>>i&1]; B=[tl[i] for i in range(n) if not mt>>i&1]
            for mc in range(1<<len(rest)):
                C1=tuple(rest[i] for i in range(len(rest)) if mc>>i&1); C2=tuple(rest[i] for i in range(len(rest)) if not mc>>i&1)
                d1=sum(A)+a+sum(C1)-len(A)-1+3
                if d1%3: continue
                g1=d1//3; g2=g-g1
                if g1<0 or g2<0: continue
                v=hod(g1,tuple(sorted(A+[a])),C1)
                if v==0: continue
                ss+=v*hod(g2,tuple(sorted(B+[b])),C2)
        s+=Fr(sg,2)*ss
    return pref*s
def partitions(n,m=None,parts=None):
    if m is None: m=n
    if n==0: yield (); return
    for k in range(min(n,m),0,-1):
        if parts is not None and k not in parts: continue
        for r in partitions(n-k,k,parts): yield (k,)+r
def mult(p):
    r=1
    for v in set(p): r*=factorial(p.count(v))
    return r
def Z(h,tau):
    # Omega = exp(-sum_{k odd}(k-1)! P_k ch_k), P_k=1+tau^-k+(-1-tau)^-k
    def P(k): return 1+tau**(-k)+(-1-tau)**(-k)
    dim=3*h-3; tot=Fr(0)
    odd=set(range(1,dim+1,2))
    for j in range(dim+1):
        for chm in partitions(j,None,odd):
            co=Fr(1)
            for k in chm: co*=-factorial(k-1)*P(k)
            co/=mult(chm)
            ch=tuple(sorted(chm))
            for pp in partitions(dim-j):
                ds=tuple(sorted(q+1 for q in pp))
                tot+=co*Fr((-1)**len(ds),mult(pp))*hod(h,ds,ch)
    return tot
def beta(n,a,b):
    # (n-1)! [t^n] 1/(4 sinh(at/2) sinh(bt/2))
    N=n+4
    def S(c): # (ct/2)/sinh(ct/2) to t^N
        inv=[Fr(0)]*(N+1)
        for i in range(0,N+1,2): inv[i]=(c/2)**i/factorial(i+1)
        o=[Fr(0)]*(N+1); o[0]=Fr(1)
        for i in range(1,N+1):
            o[i]=-sum(inv[j]*o[i-j] for j in range(1,i+1))
        return o
    A=S(a);B=S(b)
    # 1/(4 sinh sinh) = S S/(a b t^2)
    co=sum(A[i]*B[n+2-i] for i in range(n+3))/(a*b)
    return factorial(n-1)*co
if __name__=="__main__":
    print("check: lambda1^3 g=2 via ch: <;ch1,ch1,ch1>_2 =",hod(2,(),(1,1,1)),"(expected 1/2880)")
    print("check: <tau_4 lam1... > ch1 g=1: <tau_0;ch1>_1,1 =",hod(1,(0,),(1,)),"(0 by dim?) <;ch1>_{1,1} dim1:",hod(1,(0,),(1,)))
    hmax=int(sys.argv[1]) if len(sys.argv)>1 else 4
    for tau in [Fr(3,11),Fr(-7,4),Fr(5,2)]:
        sig=tau*(1+tau)
        for h in range(2,hmax+1):
            z=Z(h,tau); lhs=sig**(h-1)*z; rhs=beta(2*h-2,-tau,1+tau)
            print("tau",tau,"h",h,"sigma^(h-1)Z",lhs,"beta",rhs,"OK" if lhs==rhs else "MISMATCH", "ctl(beta w/ wrong sigma power)", lhs==rhs*sig)
            sys.stdout.flush()
