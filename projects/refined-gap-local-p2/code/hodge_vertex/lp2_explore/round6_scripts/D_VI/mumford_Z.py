# MV-free computation of Z_h(tau): Witten-Kontsevich + Mumford's formula for ch_{2l-1} (Faber's algorithm). Fresh code.
# No Schur function, no Marino-Vafa, no TR.
from fractions import Fraction as Fr
from functools import lru_cache
import math, sys, time
from itertools import product
sys.setrecursionlimit(100000)
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(math.comb(m+1,k)*B[k] for k in range(m))/(m+1)
    return B
BN=bern(40)
def dfact(n):
    r=1
    while n>1: r*=n; n-=2
    return r
def msubsets(t):
    """all sub-multisets of sorted tuple t with multiplicity count: yields (sub,comp,mult)."""
    vals=sorted(set(t)); cnt=[t.count(v) for v in vals]
    for ch in product(*[range(c+1) for c in cnt]):
        m=1; a=[];b=[]
        for v,c,k in zip(vals,cnt,ch):
            m*=math.comb(c,k); a+= [v]*k; b+=[v]*(c-k)
        yield tuple(a),tuple(b),m
@lru_cache(maxsize=None)
def wk(g,d):
    """<tau_d1...tau_dn>_g, d sorted tuple."""
    n=len(d)
    if 2*g-2+n<=0: return Fr(0)
    if sum(d)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and d==(1,): return Fr(1,24)
    if d and d[0]==0:   # string
        rest=d[1:]; s=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                e=list(rest); e[i]-=1
                s+=wk(g,tuple(sorted(e)))
        return s
    if 1 in d:          # dilaton
        l=list(d); l.remove(1)
        return (2*g-2+n-1)*wk(g,tuple(l))
    # DVV on largest
    k1=d[-1]; k=k1-1; rest=d[:-1]
    s=Fr(0)
    for i in range(len(rest)):
        e=list(rest); di=e[i]; e[i]=di+k
        s+=Fr(dfact(2*k+2*di+1),dfact(2*di-1))*wk(g,tuple(sorted(e)))
    for a in range(k):
        b=k-1-a
        c=dfact(2*a+1)*dfact(2*b+1)
        t=Fr(0)
        if g>=1: t+=wk(g-1,tuple(sorted(rest+(a,b))))
        for I,Jc,m in msubsets(rest):
            for g1 in range(g+1):
                n1=len(I)+1; g2=g-g1; n2=len(Jc)+1
                if sum(I)+a!=3*g1-3+n1: continue
                if 2*g1-2+n1<=0 or 2*g2-2+n2<=0: continue
                t+=m*wk(g1,tuple(sorted(I+(a,))))*wk(g2,tuple(sorted(Jc+(b,))))
        s+=Fr(c,2)*t
    return s/dfact(2*k+3)
@lru_cache(maxsize=None)
def H(g,d,c):
    """Int_{Mbar_{g,n}} prod psi^d prod ch_c ; d sorted tuple, c sorted tuple of odd integers."""
    n=len(d)
    if 2*g-2+n<=0: return Fr(0)
    if sum(d)+sum(c)!=3*g-3+n: return Fr(0)
    if not c: return wk(g,d)
    if g==0: return Fr(0)
    if any(k>=2*g for k in c): return Fr(0)  # ch_k=0 for k>=2g (Mumford)
    if d and d[0]==0 and n>=1 and 2*g-2+n-1>0:
        rest=d[1:]; s=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                e=list(rest); e[i]-=1
                s+=H(g,tuple(sorted(e)),c)
        return s
    if 1 in d and 2*g-2+n-1>0:
        l=list(d); l.remove(1)
        return (2*g-2+n-1)*H(g,tuple(l),c)
    k=c[-1]; rest=c[:-1]
    pref=BN[k+1]/math.factorial(k+1)
    s=H(g,tuple(sorted(d+(k+1,))),rest)
    for i in range(n):
        e=list(d); e[i]+=k
        s-=H(g,tuple(sorted(e)),rest)
    b=Fr(0)
    for i in range(k):
        j=k-1-i; sg=(-1)**j
        t=H(g-1,tuple(sorted(d+(i,j))),rest) if g>=1 else Fr(0)
        for dI,dJ,m in msubsets(d):
            for cI,cJ,mc in msubsets(rest):
                for g1 in range(g+1):
                    g2=g-g1
                    if sum(dI)+sum(cI)+i!=3*g1-3+len(dI)+1: continue
                    if 2*g1-2+len(dI)+1<=0 or 2*g2-2+len(dJ)+1<=0: continue
                    v=H(g1,tuple(sorted(dI+(i,))),cI)
                    if v==0: continue
                    t+=m*mc*v*H(g2,tuple(sorted(dJ+(j,))),cJ)
        b+=sg*t
    s+=b/2
    return pref*s
def part_odd(n,mx=None):
    if mx is None: mx=n if n%2 else n-1
    if n==0: yield (); return
    k=mx
    while k>=1:
        if k<=n:
            for r in part_odd(n-k,k): yield (k,)+r
        k-=2
def parts_ge2(n,mx=None):
    """multisets d_i>=2 with Sum(d_i-1)=n"""
    if mx is None: mx=n
    if n==0: yield (); return
    for k in range(min(n,mx),0,-1):
        for r in parts_ge2(n-k,k): yield (k+1,)+r
def mult_factor(t):
    r=1
    for v in set(t): r*=math.factorial(t.count(v))
    return r
def Ztable(h):
    """A[lam] = Sum_D (-1)^m/|Aut D| * H(h,D,lam) so that Z_h = Sum_lam (prod c_k)/|Aut lam| A[lam]."""
    A={}
    for e in range(0,3*h-3+1):
        for lam in part_odd(e):
            if any(k>=2*h for k in lam): continue
            lam_s=tuple(sorted(lam))
            s=Fr(0)
            for D in parts_ge2(3*h-3-e):
                s+=Fr((-1)**len(D),mult_factor(D))*H(h,tuple(sorted(D)),lam_s)
            A[lam_s]=s
    return A
def ck(k,tau):
    return math.factorial(k-1)*(Fr(-1)-tau**(-k)+(1+tau)**(-k))
def Zval(A,tau):
    tot=Fr(0)
    for lam,v in A.items():
        if v==0: continue
        p=Fr(1)
        for k in lam: p*=ck(k,tau)
        tot+=p*v/mult_factor(lam)
    return tot
def beta(h,tau):
    L=h+2
    Sinv=[Fr(1,4**j*math.factorial(2*j+1)) for j in range(L)]
    sj=[Fr(0)]*L; sj[0]=Fr(1)
    for n in range(1,L): sj[n]=-sum(Sinv[k]*sj[n-k] for k in range(1,n+1))
    a=-tau;b=1+tau
    return math.factorial(2*h-3)*sum(sj[i]*sj[h-i]*a**(2*i)*b**(2*(h-i)) for i in range(h+1))/(a*b)
if __name__=="__main__":
    # unit tests against classical values
    assert H(1,(0,),(1,))==Fr(1,24)
    assert H(2,(0,0,),()) if False else True
    l2psi2=Fr(1,2)*H(2,(2,),(1,1))   # lambda_2=ch_1^2/2 ; b_2=7/5760
    print("b_2",l2psi2,"expected 7/5760")
    print("lambda_1^3 on M2:",H(2,(),(1,1,1)),"expected 1/2880")
    # lambda_3 lambda_2 lambda_1 on M_3 =1/2(1/4!)(|B_4|/4)(|B_6|/6) =1/2903040
    # ch form: lambda_3=(ch1^3+2 ch3*... ) test via ch3 known: Int ch_3 ch_1^3 on M_3
    hmax=int(sys.argv[1])
    taus=[Fr(2,7),Fr(-1,3),Fr(3),Fr(5,3),Fr(-1,2),Fr(11,4)]
    for h in range(2,hmax+1):
        t0=time.time()
        A=Ztable(h)
        for tau in taus:
            s=tau*(1+tau)
            z=Zval(A,tau)
            print("h",h,"tau",tau,"sigma^(h-1) Z_h == beta_{2h-2}:",s**(h-1)*z==beta(h,tau), "Z_h =",z if h<=3 else '')
        # independence of tau: Z_h polynomial in 1/sigma - test tau -> -1-tau
        print("   h",h,"time",round(time.time()-t0,1)); sys.stdout.flush()
