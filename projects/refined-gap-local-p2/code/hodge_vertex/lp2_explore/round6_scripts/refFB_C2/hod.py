# Referee's own Hodge integrals (independent of FB_C/hodge_lib.py). Evidence only.
import math
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations

def df(n):  # double factorial, (-1)!!=1
    r=1
    while n>1: r*=n; n-=2
    return r

B=[Fr(1)]
for m in range(1,41):
    B.append(-sum(Fr(math.comb(m+1,k))*B[k] for k in range(m))/(m+1))

def splits(t):
    t=tuple(t); n=len(t)
    for mask in range(2**n):
        yield tuple(t[i] for i in range(n) if mask>>i&1), tuple(t[i] for i in range(n) if not mask>>i&1)

@lru_cache(maxsize=None)
def tau(g,d):
    """<tau_d1...tau_dn>_g, d sorted tuple."""
    n=len(d)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(d)!=3*g-3+n: return Fr(0)
    if n==3 and g==0: return Fr(1)
    if n==1 and g==1: return Fr(1,24)
    if 0 in d:
        i=d.index(0); rest=d[:i]+d[i+1:]
        s=Fr(0)
        for j in range(len(rest)):
            if rest[j]>0:
                s+=tau(g,tuple(sorted(rest[:j]+(rest[j]-1,)+rest[j+1:])))
        return s
    if 1 in d:
        i=d.index(1); rest=d[:i]+d[i+1:]
        return (2*g-2+len(rest))*tau(g,rest)
    # DVV on last element
    k=d[-1]-1; rest=d[:-1]
    s=Fr(0)
    for j in range(len(rest)):
        s+=Fr(df(2*k+2*rest[j]+1),df(2*rest[j]-1))*tau(g,tuple(sorted(rest[:j]+(rest[j]+k,)+rest[j+1:])))
    for r in range(k):
        q=k-1-r
        c=Fr(df(2*r+1)*df(2*q+1),2)
        s+=c*tau(g-1,tuple(sorted(rest+(r,q))))
        for I,J in splits(rest):
            for g1 in range(g+1):
                a=tau(g1,tuple(sorted(I+(r,))))
                if a: s+=c*a*tau(g-g1,tuple(sorted(J+(q,))))
    return s/df(2*k+3)

@lru_cache(maxsize=None)
def chint(g,a,K):
    """Int prod psi^a prod ch_k, k odd, by Mumford."""
    n=len(a)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(a)+sum(K)!=3*g-3+n: return Fr(0)
    if not K: return tau(g,a)
    if g==0: return Fr(0)
    k=K[-1]; KK=K[:-1]; l=(k+1)//2
    s=chint(g,tuple(sorted(a+(k+1,))),KK)
    for i in range(n):
        s-=chint(g,tuple(sorted(a[:i]+(a[i]+k,)+a[i+1:])),KK)
    for p in range(k):
        q=k-1-p
        sg=(-1)**p
        b=chint(g-1,tuple(sorted(a+(p,q))),KK)
        for I,J in splits(a):
            for K1,K2 in splits(KK):
                for g1 in range(g+1):
                    x=chint(g1,tuple(sorted(I+(p,))),tuple(sorted(K1)))
                    if x: b+=x*chint(g-g1,tuple(sorted(J+(q,))),tuple(sorted(K2)))
        s+=Fr(sg,2)*b
    return s*B[2*l]/math.factorial(2*l)

def oddparts(D,mx=None):
    mx=D if mx is None else mx
    if D==0: yield (); return
    p=min(D,mx)
    if p%2==0: p-=1
    while p>=1:
        for r in oddparts(D-p,p): yield (p,)+r
        p-=2

class Hodge:
    """Int_{Mbar_{g,n}} Omega_v prod psi^a, Omega_v=prod Lambda^vee(v_a)/v_a^g."""
    def __init__(s,wts):
        s.w=list(wts); s.m={}
        s.p={}
    def pw(s,k):
        if k not in s.p: s.p[k]=sum(Fr(1)/x**k for x in s.w)
        return s.p[k]
    def __call__(s,g,a):
        a=tuple(sorted(a)); key=(g,a)
        if key in s.m: return s.m[key]
        D=3*g-3+len(a)-sum(a); t=Fr(0)
        if D>=0 and 2*g-2+len(a)>0:
            for part in oddparts(D):
                c=Fr(1); mult={}
                for q in part:
                    mult[q]=mult.get(q,0)+1
                    c*=-math.factorial(q-1)*s.pw(q)
                for x in mult.values(): c/=math.factorial(x)
                t+=c*chint(g,a,tuple(sorted(part)))
        s.m[key]=t
        return t

def lam_poly_int(g,lam_exps):
    """Int prod lambda_i^{e_i} on Mbar_g from ch's, by extracting via a numeric Hodge? use Newton: lambda from ch."""
    raise NotImplementedError

if __name__=="__main__":
    print("tau",tau(2,(4,)),tau(3,(7,)),tau(2,(2,3)),tau(3,(2,6)),"expect 1/1152 1/82944 29/5760 77/414720")
    # lambda_1 = ch1. lambda_1^3 on M2, lambda_1^6 M3
    print(chint(2,(),(1,1,1)),"exp 1/2880;",chint(3,(),(1,)*6),"exp 1/90720")
    # lambda_g psi^(2g-2): expand Hodge at one weight: coefficient of u^0: Lambda^vee(u)=u^g-..+(-1)^g lambda_g
    # single weight Omega_(u): Int Omega_u psi^(2g-2)= (-1)^g u^-g Int lambda_g psi^(2g-2) on M_{g,1}
    for g,exp in ((1,Fr(1,24)),(2,Fr(7,5760)),(3,Fr(31,967680))):
        h=Hodge([Fr(2)])
        print("lam_g psi^{2g-2}",g,h(g,(2*g-2,))*Fr(2)**g*(-1)**g,exp)
    # lambda_{g-2}lambda_{g-1}lambda_g via three weights: top coefficient
    # FP: Int lambda_{g-2} lambda_{g-1} lambda_g = 1/(2(2g-2)!) |B_{2g-2}|/(2g-2) |B_{2g}|/(2g)
    for g in (2,3,4):
        exp=Fr(1,2*math.factorial(2*g-2))*abs(B[2*g-2])/(2*g-2)*abs(B[2*g])/(2*g)
        val=(-1)**(g-1)*math.factorial(2*g-3)*chint(g,(),(2*g-3,))  # (2g-1)! ch_{2g-1} = (-1)^{g-1} lambda_g lambda_{g-1} ; use lambda_{g-2} lambda_{g-1} lambda_g through ch below
        print("lam3 g",g,"FP",exp)
