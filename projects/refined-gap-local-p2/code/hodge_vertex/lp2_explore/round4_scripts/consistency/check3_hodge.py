# Genus-2 exact Hodge check of s.1.5, s.7.1, s.7.4, s.7.6: normalised Omega_w, leaves, dilaton point q_0.
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
# base values, genus 2, lambda monomials (l1 power reduced: l1^2=2 l2) labelled 'd' = (e1,e2) exponents of l1,l2
BASE={ ((0,0),(4,)):F(1,1152), ((0,0),(2,3)):F(29,5760), ((0,0),(2,2,2)):F(7,240),
       ((1,0),(3,)):F(1,480), ((1,0),(2,2)):F(5,576),
       ((0,1),(2,)):F(7,5760),
       ((1,1),()):F(1,5760)}
@lru_cache(None)
def I(mon,a):
    a=tuple(sorted(a))
    l1,l2=mon
    deg=l1+2*l2
    n=len(a)
    if sum(a)+deg!=3+n: return F(0)
    if min(a,default=9)==0:
        rest=list(a);rest.remove(0)
        t=F(0)
        for i in range(len(rest)):
            if rest[i]>=1:
                r=rest[:];r[i]-=1
                t+=I(mon,tuple(r))
        return t
    if 1 in a:
        rest=list(a);rest.remove(1)
        return (2+len(rest))*I(mon,tuple(rest))
    return BASE[(mon,a)]
assert I((1,1),(1,))==F(1,2880)
def omega(w):
    """returns dict mon->coef for prod_a (1 - l1 u_a + l2 u_a^2), reduced by l1^2=2 l2, l2^2=0"""
    r={(0,0):F(1)}
    for wa in w:
        u=1/wa
        new={}
        for (m,c) in r.items():
            for (d,cc) in (((0,0),F(1)),((1,0),-u),((0,1),u*u)):
                k=(m[0]+d[0],m[1]+d[1])
                new[k]=new.get(k,F(0))+c*cc
        r=new
    # reduce
    out={}
    for (m,c) in r.items():
        l1,l2=m
        coef=c
        while l1>=2:
            l1-=2;l2+=1;coef*=2
        if l2>=2: continue
        if l1+2*l2>3: continue
        out[(l1,l2)]=out.get((l1,l2),F(0))+coef
    return out
def Xi(R,e1,e2,scale=None,w=None,mmax=None):
    """Xi_2(R): Sum_m 1/m! Int Omega prod lam_j(psi_i) prod T(psi_k)"""
    E=e1+e2
    c=E if scale is None else scale
    if w is None: w=(E,-e1,-e2)
    Om=omega(w)
    r=len(R); mmax=3+r
    tot=F(0)
    for m in range(0,mmax+1):
        n=r+m
        dim=3+n
        # T = -Sum_{a>=2} psi^a c^{1-a}
        for dils in product(range(2,dim+1),repeat=m):
            sd=sum(dils)
            if sd>dim: continue
            cd=F(1)
            for a in dils: cd*= -c**(1-a)
            for leg in product(range(0,dim+1),repeat=r):
                sl=sum(leg)
                if sl+sd>dim: continue
                cl=F(1)
                for j,a in zip(R,leg):
                    cl*= F(j)**a * c**(-1-a)
                for mon,co in Om.items():
                    if mon[0]+2*mon[1]+sl+sd!=dim: continue
                    tot+=co*cd*cl*I(mon,tuple(leg)+tuple(dils))/factorial(m)
    return tot
def beta2(e1,e2):
    return (7*e1**4+10*e1**2*e2**2+7*e2**4)/(5760*e1*e2)
for (e1,e2) in [(F(7),F(3)),(F(1,3),F(1)),(F(1,8),F(1)),(F(4),F(9))]:
    E=e1+e2
    Z=Xi((),e1,e2)
    print("eps",e1,e2,"(-e1e2E^3)Z[2,0] =",-e1*e2*E**3*Z," beta_2 =",beta2(e1,e2)," (-3)beta_2=",-3*beta2(e1,e2))
    print("   M_irr(2) Xi({2}) =",Xi((2,),e1,e2),"  Thm C Xi({1}) =",Xi((1,),e1,e2), " Xi({1,1}) =",Xi((1,1),e1,e2))
    print("   string Xi(1,2) vs (2/E^2) Xi(2):",Xi((1,2),e1,e2),"; Xi(1,3)-(3/E^2)Xi(3) =",Xi((1,3),e1,e2)-F(3)/E**2*Xi((3,),e1,e2),"; Xi(3)=",Xi((3,),e1,e2))
    print("   L0: Xi({0}) - ((2-2h)/E^2) Z =",Xi((0,),e1,e2)-F(-2)/E**2*Z)
    # L0 with R={3}: D_R=2, h=2 : (D_R-2h+2)=0
    print("   L0: Xi({0,3}) - ((2-4+2)/E^2)Xi(3) =",Xi((0,3),e1,e2))
    # scale law (7.2): Xi(R)(e1,e2)=E^{-(3h-3+2r)} Xi^{c=1}(tau)
    tau=-e1/E
    for R in ((),(3,),(2,3)):
        c1v=Xi(R,F(1),F(1),scale=F(1),w=(F(1),tau,-1-tau))
        print("   scale law R=%s: diff"%(R,),Xi(R,e1,e2)-E**(-(3+2*len(R)))*c1v)
    # controls
    print("   control mode A (w->-w): (-e1e2E^3)Z =",-e1*e2*E**3*Xi((),e1,e2,w=(-E,e1,e2)), "; Xi({2}) =",Xi((2,),e1,e2,w=(-E,e1,e2)))
    print("   control scale wrong (leaf scale E->e1): Xi({2})=",Xi((2,),e1,e2,scale=e1))
    # s-decomposition / R2 cor 7 h=2 sign : sigma^2 Z^{c=1} = -(1)! [l^4] S(e1 l)S(e2 l), with S = ?
