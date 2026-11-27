# (a) Lemma 3.2 as a series identity with Hodge-side M and Q (g=2,3,4); controls.
# (b) (4.1) at all low n including n=0..3 (report's checks use n=2G,2G+1).
# (c) LLZ V_nu eq (3) vs hook product numerically.
# (d) Xi_h(empty) as a polynomial in 1/sigma, closed-form probe.
import math
from fractions import Fraction as Fr
from math import factorial
from hodge import *
import t2_relation as T2
from t2_relation import relation,Opts

# ---- (a)
ORD=9
def smul(a,b):
    r=[Fr(0)]*(ORD+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>ORD: break
            r[i+j]+=x*y
    return r
def sexp(c):
    return [Fr(c)**i/factorial(i) for i in range(ORD+1)]
def spow(a,n):
    r=[Fr(1)]+[Fr(0)]*ORD
    for _ in range(n): r=smul(r,a)
    return r
w=[Fr(1)]*(ORD+1)  # 1/(1-u)
winv=[Fr(1),Fr(-1)]+[Fr(0)]*(ORD-1)
def Qcoef(g,U,m,tau):
    # (1/m!) <Omega lam_U (psi^2/(1-psi))^m> , s^m coefficient
    r=len(U); tot=Fr(0)
    N=3*g-3+r
    for d in comps(r,N+m):
        co=Fr(1)
        for u,di in zip(U,d): co*=u**di
        rem=N+m-sum(d)
        for lp in partitions(rem):
            pass
    return tot
def Qall(g,U,tau):
    """dict m-> s^m coefficient, using correlators directly"""
    r=len(U); N=3*g-3+r; res={}
    for m in range(0,N+1):
        val=Fr(0)
        for d in comps(r,N+m):
            sd=sum(d)
            co=Fr(1)
            for u,di in zip(U,d): co*=u**di
            # leaves with degrees >=2, exactly m of them, ch monomial total c
            for c in range(0,N+m-sd+1):
                ldeg=N+m-sd-c
                if ldeg<2*m: continue
                if m==0 and ldeg!=0: continue
                for lp in partitions(ldeg):
                    if len(lp)!=m or (m>0 and min(lp)<2): continue
                    ls=tuple(sorted(lp))
                    cl=Fr(1)
                    for v in set(ls): cl/=factorial(ls.count(v))
                    for mon in odd_parts(c):
                        if mon and max(mon)>=2*g: continue
                        val+=co*cl*ch_weight(mon,tau)*corr(g,tuple(sorted(d+ls)),mon)
        res[m]=val
    return res
def lemma32(g,U,tau):
    r=len(U); J=sum(U); n=2*g-2+r
    Q=Qall(g,U,tau)
    # consistency Q(-1)=X
    q_m1=sum(v*(-1)**m for m,v in Q.items())
    assert q_m1==X(g,U,tau),("Q(-1)!=X",g,U)
    # RHS series in u: e^{uJ} w^n Q(w-1), w-1 = u/(1-u)
    wm1=[Fr(0)]+[Fr(1)]*ORD
    rhs=[Fr(0)]*(ORD+1)
    for m,v in Q.items():
        t=spow(wm1,m); rhs=[a+v*b for a,b in zip(rhs,t)]
    rhs=smul(smul(sexp(J),spow(w,n)),rhs)
    # LHS: Sum_a M_g(U+1^a) y^a/a!, y=u e^{-u}
    y=smul([Fr(0),Fr(1)]+[Fr(0)]*(ORD-1),sexp(-1))
    lhs=[Fr(0)]*(ORD+1)
    for a in range(0,ORD+1):
        co=hodge_M_strip(g,U,a,tau)/factorial(a)
        t=spow(y,a); lhs=[p+co*q for p,q in zip(lhs,t)]
    return lhs==rhs,sum(1 for x in lhs if x!=0)
def hodge_M_strip(g,U,a,tau):
    # direct M (no string) for small a, string for larger
    if a<=3 or len(U)+a<=4: return M(g,tuple(U)+(1,)*a,tau)
    J=sum(U); base=M(g,tuple(U)+(1,)*3,tau); r=base
    for i in range(3,a): r*=(J+i)
    return r
tau=Fr(3,5)
for g,U in [(1,(2,)),(2,(2,)),(2,(3,)),(2,(2,3)),(3,(2,)),(3,(4,)),(3,(2,5)),(4,(3,))]:
    print("Lemma 3.2 g",g,"U",U,lemma32(g,U,tau),flush=True)
# negative control: wrong w power
# ---- (b)
fails=0;tot=0
for t in [Fr(3,5),Fr(-5,3)]:
    for n in range(0,6):
        for k0,R in [(2,()),(3,()),(4,()),(5,()),(6,()),(2,(2,)),(2,(3,)),(3,(3,)),(4,(5,)),(2,(2,3)),(3,(2,2)),(8,())]:
            if k0==8 and n>3: continue
            v=relation(k0,n,R,t); tot+=1
            if v!=0:
                fails+=1; print("FAIL",t,n,k0,R,v)
print("low n relations",tot,"failures",fails,flush=True)
# ---- (c)
def V3(nu,lam):
    l=len(nu); r=1.0
    for a in range(l):
        for b in range(a+1,l):
            r*=math.sin((nu[a]-nu[b]+b-a)*lam/2)/math.sin((b-a)*lam/2)
    for i in range(1,l+1):
        for v in range(1,nu[i-1]+1):
            r/=2*math.sin((v-i+l)*lam/2)
    return r
def Vh(nu,lam):
    c=[sum(1 for x in nu if x>i) for i in range(nu[0])]
    r=1.0
    for i,ri in enumerate(nu):
        for j in range(ri):
            h=ri-j-1+c[j]-i-1+1
            r/=2*math.sin(h*lam/2)
    return r
mx=0
for n in range(1,10):
    for nu in partitions(n):
        mx=max(mx,abs(V3(nu,0.137)/Vh(nu,0.137)-1))
print("V_nu eq(3) vs hook product, max rel diff",mx,flush=True)
# ---- (d)
for h in [2,3,4]:
    # Xi_h(empty) polynomial in s=1/sigma degree<=h: interpolate at h+1 taus with exact solve, verify at another
    taus=[Fr(2,7),Fr(3,5),Fr(5,3),Fr(7,2),Fr(11,4),Fr(13,9)][:h+2]
    pts=[(1/(t*(1+t)),X(h,(),t)) for t in taus]
    # Newton/Lagrange through first h+1
    def lag(x):
        tot=Fr(0)
        for i in range(h+1):
            xi,yi=pts[i]; p=yi
            for j in range(h+1):
                if j!=i: p*=(x-pts[j][0])/(xi-pts[j][0])
            tot+=p
        return tot
    # coefficients by solving
    A=[[pts[i][0]**j for j in range(h+1)]+[pts[i][1]] for i in range(h+1)]
    for c in range(h+1):
        pv=next(r for r in range(c,h+1) if A[r][c]!=0); A[c],A[pv]=A[pv],A[c]
        A[c]=[x/A[c][c] for x in A[c]]
        for r in range(h+1):
            if r!=c: A[r]=[x-A[r][c]*y for x,y in zip(A[r],A[c])]
    print("h",h,"Xi(empty) coeffs in (1/sigma)^j:",[A[i][h+1] for i in range(h+1)],"predict check",lag(pts[h+1][0])==pts[h+1][1],flush=True)
