# Evidence only. Conifold convention w=(c,-x,-1): port of draft s.1.5, (H0), Theorem VI, M_irr.
import os
import sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "L_comp"))
from fractions import Fraction as Fr
from math import factorial, comb
from functools import lru_cache
from hodge import int_omega_poly, peval

# beta_n(a,b) = (n-1)! [t^n] 1/(4 sinh(at/2) sinh(bt/2))
def S_coefs(N):
    # (y/2)/sinh(y/2) = sum s_i y^(2i); invert sinh(y/2)/(y/2) = sum y^(2i)/(4^i (2i+1)!)
    a=[Fr(1,4**i*factorial(2*i+1)) for i in range(N+1)]
    s=[Fr(0)]*(N+1); s[0]=Fr(1)
    for n in range(1,N+1):
        s[n]=-sum(a[k]*s[n-k] for k in range(1,n+1))
    return s
def beta(n,a,b):
    h=(n+2)//2
    s=S_coefs(h)
    return factorial(n-1)*sum(s[i]*s[h-i]*a**(2*i)*b**(2*(h-i)) for i in range(h+1))/(a*b)

class Th:
    def __init__(self,x):
        self.x=Fr(x); self.c=1+self.x; self.tau=-self.x/self.c; self.sig=self.tau*(1+self.tau)
        self.z=-1/self.sig; self.memo={}; self.memoT={}
    def omega_int(self,h,ds):
        # Int_{Mbar_{h,n}} Omega_w prod psi^ds, w=(c,-x,-1): degree-i part of Omega_w scales as c^-i
        ds=tuple(sorted(ds)); k=(h,ds)
        if k in self.memo: return self.memo[k]
        n=len(ds); deg=3*h-3+n-sum(ds)
        if deg<0: v=Fr(0)
        else: v=peval(int_omega_poly(h,ds),self.z)*self.c**(-deg)
        self.memo[k]=v; return v
    def with_T(self,h,ds):
        # Sum_m (1/m!) Int Omega prod psi^ds prod T(psi), T=-psi^2/(c-psi)=-Sum_{e>=2} psi^e c^(1-e)
        ds=tuple(sorted(ds)); k=(h,ds)
        if k in self.memoT: return self.memoT[k]
        tot=Fr(0); n=len(ds)
        for m in range(0,3*h-3+n+1):
            # parts e_i>=2, dims: sum ds + sum e = 3h-3+n+m - deg Omega, deg Omega>=0
            maxsum=3*h-3+n+m-sum(ds)
            if 2*m>maxsum: break
            def rec(j,minE,rem,cur,mult,last,run):
                pass
            # enumerate nondecreasing e-tuples with multinomial
            def gen(left,mn,rem):
                if left==0:
                    yield (); return
                for e in range(mn,rem//left+1):
                    for rest in gen(left-1,e,rem-e):
                        yield (e,)+rest
            for es in gen(m,2,maxsum):
                mult=factorial(m)
                from collections import Counter
                for v in Counter(es).values(): mult//=factorial(v)
                coef=Fr(mult,factorial(m))
                for e in es: coef*=-self.c**(1-e)
                tot+=coef*self.omega_int(h,ds+es)
        self.memoT[k]=tot; return tot
    def dd(self,h,gs):
        # gs: list of psi-series as coefficient lists
        dim=3*h-3+len(gs)
        tot=Fr(0)
        for ds in itertools.product(*[range(min(len(g),dim+1)) for g in gs]):
            if sum(ds)>dim+200: continue
            cf=Fr(1)
            for g,d in zip(gs,ds):
                cf*=g[d]
                if cf==0: break
            if cf==0: continue
            tot+=cf*self.with_T(h,ds)
        return tot
    def lam(self,j,N):
        return [Fr(j)**a*self.c**(-a-1) if not (j==0 and a>0) else Fr(0) for a in range(N+1)]
    def q0(self,N):
        return [Fr(0)]+[-self.c**(1-a) for a in range(1,N+1)]  # -c psi/(c-psi)

def run():
    fails=0
    xs=[Fr(2,3),Fr(3),Fr(-5,7),Fr(1)]
    for x in xs:
        th=Th(x); c=th.c
        # Theorem VI at the conifold specialisation
        for h in (2,3,4,5):
            if h==5 and x not in (Fr(2,3),): continue
            Z=th.with_T(h,())
            lhs=(-x*c**3)**(h-1)*Z
            rhs=beta(2*h-2,x,Fr(1))
            ok=lhs==rhs
            fails+=not ok
            print("VI x=%s h=%d %s value=%s"%(x,h,ok,lhs))
        # M_irr in the c-convention, non-empty; empty and D=2h-2 controls non-zero
        for h,Rs,ctl in ((2,[(2,)],[(3,),(2,2)]),(3,[(2,),(3,),(4,),(2,2),(2,3),(2,2,2)],[(5,),(3,3),(2,2,2,2)])):
            N=3*h+6
            for R in Rs:
                v=th.dd(h,[th.lam(j,N) for j in R]); fails+=(v!=0)
                print("M_irr x=%s h=%d R=%s zero=%s"%(x,h,R,v==0))
            for R in ctl:
                v=th.dd(h,[th.lam(j,N) for j in R])
                print("control x=%s h=%d R=%s nonzero=%s"%(x,h,R,v!=0))
    # (H0): base point, and full statement on the domain, h=2, t-adic to t^3
    th=Th(Fr(2,3)); c=th.c; h=2; N=9
    q0=th.q0(N)
    g1=[Fr(1),Fr(2),Fr(-1),Fr(3)]; g2=[Fr(0),Fr(1),Fr(5)]
    for gs in ([],[g1],[g1,g2]):
        L=th.dd(h,gs+[q0]); Rr=(2-2*h-len(gs))*th.dd(h,gs)
        fails+=(L!=Rr); print("H0 base s=%d %s"%(len(gs),L==Rr))
    # full: lam=1+m1 t+m2 t^2+m3 t^3, delta=t d1+t^2 d2+t^3 d3, d's have psi-linear terms
    TO=3
    def dF(qdelta,fs):
        # qdelta: dict k->psi series, k>=1. returns t-coefficients of Sum_r (1/r!) dd(fs,delta^r), to t^TO
        out=[Fr(0)]*(TO+1)
        for r in range(0,TO+1):
            for ks in itertools.product(range(1,TO+1),repeat=r):
                if sum(ks)>TO: continue
                out[sum(ks)]+=th.dd(h,fs+[qdelta[k] for k in ks])/factorial(r)
        return out
    def smul(a,s): return [a*v for v in s]
    def sadd(u,v):
        n=max(len(u),len(v)); return [(u[i] if i<len(u) else 0)+(v[i] if i<len(v) else 0) for i in range(n)]
    d={1:[Fr(0),Fr(-1),Fr(2)],2:[Fr(1),Fr(1,2)],3:[Fr(0),Fr(3)]}   # includes psi-linear and constant terms
    mu={1:Fr(2),2:Fr(-1),3:Fr(1,3)}
    # delta' = (lam-1) q0 + lam delta
    dp={}
    for k in range(1,TO+1):
        s=smul(mu[k],q0)
        s=sadd(s,d[k])
        for a in range(1,k):
            s=sadd(s,smul(mu[a],d[k-a]))
        dp[k]=s
    def powlam(e):
        # (1+m1 t+...)^e to t^TO, e integer
        lam=[Fr(1),mu[1],mu[2],mu[3]]
        def mul(u,v): return [sum(u[i]*v[k-i] for i in range(k+1)) for k in range(TO+1)]
        inv=[Fr(1)]+[Fr(0)]*TO
        for k in range(1,TO+1): inv[k]=-sum(lam[i]*inv[k-i] for i in range(1,k+1))
        r=[Fr(1)]+[Fr(0)]*TO
        for _ in range(abs(e)): r=mul(r,lam if e>0 else inv)
        return r
    for fs in ([],[g1]):
        n=len(fs)
        L=dF(dp,fs); Rr=dF(d,fs); pw=powlam(2-2*h-n)
        Rr2=[sum(pw[i]*Rr[k-i] for i in range(k+1)) for k in range(TO+1)]
        fails+=(L!=Rr2); print("H0 full n=%d %s coeffs=%s"%(n,L==Rr2,L))
    # identity for gap: README polynomials
    x=Fr(5,3)
    print("g2",beta(2,x,Fr(1))==(7+10*x**2+7*x**4)/(5760*x))
    print("g3",beta(4,x,Fr(1))==-(1+x**2)*(31+18*x**2+31*x**4)/(161280*x))
    print("g4",beta(6,x,Fr(1))==(381*x**8+620*x**6+686*x**4+620*x**2+381)/(3870720*x))
    print("FAILS",fails)
run()
