# Extras: hook vs p_k form of s_nu(q^rho); row-hook lemma multiset; LLZ (1)-(3) literal over Q(i) vs Hodge engine;
# s.8.2 at h=1 and R empty; s.7.3 scaled string/part-0 relations; mode A control; symmetric/tau=0 consistency.
from schur_side import *
from hodge import *
from mk import X0
import sys
counts=Counter(); Fails=[]
def chk(name,cond,info=""):
    counts[name]+=1
    if not cond:
        Fails.append((name,info)); print("FAIL",name,info,flush=True)

# --- series helpers (real Laurent treated as list with offset)
def sinh_ser(h,P):
    """2 sinh(h x/2)= h x * S(hx). returns S series in x len P where S(y)=sinh(y/2)/(y/2)"""
    r=[Fr(0)]*P
    for j in range(0,P//2+1):
        if 2*j<P:
            r[2*j]=Fr(h)**(2*j)/(Fr(4)**j*factorial(2*j+1))
    return r
def smul(a,b,P):
    r=[Fr(0)]*P
    for i,x in enumerate(a):
        if x:
            for j in range(P-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def sinv(a,P):
    r=[Fr(0)]*P
    r[0]=1/a[0]
    for n in range(1,P):
        s=Fr(0)
        for i in range(1,n+1):
            s+=a[i]*r[n-i]
        r[n]=-s/a[0]
    return r
def sexp_lin(c,P):
    return [Fr(c)**i/factorial(i) for i in range(P)]

# hook form vs character form of s_nu(q^rho)
P=12
for N in range(1,8):
    for nu in partitions(N):
        # character: Sum chi_nu(mu)/z_mu prod 1/(2 sinh(mu_i x/2)), starts at x^{-l(mu)}; put all at offset x^{-N}
        acc=[Fr(0)]*P
        for mu in partitions(N):
            ch=chi(nu,mu)
            if not ch: continue
            s=[Fr(0)]*P; s[0]=Fr(1)
            for m in mu:
                s=smul(s,sinv(sinh_ser(m,P),P),P)
            pre=Fr(ch,zee(mu))
            for m in mu: pre/=m
            shift=N-len(mu)   # x^{-l}=x^{-N} x^{N-l}
            for i in range(P-shift):
                acc[i+shift]+=pre*s[i]
        # hook form e^{f2 x/2} prod 1/(2 sinh(h x/2))
        s=sexp_lin(Fr(contents_sum(nu),2),P)
        pre=Fr(1)
        for h in hooks(nu):
            s=smul(s,sinv(sinh_ser(h,P),P),P); pre/=h
        s=[pre*x for x in s]
        chk("s_nu(q^rho) hook = character form",s==acc,nu)
# row-hook multiset lemma (Prop 2.1(a))
for N in range(1,13):
    for nu in partitions(N):
        l=len(nu); bs=[nu[i]+l-1-i for i in range(l)]
        c=conj(nu)
        ok=True
        for i in range(l):
            rowh=[nu[i]-j-1+c[j]-i-1+1 for j in range(nu[i])]
            other=[bs[i]-bs[k] for k in range(i+1,l)]
            if sorted(rowh+other)!=list(range(1,bs[i]+1)): ok=False
        chk("row-hook lemma",ok,nu)

# --- literal LLZ over Q(i): series in lambda, complex coefficients as pairs
class C:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): return C(s.a+o.a,s.b+o.b)
    def __sub__(s,o): return C(s.a-o.a,s.b-o.b)
    def __mul__(s,o): return C(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def scal(s,c): return C(s.a*c,s.b*c)
    def iszero(s): return s.a==0 and s.b==0
    def __eq__(s,o): return s.a==o.a and s.b==o.b
    def __repr__(s): return f"({s.a}+{s.b}i)"
def ipow(n):
    n%=4
    return [C(1,0),C(0,1),C(-1,0),C(0,-1)][n]

def sin2_ser(n,P):
    """2 sin(n lam/2)= n lam * Sn(lam), Sn(lam)=sin(n lam/2)/(n lam/2); returns Sn series. n may be negative: even function; pre n"""
    r=[Fr(0)]*P
    for j in range(0,P//2+1):
        if 2*j<P:
            r[2*j]=(-1)**j*Fr(n)**(2*j)/(Fr(4)**j*factorial(2*j+1))
    return r
def V_literal(nu,P):
    """(LLZ-3) literally: returns (power of lambda offset, real series)"""
    l=len(nu)
    s=[Fr(0)]*P; s[0]=Fr(1); pre=Fr(1); off=0
    for a in range(l):
        for b in range(a+1,l):
            top=nu[a]-nu[b]+(b-a)      # 1-index difference same
            bot=b-a
            # sin(top lam/2)/sin(bot lam/2)
            s=smul(s,sin2_ser(top,P),P); pre*=top
            s=smul(s,sinv(sin2_ser(bot,P),P),P); pre/=bot
    for i0,ni in enumerate(nu):
        i=i0+1
        for v in range(1,ni+1):
            n=v-i+l
            s=smul(s,sinv(sin2_ser(n,P),P),P); pre/=n; off-=1
    return off,[pre*x for x in s]
def Z_LLZ(tau,W,P):
    """Z[mu]=coefficient of p_mu: list of complex coefficient with index i meaning lambda^{-|mu|+i}."""
    Z={}
    for N in range(W+1):
        for mu in partitions(N):
            Z[mu]=[C() for _ in range(P)]
    for N in range(W+1):
        for nu in partitions(N):
            off,V=V_literal(nu,P)
            assert off==-N
            kappa=2*contents_sum(nu)
            # exp(i (tau+1/2) kappa lam/2)
            c=(tau+Fr(1,2))*kappa/2
            ex=[]
            for j in range(P):
                ex.append(ipow(j).scal(c**j/factorial(j)))
            ser=[C() for _ in range(P)]
            for i in range(P):
                for j in range(P-i):
                    if V[j]:
                        ser[i+j]=ser[i+j]+ex[i].scal(V[j])
            for mu in partitions(N):
                ch=chi(nu,mu)
                if ch:
                    cc=Fr(ch,zee(mu))
                    for i in range(P):
                        Z[mu][i]=Z[mu][i]+ser[i].scal(cc)
    return Z
def mono_mul(m1,m2): return tuple(sorted(m1+m2,reverse=True))
def log_Z(Z,W,P):
    F={}
    all_m=[mu for N in range(1,W+1) for mu in partitions(N)]
    all_m.sort(key=len)
    # use weight grading by |mu|: |m| F_m = |m| Z_m - Sum_{m1 m2=m} |m1| F_m1 Z_m2
    all_m.sort(key=sum)
    from itertools import product as prod
    for m in all_m:
        w=sum(m)
        res=[x.scal(w) for x in Z[m]]
        cnt=Counter(m); vals=sorted(cnt)
        for ev in prod(*[range(cnt[v]+1) for v in vals]):
            m1=();
            for v,e in zip(vals,ev): m1+=(v,)*e
            if not m1 or sum(m1)==w: continue
            m1=tuple(sorted(m1,reverse=True))
            m2=[]
            for v,e in zip(vals,ev): m2+=[v]*(cnt[v]-e)
            m2=tuple(sorted(m2,reverse=True))
            w1=sum(m1)
            A=F[m1];B=Z[m2]
            for i in range(P):
                if A[i].iszero(): continue
                for j in range(P-i):
                    if B[j].iszero(): continue
                    res[i+j]=res[i+j]-(A[i]*B[j]).scal(w1)
        F[m]=[x.scal(Fr(1,w)) for x in res]
    return F
for tau in [Fr(2,7),Fr(-5,3)]:
    sig=tau*(1+tau); z=z_of_tau(tau)
    W=6;P=16
    Z=Z_LLZ(tau,W,P)
    F=log_Z(Z,W,P)
    for m,ser in F.items():
        l=len(m); N=sum(m)
        PiM=Fr(1)
        for j in m: PiM*=Pi(j,tau)
        for i,val in enumerate(ser):
            pw=-N+i
            gg=pw+2-l
            if gg<0 or gg%2:
                chk("LLZ literal: other powers vanish",val.iszero(),(m,pw,val)); continue
            g=gg//2
            if g>3: continue
            Ig=Fr(sum(m))**(l-3) if g==0 else (-sig)**g*peval(M_poly(g,m),z)
            exp=ipow(N+l).scal(-sig**(l-1)*PiM*Ig/aut(m))
            chk("LLZ (LLZ-def),(LLZ-1),(LLZ-3) literal vs Hodge engine",val==exp,(tau,m,g,val,exp))

# s.8.2 at h=1, and R empty
K=Fld(None)
for tau in [Fr(2,7),Fr(-5,3),Fr(13,9)]:
    z=z_of_tau(tau); sig=tau*(1+tau)
    for h,R in [(1,(2,)),(1,(3,)),(1,(2,2)),(1,(2,3)),(1,(2,2,2)),(2,(2,)),(2,(3,)),(3,(4,)),(2,(2,3))]:
        r=len(R);N=3*h-3+r;J=sum(R);Pp=2*h-2+r+2*N+J+1
        F,items=compute_F(R,N,Pp,tau,K)
        chk("s.8.2 vs definition (Hodge)",Xi_from_F(F,items,R,h,tau,K)==peval(Xi_poly(h,R),z),(tau,h,R))
    for h in (2,3,4):
        N=3*h-3;Pp=2*h-2+2*N+1
        F,items=compute_F((),N,Pp,tau,K)
        tot=laguerre(N,2*h-2,0)*peval(int_omega_poly(h,()),z)
        for a in range(1,N+1):
            tot+=(-1)**a*(1/sig)**a*laguerre(N-a,2*h-2+a,-a)*coeff(F,items,(),a,2*h-2+a)
        tot/=sig**(h-1)
        chk("s.8.2 empty profile",tot==peval(Xi_poly(h,()),z),(tau,h))
# s.7.3 scaled relations + mode A
def XiW(h,U,E,w1,w2,w3):
    """Omega with weights (w1,w2,w3) normalised; leaves at scale E"""
    def cc(l): return (Fr(-1)/w1)**l+(Fr(-1)/w2)**l+(Fr(-1)/w3)**l
    def om(chs):
        r=Fr(1)
        for l,m in Counter(chs).items():
            r*=(factorial(l-1)*cc(l))**m/factorial(m)
        return r
    r=len(U); N=3*h-3+r
    tot=Fr(0)
    for e,wt in leaf_groups(U,N).items():
        Mleft=N-sum(e)
        for M1 in range(Mleft+1):
            for par in partitions_ge(M1,M1 if M1>0 else 1):
                ds=tuple(sorted(e+tuple(p+1 for p in par)))
                w=Fr(1)
                for v,c in Counter(par).items(): w/=factorial(c)
                w*=Fr((-1)**len(par))
                ex=-sum(x+1 for x in e)-sum(par)+len(par)-len(par)
                ex=-sum(x+1 for x in e)-sum(p for p in par)
                deg=Mleft-M1
                s=Fr(0)
                for chs in odd_partitions(deg,2*h-1):
                    v=chint(h,chs,ds)
                    if v: s+=om(chs)*v
                tot+=wt*w*s*Fr(E)**ex
    return tot
for (E,w2,w3) in [(Fr(3),Fr(2),Fr(-5)),(Fr(-7,2),Fr(1,3),Fr(19,6))]:
    for h,U in [(2,(3,)),(2,(2,2)),(3,(5,)),(2,(2,3))]:
        a=XiW(h,U+(1,),E,E,w2,w3); b=Fr(sum(U))/E**2*XiW(h,U,E,E,w2,w3)
        chk("s.7.3 scaled part 1",a==b,(h,U))
        D=sum(U)-len(U)
        a=XiW(h,U+(0,),E,E,w2,w3); b=Fr(D-2*h+2)/E**2*XiW(h,U,E,E,w2,w3)
        chk("s.7.3 scaled part 0",a==b,(h,U))
    # in-range vanishing at general scale, and mode A nonvanishing
    chk("Cor L3 in-range vanishing h=2 (2)",XiW(2,(2,),E,E,w2,w3)==0)
    chk("Cor L3 in-range vanishing h=3 (2,2,2)",XiW(3,(2,2,2),E,E,w2,w3)==0)
    v=XiW(2,(2,),E,-E,-w2,-w3)
    print("mode A Xi_2(2) at",E,w2,":",v)
    chk("mode A control nonzero",v!=0)
# conifold & P2 conventions
for x in [Fr(2,5),Fr(-3,7)]:
    c=1+x
    chk("s.7.4 conifold h=2 (2)",XiW(2,(2,),c,c,-x,Fr(-1))==0)
    chk("s.7.4 conifold h=3 (3,2)",XiW(3,(3,2),c,c,-x,Fr(-1))==0)
    chk("s.7.4 exponent",XiW(2,(3,),c,c,-x,Fr(-1))==c**(-(3*2-3+2))*peval(Xi_poly(2,(3,)),z_of_tau(-x/c)))
for (e1,e2) in [(Fr(1),Fr(1)),(Fr(2),Fr(3)),(Fr(1,3),Fr(-5))]:
    E=e1+e2
    chk("s.7.5 P2 h=2 (2)",XiW(2,(2,),E,E,-e1,-e2)==0)
    chk("s.7.5 P2 h=3 (4)",XiW(3,(4,),E,E,-e1,-e2)==0)
    chk("s.7.5 P2 Cor L2 h=3 (4,1,1)",XiW(3,(4,1,1),E,E,-e1,-e2)==0)
    chk("s.7.5 out of range nonzero h=2 (3)",XiW(2,(3,),E,E,-e1,-e2)!=0)
# tau->-1-tau symmetry of Z (s.2.2 (ii)) and tau=0: from F coefficient
for tau in [Fr(2,7),Fr(-5,3)]:
    F1,it1=compute_F((3,2),3,14,tau,K); F2,it2=compute_F((3,2),3,14,-1-tau,K)
    for key in F1:
        ev,a=key
        sub=()
        for (v,_),e in zip(it1,ev): sub+=(v,)*e
        sg=1
        for j in sub+(1,)*a: sg*=(-1)**(j-1)
        chk("s.2.2 (ii) nu->nu' symmetry",[sg*x for x in F1[key]]==F2[key],key)
F0,it0=compute_F((3,2),3,14,Fr(0),K)
for key,val in F0.items():
    ev,a=key
    sub=()
    for (v,_),e in zip(it0,ev): sub+=(v,)*e
    mu=sub+(1,)*a
    if len(mu)>=2:
        chk("s.2.2 (i) tau=0 no l>=2",all(x==0 for x in val),mu)
print("---- counts")
for k,v in counts.items(): print(k,v)
print("FAILS",len(Fails))
for f in Fails[:30]: print(f)
