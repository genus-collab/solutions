# Test (v): Lemma 3.1 (series check), 3.2 (incl non-integer, part 0), 3.3, 3.4, Cor 3.5, Rem 3.6, s.7.2, s.7.3, table s.8.4, s.8.2 vs definition, s.8.5.
from hodge import *
from mk import X0
import sys
Fails=[]
counts=Counter()
def chk(name,cond,info=""):
    counts[name]+=1
    if not cond:
        Fails.append((name,info)); print("FAIL",name,info,flush=True)

ORD=40
def ser_mul(a,b,P=ORD):
    r=[Fr(0)]*P
    for i,x in enumerate(a[:P]):
        if x:
            for j in range(P-i):
                y=b[j] if j<len(b) else 0
                if y: r[i+j]+=x*y
    return r
def one_plus_s_pow(n,P=ORD):
    """(1+s)^n, n integer"""
    r=[Fr(0)]*P
    c=Fr(1)
    for i in range(P):
        r[i]=c
        c=c*(n-i)/(i+1)
    return r
U_SER=None
def u_series(P=ORD):
    r=[Fr(0)]*P
    # s/(1+s)
    for i in range(1,P):
        r[i]=Fr((-1)**(i-1))
    return r
def ser_exp(a,P=ORD):
    # a[0]=0
    e=[Fr(0)]*P; e[0]=Fr(1)
    for n in range(1,P):
        s=Fr(0)
        for i in range(1,n+1):
            if a[i]: s+=i*a[i]*e[n-i]
        e[n]=s/n
    return e
US=u_series()
def e_u(X,P=ORD):
    return ser_exp([Fr(X)*x for x in US],P)
def poly_from_Q(g,U,z):
    """coefficient list in s of Q_g(s;U) at z"""
    q=Q_poly(g,U)
    N=3*g-3+len(U)
    out=[Fr(0)]*(N+1)
    for m,p in q.items():
        out[m]=peval(p,z)
    return out

def M_a(g,U,a,z):
    mu=tuple(U)+(1,)*a
    return peval(M_poly(g,mu),z) if g>0 else Fr(sum(mu))**(len(mu)-3)

def rhs_31(g,U,z,Ord):
    """(1+s)^{-(2g-2+r)} e^{-uJ} Sum_a M_a (u e^{-u})^a/a! as series to order Ord in s"""
    r=len(U); J=sum(U,Fr(0)); n=2*g-2+r
    tot=[Fr(0)]*Ord
    for a in range(Ord):
        Ma=M_a(g,U,a,z)
        # s^a (1+s)^{-n-a} e^{-u(J+a)}
        base=ser_mul(one_plus_s_pow(-n-a,Ord),e_u(-J-a,Ord),Ord)
        for i in range(Ord-a):
            tot[i+a]+=Ma*base[i]/factorial(a)
    return tot

taus=[Fr(2,7),Fr(-5,3),Fr(-1,2),Fr(-1,3)]
# Lemma 3.2 incl non-integer leaves and parts 0, parts 1
L32cases=[(1,(1,)),(1,(2,)),(1,(3,)),(1,(2,2)),(1,(Fr(1,2),)),(1,(0,)),(1,(0,3)),(1,(1,1)),(2,()),(2,(1,)),(2,(2,)),(2,(3,)),(2,(0,)),(2,(Fr(1,2),)),(2,(Fr(1,3),2)),(2,(2,2)),(3,()),(3,(2,)),(3,(5,)),(3,(0,2)),(3,(Fr(3,2),)),(0,(2,3,4)),(0,(2,2,2,2)),(0,(1,1,1))]
for tau in taus:
    z=z_of_tau(tau)
    for g,U in L32cases:
        N=3*g-3+len(U)
        Ord=N+6
        if g==0:
            # M_a from genus 0 multinomial for r>=3; Q from definition with engine
            pass
        Qdef=poly_from_Q(g,U,z)
        # for non-integer U in M_a the genus-0 formula also holds
        rhs=rhs_31(g,U,z,Ord)
        ok=all(rhs[i]==(Qdef[i] if i<len(Qdef) else 0) for i in range(Ord))
        chk("Lemma 3.2 (3.1)",ok,(tau,g,U))
        # Laguerre X_g
        r=len(U); J=sum(U,Fr(0)); X=sum(Fr((-1)**i)*c for i,c in enumerate(Qdef))
        def L(nn,al,x):
            return sum((-1)**i*comb(nn+al,nn-i)*Fr(x)**i/factorial(i) for i in range(nn+1))
        X2=sum(Fr((-1)**a,factorial(a))*L(N-a,2*g-2+r+a,-J-a)*M_a(g,U,a,z) for a in range(N+1))
        chk("Laguerre X_g",X==X2,(tau,g,U))
        for m in range(N+1,N+4):
            v=sum(M_a(g,U,a,z)/factorial(a)*(-1)**(m-a)*L(m-a,2*g-3+r+a,-J-a) for a in range(m+1))
            chk("Laguerre [s^m]=0 for m>N",v==0,(tau,g,U,m))

# Lemma 3.3 / X_0 vs actual Q_0(-1) for r>=3
for U in [(2,3,4),(2,2,2),(1,1,1),(3,3,3,3),(2,2,2,2),(5,2,1,1,1),(2,2,2,2,2)]:
    q=Xi_poly(0,U)
    val=q[0] if q else Fr(0)
    chk("X_0 = Xi_0 r>=3",val==X0(U),(U,val,X0(U)))

# Lemma 3.4: e^{-uJ}(J+theta)^k G^(g)_U = w^{2g-2+r+k}( (delta)_k X_g + w poly(w))
def theta_coeffs(J,Ma,Ord):
    return Ma
def lemma34(g,U,z,kmax):
    r=len(U); J=sum(U); D=J-r; delta=D-2*g+2
    N=3*g-3+r if g>=1 else 0
    for k in range(kmax+1):
        n0=2*g-2+r+k
        # poly degree expected: Q degree N plus k: total degree of polynomial in w at most N+k (g>=1). For g=0, r-2+... small
        degbound=(N+k) if g>=1 else (k+max(r,3))
        Ord=degbound+6
        tot=[Fr(0)]*Ord
        for a in range(Ord):
            Ma=M_a(g,U,a,z)*Fr(J+a)**k
            base=ser_mul(one_plus_s_pow(-n0-a,Ord),e_u(-J-a,Ord),Ord)
            for i in range(Ord-a):
                tot[i+a]+=Ma*base[i]/factorial(a)
        # tot = series in s of polynomial in w=1+s; check truncation
        ok=all(tot[i]==0 for i in range(degbound+1,Ord))
        chk("Lemma 3.4 polynomiality",ok,(g,U,k))
        val=sum(Fr((-1)**i)*tot[i] for i in range(degbound+1))
        if g==0:
            Xg=X0(tuple(U))
        else:
            Xg=peval(Xi_poly(g,tuple(sorted(U))),z)
        chk("Lemma 3.4 leading coeff",val==falling(delta,k)*Xg,(g,U,k,val,falling(delta,k)*Xg))
for tau in [Fr(2,7),Fr(-1,2)]:
    z=z_of_tau(tau)
    for g,U in [(0,(1,)),(0,(2,)),(0,(3,)),(0,(1,1)),(0,(2,3)),(0,(1,1,1)),(0,(2,3,4)),(0,(1,1,1,1)),(1,(1,)),(1,(2,)),(1,(1,1)),(1,(2,3)),(2,(1,)),(2,(2,)),(2,(3,)),(2,(1,1)),(2,(2,2)),(2,(4,)),(3,(1,)),(3,(2,)),(3,(5,)),(3,(2,3))]:
        lemma34(g,U,z,6 if g<3 else 4)

# Cor 3.5, (S), Rem 3.6 raw
for tau in taus:
    z=z_of_tau(tau)
    for g in (1,2,3,4):
        for U in [(2,),(3,),(5,),(2,2),(2,3),(1,),(1,2),(1,1),(7,)]:
            if g==4 and len(U)>1: continue
            a=peval(Xi_poly(g,tuple(sorted(U+(1,)))),z)
            b=sum(U)*peval(Xi_poly(g,tuple(sorted(U))),z)
            chk("Cor 3.5 X(U+1)=J X(U)",a==b,(tau,g,U))
    for g in (2,3,4):
        for v in range(1,6):
            chk("Cor 3.5 X_g(1^v)=0",not Xi_poly(g,(1,)*v),(g,v))
    for g in (2,3):
        for Up,v in [((2,),3),((3,),2),((2,2),2),((4,3),3),((7,),4)]:
            a=peval(Xi_poly(g,tuple(sorted(Up+(1,)*v))),z)
            b=rising(sum(Up),v)*peval(Xi_poly(g,tuple(sorted(Up))),z)
            chk("(S) rising",a==b,(g,Up,v))
    # Rem 3.6
    sigma=tau*(1+tau)
    X11=peval(Xi_poly(1,(1,)),z)
    chk("Rem 3.6 X_1(1)",X11==-(1+1/sigma)/24,(tau,X11))
    for v in range(1,7):
        chk("Rem 3.6 X_0(1^v)",X0((1,)*v)==Fr(factorial(v-1),2),v)
        lc=[falling(2,k)*X0((1,)*v) for k in range(5)]
        exp=[Fr(factorial(v-1),2),Fr(factorial(v-1)),Fr(factorial(v-1)),0,0]
        chk("Rem 3.6 table g=0",lc==exp,v)
        chk("Rem 3.6 g=1",peval(Xi_poly(1,(1,)*v),z)==factorial(v-1)*X11,v)

# s.7.2 operators on Q and Xi_g(U+0)
for tau in [Fr(2,7),Fr(-1,2)]:
    z=z_of_tau(tau)
    for g,U in [(1,(2,)),(1,(3,2)),(1,(Fr(1,2),)),(2,()),(2,(2,)),(2,(3,)),(2,(Fr(1,2),2)),(3,()),(3,(2,)),(3,(4,))]:
        r=len(U);J=sum(U,Fr(0));n=2*g-2+r
        Q=poly_from_Q(g,U,z)
        dQ=[(i)*Q[i] for i in range(1,len(Q))]+[Fr(0)]
        def padd2(a,b):
            L=max(len(a),len(b));return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(L)]
        def pm(a,b):
            o=[Fr(0)]*(len(a)+len(b)-1)
            for i,x in enumerate(a):
                for j,y in enumerate(b):
                    o[i+j]+=x*y
            return o
        Q1=padd2(pm([Fr(1),Fr(2),Fr(1)],dQ),pm([J+n,Fr(n)],Q))
        Qn=poly_from_Q(g,U+(1,),z)
        L=max(len(Q1),len(Qn))
        chk("s.7.2 part 1 operator",all((Q1[i] if i<len(Q1) else 0)==(Qn[i] if i<len(Qn) else 0) for i in range(L)),(g,U))
        Q0=padd2(pm([J,Fr(n)],Q),pm([Fr(0),Fr(1),Fr(1)],dQ))
        Qn0=poly_from_Q(g,U+(0,),z)
        L=max(len(Q0),len(Qn0))
        chk("s.7.2 part 0 operator",all((Q0[i] if i<len(Q0) else 0)==(Qn0[i] if i<len(Qn0) else 0) for i in range(L)),(g,U))
        a=sum(Fr((-1)**i)*c for i,c in enumerate(Qn0)); b=(J-r-2*g+2)*sum(Fr((-1)**i)*c for i,c in enumerate(Q))
        chk("s.7.2 Xi(U+0)=(D-2g+2)Xi(U)",a==b,(g,U))

# Cor L1 specialisation at tau=-1/2 etc is automatic from polys. Cor L3: direct general-scale integral
def XiE(h,U,E,w2,w3):
    """Direct: Omega_w=prod Lambda^vee(w_a)/w_a^h; lam^E; T^E. Exponential of ch with t=-1/w_a."""
    # power sums
    def cc(l): return (Fr(-1)/E)**l+(Fr(-1)/w2)**l+(Fr(-1)/w3)**l
    def om(chs):
        r=Fr(1)
        for l,m in Counter(chs).items():
            r*=(factorial(l-1)*cc(l))**m/factorial(m)
        return r
    r=len(U); N=3*h-3+r
    tot=Fr(0)
    # leaves 1/(E-j psi)=Sum j^e psi^e/E^{e+1}; T^E=-psi^2/(E-psi)=-Sum_{d>=2} psi^d/E^{d-1}
    for e,wt in leaf_groups(U,N).items():
        # wt contains prod j^e
        Mleft=N-sum(e)
        for M1 in range(Mleft+1):
            for par in partitions_ge(M1,M1 if M1>0 else 1):
                ds=tuple(sorted(e+tuple(p+1 for p in par)))
                w=Fr(1)
                for v,c in Counter(par).items(): w/=factorial(c)
                w*=Fr((-1)**len(par))
                # E powers
                ex=-sum(x+1 for x in e)-sum(p for p in par)
                deg=Mleft-M1
                s=Fr(0)
                for chs in odd_partitions(deg,2*h-1):
                    v=chint(h,chs,ds)
                    if v: s+=om(chs)*v
                tot+=wt*w*s*Fr(E)**ex
    return tot
for (E,w2,w3) in [(Fr(3),Fr(2),Fr(-5)),(Fr(-7,2),Fr(1,3),Fr(19,6)),(Fr(1,5),Fr(1,5),Fr(-2,5))]:
    tau=w2/E; z=z_of_tau(tau)
    for h,U in [(2,(3,)),(2,(2,2)),(3,(5,)),(2,(2,)),(3,(3,)),(2,(1,3)),(2,())]:
        a=XiE(h,U,E,w2,w3)
        b=Fr(E)**(-(3*h-3+2*len(U)))*peval(Xi_poly(h,U),z)
        chk("Cor L3 scale",a==b,(E,w2,h,U,a,b))

# Table s.8.4
def QQ_series(g,R,z,Ord):
    r=len(R);J=sum(R);n=2*g-2+r
    tot=[Fr(0)]*Ord
    for a in range(Ord):
        if (not R) and a==0: continue
        Ma=M_a(g,R,a,z)
        base=ser_mul(one_plus_s_pow(-n-a,Ord),e_u(-J-a,Ord),Ord)
        for i in range(Ord-a):
            tot[i+a]+=Ma*base[i]/factorial(a)
    return tot
def ser_eq(a,b): return all(x==y for x,y in zip(a,b))
for tau in [Fr(2,7),Fr(-5,3),Fr(-1,2)]:
    z=z_of_tau(tau); sigma=tau*(1+tau); Ord=14
    for g,R in [(1,(2,)),(1,(3,2)),(2,(2,)),(2,(3,3)),(3,(2,)),(0,(2,3,4)),(0,(2,2,2,2))]:
        Q=poly_from_Q(g,R,z); Q=Q+[Fr(0)]*(Ord-len(Q))
        S=QQ_series(g,R,z,Ord)
        chk("s.8.4 row r>=1",ser_eq(Q[:Ord],S),(g,R))
    for j in (2,3,5,8):
        exp=[(j-1-1/sigma)/24,Fr(1,24)]+[Fr(0)]*(Ord-2)
        chk("s.8.4 (1,(j))",ser_eq(exp,QQ_series(1,(j,),z,Ord)),j)
        exp2=[Fr(1,j**2),Fr(1,j**2)-Fr(1,j*(j+1))]+[Fr(0)]*(Ord-2)
        chk("s.8.4 (0,(j))",ser_eq(exp2,QQ_series(0,(j,),z,Ord)),j)
        chk("s.8.4 (0,(j)) value",(exp2[0]-exp2[1])==Fr(1,j*(j+1)),j)
        for k in (2,3):
            exp3=[Fr(1,j+k)]+[Fr(0)]*(Ord-1)
            chk("s.8.4 (0,(j,k))",ser_eq(exp3,QQ_series(0,(j,k),z,Ord)),(j,k))
    # (0,empty)
    exp=[Fr(0)]*Ord
    exp[1]+=1; exp[2]+=Fr(1,4)
    ser=[Fr(0)]*Ord
    ser6=one_plus_s_pow(-1,Ord)
    for i in range(Ord-3): exp[i+3]+=ser6[i]/6
    # M_0(1^a)=a^{a-3}
    tot=[Fr(0)]*Ord
    for a in range(1,Ord):
        base=ser_mul(one_plus_s_pow(2-a,Ord),e_u(-a,Ord),Ord)
        for i in range(Ord-a):
            tot[i+a]+=Fr(a)**(a-3)*base[i]/factorial(a)
    chk("s.8.4 (0,empty)",ser_eq(exp,tot),tau)
    # (1,empty): (1/24)[log(1+s)+(1+1/sigma)(1/(1+s)-1)]
    exp=[Fr(0)]*Ord
    for i in range(1,Ord):
        exp[i]+=Fr((-1)**(i-1),i)/24
        exp[i]+=(1+1/sigma)*Fr((-1)**i)/24
    tot=QQ_series(1,(),z,Ord) if False else None
    tot=[Fr(0)]*Ord
    for a in range(1,Ord):
        Ma=M_a(1,(),a,z)
        base=ser_mul(one_plus_s_pow(-a,Ord),e_u(-a,Ord),Ord)
        for i in range(Ord-a):
            tot[i+a]+=Ma*base[i]/factorial(a)
    chk("s.8.4 (1,empty)",ser_eq(exp,tot),tau)
    # g>=2 empty
    for g in (2,3):
        Ord2=3*g-3+6
        Q=poly_from_Q(g,(),z); Q=Q+[Fr(0)]*(Ord2-len(Q))
        om=peval(int_omega_poly(g,()),z)
        pw=one_plus_s_pow(-(2*g-2),Ord2)
        exp=[Q[i]-om*pw[i] for i in range(Ord2)]
        tot=[Fr(0)]*Ord2
        for a in range(1,Ord2):
            Ma=M_a(g,(),a,z)
            base=ser_mul(one_plus_s_pow(-(2*g-2)-a,Ord2),e_u(-a,Ord2),Ord2)
            for i in range(Ord2-a):
                tot[i+a]+=Ma*base[i]/factorial(a)
        chk("s.8.4 g>=2 empty",ser_eq(exp,tot),(tau,g))
        # s.8.5 <Omega_g>=-A_g sigma^{1-g}
    A={2:Fr(1,5760),3:Fr(1,1451520),4:Fr(1,87091200)}
    for g in (2,3,4):
        chk("s.8.5 <Omega_g>",peval(int_omega_poly(g,()),z)==-A[g]*sigma**(1-g),(g,tau))
    # s.8.5 determined by MV coefficients: relation at m=N+1
    for h in (2,3):
        N=3*h-3
        m=N+1
        def L(nn,al,x):
            return sum((-1)**i*comb(nn+al,nn-i)*Fr(x)**i/factorial(i) for i in range(nn+1))
        v=sum(M_a(h,(),a,z)/factorial(a)*(-1)**(m-a)*L(m-a,2*h-3+a,-a) for a in range(m+1))
        chk("s.8.5 relation m=N+1",v==0,(h,tau))
        chk("s.8.5 a=0 coefficient",(-1)**m*L(m,2*h-3,0)*(-1)**m==comb(N+2*h-2,N+1),h)
    # s.8.5 R2 Cor 7 gap formula (cited): sigma^h Xi_h(empty) = -(2h-3)! [l^{2h}] S(e1 l)S(e2 l)
    def Scoef(j):
        # y/2 / sinh(y/2): coefficient of y^{2j}: -(2^{2j-1}-1)B_{2j}/(2^{2j-1}(2j)!)
        if j==0: return Fr(1)
        return -(Fr(2**(2*j-1)-1,2**(2*j-1)))*BERN[2*j]/factorial(2*j)
    for h in (2,3,4,5):
        e1=-tau;e2=1+tau
        cf=sum(Scoef(j)*Scoef(h-j)*e1**(2*j)*e2**(2*(h-j)) for j in range(h+1))
        chk("s.8.5 R2 Cor 7 formula",sigma**h*peval(Xi_poly(h,()),z)==-factorial(2*h-3)*cf,(h,tau))
print("---- counts")
for k,v in counts.items(): print(k,v)
print("FAILS",len(Fails))
for f in Fails[:40]: print(f)
