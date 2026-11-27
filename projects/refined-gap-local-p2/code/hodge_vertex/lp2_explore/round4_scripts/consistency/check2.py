from lser import *
import sys
def pair(e1,e2,a,b):
    E=e1+e2
    return (e1*e2*E**3/9)*a[0]*b[0]+(-3*e1*e2*E)*a[1]*b[1]+(-3*e1*e2)*(a[1]*b[2]+a[2]*b[1])
for (e1,e2) in [(F(7),F(3)),(F(4),F(9))]:
    print("=== eps",e1,e2)
    N=9
    ev=Env(e1,e2,N);E=ev.E
    print("s(t_c)",ev.s.c[:3], "Del",ev.D.c[:3], "wc",ev.wc[:3])
    v=ev.v_of(0)
    rho=-E*v/ev.x
    print("rho coeffs",rho.c[:3]," E v:",[x*E for x in v.c[:3]]," expect -1,-1/18")
    rr=rec(30);V=E*v
    fact=1
    Sinf=L(0,[F(0)],10**6);Vp=V
    for l in range(0,20):
        fact*= (l+1)
        Sinf=Sinf+Vp*(rr[l]/fact);Vp=Vp*V
    print("Z-inf S_inf + t_c:",(Sinf+ev.x).c[:8],"prec",(Sinf+ev.x).prec)
    psis=[F(1,3),F(-2,5),F(3,7),F(5,11),F(-1,13),F(2,17),F(-3,19),F(4,23),F(7,29),F(-5,31),F(9,37),F(-7,41),F(11,43)]
    sig={}
    for p in psis:
        R=ev.rows(-p)
        sig[p]=[rho*x for x in R]
    p=psis[0]
    print("sigma valuations at psi=1/3",[x.val for x in sig[p]],"prec",[x.prec for x in sig[p]])
    lam=lambda j,ps:1/(E-j*ps)
    print(" sigma^(-1)=",[x.coef(-1) for x in sig[p]],"expected",[0,E,-E*E])
    print(" sigma^(0) =",[x.coef(0) for x in sig[p]],"expected",[E*lam(1,p),E*E/9*lam(1,p),0])
    for p in psis[:3]:
        a=[rho*x for x in ev.rows(-p)];b=[rho*x for x in ev.rows(p)]
        K=pair(e1,e2,a,b)
        K=K-L(-2,[3*e1*e2*E**3],K.prec)
        print(" K=1 residual prec",K.prec,"nonzero" if K.c else "zero")
    def fit(vals,pts,K):
        A=[[lam(j,pts[i]) for j in range(1,K+1)]+[vals[i]] for i in range(K)]
        n=K
        for c in range(n):
            piv=next(r for r in range(c,n) if A[r][c])
            A[c],A[piv]=A[piv],A[c]
            A[c]=[x/A[c][c] for x in A[c]]
            for r in range(n):
                if r!=c and A[r][c]:
                    A[r]=[x-A[r][c]*y for x,y in zip(A[r],A[c])]
        sol=[A[i][n] for i in range(n)]
        ok=all(sum(sol[j]*lam(j+1,pts[i]) for j in range(n))==vals[i] for i in range(n,len(pts)))
        return sol,ok
    leaf={}
    res=[]
    for k in range(0,5):
        for a,Kb in ((0,k+1),(1,k+1),(2,k)):
            vals=[sig[p][a].coef(k) for p in psis]
            if Kb==0:
                ok=all(x==0 for x in vals);sol=[]
            else:
                sol,ok=fit(vals,psis,Kb)
            leaf[(a,k)]=sol
            res.append(ok)
    print(" H2(b) leaf support k=0..4 all OK:",all(res))
    print(" sigma_2^(1) leaf:",leaf[(2,1)])
    print(" Z-inf Sum beta/j :",[sum(leaf[(0,k)][j]/(j+1) for j in range(k+1)) for k in range(5)])
    for label,c1 in (("on-shell",ev.g),("offshell",ev.g*(1+F(1,10)*ev.D))):
        C1i=c1.inv();th=ev.theta
        z=F(2,7)
        R=ev.rows(z,c1=c1)
        Sm1=ev.S(v,z,m=1);Sm2=ev.S(v,z,m=2)
        thv=th(v)
        P1=E*C1i*(thv/v)*Sm1
        d2=E*(C1i*C1i-1)/(2*ev.D)
        lg=th(C1i*thv)/(C1i*thv)
        P2=E*c1*ev.D*thv*(E*thv*Sm2/(v*v)+(z*lg-d2)*Sm1/v)
        print(" Lemma 4.1",label,"P1 diff",(P1-R[1]).c[:3],(P1-R[1]).prec,"P2 diff",(P2-R[2]).c[:3],(P2-R[2]).prec)
    print(" theta v - (1/E - v/3):",(ev.theta(v)-(1/E-v*F(1,3))).c[:3])
    def Ehat(ps,pp,c1=None):
        A=ev.rows(-ps,c1=c1);B=ev.rows(-pp,c1=c1)
        Ad=ev.rows(pp,c1=c1)
        Pzw=-1*pair(e1,e2,A,B)
        Pww=-1*pair(e1,e2,Ad,B)
        return (Pzw-Pww)*(1/(ps+pp))
    ps,pp=F(1,3),F(2,5)
    Eh=Ehat(ps,pp)
    print(" Ehat_00 val",Eh.val,"prec",Eh.prec,"at 0:",Eh.coef(0),"formula",-(e1*e2*E**4/27)*lam(1,ps)*lam(1,pp))
    Ebar=rho*rho*ev.x*ev.x*Eh
    for n in range(2,6):
        tot=F(0)
        for i in range(1,n):
            for l in range(1,n-i+1):
                C=F(0)
                for k1 in range(0,n-1):
                    k2=n-2-k1
                    a=[ (leaf[(q,k1)][i-1] if i-1<len(leaf[(q,k1)]) else F(0)) for q in range(3)]
                    b=[ (leaf[(q,k2)][l-1] if l-1<len(leaf[(q,k2)]) else F(0)) for q in range(3)]
                    C-=pair(e1,e2,a,b)
                tot+=F(i*l)/(E*(i+l))*C*lam(i,ps)*lam(l,pp)
        print("  E00(i) n=%d: direct-theorem ="%n, Ebar.coef(n)-tot)
    if (e1,e2)==(F(7),F(3)):
        st={0:[[F(-70000,9)]],
            1:[[F(70000,9),F(-280000,81)],[F(-280000,81),F(0)]],
            2:[[F(140000,81),F(140000,81),F(-350000,243)],[F(140000,81),F(-140000,81),F(0)],[F(-350000,243),F(0),F(0)]],
            3:[[F(-700000,2187),F(1400000,729),F(-70000,243),F(-1120000,2187)],[F(1400000,729),F(0),F(-560000,729),F(0)],[F(-70000,243),F(-560000,729),F(0),F(0)],[F(-1120000,2187),F(0),F(0),F(0)]]}
        ser=L(0,[F(0)],4)
        sp=L(0,[F(1)],10**6)
        for k in range(4):
            val=sum(st[k][i][l]*lam(i+1,ps)*lam(l+1,pp) for i in range(k+1) for l in range(k+1))
            ser=ser+sp*val
            sp=sp*ev.s
        d=Eh-ser
        print(" stored a_il (s^0..s^3) vs direct Ehat_00: residual prec",d.prec,"nonzero" if d.c else "zero",d.c[:3])
    for col in (1,2):
        w=F(2,5);ps=F(1,3)
        A=[x for x in ev.rows(-ps,col=0)];B=ev.rows(-w,col=col);Ad=ev.rows(w,col=0)
        Pww=pair(e1,e2,Ad,B)
        print(" orthogonality col",col,"prec",Pww.prec,"coefs",Pww.c[:3])
        Pz=-1*pair(e1,e2,[rho*x for x in A],B)
        E0j=Pz*(1/(ps+w))
        print("  rho Ehat_0j val",E0j.val)
        for k in range(0,3):
            tot=W(0)
            for k1 in range(0,k+1):
                for l in range(1,k1+2):
                    Sl=[(leaf[(q,k1)][l-1] if l-1<len(leaf[(q,k1)]) else F(0)) for q in range(3)]
                    vv=[B[q].coef(k-k1) for q in range(3)]
                    tot=tot-l*pair(e1,e2,Sl,vv)*lam(l,ps)/(E+l*w)
            print("   M-b k=%d diff"%k, E0j.coef(k)-tot)
        vj=ev.v_of(col)
        Dj=3*e1*e2*E/(vj*vj)
        om=W(0,1) if col==1 else W(-1,-1)
        print("  Delta_j(0)-formula",Dj.coef(0)-e1*e2*E**3*(om*om)/(3*(1-om)*(1-om)))
        print("  rows at col %d, t_c=0: P0,P1,P2 ="%col,B[0].coef(0),B[1].coef(0),B[2].coef(0))
    for label,c1 in (("on",ev.g),("off",ev.g*(1+F(1,10)*ev.D))):
        z=F(2,7);C1i=c1.inv();th=ev.theta
        Y=1/v-E/3
        d2=E*(C1i*C1i-1)/(2*ev.D)
        R=ev.rows(z,c1=c1)
        P=[v*x for x in R]
        XX=th(c1)/c1
        Dc=ev.D;yy=(Dc-1)*F(1,27)
        Rd=Dc*th(XX)+Dc*XX*XX+27*yy*XX+6*yy
        lhs=C1i*(Y*P[2]+z*th(P[2]))
        rhs=(-1*E**3/27)*P[0]+(-1*c1*(E*d2+Dc*d2*d2+E*E*F(1,3)))*P[1]+(-1*C1i*(E/Dc+d2))*P[2]-z*z*c1*Rd*P[1]
        print(" Prop 3.1 row2",label,"diff",(lhs-rhs).c[:3],"prec",(lhs-rhs).prec, "Rd val",Rd.val,"XX(0)",XX.coef(0))
