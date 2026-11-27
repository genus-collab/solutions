from lser import *
def pair(e1,e2,a,b):
    E=e1+e2
    return (e1*e2*E**3/9)*a[0]*b[0]+(-3*e1*e2*E)*a[1]*b[1]+(-3*e1*e2)*(a[1]*b[2]+a[2]*b[1])
def zrows(ev,col,Z,c1=None):
    """z-series rows: lists of L for P^_0,P^_1,P^_2, z^0..z^Z"""
    E=ev.E
    if c1 is None: c1=ev.g
    C1i=c1.inv();v=ev.v_of(col);Y=1/v-E/3;th=ev.theta;Ds=ev.D
    d2=E*(C1i*C1i-1)/(2*Ds)
    rr=rec(Z+3)
    # [z^n] 1/prod_{j<=l+1}(1+jz/E) by series
    S=[L(0,[F(0)],10**6) for _ in range(Z+1)]
    for l in range(Z+1):
        ser=[F(1)]+[F(0)]*Z
        for j in range(1,l+2):
            new=[F(0)]*(Z+1)
            for n in range(Z+1):
                new[n]=ser[n]-(F(j)/E*new[n-1] if n>=1 else 0)
            ser=new
        vl=v**l
        for n in range(Z+1-l):
            S[n+l]=S[n+l]+vl*(rr[l]*ser[n])
    P0B=[v*x for x in S]
    def Th(f,extra=None):
        out=[]
        for n in range(len(f)):
            t=Y*f[n]
            if n>=1: t=t+th(f[n-1])
            if extra is not None: t=t-extra*f[n]
            out.append(t)
        return out
    Q=[C1i*x for x in Th(P0B)]
    P2=[Ds*c1*c1*x for x in Th(Q,d2)]
    return [S,[x/v for x in Q],[x/v for x in P2]]
for (e1,e2) in [(F(7),F(3)),(F(4),F(9))]:
    print("=== eps",e1,e2)
    ev=Env(e1,e2,8);E=ev.E;Z=6
    v0=ev.v_of(0);rho=-E*v0/ev.x
    R0=zrows(ev,0,Z)
    lam=lambda j,ps:1/(E-j*ps)
    for col in (1,2):
        Rj=zrows(ev,col,Z)
        # orthogonality <P^0(w),P^j(-w)>=0 ; diag
        bad=[]
        for n in range(Z+1):
            t=L(-1,[F(0)],6)
            for p in range(n+1):
                q=n-p
                sg=(-1)**q
                t=t+sg*pair(e1,e2,[R0[a][p] for a in range(3)],[Rj[a][q] for a in range(3)])
            bad.append((t.prec,bool(t.c)))
        print(" colour",col,"orthogonality z^0..z^%d (prec,nonzero):"%Z,bad)
        print("  P^_2^(j)[z^n] at t_c=0:",[Rj[2][n].coef(0) for n in range(4)])
        # (0,j) with stored values at (7,3)
        if (e1,e2)==(F(7),F(3)):
            print("  P^_0^(j)(z;0) z^0..3:",[Rj[0][n].coef(0) for n in range(4)])
            # stored chi_j(w) w^0..w^3
            Nw=4
            # chi_j(w)= -<b,v_{j,0}(-w)>/(E+w), b=(E,E^2/9,0)
            b=[E,E*E/9,F(0)]
            coef=[]
            for n in range(Nw):
                coef.append(pair(e1,e2,b,[Rj[a][n].coef(0)*((-1)**n) for a in range(3)])*(-1))
            chi=[]
            for n in range(Nw):
                x=coef[n]
                if n>=1: x=x-chi[n-1]
                chi.append(x/E)
            print("  chi_j(w) w^0..3:",chi)
        # M-b at numeric psi, series in w up to Nw
        ps=F(1,3)
        # sigma(psi) rows numeric
        sig=[rho*x for x in ev.rows(-ps,col=0)]
        Nw=4
        # N(psi,w)=-<sigma(psi),P^j(-w)>; rho Ehat_0j=N/(psi+w): series in w
        Ncoef=[-1*pair(e1,e2,sig,[Rj[a][n]*((-1)**n) for a in range(3)]) for n in range(Nw+1)]
        Ecoef=[]
        for n in range(Nw+1):
            x=Ncoef[n]
            if n>=1: x=x-Ecoef[n-1]
            Ecoef.append(x*(1/ps))
        print("  rho Ehat_0j [w^b] valuations:",[x.val for x in Ecoef],"prec",Ecoef[0].prec)
        # leaf coefficients of sigma by fit
        psis=[F(1,3),F(-2,5),F(3,7),F(5,11),F(-1,13),F(2,17),F(-3,19),F(4,23)]
        sigs={p:[rho*x for x in ev.rows(-p,col=0)] for p in psis}
        def leafs(a,k):
            K=k+1
            A=[[lam(j,psis[i]) for j in range(1,K+1)]+[sigs[psis[i]][a].coef(k)] for i in range(K)]
            for c in range(K):
                piv=next(r for r in range(c,K) if A[r][c])
                A[c],A[piv]=A[piv],A[c]
                A[c]=[x/A[c][c] for x in A[c]]
                for r in range(K):
                    if r!=c and A[r][c]:
                        A[r]=[x-A[r][c]*y for x,y in zip(A[r],A[c])]
            return [A[i][K] for i in range(K)]
        for k in range(0,3):
            for bb in range(0,3):
                # RHS [w^b] of -Sum l <S_l,v_{j,k2}(-w)> lam_l(psi)/(E+l w)
                tot=W(0)
                for k1 in range(0,k+1):
                    lf=[leafs(q,k1) for q in range(3)]
                    for l in range(1,k1+2):
                        Sl=[lf[q][l-1] for q in range(3)]
                        for n in range(bb+1):
                            vv=[Rj[a][n].coef(k-k1)*((-1)**n) for a in range(3)]
                            # [w^(bb-n)] 1/(E+l w)
                            g=(F(-l)/E)**(bb-n)/E
                            tot=tot-l*pair(e1,e2,Sl,vv)*lam(l,ps)*g
                print("   M-b k=%d b=%d diff="%(k,bb),Ecoef[bb].coef(k)-tot)
    # off-shell controls: H2(d), pole of Ehat_00
    Q=F(1,10)
    f=1+Q*ev.D
    c1=ev.g*f
    ps,pp=F(1,3),F(2,5)
    A=ev.rows(-ps,c1=c1);B=ev.rows(-pp,c1=c1);Ad=ev.rows(pp,c1=c1)
    Eh=(-1*pair(e1,e2,A,B)+pair(e1,e2,Ad,B))*(1/(ps+pp))
    print(" off-shell Ehat_00 val",Eh.val,"coef t^-1:",Eh.coef(-1),"predicted 3 e1 e2 E^2 (XX0+2/9)/t_c:",3*e1*e2*E*E*(-Q))
    # psi-linear coefficient
    th=ev.theta
    pred=E*ev.D*ev.g*ev.g*th(f)/ev.x
    psis=[F(1,3),F(-2,5),F(3,7),F(5,11),F(-1,13),F(2,17),F(-3,19),F(4,23),F(7,29),F(-5,31),F(9,37),F(-7,41)]
    sg={p:(rho*ev.rows(-p,c1=c1)[2]) for p in psis}
    for k in range(-1,3):
        K=k+3 if k>=-1 else 1
        K=max(K,1)
        # basis psi,1,lam_1..lam_K
        n=K+2
        A=[[psis[i],F(1)]+[lam(j,psis[i]) for j in range(1,K+1)]+[sg[psis[i]].coef(k)] for i in range(n)]
        for c in range(n):
            piv=next(r for r in range(c,n) if A[r][c])
            A[c],A[piv]=A[piv],A[c]
            A[c]=[x/A[c][c] for x in A[c]]
            for r in range(n):
                if r!=c and A[r][c]:
                    A[r]=[x-A[r][c]*y for x,y in zip(A[r],A[c])]
        sol=[A[i][n] for i in range(n)]
        chk=all(sol[0]*psis[i]+sol[1]+sum(sol[j+1]*lam(j,psis[i]) for j in range(1,K+1))==sg[psis[i]].coef(k) for i in range(n,len(psis)))
        print("  off-shell sigma_2 k=%d: psi-linear fit %s predicted %s  (fit consistent: %s)"%(k,sol[0],pred.coef(k),chk))
