# derived_checks.py: numerical checks of the by-hand generic formulas: alg, delta, top leaf, Frobenius values at t=0, row values at t=0, sigma^(0)
import sys,os,random
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
OUT=os.path.dirname(os.path.abspath(__file__))
fs=open(os.path.join(OUT,"derived_checks_out.txt"),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
M=10;K=7;N=K+10
random.seed(5)
vals=[random.randrange(2,P) for _ in range(4)]
Wr=dict(a1=vals[0],a2=vals[1],e1=vals[2],e2=vals[3]); Wr['a0']=(Wr['e1']+Wr['e2']-Wr['a1']-Wr['a2'])%P
for name,W in (('fibrewise',cyw(3,5,3,5)),('S1',cyw(3,4,5,11)),('S2',cyw(7,13,2,29)),('S3',cyw(17,-5,19,-8)),('random',Wr)):
    a0,a1,a2,e1,e2=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']
    E=(a0+a1+a2)%P; A=a1; B=(a0+a2)%P; C=a2; u=(E-2*a1)*(E-2*a2)%P
    rows=Rows(W,M,N=N)
    wB=[E,-e1%P,-e2%P]
    fr=Frame(rows,{1:wB,2:wB},gmax=2)
    Q=rows.Q
    tau=s_add(s_const(1,N+4),s_neg(s_inv(Q)))
    tt=(1,[1]+[0]*(N+3))
    U1=rows.U[1]
    ss=s_mul(s_add(s_const(1,N+4),s_neg(U1)),s_inv(tau))   # s* = (1-U1)/tau
    pr("===",name)
    # critical equation of the Euler phase
    crit=s_add(s_add(s_scale((a0+a1),s_mul(tau,s_mul(ss,ss))),s_neg(s_mul(s_add(s_const(E,N+4),s_scale(a0,tau)),ss))),s_const(B,N+4))
    pr("  (a0+a1) tau s*^2 - (E + a0 tau) s* + (a0+a2) = 0:",s_iszero(crit),"; s*(0)=(a0+a2)/E:",s_coef(ss,0)==B*inv(E)%P)
    # alg^2 closed form
    oneminus=s_add(s_const(1,N+4),s_neg(ss))
    alg2=s_mul(U1,s_add(s_scale(A,s_mul(s_add(s_const(1,N+4),s_neg(tau)),s_inv(s_mul(oneminus,oneminus)))),s_scale(B,s_inv(s_mul(ss,ss)))))
    alg2=s_scale(A*B*inv(E*E*E%P),alg2)
    lhs=s_mul(s_mul(tau,tau),s_scale(inv(-E*E*E*e1*e2),rows.Delta[1]))
    pr("  Delta_1 tau^2/(-E^3 e1 e2) = (A B/E^3) U1 [A(1-tau)/(1-s*)^2 + B/s*^2] =: alg^2:",s_iszero(s_add(lhs,s_neg(alg2))),"; alg^2(0)=1:",s_coef(alg2,0)==1)
    # delta
    pp=s_mul(s_scale(C,ss),s_inv(s_add(s_const(1,N+4),s_neg(s_mul(tau,ss)))))  # Sph_tau
    Ptr=s_add(s_add(s_scale(2,s_mul(s_mul(tau,s_add(s_const(1,N+4),s_neg(tau))),pp)),s_const(E,N+4)),s_scale(-(B+C),tau))
    d2=s_add(s_add(s_scale(E*E,s_add(s_const(1,N+4),s_neg(tau))),s_scale(-u,tau)),s_scale(a0*a0,s_mul(tau,tau)))
    pr("  eikonal tau(1-tau)p^2+(E-(a0+2a2... )tau)p - BC=0 at p=C s*/(1-tau s*):",s_iszero(s_add(s_add(s_mul(s_mul(tau,s_add(s_const(1,N+4),s_neg(tau))),s_mul(pp,pp)),s_mul(s_add(s_const(E,N+4),s_scale(-(B+C),tau)),pp)),s_const(-B*C,N+4))))
    pr("  P_tr^2 = delta^2 = E^2(1-tau) - u tau + a0^2 tau^2:",s_iszero(s_add(s_mul(Ptr,Ptr),s_neg(d2))),"; P_tr(0)=E:",s_coef(Ptr,0)==E)
    # t=0 values
    pr("  L_2(0)=a1 a2/E:",s_coef(rows.L[2],0)==a1*a2*inv(E)%P,"; Delta_2(0)=e1 e2 a1 a2 (a0+a1)(a0+a2)/E:",s_coef(rows.Delta[2],0)==e1*e2*a1*a2*(a0+a1)*(a0+a2)*inv(E)%P)
    # top leaf: from the closed formula against the fit
    rm2=(rows.Delta[1][0]+2,[x*inv(-E*E*E*e1*e2)%P for x in rows.Delta[1][1]])
    def s_sqrt(s):
        c=s[1]; n=len(c); r=[1]+[0]*(n-1); h=inv(2)
        for k in range(1,n):
            t_=c[k]
            for j in range(1,k): t_-=r[j]*r[k-j]
            r[k]=t_%P*h%P
        return (0,r)
    rho=s_inv(s_sqrt(rm2))
    G=[s_mul(rho,fr.r[1][0][m]) for m in range(M+1)]
    ok=True
    for k in range(1,6):
        # solve for leaves via Vandermonde on z^1..z^(k+1): use general elimination
        import fractions
        n=k+1
        Amat=[[pow(j,m,P)*pow(-1,m,P)%P*pow(inv(E),m+1,P)%P for j in range(1,n+1)]+[0 if m==0 else s_coef(G[m-1],k)] for m in range(n)]
        # Gauss
        for c in range(n):
            p=next(i for i in range(c,n) if Amat[i][c]%P)
            Amat[c],Amat[p]=Amat[p],Amat[c]
            iv=inv(Amat[c][c]); Amat[c]=[x*iv%P for x in Amat[c]]
            for i in range(n):
                if i!=c and Amat[i][c]%P:
                    f=Amat[i][c]; Amat[i]=[(x-f*y)%P for x,y in zip(Amat[i],Amat[c])]
        top=Amat[k][n]
        Pk=1
        for j in range(1,k+1):
            Pk=Pk*(((k+1)*C-j*E)%P)%P*(((k+1)*B-j*E)%P)%P
        pred=pow(-1,k+1,P)*Pk%P*inv((k+1)*math.factorial(k)**2%P)%P*pow(inv(E),2*k-2,P)%P
        if pred!=top: ok=False; pr("   top mismatch k",k,fmt(top),fmt(pred))
    pr("  top leaf a_{k+1}^(k) = (-1)^(k+1) P_k/((k+1) k!^2 E^(2k-2)), P_k = prod_{j<=k} ((k+1)a2 - jE)((k+1)(a0+a2) - jE), k<=5:",ok)
    pr("  a^(1) = (u/2)(1,-1), rho^-1 = 1 - u t/(2E^2) + ...: ",fmt(s_coef(s_inv(rho),1)),fmt(-u*inv(2*E*E%P)%P))
    # colour 2 vertex values: r'_{1,2}(z;0)=L2(0) r'_{0,2}(z;0) (the analogue of Prop 8.6(b))
    ok=all(s_coef(fr.r[2][1][m],0)==a1*a2*inv(E)%P*s_coef(fr.r[2][0][m],0)%P for m in range(M+1))
    pr("  v_2(z;0) proportional to (1, L_2(0)):",ok)
    # e1,e2 independence of the rows
    W3=dict(W);
    shift=(e1*3+7)%P
    W3['e1']=(e1+shift)%P; W3['e2']=(e2-shift)%P
    rows3=Rows(W3,M,N=N);
    fr3=Frame(rows3,{1:[E,-W3['e1']%P,-W3['e2']%P],2:[E,-W3['e1']%P,-W3['e2']%P]},gmax=2)
    same=all(all(s_iszero(s_add(fr.r[c][a][m],s_neg(fr3.r[c][a][m]))) for m in range(M+1)) for c in (1,2) for a in (0,1))
    pr("  frame-B rows do not depend on (e1,e2) at fixed (a0,a1,a2,E):",same)
