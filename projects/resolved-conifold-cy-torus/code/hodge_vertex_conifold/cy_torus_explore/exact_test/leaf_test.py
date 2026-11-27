# leaf_test.py: F1, F2, F3, L1, L2, L3, QDE, symplecticity, E11 flow, closed form of the singular row at generic CY weights,
# directly from the Wick rows (nothing is assumed about the rows beyond their definition).
import sys,os,time,random
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
OUT=os.path.dirname(os.path.abspath(__file__))
fs=open(os.path.join(OUT,"leaf_test_out.txt"),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
M=int(sys.argv[1]) if len(sys.argv)>1 else 12
K=int(sys.argv[2]) if len(sys.argv)>2 else 6
N=K+8
def solve(A,b):
    """Gaussian elimination mod P on a rectangular system; returns (solution, rank) or None if inconsistent"""
    m=len(A); n=len(A[0]) if A else 0
    R=[row[:]+[b[i]] for i,row in enumerate(A)]
    piv=[];r=0
    for c in range(n):
        p=None
        for i in range(r,m):
            if R[i][c]%P: p=i;break
        if p is None: continue
        R[r],R[p]=R[p],R[r]
        iv=inv(R[r][c])
        R[r]=[x*iv%P for x in R[r]]
        for i in range(m):
            if i!=r and R[i][c]%P:
                f=R[i][c]
                R[i]=[(x-f*y)%P for x,y in zip(R[i],R[r])]
        piv.append(c);r+=1
    for i in range(r,m):
        if R[i][n]%P: return None
    sol=[0]*n
    for i,c in enumerate(piv): sol[c]=R[i][n]
    return sol,len(piv)
def check_sample(name,W):
    t0=time.time()
    a0,a1,a2,e1,e2=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']
    E=(a0+a1+a2)%P; A=a1; B=(a0+a2)%P; C=a2
    u=(E-2*a1)*(E-2*a2)%P
    rows=Rows(W,M,N=N)
    wB=[E,-e1%P,-e2%P]
    fr=Frame(rows,{1:wB,2:wB},gmax=2)
    pr("\n=== ",name," rows+frame",round(time.time()-t0,1),"s; M=",M,"K=",K)
    tau=s_add(s_const(1,N),s_neg(s_inv(rows.Q)))  # tau=1-1/Q
    taui=s_inv(tau)
    L1,L2=rows.L[1],rows.L[2]
    print_ok=lambda lab,x: pr("  ",lab,":", x)
    print_ok("L1+L2 = E/tau - a0",s_iszero(s_add(s_add(L1,L2),s_add(s_scale(-E,taui),s_const(a0,N)))))
    print_ok("L1 L2 = a1 a2/tau",s_iszero(s_add(s_mul(L1,L2),s_scale(-a1*a2,taui))))
    print_ok("eta^-1 = e1 e2 {{a0 a1 a2, a1 a2},{a1 a2, -E}}",rows.etainv==[[e1*e2*a0*a1*a2%P,e1*e2*a1*a2%P],[e1*e2*a1*a2%P,-e1*e2*E%P]])
    D1=rows.Delta[1]; D2=rows.Delta[2]
    c0=(-E*E*E*e1*e2)%P
    print_ok("F1: t^2 Delta_1 -> -E^3 e1 e2",D1[0]==-2 and D1[1][0]==c0)
    rm2=(D1[0]+2,[x*inv(c0)%P for x in D1[1]])
    def s_sqrt(s):
        assert s[0]==0 and s[1][0]==1
        c=s[1]; n=len(c); r=[1]+[0]*(n-1); h=inv(2)
        for k in range(1,n):
            tot=c[k]
            for j in range(1,k): tot-=r[j]*r[k-j]
            r[k]=tot%P*h%P
        return (0,r)
    rhoi=s_sqrt(rm2); rho=s_inv(rhoi)
    pr("   rho^-1 = 1 + (%s) t, -u/(2E^2) = %s"%(fmt(rhoi[1][1]),fmt(-u*inv(2*E*E%P)%P)))
    for c in (1,2):
        U=rows.U[c]; Q=rows.Q
        QU=s_mul(Q,U); X1=s_mul(QU,s_inv(s_add(s_const(1,N+4),s_neg(QU)))); X2=s_mul(U,s_inv(s_add(s_const(1,N+4),s_neg(U))))
        W2=s_add(s_scale(-A,s_add(X1,s_mul(X1,X1))),s_scale(-B,s_add(X2,s_mul(X2,X2))))
        print_ok("Delta_%d = e1 e2 a1 (a0+a2) W_2(U_%d) ?"%(c,c),s_iszero(s_add(rows.Delta[c],s_scale(-e1*e2*A*B,W2))))
    pr("   Delta_2(0) =",fmt(s_coef(D2,0)),"; Delta_2'(0) =",fmt(s_coef(D2,1)))
    Th=s_scale(e1*e2,s_add(s_const(a1*a2,N),s_scale(-E,L2)))
    print_ok("Delta_2 = (L2-L1) Theta, Theta=e1 e2(a1 a2-E L2)",s_iszero(s_add(D2,s_neg(s_mul(s_add(L2,s_neg(L1)),Th)))))
    print_ok("Theta(0)=0",s_coef(Th,0)==0)
    pr("   Theta_1 =",fmt(s_coef(Th,1)))
    ei=rows.etainv
    r=fr.r
    def pair(i,j,p,q):
        tot=(10**6,[])
        for a in (0,1):
            for b in (0,1):
                if ei[a][b]%P: tot=s_add(tot,s_scale(ei[a][b],s_mul(r[i][a][p],r[j][b][q])))
        return tot
    ok=True
    for i in (1,2):
        for j in (1,2):
            for n in range(0,M+1):
                tot=(10**6,[])
                for p in range(n+1):
                    tot=s_add(tot,s_scale((-1)**(n-p),pair(i,j,p,n-p)))
                if n==0 and i==j: tot=s_add(tot,s_neg(rows.Delta[i]))
                if not s_iszero(tot): ok=False; pr("   symplecticity fails",i,j,n,tot[0])
    print_ok("<v_i(w),v_j(-w)> = delta_ij Delta_i to w^%d"%M,ok)
    ok=True
    for c in (1,2):
        dlog=s_mul(s_der(rows.Delta[c]),s_inv(s_scale(2,rows.Delta[c])))
        for m in range(0,M):
            lhs=(10**6,[])
            if m>=1:
                lhs=s_add(s_der(r[c][1][m-1]),s_neg(s_mul(r[c][1][m-1],dlog)))
            rhs=s_add(s_add(s_scale(-a1*a2,s_mul(taui,r[c][0][m])),s_mul(s_add(s_scale(E,taui),s_const(-a0,N)),r[c][1][m])),s_neg(s_mul(rows.L[c],r[c][1][m])))
            if not s_iszero(s_add(lhs,s_neg(rhs))): ok=False; pr("   QDE row2 fails",c,m)
    print_ok("QDE row 2 with A=[[0,1],[-a1a2/tau,E/tau-a0]] at both colours, z^0..z^%d"%(M-1),ok)
    vals={}
    for c in (1,2):
        vals[c]=min(s[0] for s in fr.dil[c].values())
    pr("   min valuation of Dil_1, Dil_2:",vals[1],vals[2])
    pr("   min valuation of E_11,E_12,E_21,E_22:",[min(s[0] for s in fr.ed[ij].values()) for ij in ((1,1),(1,2),(2,1),(2,2))],"; E_21(b,a)=E_12(a,b):",all(s_iszero(s_add(fr.ed[(2,1)][(n,m)],s_neg(fr.ed[(1,2)][(m,n)]))) for (m,n) in fr.ed[(1,2)]))
    ok=all(s_coef(r[1][0][m],0)==pow(-1,m,P)*pow(inv(E),m,P)%P for m in range(M+1))
    print_ok("F2: r'_{0,1}(z;0)=E/(E+z) to z^%d"%M,ok)
    r11=[s_mul(rho,r[1][1][m]) for m in range(M+1)]
    print_ok("rho r_{1,1}: t^-1 coefficient = E delta_{m0}",all((s_coef(r11[m],-1)==(E if m==0 else 0)) for m in range(M+1)))
    X=s_coef(r11[0],0)*E%P
    ok=all(s_coef(r11[m],0)==X*pow(-1,m,P)%P*pow(inv(E),m+1,P)%P for m in range(M+1))
    pr("   F2: [t^0] rho r_{1,1}(z) = X/(E+z) to z^%d: %s; X = %s; a1 a2 = %s, a1 a2 + a0 E = %s"%(M,ok,fmt(X),fmt(a1*a2%P),fmt((a1*a2+a0*E)%P)))
    G=[s_mul(rho,r[1][0][m]) for m in range(M+1)]
    for k in range(0,K+1):
        ser={m+1:s_coef(G[m],k) for m in range(M+1)}
        ser[0]=0
        if k==0:
            mat=[[pow(j,m,P)*pow(-1,m,P)%P*pow(inv(E),m+1,P)%P for j in (1,)]+[1 if m==0 else 0] for m in range(M+2)]
            sol=solve(mat,[ser[m] for m in range(M+2)])
            pr("   L1 k=0: fit",None if sol is None else (fmt(sol[0][0]),fmt(sol[0][1])),"expect -E^2, E:",fmt(-E*E%P),fmt(E))
            continue
        mat=[[pow(j,m,P)*pow(-1,m,P)%P*pow(inv(E),m+1,P)%P for j in range(1,k+2)] for m in range(M+2)]
        sol=solve(mat,[ser[m] for m in range(M+2)])
        if sol is None: pr("   L1 k=%d: NOT in Lam_{k+1}"%k); continue
        mat2=[row[:-1] for row in mat]
        sol2=solve(mat2,[ser[m] for m in range(M+2)])
        pr("   L1 k=%d: in Lam_%d (rank %d, %d equations): True ; top coefficient non-zero: %s ; in Lam_%d: %s"%(k,k+1,sol[1],M+2,sol[0][k]%P!=0,k,sol2 is not None))
        if k<=2: pr("      a_j^(%d) ="%k,[fmt(x) for x in sol[0]],"; u/2 =",fmt(u*inv(2)%P),"; sum =",fmt(sum(sol[0])%P))
    r2=s_mul(rho,rho)
    for k in range(0,min(K,5)+1):
        idx=[(i,l) for i in range(1,k+2) for l in range(1,k+2) if (i-1)+(l-1)<=k]
        mat=[];b=[]
        for (m,n),s in fr.ed[(1,1)].items():
            x=s_mul(r2,s)
            b.append(s_coef(x,k))
            mat.append([pow(i,m,P)*pow(l,n,P)%P*pow(inv(E),m+n+2,P)%P for i,l in idx])
        sol=solve(mat,b)
        if sol is None: pr("   L3 k=%d: FAILS"%k); continue
        idx2=[j for j,(i,l) in enumerate(idx) if (i-1)+(l-1)<k]
        sol2=solve([[row[j] for j in idx2] for row in mat],b) if idx2 else None
        pr("   L3 k=%d: in span (%d unknowns, %d equations, rank %d): True; smaller span (no antidiagonal): %s"%(k,len(idx),len(mat),sol[1],sol2 is not None))
    for j in (1,2):
        res=[]
        for b_ in range(0,4):
            for k in range(0,min(K,5)+1):
                mat=[];rhs=[]
                for m in range(0,M-b_):
                    x=s_mul(rho,fr.ed[(1,j)][(m,b_)])
                    rhs.append(s_coef(x,k))
                    mat.append([pow(i,m,P)*pow(inv(E),m+1,P)%P for i in range(1,k+2)])
                sol=solve(mat,rhs)
                res.append(sol is not None)
        pr("   L2 j=%d, b=0..3, k=0..%d, in Lam_{k+1}:"%(j,min(K,5)),all(res),"(%d cases)"%len(res))
    ok=True
    Dg=[]
    for m in range(M+1):
        tg=s_mul((1,[1]+[0]*(N-1)),s_scale((-1)**m,G[m]))
        Dg.append(s_der(tg))
    for (m,n),s in fr.ed[(1,1)].items():
        tt=(2,[1]+[0]*(N-1))
        lhs=s_der(s_mul(tt,s_mul(r2,s)))
        rhs=s_mul(Th,s_mul(Dg[m],Dg[n]))
        if not s_iszero(s_add(lhs,s_neg(rhs))): ok=False
    print_ok("E11 flow d/dt[t^2 rho^2 E11]=Theta Dg Dg' for m+n<=%d"%(M-1),ok)
    r02=[s_coef(r[2][0][m],0) for m in range(M+1)]
    def Tseries(lst):
        ex=[0]*(M+1);k=1
        while 2*k-1<=M:
            c=fr2p(bernoulli(2*k)/(2*k*(2*k-1)))
            tot=0
            for (v,sg) in lst: tot+=sg*pow(inv(v),2*k-1,P)
            ex[2*k-1]=c*tot%P
            k+=1
        return series_exp(ex,M)
    cand={}
    cand['T(E)^-2 T(a1) T(a0+a2) T(a2) T(a0+a1), T(v)=Exp[Sum b_k (z/v)^(2k-1)]']=Tseries([(a1,1),(B,1),(a2,1),((a0+a1)%P,1),(E,-2)])
    cand['the inverse of that']=Tseries([(a1,-1),(B,-1),(a2,-1),((a0+a1)%P,-1),(E,2)])
    for n_,v in cand.items():
        pr("   r'_{0,2}(z;0) =",n_,":",v==r02)
    closed_form_test(W,fr,rho,G,E,A,B,C,a0,a1,a2)
def closed_form_test(W,fr,rho,G,E,A,B,C,a0,a1,a2):
    NT=K+3
    s=[B*inv(E)%P]+[0]*(NT-1)
    def pmul(x,y):
        out=[0]*NT
        for i,a in enumerate(x):
            if a:
                for j in range(NT-i): out[i+j]+=a*y[j]
        return [v%P for v in out]
    tau=[0,1]+[0]*(NT-2)
    for it in range(NT):
        s2=pmul(s,s)
        ts2=pmul(tau,s2); ts=pmul(tau,s)
        s=[((B if i==0 else 0)+(a0+a1)*ts2[i]-a0*ts[i])%P*inv(E)%P for i in range(NT)]
    den=[(-x)%P for x in pmul(tau,s)]; den[0]=(den[0]+1)%P
    di=[0]*NT; di[0]=inv(den[0])
    for k in range(1,NT):
        di[k]=(-sum(den[j]*di[k-j] for j in range(1,k+1)))%P*di[0]%P
    dphip=[C*x%P for x in pmul(s,di)]
    dphi=[0]+[dphip[i]*inv(i+1)%P for i in range(NT-1)]
    LO=NT+2; HI=M+2*NT+3
    size=LO+HI+1
    def lz_zero(): return [0]*size
    def lz_mul(x,y):
        out=[0]*size
        for i,a in enumerate(x):
            if a:
                pa=i-LO
                for j,b in enumerate(y):
                    if b:
                        pw=pa+j-LO
                        if -LO<=pw<=HI: out[pw+LO]+=a*b
        return [v%P for v in out]
    def lz_const(c):
        o=lz_zero(); o[LO]=c%P; return o
    def lz_lin_inv(c0,c1):
        o=lz_zero(); ic=inv(c0)
        cur=ic
        for p in range(0,HI+1):
            o[LO+p]=cur
            cur=cur*(-c1)%P*ic%P
        return o
    def lz_lin(c0,c1):
        o=lz_zero(); o[LO]=c0%P; o[LO+1]=c1%P; return o
    def lz_zpow(k):
        o=lz_zero(); o[LO+k]=1; return o
    Xs=[lz_const(1)]+[None]*(NT-1)
    zi=lz_zpow(-1)
    for m in range(1,NT):
        tot=lz_zero()
        for l in range(1,m+1):
            if dphi[l]:
                t=lz_mul(lz_mul(lz_const(-l*dphi[l]),zi),Xs[m-l])
                tot=[(a+b)%P for a,b in zip(tot,t)]
        Xs[m]=[v*inv(m)%P for v in tot]
    Fs=[]
    for n in range(NT):
        cur=lz_const(inv(math.factorial(n)))
        for j in range(1,n+1):
            cur=lz_mul(cur,lz_lin(C,j)); cur=lz_mul(cur,lz_lin(B,j)); cur=lz_mul(cur,zi)
        for j in range(2,n+2):
            cur=lz_mul(cur,lz_lin_inv(E,j))
        Fs.append(cur)
    Hs=[]
    for n in range(NT):
        tot=lz_zero()
        for n1 in range(n+1):
            t=lz_mul(Fs[n1],Xs[n-n1])
            tot=[(a+b)%P for a,b in zip(tot,t)]
        Hs.append(tot)
    regular=all(all(H[LO+p]==0 for p in range(-LO,0)) for H in Hs)
    pr("   H_n=[tau^n]F Ephi regular at z=0, n<%d: %s"%(NT,regular))
    taut=[0]+[(-(-1)**i*inv(math.factorial(i)))%P for i in range(1,NT+3)]
    def tmul(x,y):
        out=[0]*(NT+3)
        for i,a in enumerate(x):
            if a:
                for j in range(NT+3-i): out[i+j]+=a*y[j]
        return [v%P for v in out]
    tauot=taut[1:]+[0]
    powers=[[1]+[0]*(NT+2)]
    for n in range(1,NT): powers.append(tmul(powers[-1],taut))
    pref=lz_mul(lz_const(E),lz_lin_inv(E,1))
    ok=True; checked=0
    for k in range(0,K+1):
        tot=lz_zero()
        for n in range(0,min(k,NT-1)+1):
            d=tmul(tauot,powers[n])[k]
            tot=[(a+d*b)%P for a,b in zip(tot,Hs[n])]
        tot=lz_mul(pref,tot)
        for m in range(0,M+1):
            lhs=s_coef(G[m],k)
            if lhs!=tot[LO+m]: ok=False
            checked+=1
    pr("   KEY FORMULA rho r'_{0,1}(z;t) = (E/(E+z))(tau/t) 2F1(a2/z+1,(a0+a2)/z+1;E/z+2;tau) Exp[-dphi/z]: %s (%d coefficients, t^0..t^%d, z^0..z^%d)"%(ok,checked,K,M))
random.seed(7)
vals=[random.randrange(2,P) for _ in range(4)]
Wr=dict(a1=vals[0],a2=vals[1],e1=vals[2],e2=vals[3]); Wr['a0']=(Wr['e1']+Wr['e2']-Wr['a1']-Wr['a2'])%P
for nm,W in (('fibrewise (3,5,3,5)',cyw(3,5,3,5)),('S1 (3,4,5,11)',cyw(3,4,5,11)),('S2 (7,13,2,29)',cyw(7,13,2,29)),('S3 (17,-5,19,-8)',cyw(17,-5,19,-8)),('random F_p',Wr)):
    check_sample(nm,W)
