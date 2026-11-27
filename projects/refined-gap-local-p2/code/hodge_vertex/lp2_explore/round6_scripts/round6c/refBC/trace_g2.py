# Evidence only. Fresh code, arithmetic mod p = 2^61-1 (p = 1 mod 3), truncated in a parameter x (x^K = 0).
# Mode diag: alpha = E/3 + x*d, fixed y = t^3. Mode yadic: generic alpha, y = x.
# All rows come from the Wick expansion Phi(alpha,Y;z) at the root Y_j, never from the closed form.
import sys, itertools, math
from fractions import Fraction as Fr
p=2**61-1
def fr(a):
    a=Fr(a); return a.numerator%p*pow(a.denominator,-1,p)%p
g=3
omega=pow(g,(p-1)//3,p)
while omega==1:
    g+=1; omega=pow(g,(p-1)//3,p)
K=2
def sc(a): return [fr(a)]+[0]*(K-1)
def szero(): return [0]*K
def sadd(a,b): return [(x+y)%p for x,y in zip(a,b)]
def ssub(a,b): return [(x-y)%p for x,y in zip(a,b)]
def smul(a,b):
    r=[0]*K
    for i,x in enumerate(a):
        if x:
            for j in range(K-i): r[i+j]=(r[i+j]+x*b[j])%p
    return r
def sscal(c,a): return [c*x%p for x in a]
def sinv(a):
    i0=pow(a[0],-1,p); r=[0]*K; r[0]=i0
    for n in range(1,K):
        s=0
        for j in range(1,n+1): s=(s+a[j]*r[n-j])%p
        r[n]=(-s*i0)%p
    return r
DN=3
def dc(s): return [s]+[szero() for _ in range(DN-1)]
def dadd(a,b): return [sadd(x,y) for x,y in zip(a,b)]
def dsub(a,b): return [ssub(x,y) for x,y in zip(a,b)]
def dmul(a,b):
    r=[szero() for _ in range(DN)]
    for i in range(DN):
        for j in range(DN-i):
            r[i+j]=sadd(r[i+j],smul(a[i],b[j]))
    return r
def dscal(c,a): return [sscal(c,x) for x in a]
def dinv(a):
    i0=sinv(a[0]); r=[szero() for _ in range(DN)]; r[0]=i0
    for n in range(1,DN):
        s=szero()
        for j in range(1,n+1): s=sadd(s,smul(a[j],r[n-j]))
        r[n]=sscal(p-1,smul(s,i0))
    return r
def ddelta(c):
    return [sscal(i+1,c[i+1]) if i+1<DN else szero() for i in range(DN)]
BERN={2:Fr(1,6),4:Fr(-1,30),6:Fr(1,42),8:Fr(-1,30)}
def theta_ser(c,N):
    l=[Fr(0)]*(N+1)
    for m in range(1,5):
        if 2*m-1<=N: l[2*m-1]=BERN[2*m]/(2*m*(2*m-1))*c**(2*m-1)
    r=[Fr(0)]*(N+1); r[0]=Fr(1)
    for n in range(1,N+1): r[n]=sum(k*l[k]*r[n-k] for k in range(1,n+1))/n
    return r
def phi_rows(Yv,al,N):
    ZERO=dc(szero()); ONE=dc(sc(1))
    Y=[Yv,sc(1),szero()]
    q=[dadd(Y,dc(a)) for a in al]
    third=fr(Fr(1,3))
    Y3=dscal(third,Y)
    H=[[dsub(q[i],Y3) if i==j else dsub(ZERO,Y3) for j in range(3)] for i in range(3)]
    def m2(a,b,c,d): return dsub(dmul(a,d),dmul(b,c))
    adj=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r=[x for x in range(3) if x!=i]; c=[x for x in range(3) if x!=j]
            mi=m2(H[r[0]][c[0]],H[r[0]][c[1]],H[r[1]][c[0]],H[r[1]][c[1]])
            adj[j][i]=dscal(p-1,mi) if (i+j)%2 else mi
    det=dadd(dadd(dmul(H[0][0],adj[0][0]),dmul(H[0][1],adj[1][0])),dmul(H[0][2],adj[2][0]))
    idet=dinv(det)
    G=[[dscal(p-1,dmul(adj[i][j],idet)) for j in range(3)] for i in range(3)]
    Mmax=2*N+2
    T={}
    for m in range(3,Mmax+1):
        tm={}
        fm=fr(Fr(1,math.factorial(m)))
        for j in range(3):
            e=[0,0,0]; e[j]=m
            tm[tuple(e)]=dscal(fm,q[j])
        c3=fr(Fr(-3,3**m))
        for a in range(m+1):
            for b in range(m+1-a):
                e=(a,b,m-a-b)
                co=fr(Fr(1,math.factorial(a)*math.factorial(b)*math.factorial(m-a-b)))
                v=dscal(c3*co%p,Y)
                tm[e]=dadd(tm[e],v) if e in tm else v
        T[m]=tm
    def pmul(f,gg):
        r={}
        for e1,x in f.items():
            for e2,y in gg.items():
                e=(e1[0]+e2[0],e1[1]+e2[1],e1[2]+e2[2])
                v=dmul(x,y)
                r[e]=dadd(r[e],v) if e in r else v
        return r
    def lap(f):
        r={}
        for e,x in f.items():
            for a in range(3):
                for b in range(3):
                    ee=list(e)
                    if a==b:
                        if ee[a]<2: continue
                        co=ee[a]*(ee[a]-1); ee[a]-=2
                    else:
                        if ee[a]<1 or ee[b]<1: continue
                        co=ee[a]*ee[b]; ee[a]-=1; ee[b]-=1
                    ee=tuple(ee)
                    v=dmul(dscal(co%p,x),G[a][b])
                    r[ee]=dadd(r[ee],v) if ee in r else v
        return r
    A={(0,0):{(0,0,0):ONE}}
    Phi=[ZERO for _ in range(N+1)]
    Phi[0]=ONE
    for k in range(1,2*N+1):
        ik=fr(Fr(1,k))
        for D_ in range(3*k,2*(N+k)+1):
            tot={}
            for m in range(3,Mmax+1):
                prev=A.get((k-1,D_-m))
                if prev is None: continue
                pr=pmul(T[m],prev)
                for e,x in pr.items():
                    tot[e]=dadd(tot[e],x) if e in tot else x
            if tot:
                A[(k,D_)]={e:dscal(ik,x) for e,x in tot.items()}
        for n in range(1,N+1):
            D_=2*(n+k)
            if (k,D_) in A:
                L=n+k
                f=A[(k,D_)]
                for _ in range(L): f=lap(f)
                val=f.get((0,0,0),ZERO)
                co=fr(Fr(1,2**L*math.factorial(L)))
                Phi[n]=dadd(Phi[n],dscal(co,val))
    return Phi,Y,q,det
def build(e1_,e2_,alpha,roots,N=3):
    ZERO=dc(szero())
    ee=fr(e1_*e2_)
    th1=theta_ser(1/Fr(e1_),N); th2=theta_ser(1/Fr(e2_),N)
    tt=[fr(sum(th1[i]*th2[n-i] for i in range(n+1))) for n in range(N+1)]
    cols=[]
    for Yv in roots:
        Phi,Y,q,det3=phi_rows(Yv,alpha,N)
        P=[ZERO for _ in range(N+1)]
        for n in range(N+1):
            for i in range(n+1):
                P[n]=dadd(P[n],dscal(tt[i],Phi[n-i]))
        W=dscal(3,det3)
        iW=dinv(W)
        prodq=dmul(dmul(q[0],q[1]),q[2])
        gg=dscal(p-1,dmul(dmul(Y,prodq),iW))
        tau=dmul(ddelta(W),dmul(gg,iW))
        half=fr(Fr(1,2))
        def theta(f): return dmul(ddelta(f),gg)
        def Top(f):
            out=[ZERO for _ in range(N+1)]
            for n in range(N+1):
                out[n]=dmul(Y,f[n])
                if n>=1:
                    out[n]=dadd(out[n],theta(f[n-1]))
                    out[n]=dsub(out[n],dscal(half,dmul(tau,f[n-1])))
            return out
        Q=Top(P); TQ=Top(Q)
        Delta=sscal(3*ee%p,W[0])
        cols.append(dict(P=[x[0] for x in P],Q=[x[0] for x in Q],TQ=[x[0] for x in TQ],D=Delta,Y=Yv))
    return cols
def wk(ds,g_):
    ds=sorted(ds)
    n=len(ds)
    if sum(ds)!=3*g_-3+n: return Fr(0)
    if g_==0:
        r=Fr(math.factorial(n-3))
        for d in ds: r/=math.factorial(d)
        return r
    if ds and ds[0]==0:
        rest=ds[1:]; tot=Fr(0)
        for i in range(len(rest)):
            if rest[i]>0:
                rr=rest[:]; rr[i]-=1; tot+=wk(rr,g_)
        return tot
    if 1 in ds:
        rr=ds[:]; rr.remove(1)
        if g_==1 and not rr: return Fr(1,24)
        return (2*g_-2+len(rr))*wk(rr,g_)
    base={(4,):Fr(1,1152),(2,3):Fr(29,5760),(2,2,2):Fr(7,240)}
    return base[tuple(ds)]
def total(cols,al,y,ee_,Fv):
    e1a=sadd(sadd(al[0],al[1]),al[2])
    e3a=smul(smul(al[0],al[1]),al[2])
    ee=fr(ee_)
    Dc=sadd(sc(1),sscal(27,y))
    eta00=sscal(3*ee%p,e3a)
    def pr(i,j,M,b):
        a=cols[i];c=cols[j]
        r=smul(eta00,smul(a['P'][M],c['P'][b]))
        r=ssub(r,sscal(3*ee%p,smul(e1a,smul(a['Q'][M],c['Q'][b]))))
        r=ssub(r,sscal(3*ee%p,smul(Dc,sadd(smul(a['Q'][M],c['TQ'][b]),smul(a['TQ'][M],c['Q'][b])))))
        s=szero()
        if M>=1: s=sadd(s,smul(a['Q'][M-1],c['Q'][b]))
        if b>=1: s=sadd(s,smul(a['Q'][M],c['Q'][b-1]))
        r=sadd(r,sscal(3*ee%p,smul(Dc,smul(Fv,s))))
        return r
    def Ed(i,j,a,b):
        r=szero()
        for k in range(b+1):
            sg=(-1)**k*(-1)**(a+b+1)
            r=ssub(r,sscal(sg%p,pr(i,j,a+k+1,b-k)))
        return r
    def Dil(i,d): return sscal(((-1)**d)%p,cols[i]['P'][d-1])
    def V(i,g_,exps):
        n=len(exps); dim=3*g_-3+n
        tot=szero()
        for m in range(0,dim+1):
            for dd in itertools.product(range(2,dim+m+1),repeat=m):
                if sum(exps)+sum(dd)!=dim+m: continue
                c=wk(list(exps)+list(dd),g_)/math.factorial(m)
                if c==0: continue
                t=sc(c)
                for d in dd: t=smul(t,Dil(i,d))
                tot=sadd(tot,t)
        Dp=sc(1)
        if g_==2: Dp=cols[i]['D']
        if g_==0: Dp=sinv(cols[i]['D'])
        return smul(tot,Dp)
    C=range(3)
    h8=fr(Fr(1,8));h2=fr(Fr(1,2));h12=fr(Fr(1,12))
    t1=szero()
    for i in C: t1=sadd(t1,V(i,2,[]))
    t2=szero()
    for i in C:
        for j in C:
            for a in range(2):
                for b in range(2):
                    t2=sadd(t2,smul(smul(V(i,1,[a]),V(j,1,[b])),Ed(i,j,a,b)))
    t3=szero()
    for i in C:
        for a in range(3):
            for b in range(3-a):
                t3=sadd(t3,smul(V(i,1,[a,b]),Ed(i,i,a,b)))
    t4=szero()
    for i in C:
        for ex in itertools.product(range(2),repeat=4):
            if sum(ex)>1: continue
            t4=sadd(t4,smul(V(i,0,list(ex)),smul(Ed(i,i,ex[0],ex[1]),Ed(i,i,ex[2],ex[3]))))
    t5=szero()
    for i in C:
        for j in C:
            for x in range(2):
                t5=sadd(t5,smul(smul(V(i,1,[x]),V(j,0,[0,0,0])),smul(Ed(i,j,x,0),Ed(j,j,0,0))))
    t6=szero()
    t7=szero()
    for i in C:
        for j in C:
            vv=smul(V(i,0,[0,0,0]),V(j,0,[0,0,0]))
            t6=sadd(t6,smul(vv,smul(Ed(i,i,0,0),smul(Ed(i,j,0,0),Ed(j,j,0,0)))))
            e=Ed(i,j,0,0)
            t7=sadd(t7,smul(vv,smul(e,smul(e,e))))
    return sadd(t1,sadd(sscal(h2,t2),sadd(sscal(h2,t3),sadd(sscal(h8,t4),sadd(sscal(h2,t5),sadd(sscal(h8,t6),sscal(h12,t7)))))))
def newton(alpha,y,Y0):
    Y=Y0
    for _ in range(K+2):
        q=[sadd(Y,a) for a in alpha]
        L=sadd(smul(smul(q[0],q[1]),q[2]),sscal(27,smul(y,smul(Y,smul(Y,Y)))))
        dL=sadd(sadd(smul(q[1],q[2]),smul(q[0],q[2])),smul(q[0],q[1]))
        dL=sadd(dL,sscal(81,smul(y,smul(Y,Y))))
        Y=ssub(Y,smul(L,sinv(dL)))
    return Y
def F20(e1,e2,X,F):
    E=e1+e2; ee=e1*e2
    iX=pow(X,-1,p)
    f=fr
    F02=(f(Fr(5,8))*F*F*F*iX+(f(Fr(5,8))*iX-f(Fr(1,2)))*F*F+f(Fr(1,96))*(13*X+20*iX-32)*F-f(Fr(1,2160))*(27*X*X*X-98*X*X+121*X-50)*iX)%p
    F11=(f(Fr(1,16))*F*F+f(Fr(1,96))*(8-7*X)*F+f(Fr(1,5760))*(84*X*X-196*X+121))%p
    F20_=(f(Fr(1,384))*X*F-f(Fr(1,5760))*(21*X*X-29*X+8))%p
    return (fr(E**4/ee)*F20_+fr(E**2)*F11+fr(ee)*F02)%p
def mode_diag(e1,e2,t,d):
    global K
    K=2
    e1,e2=Fr(e1),Fr(e2);E=e1+e2
    tm=fr(t); y=[pow(tm,3,p)]+[0]*(K-1)
    alpha=[[fr(E/3),fr(di)] for di in d]
    roots=[]
    for j in range(3):
        tj=tm*pow(omega,j,p)%p
        Y0=(-fr(E/3)*pow((1+3*tj)%p,-1,p))%p
        roots.append(newton(alpha,y,[Y0]+[0]*(K-1)))
    cols=build(e1,e2,alpha,roots)
    res=[]
    for Ff in [Fr(0),Fr(1),Fr(-2,3),Fr(5,3),Fr(11,7)]:
        v=total(cols,alpha,y,e1*e2,sc(Ff))
        X=pow((1+27*y[0])%p,-1,p)
        res.append((Ff,v,F20(e1,e2,X,fr(Ff))))
    return res
def Fseries(K_):
    # F = theta log A, as a list of mod p values in y, length K_
    a=[]
    for d in range(K_):
        a.append(fr(Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d))
    # theta A = A F
    F=[0]*K_
    ia0=pow(a[0],-1,p)
    for n in range(K_):
        s=n*a[n]%p
        for i in range(1,n+1): s=(s-a[i]*F[n-i])%p
        F[n]=s*ia0%p
    return F
def mode_yadic(e1,e2,alpha_vals,K_=3):
    # generic alpha: y-adic, returns the total mod y^K_, and the mirror map data
    global K
    K=K_
    e1,e2=Fr(e1),Fr(e2)
    y=[0,1]+[0]*(K-2)
    alpha=[[fr(a)]+[0]*(K-1) for a in alpha_vals]
    roots=[newton(alpha,y,[(-fr(a))%p]+[0]*(K-1)) for a in alpha_vals]
    cols=build(e1,e2,alpha,roots)
    Fv=Fseries(K)
    return total(cols,alpha,y,e1*e2,Fv)
if __name__=='__main__':
    for (e1,e2,t) in [(3,7,Fr(1,5)),(1,1,Fr(-2,7)),(Fr(5,2),Fr(-3,7),Fr(3,11))]:
        R=mode_diag(e1,e2,t,[1,2,5])
        for Ff,v,bs in R:
            print('eps',e1,e2,'t',t,'F',Ff,'s^0 equals eq:F2ref:',v[0]==bs,'s^1 nonzero:',v[1]!=0,flush=True)
    # y-adic along a line through alpha*; polynomial interpolation in r to r = 0
    for (e1,e2,dval1) in [(3,7,Fr(48241,1680)),(1,1,Fr(81,80))]:
        E=Fr(e1)+Fr(e2)
        d=(Fr(1),Fr(2,3),Fr(-5,7))
        pts=[]
        for r in range(1,9):
            al=[E/3+r*x for x in d]
            v=mode_yadic(e1,e2,al,K_=2)
            pts.append((r,v))
            print('line',e1,e2,'r',r,'done',flush=True)
        # interpolate e3*eps1eps2*GW0 (degree <= 7) and eps1eps2 GW1 (degree <= 4) to r=0
        def lag(xy,x0):
            tot=0
            for i,(xi,yi) in enumerate(xy):
                c=yi
                for j,(xj,_) in enumerate(xy):
                    if i!=j: c=c*((x0-xj)%p)%p*pow((xi-xj)%p,-1,p)%p
                tot=(tot+c)%p
            return tot
        gw1=[(r,v[1]) for r,v in pts]
        print('y^1 interpolated at r=0 from 6 points equals refined d=1:',lag(gw1[:6],0)==fr(dval1),' from last 6:',lag(gw1[2:],0)==fr(dval1))
        e3=[(r,fr((E/3+r*d[0])*(E/3+r*d[1])*(E/3+r*d[2]))) for r,_ in pts]
        gw0=[(r,v[0]*e[1]%p) for (r,v),e in zip(pts,e3)]
        print('y^0: e3 GW0 interpolated at r=0 / e3(alpha*) equals E^2/640:',lag(gw0,0)*pow(fr((E/3)**3),-1,p)%p==fr(E**2/640))
