# Independent referee recomputation for H2 (round 3). stdlib fractions only.
# usage: py ref.py e1 e2 NV NZ [Q] [control]
#  Q: off-shell section c1 = C1cf*(1+Q*Del); control in {none, r2, eta, sign}
import sys, math
from fractions import Fraction as F
e1,e2=F(sys.argv[1]),F(sys.argv[2]); E=e1+e2
NV=int(sys.argv[3]); NZ=int(sys.argv[4]); Q=F(sys.argv[5]) if len(sys.argv)>5 else F(0)
CTRL=sys.argv[6] if len(sys.argv)>6 else 'none'
print("e1,e2,E=",e1,e2,E,"NV,NZ=",NV,NZ,"Q=",Q,"ctrl=",CTRL)
def Z0(): return [F(0)]*(NZ+1)
CAP=NV+2
def S(d): return {k:v for k,v in d.items() if k<CAP and any(v)}
def zm(a,b):
    r=Z0()
    for i,x in enumerate(a):
        if x:
            for j in range(NZ+1-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def add(*fs):
    r={}
    for f in fs:
        for k,v in f.items():
            r[k]=[a+b for a,b in zip(r.get(k,Z0()),v)]
    return S(r)
def sc(c,f): return S({k:[c*x for x in v] for k,v in f.items()})
def mul(f,g):
    r={}
    for a,u in f.items():
        for b,w in g.items():
            if a+b<CAP:
                p=zm(u,w); r[a+b]=[x+y for x,y in zip(r.get(a+b,Z0()),p)]
    return S(r)
def cst(c,e=0):
    z=Z0(); z[0]=F(c); return {e:z}
def zf(cs,val=0):
    return S({val+i:[F(c)]+[F(0)]*NZ for i,c in enumerate(cs) if c})
def d(f): return S({k-1:[k*x for x in v] for k,v in f.items() if k})
def zz(f): return S({k:[F(0)]+v[:NZ] for k,v in f.items()})
def flip(f): return S({k:[x if n%2==0 else -x for n,x in enumerate(v)] for k,v in f.items()})
def inv(f):
    m=min(k for k in f); a=[f.get(m+i,Z0())[0] for i in range(CAP-m+3)]
    b=[1/a[0]]
    for n in range(1,CAP+3-m):
        b.append(-sum(a[j]*b[n-j] for j in range(1,min(n,len(a)-1)+1))/a[0])
    return zf(b,-m)
# ---- r_l: coefficient extraction of e^{z d1 d2/3}[e^{s1^3/z}(1-u)^{-1/2}] at s1=s2=0 ----
L=NV+3; DD=3*L+3
def pmul(A,B):
    R={}
    for (a,b),x in A.items():
        for (c,dd),y in B.items():
            if a+c<=DD and b+dd<=DD: R[(a+c,b+dd)]=R.get((a+c,b+dd),0)+x*y
    return R
u={(1,1):F(18,27),(3,0):F(-4,27),(0,3):F(-4,27),(2,2):F(1,27)}
G={(0,0):F(1)}; Pk={(0,0):F(1)}; c=F(1)
for k in range(1,2*DD+2):
    Pk=pmul(Pk,u); c=c*(2*k-1)/(2*k)
    if not Pk: break
    for key,x in Pk.items(): G[key]=G.get(key,0)+c*x
r=[F(0)]*(L+1)
for m in range(DD+1):
    for j in range(0,m//3+1):
        l=m-j
        if l<=L: r[l]+=F(math.factorial(m),3**m*math.factorial(j))*G.get((m-3*j,m),0)
print("r_l",[str(x) for x in r[:6]])
if CTRL=='r2': r[2]+=F(1,100)
def iprod(n):
    c=Z0(); c[0]=F(1)
    for j in range(1,n+1): c=zm(c,[(-F(j)/E)**k for k in range(NZ+1)])
    return c
P0={}
for l in range(0,CAP-1):
    zl=Z0(); ip=iprod(l+1)
    for k in range(NZ+1-l): zl[k+l]=ip[k]*r[l]
    P0[l+1]=zl
P0=S(P0)
thv=zf([1/E,F(-1,3)])
def th(f): return mul(thv,d(f))
Y=zf([1,-E/3],-1)
Dl=[F(0)]+[-F((n+1)*(n+2),2)*(E/3)**n for n in range(1,CAP+3)]
Del=zf(Dl)
assert all(k>=CAP-1 for k in add(th(Del),sc(-1,Del),cst(1)))
N2=CAP+5
def xd(a): return [(i+1)*a[i+1] for i in range(len(a)-1)]+[F(0)]
def xm(a,b):
    r_=[F(0)]*len(a)
    for i,x in enumerate(a):
        for j in range(len(a)-i): r_[i+j]+=x*b[j]
    return r_
xm1=[F(-1),F(1)]+[F(0)]*(N2-2)
def tx(a): return xm(xm1,xd(a))
def res(w):
    t1=tx(w); t2=tx(t1)
    inner=[9*t2[k]+9*t1[k]+2*w[k] for k in range(N2)]
    o=xm(xm1,inner)
    return [t2[i]+o[i]/9 for i in range(N2)]
w=[F(1)]+[F(0)]*(N2-1)
for n in range(1,N2-1):
    r0=res(w)[n-1]; w2=list(w); w2[n]=F(1); r1=res(w2)[n-1]-r0
    w[n]=-r0/r1
assert all(x==0 for x in res(w)[:N2-2])
gg=[sum(w[:n+1]) for n in range(N2)]
tcx=[F(0)]+[gg[n-1]/n for n in range(1,N2)]
print("wc:",[str(x) for x in w[:4]],"t_c(Del):",[str(x) for x in tcx[:4]])
def comp(cs,X):
    acc={}; P=cst(1)
    for cc in cs:
        if cc: acc=add(acc,sc(cc,P))
        P=mul(P,X)
    return acc
wc=comp(w[:CAP+2],Del); tcv=comp(tcx[:CAP+2],Del)
C1cf=sc(-1,wc)
assert all(k>=CAP-1 for k in add(th(tcv),sc(-1,C1cf)))
c1=mul(C1cf,add(cst(1),sc(Q,Del)))
C1i=inv(c1)
d2=mul(sc(E/2,add(mul(C1i,C1i),cst(-1))),inv(Del))
def YzT(f): return add(mul(Y,f),zz(th(f)))
P1=mul(C1i,YzT(P0))
P2=mul(mul(Del,mul(c1,c1)),add(YzT(P1),sc(-1,mul(d2,P1))))
a=[tcv.get(k,Z0())[0] for k in range(0,CAP+2)]
Vc=[F(0),1/a[1]]
for n in range(2,CAP+1):
    Vs=zf(Vc+[F(0)]); acc={}; P=cst(1)
    for k in range(1,n+1):
        P=mul(P,Vs); acc=add(acc,sc(a[k],P))
    Vc.append(-acc.get(n,Z0())[0]/a[1])
V=zf(Vc)
print("v(t_c) coeffs:",[str(x) for x in Vc[1:4]])
Vinv=inv(V)
Vp={0:cst(1)}
for k in range(1,CAP+1): Vp[k]=mul(Vp[k-1],V)
for k in range(1,3): Vp[-k]=mul(Vp[-k+1],Vinv)
def tov(f):
    acc={}
    for k,zl in f.items(): acc=add(acc,mul({0:zl},Vp[k]))
    return acc
tci=zf([1],-1)
sig=[sc(-E,mul(tci,tov(P))) for P in (P0,P1,P2)]
if CTRL=='sign': sig[2]=sc(-1,sig[2])
KT=NV-2
def fit(cz,Jmax):
    for J in range(0,Jmax+1):
        Pz=Z0(); Pz[0]=F(1)
        for j in range(1,J+1): Pz=zm(Pz,[E,F(j)]+[F(0)]*(NZ-1))
        Nn=zm(cz,Pz)
        if any(Nn[J+2:]): continue
        num=Nn[:J+2]
        A1=num[J+1]/Pz[J]
        rem=[num[i]-A1*(Pz[i-1] if i else 0) for i in range(J+2)]
        A0=rem[J]/Pz[J]; rem=[rem[i]-A0*Pz[i] for i in range(J+2)]
        beta=[]
        for j in range(1,J+1):
            z0=-E/j; val=sum(rem[i]*z0**i for i in range(J+2)); den=F(1)
            for i in range(1,J+1):
                if i!=j: den*=(E+i*z0)
            beta.append(val/den)
        return dict(J=J,psi=-A1,c=A0,beta=beta,spare=NZ-J-1)
    return None
viol=[]; fits={}
for A in range(3):
    for k in range(-2,KT):
        cz=sig[A].get(k,Z0()); f=fit(cz,NZ-4); fits[(A,k)]=f
        if f is None: viol.append((A,k,'nofit')); continue
        top=max([j+1 for j,b in enumerate(f['beta']) if b],default=0)
        bnd=(k+1) if A<2 else k
        if f['psi']!=0: viol.append((A,k,'psi',str(f['psi'])))
        if k>=0 and f['c']!=0: viol.append((A,k,'const',str(f['c'])))
        if top>max(bnd,0): viol.append((A,k,'top',top))
        if k<=2: print(f"sig_{A}^({k}): psi={f['psi']} c={f['c']} beta={[str(b) for b in f['beta']]} spare={f['spare']}")
print("violations:",viol[:12])
print("Zoo sums:",[str(sum(b/(j+1) for j,b in enumerate(fits[(0,k)]['beta']))) if fits[(0,k)] else None for k in range(0,KT)])
et00=e1*e2*E**3/9; et11=-3*E*e1*e2; et12=-3*e1*e2
if CTRL=='eta': et12=-et12
sm=[flip(x) for x in sig]
K=add(sc(et00,mul(sig[0],sm[0])),sc(et11,mul(sig[1],sm[1])),sc(et12,add(mul(sig[1],sm[2]),mul(sig[2],sm[1]))))
Kd=add(K,cst(-3*e1*e2*E**3,-2))
print("K=1 residual nonzero orders (<KT-1):",[k for k in sorted(Kd) if k<KT-1])
L1=tov(mul(C1i,Y)); kap=tov(mul(mul(mul(c1,c1),c1),Del)); mu=tov(sc(E/2,add(mul(c1,c1),cst(-1))))
def Dc(f): return add(mul(L1,f),zz(mul(tci,f)),zz(d(f)))
I1=add(sig[1],sc(-1,Dc(sig[0]))); I2=add(sig[2],sc(-1,add(mul(kap,Dc(sig[1])),mul(mu,sig[1]))))
print("I1 residual orders:",[k for k in sorted(I1) if k<KT-1],"I2 residual orders:",[k for k in sorted(I2) if k<KT-2])
print("L1:",[str(L1.get(k,Z0())[0]) for k in range(-1,3)],"kappa:",[str(kap.get(k,Z0())[0]) for k in range(0,3)],"mu:",[str(mu.get(k,Z0())[0]) for k in range(0,3)])
Yt=[mul(zf([1],1),s) for s in sig]
NM=NZ//2
def coefpsi(f,k):
    zl=f.get(k,Z0()); return [zl[m]*(-1)**m for m in range(NZ+1)]
eta={(0,0):et00,(1,1):et11,(1,2):et12,(2,1):et12}
for n in range(0,min(KT,8)):
    Nmat=[[F(0)]*(NM+2) for _ in range(NM+2)]
    for (aa,bb),et in eta.items():
        for p in range(0,n+1):
            A_=coefpsi(Yt[aa],p); B_=coefpsi(Yt[bb],n-p)
            for m in range(NM+2):
                for q in range(NM+2): Nmat[m][q]-=et*A_[m]*B_[q]
    if n==0: Nmat[0][0]+=3*e1*e2*E**3
    Qm=[[sum((-1)**k*Nmat[m+1+k][q-k] for k in range(q+1)) for q in range(NM-m)] for m in range(NM)]
    diag=[sum((-1)**q*Nmat[m-q][q] for q in range(m+1)) for m in range(NM+1)]
    Cil={}
    ok=None
    if all(fits.get((x,k)) for x in range(3) for k in range(0,max(n-1,0))):
        for (aa,bb),et in eta.items():
            for k1 in range(0,n-1):
                k2=n-2-k1
                fa=fits[(aa,k1)]; fb=fits[(bb,k2)]
                for i,bi in enumerate(fa['beta']):
                    for l,bl in enumerate(fb['beta']):
                        Cil[(i+1,l+1)]=Cil.get((i+1,l+1),0)-et*bi*bl
        pred=[[F(0)]*NM for _ in range(NM)]
        for (i,l),C in Cil.items():
            A_=F(i*l)/(E*(i+l))*C
            for m in range(NM):
                for q in range(NM-m):
                    pred[m][q]+=A_*F(i)**m/E**(m+1)*F(l)**q/E**(q+1)
        ok=all(Qm[m][q]==pred[m][q] for m in range(NM) for q in range(NM-m))
    zero=all(Qm[m][q]==0 for m in range(NM) for q in range(NM-m))
    sup=max([(i-1)+(l-1) for (i,l),C in Cil.items() if C],default=None)
    print(f"Ebar t_c^{n}: N(psi,-psi)=0:{all(x==0 for x in diag)} residue formula:{ok} Ebar zero:{zero} max(i-1)+(l-1)={sup} C11/(2E)={Cil.get((1,1),0)/(2*E)}")
print("target -e1e2E^4/27 =",-e1*e2*E**4/27)
pr=mul(tci,tov(sc(E*Q,mul(mul(Del,mul(C1cf,C1cf)),add(Del,cst(-1))))))
print("(d) psi coeff of sigma_2 vs prediction:",[(str(fits[(2,k)]['psi']) if fits[(2,k)] else None, str(pr.get(k,Z0())[0])) for k in range(-1,min(KT,6))])
# appended: Lemma 1 closed forms check
def Sm(mm):
    out={}
    for l in range(0,CAP-1):
        zl=Z0(); ip=iprod(max(l+1-mm,0))
        for k in range(NZ+1-l): zl[k+l]=ip[k]*r[l]
        out[l]=zl
    return S(out)
S1=Sm(1); S2_=Sm(1); iv=zf([1],-1)
Ph1=mul(iv,P1); Ph2=mul(iv,P2)
cf1=sc(E,mul(mul(C1i,mul(thv,iv)),S1))
A=mul(C1i,thv)
tl=mul(th(A),inv(A))  # theta log(C1i theta v)
inner=add(sc(E,mul(thv,mul(mul(iv,iv),S2_))), mul(add(zz(tl),sc(F(-1,2),d2)),mul(S1,iv)))
cf2=sc(E,mul(mul(Del,mul(c1,thv)),inner))
D1=add(Ph1,sc(-1,cf1)); D2=add(Ph2,sc(-1,cf2))
print("Lemma1 row1 residual v-orders:",[k for k in sorted(D1) if k<CAP-3],"row2:",[k for k in sorted(D2) if k<CAP-4])
