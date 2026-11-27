# Direct test of (*) from Schur data, Laurent-in-w structure, order eps^0, eps^1 groups, (E1), Bor dictionary.
from fractions import Fraction as Fr
import math, sys, time
from schur_slice import build, bern, smul
tau=Fr(sys.argv[1]); A=int(sys.argv[2]); hmax=int(sys.argv[3])
sig=tau*(1+tau); sp=1+2*tau; a2=1+tau; a1=-tau
t0=time.time()
FF,FT=build(tau,A,hmax)
print("tau",tau,"A",A,"hmax",hmax,"Schur time",round(time.time()-t0,1)); sys.stdout.flush()
J=2*hmax-1
# ---------- 1. direct test of (*) in y^a x^j ----------
def ser_exp(c,top):  # exp(c x) coefficients x^0..x^top
    return [c**n/math.factorial(n) for n in range(top+1)]
def direct(sign_e=1,flip_D2=1):
    """sign_e=-1 uses e^{-beta}; flip_D2=-1 flips the sign of the first term (D1^2+D2 in place of D1^2-D2)."""
    top=J+3
    Fx={a:{2*h-2:FF[a][h] for h in FF[a]} for a in FF}
    Tx={a:{2*h-2:FT[a][h] for h in FT[a]} for a in FT}
    res={}
    def add(a,j,v):
        if j<=J: res[(a,j)]=res.get((a,j),Fr(0))+v
    E={a:ser_exp(sign_e*a*tau,top+4) for a in range(0,A+1)}
    for a in range(1,A+1):
        # term1 (2/x)(1-e^{a b x}) Ft
        for j,v in Tx[a].items():
            for n in range(1,top+4):
                add(a,j+n-1,flip_D2*(-2)*E[a][n]*v)
        # term2 (a^2-a)(1+e) F
        for j,v in Fx[a].items():
            add(a,j,(a*a-a)*v)
            for n in range(0,top+4):
                add(a,j+n,(a*a-a)*E[a][n]*v)
    # term3
    G={a:{} for a in range(1,A+1)}
    for a in range(1,A+1):
        for j,v in Fx[a].items():
            for n in range(1,top+4):
                G[a][j+n]=G[a].get(j+n,Fr(0))-E[a][n]*a*v
    for a in range(1,A+1):
        for b in range(1,A+1-a):
            for j,v in G[a].items():
                for k,u in G[b].items():
                    if j+k<=J: add(a+b,j+k,v*u)
    bad=[(k,v) for k,v in res.items() if v!=0 and k[1]<=J]
    return res,bad
res,bad=direct()
print("(*) direct: monomials tested",len(res)," nonzero:",len(bad), bad[:3])
for name,kw in (("control e^{-beta}",dict(sign_e=-1)),("control D1^2+D2",dict(flip_D2=-1))):
    r,b=direct(**kw); print(name," nonzero of",len(r),":",len(b))
sys.stdout.flush()
# ---------- 2. fits to Laurent polynomials in w ----------
# y(s)=u e^{-u}, u=s/(1+s). s-series to order A.
Lf=A+1
u=[Fr(0)]*Lf
for n in range(1,Lf): u[n]=Fr((-1)**(n-1))
def sexp0(lg,L):
    e=[Fr(0)]*L;e[0]=Fr(1)
    for n in range(1,L):
        e[n]=sum(k*lg[k]*e[n-k] for k in range(1,n+1))/n
    return e
ys=smul(u,sexp0([-v for v in u],Lf),Lf)
ypow=[[Fr(1)]+[Fr(0)]*(Lf-1)]
for a in range(1,A+1): ypow.append(smul(ypow[-1],ys,Lf))
def binom_series(p,L):  # (1+s)^p
    r=[Fr(1)]
    for i in range(1,L): r.append(r[-1]*Fr(p-i+1,i))
    return r
def fit(g,bmin,bmax):
    """g: dict a->coef, a>=1. Finds c with w^{-bmin}(g+c) polynomial of degree <= bmax-bmin in s=w-1; returns dict b->coef, checks count."""
    d=bmax-bmin
    Pg=[Fr(0)]*Lf
    for a,v in g.items():
        for i in range(Lf): Pg[i]+=v*ypow[a][i]
    Pg=smul(Pg,binom_series(-bmin,Lf),Lf)
    P1=binom_series(-bmin,Lf)
    c=Fr(0)
    if not (0<=-bmin<=d):
        c=-Pg[d+1]/P1[d+1]
    P=[Pg[i]+c*P1[i] for i in range(Lf)]
    checks=[P[i] for i in range(d+1,Lf)]
    assert all(v==0 for v in checks),("fit fails",bmin,bmax,[str(v) for v in checks][:3])
    # to powers of w
    out={}
    for i in range(d+1):
        # s^i=(w-1)^i
        for k in range(i+1):
            out[k+bmin]=out.get(k+bmin,Fr(0))+P[i]*math.comb(i,k)*(-1)**(i-k)
    return {b:v for b,v in out.items() if v!=0},c,len(checks)
def D(p):
    o={}
    for b,v in p.items():
        if b==0: continue
        o[b+2]=o.get(b+2,Fr(0))+b*v; o[b+1]=o.get(b+1,Fr(0))-b*v
    return {k:v for k,v in o.items() if v!=0}
G={};Z={};Zdot={};dFl={};Ftl={}
minchecks=99
for h in range(0,hmax+1):
    g={a:FF[a][h] for a in FF}    # sigma^{h-1} G_h, a>=1
    gD={a:a*FF[a][h] for a in FF}
    rng={0:(-3,0),1:(-1,0)}.get(h,(2*h-2,5*h-5))
    if h==1:
        # log v not Laurent: fit only D G_1
        pd,_,nc=fit(gD,0,2); G[h]=None
        dFl[h]=pd
        exp={2:Fr(1,24),1:Fr(1,24)*(-2-1/sig),0:Fr(1,24)*(1+1/sig)}
        assert {k:v for k,v in exp.items() if v!=0}==pd,("DG1",pd)
    else:
        pg,c,nc=fit(g,*rng)
        minchecks=min(minchecks,nc)
        G[h]=pg; dFl[h]=D(pg)
        if h==0:
            exp={-3:-Fr(1,6)/sig,-2:-Fr(1,4)/sig}
            pgg={k:v for k,v in pg.items() if k!=0}
            assert pgg==exp,("G0",pg)
        else:
            Z[h]=pg[2*h-2]/sig**(h-1)
    # F_tau: subtract (sp/sig) D part, fit d_tau|_y part
    gt={a:FT[a][h]-(sp/sig)*gD[a] for a in FT}
    if h==1: rt=(-1,0)
    elif h==0: rt=(-3,0)
    else: rt=(2*h-2,5*h-5)
    pt,c,nc=fit(gt,*rt)
    if h>=2:
        Zdot[h]=pt.get(2*h-2,Fr(0)); minchecks=min(minchecks,nc)
    if h==1:
        exp=-sp/(24*sig**2)
        ptt={k:v for k,v in pt.items() if k!=0}
        assert ptt=={-1:exp} or (exp==0 and ptt=={}),("Ft1",pt)
    if h==0:
        ptt={k:v for k,v in pt.items() if k!=0}
        ex0={-3:sp/(6*sig**2),-2:sp/(4*sig**2)}
        assert ptt=={k:v for k,v in ex0.items() if v!=0},("Ft0",pt)
    o=dict(pt)
    for b,v in dFl[h].items(): o[b]=o.get(b,Fr(0))+(sp/sig)*v
    Ftl[h]={k:v for k,v in o.items() if v!=0}
print("Step 3 structure: all fits exact (G_0, DG_1 and d_tau pieces equal the closed forms); min number of overdetermined checks",minchecks)
# beta
L=hmax+2
Sinv=[Fr(1,4**j*math.factorial(2*j+1)) for j in range(L)]
sj=[Fr(0)]*L; sj[0]=Fr(1)
for n in range(1,L): sj[n]=-sum(Sinv[k]*sj[n-k] for k in range(1,n+1))
def betaval(h,t):   # beta_{2h-2}(-t,1+t)=(2h-3)! [t^{2h-2}] f_* , f_*=(1/(a1 a2 l^2)) S(a1 l) S(a2 l)
    a=-t;b=1+t
    return math.factorial(2*h-3)*sum(sj[i]*sj[h-i]*a**(2*i)*b**(2*(h-i)) for i in range(h+1))/(a*b)
def betadot(h,t):
    a=-t;b=1+t; S=0
    # d/dt of a^{2i-1} b^{2j-1}, a'=-1,b'=1
    for i in range(h+1):
        j=h-i
        S+=sj[i]*sj[j]*((2*i-1)*(-1)*a**(2*i-2)*b**(2*j-1)+(2*j-1)*a**(2*i-1)*b**(2*j-2))
    return math.factorial(2*h-3)*S
for h in range(2,hmax+1):
    ok1=sig**(h-1)*Z[h]==betaval(h,tau); ok2=Zdot[h]==betadot(h,tau)
    print(f"h={h}: sigma^(h-1) Z_h == beta_(2h-2): {ok1};  d_tau(sigma^(h-1)Z_h) == d_tau beta: {ok2};  Z_h={Z[h]}")
sys.stdout.flush()
# ---------- 3. (*) as identity of x^a w^b ----------
def padd(d,k,v):
    d[k]=d.get(k,Fr(0))+v
def expbD_terms(series,start):
    """series: dict a->Laurent; returns sum over n>=start of (tau x)^n/n! D^n, truncated at x^(J+3)."""
    out={}
    for a,p in series.items():
        q=dict(p); n=0; fac=Fr(1)
        while a+n<=J+3:
            if n>=start:
                dest=out.setdefault(a+n,{})
                for b,v in q.items(): padd(dest,b,v*fac)
            n+=1; q=D(q); fac=fac*tau/n
    return out
dF={2*h-2:dFl[h] for h in dFl}; Ft={2*h-2:Ftl[h] for h in Ftl}
# term 1 = (2/x)(1-e^{bD}) Ft
T1={}
for a,p in expbD_terms(Ft,1).items():
    d=T1.setdefault(a-1,{})
    for b,v in p.items(): padd(d,b,-2*v)
# term 2 = (D-1)(1+e^{bD}) DF
T2={}
ee=expbD_terms(dF,0)
tmp={}
for a,p in dF.items():
    d=tmp.setdefault(a,{})
    for b,v in p.items(): padd(d,b,v)
for a,p in ee.items():
    d=tmp.setdefault(a,{})
    for b,v in p.items(): padd(d,b,v)
for a,p in tmp.items():
    d=T2.setdefault(a,{})
    for b,v in D(p).items(): padd(d,b,v)
    for b,v in p.items(): padd(d,b,-v)
# term 3 = ((1-e^{bD})DF)^2
hh=expbD_terms(dF,1)
T3={}
for a,p in hh.items():
    for c,q in hh.items():
        if a+c>J: continue
        d=T3.setdefault(a+c,{})
        for b1,v1 in p.items():
            for b2,v2 in q.items():
                padd(d,b1+b2,v1*v2)
tot=0;nz=0
for a in range(-3,J+1):
    keys=set(T1.get(a,{}))|set(T2.get(a,{}))|set(T3.get(a,{}))
    for b in keys:
        v=T1.get(a,{}).get(b,0)+T2.get(a,{}).get(b,0)+T3.get(a,{}).get(b,0)
        tot+=1
        if v!=0: nz+=1
print("(*) as Laurent identity in x^a w^b, a<=",J,": monomials",tot,"nonzero",nz)
minorder=min(b-a for a in T1 if a<=J for b in T1[a] if T1[a][b]!=0)
print("min eps-order in term1 (expect 0)",minorder)
sys.stdout.flush()
# ---------- 4. N side ----------
K=J+6  # truncate at N^{-K}
def Nmul(p,q):
    o={}
    for i,v in p.items():
        for j,u in q.items():
            if i+j>=-K: o[i+j]=o.get(i+j,Fr(0))+v*u
    return o
def Nadd(*ps):
    o={}
    for p in ps:
        for i,v in p.items(): o[i]=o.get(i,Fr(0))+v
    return {k:v for k,v in o.items() if v!=0}
def Nscale(c,p): return {k:c*v for k,v in p.items()}
def Nd(p): return {k-1:k*v for k,v in p.items() if k!=0 and k-1>=-K}
def gbinom(p,i):
    r=Fr(1)
    for m in range(i): r=r*Fr(p-m,m+1)
    return r
def Sh(p):
    o={}
    for k,v in p.items():
        i=0
        while True:
            if k>=0 and i>k: break
            if k-i<-K: break
            o[k-i]=o.get(k-i,Fr(0))+v*gbinom(k,i)*tau**i
            i+=1
    return o
def Nabla(p): return Nadd(Sh(p),Nscale(-1,p))
ell={-i:Fr((-1)**(i-1),i)*tau**i for i in range(1,K+1)}
N1={1:Fr(1)}; Ninv={-1:Fr(1)}
k0=(1+sig)/(24*sig); k0t=-sp/(24*sig**2)
Aa={3:-1/(6*sig),1:k0}; At={3:sp/(6*sig**2),1:k0t}
def S1(p): return Nscale(-1,Nmul(ell,Sh(Nd(p))))   # S_1 f = - ell f'(M)
def build_B(Zv,Zd):
    Bs={2-2*h:sig**(h-1)*Zv[h] for h in Zv}
    Bsd={2-2*h:Zd[h] for h in Zd}
    Bp=Nadd({1:-1/(2*sig)},{-1:-Fr(1,24)},Nd(Bs))
    Bt=Nadd({2:sp/(4*sig**2)},Bsd)
    return Bs,Bsd,Bp,Bt
Bs,Bsd,Bp,Bt=build_B(Z,Zdot)
Ap=Nd(Aa); App=Nd(Ap)
def onepS(p): return Nadd(p,Sh(p))
E0a=Nscale(-2,Nabla(At))
E0b=Nscale(-1,onepS(Ap))
na=Nabla(Ap)
E0c=Nmul(na,na)
E1a=Nadd(Nscale(-2,Nabla(Nadd(Bt,Nscale(sp/sig,Ap)))),Nscale(-2,S1(At)))
E1b=Nadd(onepS(App),Nmul(Ninv,onepS(Ap)),Nscale(-1,onepS(Bp)),Nscale(-1,Nd(S1(Aa))))
E1c=Nscale(2,Nmul(na,Nadd(Nabla(Nadd(Bp,Nscale(-1,Nmul(Ap,Ninv)))),S1(Ap))))
def compare(name,Tx,Ex,eord,kmax):
    good=0;bad=[]
    for k in range(-5,kmax+1):
        v=Tx.get(k-eord,{}).get(k,Fr(0))  # x^a w^b with b=k, b-a=eord
        u=Ex.get(-k,Fr(0))
        if v==u: good+=1
        else: bad.append((k,str(v),str(u)))
    print(f"  {name}: N^-k, k=-5..{kmax}: agree {good}, disagree {len(bad)} {bad[:2]}")
print("order eps^0, three groups:")
compare("term1",T1,E0a,0,J); compare("term2",T2,E0b,0,J); compare("term3",T3,E0c,0,J)
print("sum eps^0 on N side:",Nadd(E0a,E0b,E0c))
print("order eps^1, three groups (x^{k-1} w^k against N^{-k}):")
compare("term1",T1,E1a,1,J+1); compare("term2",T2,E1b,1,J+1); compare("term3",T3,E1c,1,J+1)
S=Nadd(E1a,E1b,E1c)
print("sum eps^1 on N side, powers >= -%d nonzero:"%(J+1),{k:str(v) for k,v in S.items() if k>=-(J+1)})
# ---------- 5. (E1) ----------
def Lop(Bval,Bdot):
    Bp=Nd(Bval)
    MpN=Nadd({1:Fr(2)},{0:tau})
    return Nadd(Nscale(-2,Nabla(Bdot)),Nscale(-1,onepS(Bp)),Nscale(-1/a2,Nmul(MpN,Nabla(Bp))))
MinvNinv=onepS(Ninv)
MN=Nadd({2:Fr(1)},{1:tau})
MpN=Nadd({1:Fr(2)},{0:tau})
rhs=Nscale(1/sig**2,Nadd(Nscale(-1,Nmul(ell,Nadd(MN,{0:-sp/12}))),Nscale(tau/2,MpN),Nscale(-tau*(1+2*sig)/24,MinvNinv)))
print("rhs of (E1): nonneg powers",{k:str(v) for k,v in rhs.items() if k>=0})
L1=Lop(Bs,Bsd)
Zb={h:betaval(h,tau)/sig**(h-1) for h in Z}; Zbd={h:betadot(h,tau) for h in Z}
L2=Lop({2-2*h:sig**(h-1)*Zb[h] for h in Zb},{2-2*h:Zbd[h] for h in Zbd})
kk=2*hmax
print("(E1) Schur Z_h: agree at N^-k, k<=%d:"%kk,all(L1.get(-k,0)==rhs.get(-k,0) for k in range(-3,kk+1)))
print("(E1) beta values: agree at N^-k, k<=%d:"%kk,all(L2.get(-k,0)==rhs.get(-k,0) for k in range(-3,kk+1)))
for k in range(1,7):
    print("   k",k,"rhs",rhs.get(-k,0),"Lop",L1.get(-k,0))
# controls
Zc=dict(Z); Zc[3]=Z[3]+1
Lc=Lop({2-2*h:sig**(h-1)*Zc[h] for h in Zc},Bsd)
print("control Z_3 -> Z_3+1 : first failing k",[k for k in range(1,kk+1) if Lc.get(-k,0)!=rhs.get(-k,0)][:3])
rhs2=Nscale(1/sig**2,Nadd(Nscale(-1,Nmul(ell,Nadd(MN,{0:-sp/12}))),Nscale(tau/2,MpN),Nscale(-tau*(1+2*sig)/24*Fr(3,2),MinvNinv)))
print("control wrong kappa_0 piece: first failing k",[k for k in range(1,kk+1) if L1.get(-k,0)!=rhs2.get(-k,0)][:3])
# ---------- 6. Borel dictionary on the actual series ----------
def Bor(ser):  # ser: dict power of t -> coef, powers>=0
    return {-(k+1):math.factorial(k)*v for k,v in ser.items()}
Ltop=K+4
Pser={k:-(-tau)**k/math.factorial(k) for k in range(1,Ltop+6)}
two_minus_P={k:-v for k,v in Pser.items()}; two_minus_P[0]=Fr(2)
def Tmul(p,q,top):
    o={}
    for i,v in p.items():
        for j,u in q.items():
            if i+j<=top: o[i+j]=o.get(i+j,Fr(0))+v*u
    return o
def Lfrak(f,fdot,top):
    """f,fdot Laurent dicts in t. returns 2P(f_tau/t - f_t/a2)+(2-P) f/a2"""
    inner={}
    for k,v in fdot.items(): inner[k-1]=inner.get(k-1,Fr(0))+v
    for k,v in f.items():
        if k!=0: inner[k-1]=inner.get(k-1,Fr(0))-k*v/a2
    r1=Tmul(Pser,inner,top); r2=Tmul(two_minus_P,f,top)
    o={}
    for k,v in r1.items(): o[k]=o.get(k,Fr(0))+2*v
    for k,v in r2.items(): o[k]=o.get(k,Fr(0))+v/a2
    return {k:v for k,v in o.items() if v!=0}
fs={2*h-2:sig**(h-1)*Z[h]/math.factorial(2*h-3) for h in Z}
fsd={2*h-2:Zdot[h]/math.factorial(2*h-3) for h in Z}
top=2*hmax-2
ls=Lfrak(fs,fsd,top)
BL=Bor(ls)
print("L[Bor f_s] = Bor[Lfrak f_s] through N^-%d:"%(top+1),all(BL.get(-k,0)==L1.get(-k,0) for k in range(1,top+2)))
d0=(1+2*sig)/(24*sig); d0t=-sp/(24*sig**2)
fsing={-2:-1/sig,0:d0}; fsingd={-2:sp/sig**2,0:d0t}
lsing=Lfrak(fsing,fsingd,top)
print("Lfrak f_sing regular:",min(lsing)>=0,"; Bor[Lfrak f_sing] = -rhs(E1) through N^-%d:"%(top+1),all(Bor(lsing).get(-k,0)==-rhs.get(-k,0) for k in range(1,top+2)))
tt={}
def tadd(k,v):
    tt[k]=tt.get(k,Fr(0))+v
for k,v in Pser.items():
    tadd(k-3,2*v); tadd(k-1,-(sp/12)*v)
for k,v in two_minus_P.items():
    tadd(k-2,-tau*v); tadd(k,tau*(1+2*sig)/24*v)
tt={k:v/sig**2 for k,v in tt.items() if v!=0 and k<=top}
print("displayed formula for Lfrak f_sing correct:",tt=={k:v for k,v in lsing.items()})
fstar=dict(fsing); fstard=dict(fsingd)
for h in range(2,hmax+1):
    fstar[2*h-2]=betaval(h,tau)/math.factorial(2*h-3); fstard[2*h-2]=betadot(h,tau)/math.factorial(2*h-3)
ll=Lfrak(fstar,fstard,top)
print("Lfrak f_* = 0 through t^%d:"%top,{k:str(v) for k,v in ll.items() if k<=top}=={})
print("total time",round(time.time()-t0,1))
