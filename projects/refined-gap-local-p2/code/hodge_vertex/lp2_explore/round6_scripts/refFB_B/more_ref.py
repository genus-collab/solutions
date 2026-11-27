# Referee's checks (evidence only).
import sys,os,math,time,itertools
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from g3_ref import *
import g3_ref
from wick_ref import Phi,Theta,mul

# ---- (1) scaling of the genus-3 residual with t and sign controls
print("(1) genus 3 residual scaling and controls",flush=True)
gl=graphs(3)
for e1,e2 in [(Fr(3),Fr(7)),(Fr(5,2),Fr(-3,7))]:
    c0=hodge_const(3,e1,e2); rf=refined(3,e1,e2,3)
    for t in [Fr(1,10**3),Fr(1,10**5),Fr(1,10**7)]:
        y=t**3; Q=Qval(y)
        val=Total(e1,e2,t,Fval(y),8).total(3,gl).a
        r=(val-c0-Q*rf[1]-Q**2*rf[2]-Q**3*rf[3])/Q**3
        print("  eps",e1,e2,"t",t,": ((total-GW0-Q GW1-Q^2 GW2-Q^3 GW3)/Q^3) = %.5e  ; /y = %.5e"%(float(r),float(r/y)),flush=True)
# control: wrong h_w
class RowsC(g3_ref.Rows):
    pass
orig=g3_ref.Theta
def run_control(sE,s1):
    def Th2(c,N):
        return orig(c,N)
    base=g3_ref.Rows
    class R2(base):
        def __init__(s,e1,e2,Z):
            base.__init__(s,e1,e2,Z)
            E=e1+e2
            # multiply rows by correction h'/h
            hnew=mul(mul(orig(sE*(-1)/E,Z),orig(s1/e1,Z),Z),orig(s1/e2,Z),Z)
            hold=mul(mul(orig(-1/E,Z),orig(1/e1,Z),Z),orig(1/e2,Z),Z)
            # ratio
            inv=[Fr(0)]*(Z+1);inv[0]=Fr(1)
            for n in range(1,Z+1):
                inv[n]=-sum(hold[k]*inv[n-k] for k in range(1,n+1))
            rat=mul(hnew,inv,Z)
            s.P0=s.zmulc(rat,s.P0)
            s.Q=s.Top(s.P0)
            s.TQ=s.Top(s.Q)
    save=g3_ref.Rows
    g3_ref.Rows=R2
    e1,e2=Fr(3),Fr(7)
    t=Fr(1,10**4);y=t**3
    v=g3_ref.Total(e1,e2,t,Fval(y),8).total(3,gl)
    g3_ref.Rows=save
    return float(v.a-hodge_const(3,e1,e2))
print("  control flip sign Theta(-z/E): total-GW0=%.4e ; flip sign Theta(z/eps): total-GW0=%.4e (correct rows give 5.37e-11)"%(run_control(-1,1),run_control(1,-1)),flush=True)

# ---- (2) rank-one check of Hodge_v = h_v.1: genus 2,3, graph sum vs lambda integrals
print("(2) Hodge_v = h_v.1 with h_v = prod Theta(-z/v_a), graphs g=2,3",flush=True)
def Hodge_int(g,vs):
    # Int_{Mbar_g} prod_a Sum_i (-1)^i l_i v_a^{-i}
    def Lv(u):
        d={}
        for i in range(g+1):
            e=[0]*g
            if i>0: e[i-1]=1
            d[tuple(e)]=Fr((-1)**i)/u**i
        return d
    r={tuple([0]*g):Fr(1)}
    for u in vs:
        nx={}
        for m,c in r.items():
            for n,x in Lv(u).items():
                k=tuple(p+q for p,q in zip(m,n))
                nx[k]=nx.get(k,0)+c*x
        r=nx
    tot=Fr(0)
    for m,c in r.items():
        if sum((i+1)*x for i,x in enumerate(m))==3*g-3:
            tot+=c*mono(g,m)
    return tot
def rank1(g,s,Z):
    # s list z^n: row. Edge E(a,b)=[psi^a w^b] (1-s(-psi)s(-w))/(psi+w) telescoped; Dil(d)=(-1)^d s[d-1]; Delta=1
    def n(M,b): return Fr((-1)**(M+b))*s[M]*s[b]
    def edge(a,b):
        return -sum(Fr((-1)**k)*n(a+k+1,b-k) for k in range(b+1))
    def dil(d): return Fr((-1)**d)*s[d-1]
    memo={}
    def vert(gv,hs):
        hs=tuple(sorted(hs));key=(gv,hs)
        if key in memo: return memo[key]
        nn=len(hs);dim=3*gv-3+nn-sum(hs);tot=Fr(0)
        if dim>=0:
            def parts(rem,mx):
                if rem==0:
                    yield ();return
                for p in range(min(rem,mx),0,-1):
                    for r in parts(rem-p,p): yield (p,)+r
            for pt in parts(dim,dim):
                ds=tuple(p+1 for p in pt)
                val=wk(gv,tuple(sorted(hs+ds)))
                if val==0: continue
                mult=1
                for x in set(ds): mult*=math.factorial(ds.count(x))
                term=val/mult
                for d in ds: term*=dil(d)
                tot+=term
        memo[key]=tot
        return tot
    tot=Fr(0)
    for gs,edges,aut in graphs(g):
        he=[]
        for (u,v) in edges:
            he+= [u,v]
        nv=len(gs);val=[he.count(i) for i in range(nv)];dims=[3*gs[i]-3+val[i] for i in range(nv)]
        def rec_(h,exps,used):
            if h==len(he):
                p=Fr(1)
                for i,(u,v) in enumerate(edges):
                    p*=edge(exps[2*i],exps[2*i+1])
                    if p==0: return Fr(0)
                for i in range(nv):
                    p*=vert(gs[i],tuple(exps[q] for q in range(len(he)) if he[q]==i))
                    if p==0: return Fr(0)
                return p
            res=Fr(0);vx=he[h]
            for a in range(dims[vx]-used[vx]+1):
                used[vx]+=a
                res+=rec_(h+1,exps+[a],used)
                used[vx]-=a
            return res
        tot+=rec_(0,[],[0]*nv)/aut
    return tot
for g in (2,3):
    for vs in [(Fr(3),Fr(-7),Fr(5,2)),(Fr(1),Fr(2),Fr(5),Fr(-3),Fr(-7)),(Fr(10),Fr(-3),Fr(-7))]:
        Z=3*g-1
        s=[Fr(1)]+[Fr(0)]*Z
        sc=[Fr(1)]+[Fr(0)]*Z
        for u in vs:
            s=mul(s,Theta(-1/u,Z),Z); sc=mul(sc,Theta(1/u,Z),Z)
        print("  g",g,"v",[str(x) for x in vs],"graph==Int Hodge:",rank1(g,s,Z)==Hodge_int(g,vs),"; control Theta(+z/v):",rank1(g,sc,Z)==Hodge_int(g,vs),flush=True)

# ---- (3) d=0 at generic alpha vs T_F formula at genus 3, interpolating 3 e3 k GW along a line through alpha*
print("(3) GW^T_{g,0}=Sum_i W_i^{g-1} Int Hodge_{wt_i}: polynomiality of e3*k*GW along a line and value at alpha*",flush=True)
def GWd0(g,al,e1,e2):
    tot=Fr(0)
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        wt=(3*al[i],al[j]-al[i],al[k]-al[i],-e1,-e2)
        w=Fr(1)
        for x in wt: w*=x
        tot+=w**(g-1)*Hodge_int(g,wt)
    return tot
for g in (2,3):
  for e1,e2 in [(Fr(3),Fr(7)),(Fr(5,2),Fr(-3,7))]:
    E=e1+e2;k=e1*e2;d=(Fr(1),Fr(2,3),Fr(-5,7))
    xs=[Fr(j) for j in range(1,20)]
    ys=[]
    for s in xs:
        al=[E/3+s*x for x in d]
        ys.append(k*al[0]*al[1]*al[2]*GWd0(g,al,e1,e2))
    v=list(ys);order=0
    while any(x!=0 for x in v) and len(v)>1:
        v=[v[i+1]-v[i] for i in range(len(v)-1)];order+=1
    def lag(xs,ys):
        tot=Fr(0)
        for i,xi in enumerate(xs):
            t=ys[i]
            for j,xj in enumerate(xs):
                if j!=i: t*=(0-xj)/(xi-xj)
            tot+=t
        return tot
    dg=order-1
    val=lag(xs[:dg+1],ys[:dg+1])/(k*(E/3)**3)
    print("  g",g,"eps",e1,e2,": degree in s of e3 k GW =",dg,"; value at s=0:",val,"; T_F Hodge formula:",hodge_const(g,e1,e2),"; equal:",val==hodge_const(g,e1,e2),flush=True)

# ---- (4) Phi vs I-function at generic alpha, tiny y, exact rational (Y rational, y derived)
print("(4) generic alpha: Phi(alpha,Y;z)/Phi(alpha,-alpha_j;z) vs stripped I-function f(y,z)/f(y,0)",flush=True)
def series_Y(al,j,K):
    # root of (1+27y)Y^3+e1Y^2+e2Y+e3 with Y(0)=-al[j], Newton on truncated series
    e1=sum(al);e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2];e3=al[0]*al[1]*al[2]
    def smul(a,b):
        r=[Fr(0)]*(K+1)
        for i,x in enumerate(a):
            for jj,y in enumerate(b):
                if i+jj<=K: r[i+jj]+=x*y
        return r
    Y=[-al[j]]+[Fr(0)]*K
    for n in range(1,K+1):
        Y2=smul(Y,Y);Y3=smul(Y2,Y)
        L=[Y3[i]+(27*Y3[i-1] if i>=1 else 0)+e1*Y2[i]+e2*Y[i]+(e3 if i==0 else 0) for i in range(K+1)]
        Lp0=3*al[j]**2-2*e1*al[j]+e2
        Y[n]=-L[n]/Lp0
    return Y
def Cd(al,j,d,N,KM):
    # C_d(z)=prod_{m=0}^{3d-1}(3a_j-m z)/prod_l prod_{m=1}^d (a_l-a_j+m z) as Laurent series: dict power->coef (powers from -d to N)
    # l=j factor gives (d! z^d)^-1
    Lz=N+d+1
    ser=[Fr(1)]+[Fr(0)]*Lz
    for m in range(0,3*d):
        f=[3*al[j],Fr(-m)]+[Fr(0)]*(Lz-1)
        ser=mul(ser,f,Lz)
    for l in range(3):
        if l==j: continue
        for m in range(1,d+1):
            a=al[l]-al[j]; c=Fr(m)/a
            t=[Fr(0)]*(Lz+1);t[0]=ser[0]/a
            for n in range(1,Lz+1): t[n]=(ser[n]/a)-c*t[n-1]
            ser=t
    fac=Fr(1,math.factorial(d))
    return dict((n-d,ser[n]*fac) for n in range(Lz+1))
def fstripped(al,j,K,N):
    # returns list over k=0..K of dict zpower->coef for f(y,z)=Sum_d y^d C_d(z) exp(-ut(y)/z), theta ut = Y_j+a_j
    Y=series_Y(al,j,K)
    ut=[Fr(0)]*(K+1)
    for n in range(1,K+1): ut[n]=Y[n]/n
    # exp(-ut/z) in y-series with coefficient dict in 1/z
    E_=[{0:Fr(1)}]
    for n in range(1,K+1):
        cur={}
        for m in range(1,n+1):
            co=Fr(m,n)*(-ut[m])
            for p,c in E_[n-m].items():
                cur[p-1]=cur.get(p-1,0)+co*c
        E_.append(cur)
    out=[]
    for k in range(K+1):
        cur={}
        for d in range(k+1):
            cd=Cd(al,j,d,N+K,K) if d>0 else {0:Fr(1)}
            for p,c in cd.items():
                for q,x in E_[k-d].items():
                    if p+q<=N: cur[p+q]=cur.get(p+q,0)+c*x
        out.append(cur)
    return out
al=[Fr(1),Fr(2),Fr(5)]
N=5;K=3
for j in range(3):
    fs=fstripped(al,j,K,N)
    neg=max(abs(c) for k in fs for p,c in k.items() if p<0) if any(p<0 for k in fs for p in k) else 0
    print("  colour",j,"poles in z of stripped coefficient (max abs coefficient, must be 0):",neg,flush=True)
    dl=Fr(1,10**7)
    Y0=-al[j]+dl
    e=[Y0+x for x in al]
    y=-(e[0]*e[1]*e[2])/(27*Y0**3)
    P=Phi(al,Y0,N); P00=Phi(al,-al[j],N)
    # invert P00
    inv=[Fr(0)]*(N+1);inv[0]=Fr(1)
    for n in range(1,N+1):
        inv[n]=-sum(P00[k]*inv[n-k] for k in range(1,n+1))
    lhs=mul(P,inv,N)
    # rhs: f(y,z)/f(y,0)
    f=[sum(fs[k].get(n,Fr(0))*y**k for k in range(K+1)) for n in range(N+1)]
    rhs=[x/f[0] for x in f]
    print("   y=%.3e; differences Phi-side minus I-side at z^0..z^%d:"%(float(y),N),["%.2e"%float(lhs[n]-rhs[n]) for n in range(N+1)],flush=True)
