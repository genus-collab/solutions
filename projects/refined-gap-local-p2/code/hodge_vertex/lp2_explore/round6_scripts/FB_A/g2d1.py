# Evidence only. Genus-2 graph sum with BS canonical rows at generic alpha versus localisation (d=0,1).
# Exact fractions. Usage: py -I g2d1.py
from fractions import Fraction as Fr
from itertools import product
from functools import lru_cache
from collections import Counter
import sys, math

NY = 2      # y-truncation of rows (internal regularity tested to y^NY)
ZMAX = 6    # rows known through z^ZMAX

def bern(nmax):
    B=[Fr(0)]*(nmax+1); B[0]=Fr(1)
    for m in range(1,nmax+1):
        B[m]=-sum(Fr(math.comb(m+1,k))*B[k] for k in range(m))/(m+1)
    return B
BN=bern(30)

# ---- y-series: list of Fractions of length NY+1
def s0(c=0):
    r=[Fr(0)]*(NY+1); r[0]=Fr(c); return r
def sadd(a,b): return [x+y for x,y in zip(a,b)]
def ssub(a,b): return [x-y for x,y in zip(a,b)]
def smul(a,b):
    r=[Fr(0)]*(NY+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>NY: break
            r[i+j]+=x*y
    return r
def sscal(c,a): return [c*x for x in a]
def sinv(a):
    r=[Fr(0)]*(NY+1); r[0]=1/a[0]
    for n in range(1,NY+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def stheta(a): return [k*x for k,x in enumerate(a)]
def ssqrt1(a):  # sqrt of series with a[0]=1
    assert a[0]==1
    r=[Fr(0)]*(NY+1); r[0]=Fr(1)
    for n in range(1,NY+1):
        r[n]=(a[n]-sum(r[k]*r[n-k] for k in range(1,n)))/2
    return r
def spow(a,n):
    if n<0: return spow(sinv(a),-n)
    r=s0(1)
    for _ in range(n): r=smul(r,a)
    return r

# ---- z-series (Laurent): dict power->Fraction, truncated at zt
def lmul(a,b,zt):
    r={}
    for i,x in a.items():
        for j,y in b.items():
            if i+j<=zt: r[i+j]=r.get(i+j,Fr(0))+x*y
    return {k:v for k,v in r.items() if v!=0}
def linv_unit(a,zt):
    # a has constant term nonzero and lowest power 0
    a0=a[0]; r={0:1/a0}
    for n in range(1,zt+1):
        r[n]=-sum(a.get(k,Fr(0))*r[n-k] for k in range(1,n+1))/a0
    return r
def lexp(l,zt):
    r={0:Fr(1)}
    for n in range(1,zt+1):
        r[n]=sum(k*l.get(k,Fr(0))*r.get(n-k,Fr(0)) for k in range(1,n+1))/n
    return {k:v for k,v in r.items() if v!=0}
def theta_exp(sign,v,zt):
    # Theta(sign*z/v) as z-series
    logs={}
    m=1
    while 2*m-1<=zt:
        k=2*m-1
        logs[k]=BN[2*m]/(2*m*(2*m-1))*(Fr(sign)/v)**k
        m+=1
    return lexp(logs,zt)

# ---- WK numbers g<=2
@lru_cache(None)
def wk(g,ds):
    ds=tuple(sorted(ds)); n=len(ds)
    if 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and ds==(0,0,0): return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds[0]==0:
        rest=ds[1:]
        tot=Fr(0)
        for i,d in enumerate(rest):
            if d>0:
                new=list(rest); new[i]=d-1
                tot+=wk(g,tuple(sorted(new)))
        return tot
    if 1 in ds:
        i=ds.index(1); new=ds[:i]+ds[i+1:]
        return (2*g-2+n-1)*wk(g,new)
    table={(2,(4,)):Fr(1,1152),(2,(2,3)):Fr(29,5760),(2,(2,2,2)):Fr(7,240)}
    return table[(g,ds)]

class Theory:
    def __init__(self,al,e1,e2,csign=1):
        self.al=[Fr(x) for x in al]; self.e1=Fr(e1); self.e2=Fr(e2)
        a=self.al
        self.E1=sum(a); self.E2=a[0]*a[1]+a[1]*a[2]+a[0]*a[2]; self.E3=a[0]*a[1]*a[2]
        self.csign=csign
        self.wt=[self.weights(i) for i in range(3)]
        self.frob()
        self.rows()
    def weights(self,i):
        a=self.al
        return [3*a[i]]+[a[j]-a[i] for j in range(3) if j!=i]+[-self.e1,-self.e2]
    def frob(self):
        A=s0(0)
        for d in range(NY+1):
            A[d]=Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d
        self.A=A
        X=s0(0)
        for d in range(NY+1): X[d]=Fr((-27)**d)
        self.X=X
        self.D=smul(X,sinv(smul(A,A)))
        self.Ey=sscal(Fr(1,2),ssub(self.D,X))
    def Yroot(self,i):
        a=self.al
        Y=s0(-a[i])
        ysr=s0(0)
        if NY>=1: ysr[1]=Fr(1)
        for _ in range(NY+2):
            L=smul(smul(sadd(Y,s0(a[0])),sadd(Y,s0(a[1]))),sadd(Y,s0(a[2])))
            Y3=smul(smul(Y,Y),Y)
            L=sadd(L,sscal(Fr(27),smul(ysr,Y3)))
            dL=sadd(sadd(smul(sadd(Y,s0(a[1])),sadd(Y,s0(a[2]))),smul(sadd(Y,s0(a[0])),sadd(Y,s0(a[2])))),smul(sadd(Y,s0(a[0])),sadd(Y,s0(a[1]))))
            dL=sadd(dL,sscal(Fr(81),smul(ysr,smul(Y,Y))))
            Y=ssub(Y,smul(L,sinv(dL)))
        return Y
    def Wfn(self,Y):
        return sadd(sadd(sscal(self.E1,smul(Y,Y)),sscal(2*self.E2,Y)),s0(3*self.E3))
    def rho(self,i):
        # normalised row 0: exp(-u/z) * sum_d c_d(z) y^d, returned as list over y of z-dicts
        a=self.al; zt=ZMAX+NY+2
        cs=[]
        for d in range(NY+1):
            num={0:Fr(1)}
            for m in range(3*d):
                num=lmul(num,{0:3*a[i],1:Fr(-m)},zt+NY+1)
            den={0:Fr(1)}
            for j in range(3):
                if j==i: continue
                for m in range(1,d+1):
                    den=lmul(den,{0:a[j]-a[i],1:Fr(m)},zt+NY+1)
            q=lmul(num,linv_unit(den,zt+NY+1),zt+NY+1)
            q={k-d:v/math.factorial(d) for k,v in q.items()}
            cs.append(q)
        Y=self.Yroot(i)
        u=[Fr(0)]*(NY+1)
        for d in range(1,NY+1): u[d]=Y[d]/d
        ex=[{} for _ in range(NY+1)]; ex[0]={0:Fr(1)}
        for n in range(1,NY+1):
            tot={}
            for k in range(1,n+1):
                for p,v in ex[n-k].items():
                    tot[p-1]=tot.get(p-1,Fr(0))+Fr(k,n)*(-u[k])*v
            ex[n]=tot
        res=[]
        for n in range(NY+1):
            tot={}
            for k in range(n+1):
                for p,v in lmul(cs[k],ex[n-k],zt).items():
                    tot[p]=tot.get(p,Fr(0))+v
            tot={k:v for k,v in tot.items() if v!=0}
            if tot and min(tot)<0:
                raise Exception("pole survives at y^%d colour %d: %s"%(n,i,min(tot)))
            res.append(tot)
        return res,Y
    def rows(self):
        self.Y=[];self.Del=[];self.P=[]
        for i in range(3):
            rh,Y=self.rho(i)
            W=self.Wfn(Y)
            Del=sscal(3*self.e1*self.e2,W)
            wprod=Fr(1)
            for w in self.wt[i]: wprod*=w
            assert Del[0]==wprod
            sq=ssqrt1(sscal(1/wprod,Del))
            h={0:Fr(1)}
            for w in self.wt[i]:
                h=lmul(h,theta_exp(-self.csign,w,ZMAX+2),ZMAX+2)
            P0={}
            for n in range(NY+1):
                for p,v in rh[n].items():
                    if p<=ZMAX+2:
                        P0.setdefault(p,s0(0))
                        P0[p][n]+=v
            P0={p:smul(sq,c) for p,c in P0.items()}
            Ph={}
            for p,c in P0.items():
                for q,hv in h.items():
                    if p+q<=ZMAX+2:
                        Ph.setdefault(p+q,s0(0))
                        Ph[p+q]=sadd(Ph[p+q],sscal(hv,c))
            P0=Ph
            lg=smul(stheta(Del),sinv(Del))
            def Th(F,shift):
                out={}
                for p,c in F.items():
                    out.setdefault(p,s0(0))
                    out[p]=sadd(out[p],smul(ssub(Y,shift),c))
                    out.setdefault(p+1,s0(0))
                    out[p+1]=sadd(out[p+1],ssub(stheta(c),sscal(Fr(1,2),smul(lg,c))))
                return {p:c for p,c in out.items() if p<=ZMAX+2}
            P1=Th(P0,s0(0)); iA=sinv(self.A)
            P1={p:smul(iA,c) for p,c in P1.items()}
            shift=sscal(self.E1,self.Ey)
            P2=Th(P1,shift); iD=sinv(self.D)
            P2={p:smul(iD,c) for p,c in P2.items()}
            self.Y.append(Y);self.Del.append(Del);self.P.append([P0,P1,P2])
        k=3*self.e1*self.e2
        self.etainv={(0,0):k*self.E3,(1,1):-k*self.E1,(1,2):-k,(2,1):-k}
    def coef(self,i,a,p):
        return self.P[i][a].get(p,s0(0))
    def pair_num(self,i,j,sgn1,sgn2,k,l):
        tot=s0(0)
        for (a,b),e in self.etainv.items():
            tot=sadd(tot,sscal(e*sgn1**k*sgn2**l,smul(self.coef(i,a,k),self.coef(j,b,l))))
        return tot
    def check_sympl(self,zmax=5):
        bad=0
        for i in range(3):
            for j in range(3):
                for n in range(zmax+1):
                    tot=s0(0)
                    for k in range(n+1):
                        tot=sadd(tot,self.pair_num(i,j,1,-1,k,n-k))
                    target=self.Del[i] if (i==j and n==0) else s0(0)
                    if tot!=target: bad+=1
        return bad
    def edge(self,i,j,a,b):
        cache=self.__dict__.setdefault('_ec',{})
        key=(i,j,a,b)
        if key in cache: return cache[key]
        def n(k,l):
            r=sscal(Fr(-1),self.pair_num(i,j,-1,-1,k,l))
            if k==0 and l==0 and i==j: r=sadd(r,self.Del[i])
            return r
        def q(k,l):
            if l==0: return n(k+1,0)
            return ssub(n(k+1,l),q(k+1,l-1))
        r=q(a,b); cache[key]=r
        return r
    def dil(self,i,d):
        return sscal(Fr((-1)**d),self.coef(i,0,d-1))
    def vert(self,i,g,hs):
        cache=self.__dict__.setdefault('_vc',{})
        hs=tuple(sorted(hs)); key=(i,g,hs)
        if key in cache: return cache[key]
        n=len(hs); dim=3*g-3+n-sum(hs)
        tot=s0(0)
        if dim>=0:
            def parts(rem,maxp):
                if rem==0:
                    yield ()
                    return
                for p in range(min(rem,maxp),0,-1):
                    for r in parts(rem-p,p): yield (p,)+r
            for pt in parts(dim,dim):
                ds=tuple(p+1 for p in pt)
                val=wk(g,tuple(sorted(hs+ds)))
                if val==0: continue
                mult=1
                for c in Counter(ds).values(): mult*=math.factorial(c)
                t=s0(val/mult)
                for d in ds: t=smul(t,self.dil(i,d))
                tot=sadd(tot,t)
        res=smul(spow(self.Del[i],g-1),tot) if g!=1 else tot
        cache[key]=res
        return res
    def F2(self):
        tot=s0(0)
        col=range(3)
        V=self.vert;E=self.edge
        for i in col: tot=sadd(tot,V(i,2,()))
        for i in col:
            for a in range(3):
                for b in range(3-a):
                    tot=sadd(tot,sscal(Fr(1,2),smul(V(i,1,(a,b)),E(i,i,a,b))))
        for i in col:
            for a,b,c,d in product(range(2),repeat=4):
                if a+b+c+d>1: continue
                t=smul(V(i,0,(a,b,c,d)),smul(E(i,i,a,b),E(i,i,c,d)))
                tot=sadd(tot,sscal(Fr(1,8),t))
        for i in col:
            for j in col:
                for a in range(2):
                    for b in range(2):
                        t=smul(smul(V(i,1,(a,)),E(i,j,a,b)),V(j,1,(b,)))
                        tot=sadd(tot,sscal(Fr(1,2),t))
        for i in col:
            for j in col:
                for a in range(2):
                    t=smul(smul(V(i,1,(a,)),E(i,j,a,0)),smul(V(j,0,(0,0,0)),E(j,j,0,0)))
                    tot=sadd(tot,sscal(Fr(1,2),t))
        for i in col:
            for j in col:
                t=smul(smul(V(i,0,(0,0,0)),E(i,i,0,0)),smul(E(i,j,0,0),smul(V(j,0,(0,0,0)),E(j,j,0,0))))
                tot=sadd(tot,sscal(Fr(1,8),t))
        for i in col:
            for j in col:
                e=E(i,j,0,0)
                t=smul(smul(V(i,0,(0,0,0)),V(j,0,(0,0,0))),smul(e,smul(e,e)))
                tot=sadd(tot,sscal(Fr(1,12),t))
        return tot

# ---- Localisation, genus 2
def hodge_prod(weights,g):
    # product of Lambda^vee_g(w); returns {(i,j):coef} in monomials lam1^i lam2^j
    # g=1: w-lam1, lam1^2=0.  g=2: w^2-w lam1+lam2, lam1^2=2 lam2, lam2^2=0.
    r={(0,0):Fr(1)}
    for w in weights:
        terms={(0,0):w,(1,0):Fr(-1)} if g==1 else {(0,0):w*w,(1,0):-w,(0,1):Fr(1)}
        new={}
        for m,c in r.items():
            for t,d in terms.items():
                a=m[0]+t[0];b=m[1]+t[1];cc=c*d
                if g==1 and a>=2: continue
                while g==2 and a>=2:
                    a-=2;b+=1;cc*=2
                if g==2 and b>=2: continue
                new[(a,b)]=new.get((a,b),Fr(0))+cc
        r=new
    return r
def int_g2_n0(m):
    return {(1,1):Fr(1,5760)}.get(m,Fr(0))
def int_g2_n1(m,psi):
    tb={((0,0),4):Fr(1,1152),((1,0),3):Fr(1,480),((0,1),2):Fr(7,5760),((1,1),1):Fr(1,2880)}
    return tb.get((m,psi),Fr(0))
def int_g1_n1(m,psi):
    return {((0,0),1):Fr(1,24),((1,0),0):Fr(1,24)}.get((m,psi),Fr(0))
def loc_d0(T):
    tot=Fr(0)
    for i in range(3):
        wp=Fr(1)
        for w in T.wt[i]: wp*=w
        hp=hodge_prod(T.wt[i],2)
        tot+=sum(c*int_g2_n0(m) for m,c in hp.items())/wp
    return tot
def vert_int(weights,g,tau):
    hp=hodge_prod(weights,g)
    tot=Fr(0)
    for m,c in hp.items():
        for psi in range(0,6):
            val=int_g2_n1(m,psi) if g==2 else int_g1_n1(m,psi)
            if val!=0:
                tot+=c*val/tau**(psi+1)
    return tot
def loc_d1(T):
    a=T.al;tot=Fr(0)
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j
        tau=a[j]-a[i]
        edgefac=(2*a[i]+a[j])*(a[i]+2*a[j])/(-tau*tau*(a[k]-a[i])*(a[k]-a[j])*T.e1*T.e2)
        wi=T.wt[i];wj=T.wt[j]
        tot+=vert_int(wi,2,tau)*(-tau)*edgefac
        tot+=vert_int(wj,2,-tau)*(tau)*edgefac
        tot+=vert_int(wi,1,tau)*vert_int(wj,1,-tau)*edgefac
    return tot
def refined_d1(e1,e2):
    E=e1+e2
    return (3*(Fr(7,5760)*(e1**4+e2**4)+e1**2*e2**2/576)-E**2*(e1**2+e2**2)/24+E**4/12)/(e1*e2)

if __name__=="__main__":
    samples=[((1,2,5),3,7),((2,-3,7),1,4),((1,3,8),Fr(5,2),Fr(-3,7))]
    for al,e1,e2 in samples:
        T=Theory(al,e1,e2)
        print("alpha",al,"eps",e1,e2)
        print(" symplectic defects (all i,j,z^<=5, to y^%d):"%NY,T.check_sympl())
        F=T.F2()
        print(" graph y^0,y^1:",F[0],F[1])
        print(" local d=0,d=1:",loc_d0(T),loc_d1(T))
        print(" match d0,d1:",F[0]==loc_d0(T),F[1]==loc_d1(T))
        Tc=Theory(al,e1,e2,csign=-1)
        Fc=Tc.F2()
        print(" CONTROL opposite sign of c_j: match d0,d1:",Fc[0]==loc_d0(T),Fc[1]==loc_d1(T))
    for e1,e2 in [(Fr(3),Fr(7)),(Fr(1),Fr(1))]:
        E=e1+e2
        print("T_F limit eps",e1,e2,"refined BPS value",float(refined_d1(e1,e2)),refined_d1(e1,e2))
        for k in [2,4,6,8]:
            dl=Fr(1,10**k)
            al=(E/3+dl,E/3+2*dl-3*dl*dl,E/3-3*dl+dl*dl)
            class P: pass
            P.al=[Fr(x) for x in al];P.e1=e1;P.e2=e2
            P.wt=[[3*P.al[i]]+[P.al[j]-P.al[i] for j in range(3) if j!=i]+[-e1,-e2] for i in range(3)]
            print("  delta 1e-%d loc_d1 ="%k,float(loc_d1(P)))
