# Evidence only (exact fractions). Genus-2 colour sum as Tr(x)Tr of an element of C(x)C[F], where
# C = K[Y]/(prod(Y+alpha_j)+27 y Y^3), at rational (alpha,eps,y) -- including alpha on the diagonal, where C is a field.
# Row 0 comes from the algebraic Wick expansion Phi(alpha,Y;z) (no y-series, no I-function), P_0 = Theta(z/e1)Theta(z/e2) Phi.
import sys, math, time
from fractions import Fraction as Fr
from functools import lru_cache
from collections import Counter
from itertools import product

def bern(nmax):
    B=[Fr(0)]*(nmax+1); B[0]=Fr(1)
    for m in range(1,nmax+1):
        B[m]=-sum(Fr(math.comb(m+1,k))*B[k] for k in range(m))/(m+1)
    return B
BN=bern(30)

@lru_cache(None)
def wk(g,ds):
    ds=tuple(sorted(ds)); n=len(ds)
    if 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and ds==(0,0,0): return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds[0]==0:
        rest=ds[1:]; tot=Fr(0)
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

def zexp(l,zt):
    r=[Fr(0)]*(zt+1); r[0]=Fr(1)
    for n in range(1,zt+1):
        r[n]=sum(k*l[k]*r[n-k] for k in range(1,n+1))/n
    return r
def theta_series(c,zt):
    # Theta(c z) as list in z; Theta(x)=exp sum B_{2m} x^{2m-1}/(2m(2m-1))
    l=[Fr(0)]*(zt+1)
    m=1
    while 2*m-1<=zt:
        l[2*m-1]=BN[2*m]/(2*m*(2*m-1))*c**(2*m-1); m+=1
    return zexp(l,zt)
def zmul(a,b,zt):
    r=[Fr(0)]*(zt+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=zt: r[i+j]+=x*y
    return r

class Ring:
    """C=K[Y]/L; elements are 3-lists; L = D Y^3+e1 Y^2+e2 Y+e3, D=1+27y."""
    def __init__(self,al,y):
        self.al=[Fr(x) for x in al]; self.y=Fr(y)
        a=self.al
        self.e1=sum(a); self.e2=a[0]*a[1]+a[0]*a[2]+a[1]*a[2]; self.e3=a[0]*a[1]*a[2]
        self.Dc=1+27*self.y
        r=[-self.e3/self.Dc,-self.e2/self.Dc,-self.e1/self.Dc]  # Y^3
        self.r3=r
        self.r4=[r[2]*r[0], r[0]+r[2]*r[1], r[1]+r[2]*r[2]]
        self.one=[Fr(1),Fr(0),Fr(0)]; self.zero=[Fr(0)]*3; self.Y=[Fr(0),Fr(1),Fr(0)]
        self.tau=[self.tr_pow(k) for k in range(3)]
    def mul(self,x,y):
        c=[x[0]*y[0], x[0]*y[1]+x[1]*y[0], x[0]*y[2]+x[1]*y[1]+x[2]*y[0], x[1]*y[2]+x[2]*y[1], x[2]*y[2]]
        r3,r4=self.r3,self.r4
        return [c[0]+c[3]*r3[0]+c[4]*r4[0], c[1]+c[3]*r3[1]+c[4]*r4[1], c[2]+c[3]*r3[2]+c[4]*r4[2]]
    def add(self,x,y): return [a+b for a,b in zip(x,y)]
    def sub(self,x,y): return [a-b for a,b in zip(x,y)]
    def scal(self,c,x): return [c*a for a in x]
    def matrix(self,x):
        cols=[x,self.mul(x,self.Y)]; cols.append(self.mul(cols[1],self.Y))
        return [[cols[j][i] for j in range(3)] for i in range(3)]
    def inv(self,x):
        M=self.matrix(x)
        A=[M[i]+[Fr(1) if i==0 else Fr(0)] for i in range(3)]
        for c in range(3):
            p=next(r for r in range(c,3) if A[r][c]!=0)
            A[c],A[p]=A[p],A[c]
            iv=1/A[c][c]; A[c]=[v*iv for v in A[c]]
            for r in range(3):
                if r!=c and A[r][c]!=0:
                    f=A[r][c]; A[r]=[u-f*v for u,v in zip(A[r],A[c])]
        return [A[i][3] for i in range(3)]
    def tr_pow(self,k):
        v=self.one
        for _ in range(k): v=self.mul(v,self.Y)
        M=self.matrix(v)
        return M[0][0]+M[1][1]+M[2][2]
    def tr(self,x): return sum(a*b for a,b in zip(x,self.tau))

class Jet:
    """C[eps]/eps^3, flow along theta: f(Y(s+eps)) = f + eps theta f + eps^2/2 theta^2 f"""
    def __init__(self,R): self.R=R
    def mul(self,a,b):
        R=self.R
        c0=R.mul(a[0],b[0])
        c1=R.add(R.mul(a[0],b[1]),R.mul(a[1],b[0]))
        c2=R.add(R.add(R.mul(a[0],b[2]),R.mul(a[1],b[1])),R.mul(a[2],b[0]))
        return [c0,c1,c2]
    def add(self,a,b):
        R=self.R; return [R.add(a[i],b[i]) for i in range(3)]
    def scal(self,c,a):
        R=self.R; return [R.scal(c,a[i]) for i in range(3)]
    def const(self,x):
        return [x,self.R.zero,self.R.zero]
    def inv(self,a):
        R=self.R
        i0=R.inv(a[0])
        i1=R.scal(Fr(-1),R.mul(i0,R.mul(a[1],i0)))
        i2=R.scal(Fr(-1),R.mul(i0,R.add(R.mul(a[2],i0),R.mul(a[1],i1))))
        return [i0,i1,i2]

def yjets(R):
    # theta Y = g(Y) with g = -Y prod(Y+al)/W; theta^2 Y = g'(Y) g(Y); computed in C using explicit derivative
    a=R.al; Y=R.Y
    q=[R.add(Y,R.scal(a[j],R.one)) for j in range(3)]
    Yp=R.mul(Y,Y)
    W=R.add(R.add(R.scal(R.e1,Yp),R.scal(2*R.e2,Y)),R.scal(3*R.e3,R.one))
    Wp=R.add(R.scal(2*R.e1,Y),R.scal(2*R.e2,R.one))
    Pi=R.mul(R.mul(q[0],q[1]),q[2])
    Pip=R.add(R.add(R.mul(q[1],q[2]),R.mul(q[0],q[2])),R.mul(q[0],q[1]))
    iW=R.inv(W)
    g=R.scal(Fr(-1),R.mul(R.mul(Y,Pi),iW))
    # g' = -(Pi + Y Pi')/W + Y Pi W'/W^2
    gp=R.add(R.scal(Fr(-1),R.mul(R.add(Pi,R.mul(Y,Pip)),iW)),R.mul(R.mul(R.mul(Y,Pi),Wp),R.mul(iW,iW)))
    return Y,g,R.mul(gp,g),W,iW

def wick_rows(R,ZTOP):
    """Phi_n in Jet ring, n=0..ZTOP. Formula: exp(z/2 G dd) exp(sum_{m>=3}T_m/z)|_0 with G=-H^{-1}, jets polynomial in (alpha,Y)."""
    J=Jet(R); a=R.al
    Y0,Y1,Y2,W0,iW0=yjets(R)
    YJ=[Y0,Y1,R.scal(Fr(1,2),Y2)]
    q=[J.add(YJ,J.const(R.scal(a[j],R.one))) for j in range(3)]
    P=J.scal(Fr(-3),YJ)
    c=J.scal(Fr(1,9),P)  # H=diag(q)+c 11^T
    def adj(i,j):
        if i==j:
            o=[k for k in range(3) if k!=i]
            return J.add(J.mul(q[o[0]],q[o[1]]),J.mul(c,J.add(q[o[0]],q[o[1]])))
        k=3-i-j
        return J.scal(Fr(-1),J.mul(c,q[k]))
    det=J.add(J.mul(J.mul(q[0],q[1]),q[2]),J.mul(c,J.add(J.add(J.mul(q[0],q[1]),J.mul(q[0],q[2])),J.mul(q[1],q[2]))))
    idet=J.inv(det)
    G=[[J.scal(Fr(-1),J.mul(adj(i,j),idet)) for j in range(3)] for i in range(3)]
    def padd(A,B,cc=Fr(1)):
        r=dict(A)
        for m,x in B.items():
            y=J.scal(cc,x)
            r[m]=J.add(r[m],y) if m in r else y
        return r
    def pmul(A,B):
        r={}
        for ma,ca in A.items():
            for mb,cb in B.items():
                m=(ma[0]+mb[0],ma[1]+mb[1],ma[2]+mb[2])
                t=J.mul(ca,cb)
                r[m]=J.add(r[m],t) if m in r else t
        return r
    def T_m(m):
        out={}
        for j in range(3):
            e=[0,0,0]; e[j]=m
            out[tuple(e)]=J.scal(Fr(1,math.factorial(m)),q[j])
        for e0 in range(m+1):
            for e1_ in range(m+1-e0):
                e2_=m-e0-e1_
                mult=Fr(math.factorial(m),math.factorial(e0)*math.factorial(e1_)*math.factorial(e2_))
                cc=J.scal(mult*Fr(1,3**m)/math.factorial(m),P)
                k=(e0,e1_,e2_)
                out[k]=J.add(out[k],cc) if k in out else cc
        return out
    Ts={m:T_m(m) for m in range(3,2*ZTOP+3)}
    Es={0:{(0,0,0):J.const(R.one)}}
    for n in range(1,2*ZTOP+1):
        cur={}
        for k in range(1,n+1):
            cur=padd(cur,pmul(Ts[k+2],Es[n-k]),Fr(k,n))
        Es[n]=cur
    memo={}
    def mom(b):
        if b in memo: return memo[b]
        if sum(b)==0:
            r=J.const(R.one)
        elif sum(b)%2:
            r=None
        else:
            cc=next(j for j in range(3) if b[j]>0)
            r=None
            for d in range(3):
                nb=list(b); nb[cc]-=1
                if nb[d]==0: continue
                mult=nb[d]; nb[d]-=1
                sub=mom(tuple(nb))
                if sub is None: continue
                t=J.scal(Fr(mult),J.mul(G[cc][d],sub))
                r=t if r is None else J.add(r,t)
        memo[b]=r
        return r
    Phi=[]
    for n in range(ZTOP+1):
        tot=[R.zero,R.zero,R.zero]
        for m,cf in Es[2*n].items():
            mm=mom(m)
            if mm is not None: tot=J.add(tot,J.mul(cf,mm))
        Phi.append(tot)
    return J,Phi,W0,iW0

class Total:
    def __init__(self,al,e1,e2,y,ZTOP=3,hsign=1,wick=None):
        self.R=R=Ring(al,y); self.e1=Fr(e1); self.e2=Fr(e2); self.ZTOP=ZTOP
        J,Phi,W0,iW0=wick_rows(R,ZTOP)
        self.J=J; self.Phi=Phi
        h=zmul(theta_series(hsign/self.e1,ZTOP),theta_series(hsign/self.e2,ZTOP),ZTOP)
        P0=[]
        for n in range(ZTOP+1):
            tot=[R.zero,R.zero,R.zero]
            for k in range(n+1):
                tot=J.add(tot,J.scal(h[k],Phi[n-k]))
            P0.append(tot)
        self.P0=[p[0] for p in P0]
        kap=self.e1*self.e2
        self.Del=R.scal(3*kap,W0)
        self.iDel=R.scal(1/(3*kap),iW0)
        # tau=theta log W to first order in eps
        Y0,Y1,Y2,_,_=yjets(R)
        YJ=[Y0,Y1,R.scal(Fr(1,2),Y2)]
        WJ=J.add(J.add(J.scal(R.e1,J.mul(YJ,YJ)),J.scal(2*R.e2,YJ)),J.const(R.scal(3*R.e3,R.one)))
        dWJ=[WJ[1],R.scal(Fr(2),WJ[2]),R.zero]
        tauJ=J.mul(dWJ,J.inv(WJ))
        tau0,tau1=tauJ[0],tauJ[1]
        Y=R.Y
        def d(f): return [f[1],R.scal(Fr(2),f[2])]   # D f, valid mod eps^2
        Q=[]
        for n in range(ZTOP+1):
            f=P0[n]
            q0=R.mul(Y,f[0]); q1=R.add(R.mul(Y,f[1]),R.mul(Y1,f[0]))
            if n>0:
                g=P0[n-1]; dg=d(g)
                q0=R.add(q0,R.sub(dg[0],R.scal(Fr(1,2),R.mul(tau0,g[0]))))
                q1=R.add(q1,R.sub(dg[1],R.scal(Fr(1,2),R.add(R.mul(tau0,g[1]),R.mul(tau1,g[0])))))
            Q.append([q0,q1])
        self.Q=[q[0] for q in Q]
        TQ=[]
        for n in range(ZTOP+1):
            x=R.mul(Y,Q[n][0])
            if n>0:
                x=R.add(x,R.sub(Q[n-1][1],R.scal(Fr(1,2),R.mul(tau0,Q[n-1][0]))))
            TQ.append(x)
        self.TQ=TQ
        self.eta00=3*R.e3*kap
        self.etaQ=-3*R.e1*kap
        self.eta12=-3*kap
        self._e={}
        self._pw=[R.one]
        for k in range(1,5): self._pw.append(R.mul(self._pw[-1],R.Y))
    # C (x) C tensors: 3x3 lists
    def tens(self,x,y): return [[a*b for b in y] for a in x]
    def tadd(self,A,B,c=Fr(1)): return [[A[i][j]+c*B[i][j] for j in range(3)] for i in range(3)]
    def tzero(self): return [[Fr(0)]*3 for _ in range(3)]
    def tmul(self,A,B):
        pw=self._pw
        r=self.tzero()
        for k in range(3):
            for l in range(3):
                a=A[k][l]
                if a==0: continue
                for m in range(3):
                    for n in range(3):
                        b=B[m][n]
                        if b==0: continue
                        c=a*b; u=pw[k+m]; v=pw[l+n]
                        for i in range(3):
                            if u[i]==0: continue
                            for j in range(3):
                                if v[j]!=0: r[i][j]+=c*u[i]*v[j]
        return r
    def mu(self,A):
        pw=self._pw
        r=[Fr(0)]*3
        for k in range(3):
            for l in range(3):
                if A[k][l]!=0:
                    r=[r[i]+A[k][l]*pw[k+l][i] for i in range(3)]
        return r
    def trtr(self,A):
        t=self.R.tau
        return sum(A[k][l]*t[k]*t[l] for k in range(3) for l in range(3))
    def pairing(self,M,b,F):
        # [z^M w^b] of <P(z),P(w)> with F numeric; A has already cancelled
        R=self.R; kap=self.e1*self.e2; Dc=R.Dc
        T=self.tens(self.P0[M],self.P0[b]); T=[[self.eta00*x for x in r] for r in T]
        T=self.tadd(T,self.tens(self.Q[M],self.Q[b]),self.etaQ)
        T=self.tadd(T,self.tens(self.Q[M],self.TQ[b]),self.eta12*Dc)
        T=self.tadd(T,self.tens(self.TQ[M],self.Q[b]),self.eta12*Dc)
        if M>=1: T=self.tadd(T,self.tens(self.Q[M-1],self.Q[b]),3*kap*Dc*F)
        if b>=1: T=self.tadd(T,self.tens(self.Q[M],self.Q[b-1]),3*kap*Dc*F)
        return T
    def pairing_explicit(self,M,b,F,A):
        # same pairing from the explicit rows P1=Q/A, P2=Dc A [TQ - zF Q - e1 X (A^-2-1)/2 Q]; tests the cancellation of A
        R=self.R; kap=self.e1*self.e2; Dc=R.Dc; X=1/Dc; e1=R.e1
        ee=e1*X*(1/A**2-1)/2
        def P1(n): return R.scal(1/A,self.Q[n])
        def P2(n):
            x=R.sub(self.TQ[n],R.scal(ee,self.Q[n]))
            if n>=1: x=R.sub(x,R.scal(F,self.Q[n-1]))
            return R.scal(Dc*A,x)
        T=self.tens(self.P0[M],self.P0[b]); T=[[self.eta00*x for x in r] for r in T]
        T=self.tadd(T,self.tens(P1(M),P1(b)),-3*e1*kap)
        T=self.tadd(T,self.tens(P1(M),P2(b)),-3*kap)
        T=self.tadd(T,self.tens(P2(M),P1(b)),-3*kap)
        return T
    def sympl_defect(self,n,F):
        T=self.tzero()
        for M in range(n+1):
            T=self.tadd(T,self.pairing(M,n-M,F),Fr((-1)**(n-M)))
        return any(x!=0 for r in T for x in r)
    def edge(self,a,b,F):
        key=(a,b,F)
        if key in self._e: return self._e[key]
        T=self.tzero()
        for k in range(b+1):
            M=a+k+1; bb=b-k
            T=self.tadd(T,self.pairing(M,bb,F),Fr(-1)*Fr((-1)**k)*Fr((-1)**(M+bb)))
        self._e[key]=T
        return T
    def dil(self,d):
        return self.R.scal(Fr((-1)**d),self.P0[d-1])
    def vert(self,g,hs):
        R=self.R
        hs=tuple(sorted(hs)); n=len(hs); dim=3*g-3+n-sum(hs)
        tot=R.zero
        if dim>=0:
            def parts(rem,maxp):
                if rem==0:
                    yield (); return
                for p in range(min(rem,maxp),0,-1):
                    for r in parts(rem-p,p): yield (p,)+r
            for pt in parts(dim,dim):
                ds=tuple(p+1 for p in pt)
                val=wk(g,tuple(sorted(hs+ds)))
                if val==0: continue
                mult=1
                for c in Counter(ds).values(): mult*=math.factorial(c)
                t=R.scal(val/mult,R.one)
                for d in ds: t=R.mul(t,self.dil(d))
                tot=R.add(tot,t)
        if g==2: tot=R.mul(self.Del,tot)
        if g==0: tot=R.mul(self.iDel,tot)
        return tot
    def F2(self,F):
        R=self.R; V=self.vert
        E=lambda a,b: self.edge(a,b,F)
        tot=R.tr(V(2,()))
        for a in range(3):
            for b in range(3-a):
                tot+=Fr(1,2)*R.tr(R.mul(V(1,(a,b)),self.mu(E(a,b))))
        for a,b,c,d in product(range(2),repeat=4):
            if a+b+c+d>1: continue
            tot+=Fr(1,8)*R.tr(R.mul(V(0,(a,b,c,d)),R.mul(self.mu(E(a,b)),self.mu(E(c,d)))))
        for a in range(2):
            for b in range(2):
                t=self.tmul(self.tens(V(1,(a,)),V(1,(b,))),E(a,b))
                tot+=Fr(1,2)*self.trtr(t)
        V3=V(0,(0,0,0)); L00=R.mul(V3,self.mu(E(0,0)))
        for a in range(2):
            t=self.tmul(self.tens(V(1,(a,)),L00),E(a,0))
            tot+=Fr(1,2)*self.trtr(t)
        t=self.tmul(self.tens(L00,L00),E(0,0))
        tot+=Fr(1,8)*self.trtr(t)
        e=E(0,0)
        t=self.tmul(self.tmul(self.tmul(self.tens(V3,V3),e),e),e)
        tot+=Fr(1,12)*self.trtr(t)
        return tot

def bs_F2ref(e1,e2,y,F):
    e1=Fr(e1); e2=Fr(e2)
    X=1/(1+27*Fr(y)); F=Fr(F); E=e1+e2; k=e1*e2
    f02=5*F**3/(8*X)+(5/(8*X)-Fr(1,2))*F**2+Fr(1,96)*(13*X+20/X-32)*F-(27*X**3-98*X**2+121*X-50)/(2160*X)
    f11=F**2/16+Fr(1,96)*(8-7*X)*F+(84*X**2-196*X+121)/5760
    f20=X/384*F-(21*X**2-29*X+8)/5760
    return E**4/k*f20+E**2*f11+k*f02
