# Rows of local P2 x C^2 at generic alpha, from the I-function restricted to P_j (A2.1-A2.6). Fresh code.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from fractions import Fraction as Fr
import math

def ser_mul(a,b,N,zero):
    r=[zero]*(N+1)
    for i,x in enumerate(a):
        if x==zero: continue
        for j,y in enumerate(b):
            if i+j>N: break
            r[i+j]=r[i+j]+x*y
    return r
def ser_inv(a,N,F):
    r=[F(0)]*(N+1); r[0]=F(1)/a[0]
    for n in range(1,N+1):
        s=F(0)
        for k in range(1,min(n,len(a)-1)+1):
            s=s+a[k]*r[n-k]
        r[n]=-s*r[0]
    return r

class P2Data:
    """Frobenius data and rows.  alpha, eps: field elements. D: y-order. K: z-order."""
    def __init__(s,F,alpha,eps,D,K):
        s.F=F;s.al=[F(a) for a in alpha];s.ep=[F(e) for e in eps];s.D=D;s.K=K
        s.t=Ctx(F,D,K+D+2)   # z-room for the Laurent stage
        s.tt=Ctx(F,D,K)
        s.build()
    def wt(s,i):
        a=s.al
        j,k=[l for l in range(3) if l!=i]
        return [3*a[i],a[j]-a[i],a[k]-a[i],-s.ep[0],-s.ep[1]]
    def W0(s,i):
        r=s.F(1)
        for w in s.wt(i): r=r*w
        return r
    def h(s,i,sign=1):
        return s.hv(s.wt(i),sign)
    def hv(s,v,sign=1):
        """h_v(z)^sign = exp(-sign Sum_k B_2k/(2k(2k-1)) p_{2k-1}(1/v) z^{2k-1})"""
        F=s.F;t=s.tt
        ex=Bi(t)
        k=1
        while 2*k-1<=t.K:
            b=F(bern(2*k)/(2*k*(2*k-1)))
            p=F(0)
            for a in v: p=p+F(1)/a**(2*k-1)
            ex=ex+Bi(t,{(0,2*k-1):-b*p*sign})
            k+=1
        return bexp(ex)
    def c_at_zero(s,d,j,N):
        """z^d c_{d,j}(z) as a power series in z through z^N (list); c_{d,j} the y^d coefficient of S(1)|P_j stripped of y^(-alpha_j/z)."""
        F=s.F;a=s.al
        num=[F(1)]
        for m in range(0,3*d):
            num=ser_mul(num,[3*a[j],F(-m)],N,F(0))
        den=[F(1)]
        for l in range(3):
            for m in range(1,d+1):
                if l==j:
                    continue
                den=ser_mul(den,[a[l]-a[j],F(m)],N,F(0))
        res=ser_mul(num,ser_inv(den,N,F),N,F(0))
        fac=F(1)/F(math.factorial(d))
        return [x*fac for x in res]
    def c_at_inf(s,d,j,N):
        """c_{d,j}(1/s) as a power series in s through s^N."""
        F=s.F;a=s.al
        num=[F(0),3*a[j]]
        for m in range(1,3*d):
            num=ser_mul(num,[F(-m),3*a[j]],N,F(0))
        den=[F(1)]
        for l in range(3):
            for m in range(1,d+1):
                den=ser_mul(den,[F(m),a[l]-a[j]],N,F(0))
        return ser_mul(num,ser_inv(den,N,F),N,F(0))
    def build(s):
        F=s.F;D=s.D;K=s.K;t=s.t;tt=s.tt;a=s.al
        s.Ay=[F(math.factorial(3*d)//math.factorial(d)**3*(-1)**d) for d in range(D+1)]
        # Frobenius data from the z=infinity expansion of S(1), S(H), S(H^2): independent of BS's D, E.
        Ns=6
        A=s.Ay
        Ainv=ser_inv(A,D,F)
        def smul_y(u,v):
            r=[[F(0)]*(Ns+1) for _ in range(D+1)]
            for d1,x in enumerate(u):
                for d2,y in enumerate(v):
                    if d1+d2>D: break
                    for n1,p in enumerate(x):
                        if p==0: continue
                        for n2,q in enumerate(y):
                            if n1+n2>Ns: break
                            r[d1+d2][n1+n2]=r[d1+d2][n1+n2]+p*q
            return r
        def yser(lst):
            return [[x]+[F(0)]*Ns for x in lst]
        def op(u,j):
            # (z theta_y - alpha_j) u, for u in y,s with z=1/s. Requires u[d][0]=0 for d>0.
            r=[[F(0)]*(Ns+1) for _ in range(D+1)]
            for d in range(D+1):
                for n in range(Ns+1):
                    r[d][n]=r[d][n]-a[j]*u[d][n]
                    if n+1<=Ns and d>0:
                        r[d][n]=r[d][n]+u[d][n+1]*d
                assert d==0 or u[d][0]==0
            return r
        S0={}
        for j in range(3):
            u=[[F(0)]*(Ns+1) for _ in range(D+1)]
            u[0][0]=F(1)
            for d in range(1,D+1):
                u[d]=s.c_at_inf(d,j,Ns)
            S0[j]=u
        S1={j:smul_y(yser(Ainv),op(S0[j],j)) for j in range(3)}
        for j in range(3):
            for d in range(D+1):
                exp=-a[j] if d==0 else F(0)
                assert S1[j][d][0]==exp,'S(H)=H+O(1/z) fails'
        # Note: S1's s-truncation drops one order at each op; only low orders are read below.
        nodes=[-x for x in a]
        V=[[nodes[j]**c for c in range(3)] for j in range(3)]
        Vi=inv_field(V,F(1))
        def interp(vals):
            return [[sum((Vi[c][j]*vals[j][d] for j in range(3)),F(0)) for d in range(D+1)] for c in range(3)]
        T1={j:smul_y(yser(Ainv),op(S1[j],j)) for j in range(3)}
        C11=interp({j:[T1[j][d][0] for d in range(D+1)] for j in range(3)})
        s.C11=C11
        C2i=ser_inv(C11[2],D,F)
        S2={}
        for j in range(3):
            u=[[T1[j][d][n]-sum((C11[1][e]*S1[j][d-e][n]+C11[0][e]*S0[j][d-e][n] for e in range(d+1)),F(0)) for n in range(Ns+1)] for d in range(D+1)]
            S2[j]=smul_y(yser(C2i),u)
            for d in range(D+1):
                exp=nodes[j]**2 if d==0 else F(0)
                assert S2[j][d][0]==exp,'S(H^2)=H^2+O(1/z) fails'
        T2={j:smul_y(yser(Ainv),op(S2[j],j)) for j in range(3)}
        C12=interp({j:[T2[j][d][0] for d in range(D+1)] for j in range(3)})
        s.C12=C12
        def yb(lst,tc): return Bi(tc,{(d,0):x for d,x in enumerate(lst) if x!=0})
        s.yb=yb
        def Cmat(tc):
            return [[const(tc,0),const(tc,1),const(tc,0)],
                    [yb(C11[0],tc),yb(C11[1],tc),yb(C11[2],tc)],
                    [yb(C12[0],tc),yb(C12[1],tc),yb(C12[2],tc)]]
        s.Cm=Cmat(tt)          # Cm[a][c]=C^c_{1a}
        s.Hs=mT(s.Cm)          # H* in flat coordinates: columns = images of H^a
        eta=[[sum((nodes[i]**(p+q)/s.W0(i) for i in range(3)),F(0)) for q in range(3)] for p in range(3)]
        s.eta_f=eta
        s.etai_f=inv_field(eta,F(1))
        s.eta=mconst(tt,eta)
        # Y_j, U_j
        s.Y={};s.U={}
        for j in range(3):
            Y=[F(0)]*(D+1); Y[0]=-a[j]
            def Lhat(Y):
                pr=[F(1)]+[F(0)]*D
                for m in range(3):
                    f=list(Y); f[0]=f[0]+a[m]
                    pr=ser_mul(pr,f,D,F(0))
                y3=ser_mul(ser_mul(Y,Y,D,F(0)),Y,D,F(0))
                return [pr[d]+(F(27)*y3[d-1] if d>0 else F(0)) for d in range(D+1)]
            dL0=F(1)
            for m in range(3):
                if m!=j: dL0=dL0*(a[m]-a[j])
            for d in range(1,D+1):
                L=Lhat(Y)
                Y[d]=-L[d]/dL0
            assert all(x==0 for x in Lhat(Y)),'Y root'
            s.Y[j]=Y
            U=[F(0)]*(D+1)
            for d in range(1,D+1): U[d]=Y[d]/F(d)
            s.U[j]=U
        # rR_{0j}=sqrt(W_j) R~_{0j}, from the I-function at z=0 with e^{-U/z}
        KK=t.K
        s.rR={}
        s.rawR0={}
        for j in range(3):
            c=Bi(t,{(0,0):F(1)})
            for d in range(1,D+1):
                lst=s.c_at_zero(d,j,KK+d)
                dd={(d,n-d):x for n,x in enumerate(lst) if x!=0 and n-d<=KK}
                c=c+Bi(t,dd)
            ex=Bi(t,{(d,-1):-s.U[j][d] for d in range(1,D+1) if s.U[j][d]!=0})
            r0=bexp(ex)*c
            r0=Bi(t,{k:v for k,v in r0.c.items() if k[1]<=K+1})
            assert all(n>=0 for (d,n) in r0.c),'z-pole in R~_0 at colour %d'%j
            r0=Bi(tt,{k:v for k,v in r0.c.items() if k[1]<=K})
            Ainvb=yb(Ainv,tt);Yb=yb(s.Y[j],tt)
            def raise1(f):
                return f.theta().zshift(1)+Yb*f
            r1=Ainvb*raise1(r0)
            De=yb(ser_mul(A,C11[2],D,F(0)),tt)      # D(y)=A C^2_11
            eE=yb(ser_mul(A,C11[1],D,F(0)),tt)      # e_1 E(y)=A C^1_11
            r2=binv(De)*(raise1(r1)-eE*r1)
            s.rR[j]=[r0,r1,r2]
        e1=a[0]+a[1]+a[2];e2=a[0]*a[1]+a[0]*a[2]+a[1]*a[2];e3=a[0]*a[1]*a[2]
        s.e=(e1,e2,e3)
        s.Delta={}
        s.DeltaHalfOverS={}
        for j in range(3):
            Yb=yb(s.Y[j],tt)
            Wb=const(tt,e1)*Yb*Yb+const(tt,2*e2)*Yb+const(tt,3*e3)
            Dl=const(tt,3*s.ep[0]*s.ep[1])*Wb
            s.Delta[j]=Dl
            assert Dl.get(0,0)==s.W0(j),'Delta(0)=W0'
            s.DeltaHalfOverS[j]=bsqrt1(Dl*(F(1)/s.W0(j)))
        s.A=yb(A,tt);s.Ainv=yb(Ainv,tt)
    def thetaQ(s,a):
        return s.Ainv*a.theta()
    def rows_rational(s,hsign=1):
        """rR_{aj} h_j^hsign = sqrt(W_j) R_{aj}; matrix [a][j]"""
        out=[[None]*3 for _ in range(3)]
        for j in range(3):
            hh=s.h(j,hsign) if hsign else const(s.tt,1)
            for a in range(3):
                out[a][j]=s.rR[j][a]*hh
        return out
    def Prows(s,hsign=1):
        """P_a^{(j)}(z)=Delta_j^{1/2} R_{aj}(z), rational: (Delta/W0)^{1/2} rR h"""
        rr=s.rows_rational(hsign)
        return [[rr[a][j]*s.DeltaHalfOverS[j] for j in range(3)] for a in range(3)]
