# Evidence only (exact Fractions). Tests of FB:A3.5, draft Lemma 2.3, Lemma 2.5(i) at finite order.
import sys
from fractions import Fraction as Fr
from math import factorial as fa
from itertools import product
N=3      # y-order
NZ=int(sys.argv[1]) if len(sys.argv)>1 else 3  # z-order
class S:
    __slots__=('c',)
    def __init__(s,c=None):
        s.c=c if c is not None else [Fr(0)]*(N+1)
    @staticmethod
    def k(x):
        r=S(); r.c=[Fr(x)]+[Fr(0)]*N; return r
    def __add__(s,o):
        o=o if isinstance(o,S) else S.k(o)
        return S([a+b for a,b in zip(s.c,o.c)])
    __radd__=__add__
    def __neg__(s): return S([-a for a in s.c])
    def __sub__(s,o):
        o=o if isinstance(o,S) else S.k(o)
        return S([a-b for a,b in zip(s.c,o.c)])
    def __rsub__(s,o): return S.k(o)-s
    def __mul__(s,o):
        if not isinstance(o,S):
            o=Fr(o); return S([a*o for a in s.c])
        r=[Fr(0)]*(N+1)
        for i,a in enumerate(s.c):
            if a==0: continue
            for j,b in enumerate(o.c):
                if i+j>N: break
                if b!=0: r[i+j]+=a*b
        return S(r)
    __rmul__=__mul__
    def inv(s):
        r=[Fr(0)]*(N+1); r[0]=1/s.c[0]
        for n in range(1,N+1):
            r[n]=-sum(s.c[k]*r[n-k] for k in range(1,n+1))/s.c[0]
        return S(r)
    def th(s): return S([i*a for i,a in enumerate(s.c)])
    def sqrt1(s): # sqrt of 1+O(y)
        assert s.c[0]==1
        r=[Fr(0)]*(N+1); r[0]=Fr(1)
        for n in range(1,N+1):
            t=s.c[n]-sum(r[k]*r[n-k] for k in range(1,n))
            r[n]=t/2
        return S(r)
Yv=S(); Yv.c=[Fr(0),Fr(1)]+[Fr(0)]*(N-1)  # y

def root(al,i):
    # Newton in y-series for Lhat(Y)=prod(Y+al_j)+27 y Y^3, Y(0)=-al_i
    Y=S.k(-al[i])
    for _ in range(N+2):
        L=(Y+al[0])*(Y+al[1])*(Y+al[2])+27*Yv*Y*Y*Y
        e1=al[0]+al[1]+al[2];e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2]
        dL=3*Y*Y+2*e1*Y+e2+81*Yv*Y*Y
        Y=Y-L*dL.inv()
    return Y
def Wf(al,Y):
    e1=al[0]+al[1]+al[2];e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2];e3=al[0]*al[1]*al[2]
    return e1*Y*Y+2*e2*Y+3*e3
# polynomials in xi: dict exps->S, stored by degree
def pmul(p,q):
    r={}
    for a,x in p.items():
        for b,y in q.items():
            k=(a[0]+b[0],a[1]+b[1],a[2]+b[2])
            r[k]=r.get(k,S())+x*y
    return r
def multinom(m,a):
    return Fr(fa(m),fa(a[0])*fa(a[1])*fa(a[2]))
def Tr(al,Y,r):
    d={}
    for j in range(3):
        e=[0,0,0];e[j]=r
        d[tuple(e)]=d.get(tuple(e),S())+(Y+al[j])*Fr(1,fa(r))
    c=(-3*Y)*Fr(1,fa(r))*Fr(1,3**r)
    for a in product(range(r+1),repeat=3):
        if sum(a)==r:
            d[a]=d.get(a,S())+c*multinom(r,a)
    return d
def lap(p,G):
    # (1/2) G^{ab} d_a d_b
    r={}
    for e,x in p.items():
        for a in range(3):
            for b in range(3):
                f=list(e); co=f[a]
                if co==0: continue
                f[a]-=1
                co2=f[b]
                if co2==0: continue
                f[b]-=1
                k=tuple(f)
                r[k]=r.get(k,S())+x*G[a][b]*Fr(co*co2,2)
    return r
def Phi(al,Y,nz):
    # H, G=-H^{-1}
    q=[Y+a for a in al]
    H=[[ (q[a] if a==b else S.k(0)) - Y*Fr(1,3) for b in range(3)] for a in range(3)]
    def cof(a,b):
        m=[[H[i][j] for j in range(3) if j!=b] for i in range(3) if i!=a]
        v=m[0][0]*m[1][1]-m[0][1]*m[1][0]
        return v if (a+b)%2==0 else -v
    det=H[0][0]*cof(0,0)+H[0][1]*cof(0,1)+H[0][2]*cof(0,2)
    W=Wf(al,Y)
    d=det-W*Fr(1,3)
    assert all(x==0 for x in d.c), "det H != W/3"
    di=det.inv()
    G=[[-(cof(b,a)*di) for b in range(3)] for a in range(3)]
    Ts={r:Tr(al,Y,r) for r in range(3,2*nz+4)}
    # powers of T by degree: pw[k][deg]=poly
    out=[S.k(1)]
    pw={0:{0:{(0,0,0):S.k(1)}}}
    maxk=2*nz
    for k in range(1,maxk+1):
        pw[k]={}
        for deg,pl in pw[k-1].items():
            for r,tr in Ts.items():
                D=deg+r
                if D>2*nz+2*k: continue
                t=pw[k].setdefault(D,{})
                for e,x in pmul(pl,tr).items():
                    t[e]=t.get(e,S())+x
    for n in range(1,nz+1):
        tot=S()
        for k in range(0,2*n+1):
            m=n+k
            D=2*m
            if D not in pw[k]: continue
            p=pw[k][D]
            for _ in range(m): p=lap(p,G)
            v=p.get((0,0,0),S())
            tot=tot+v*Fr(1,fa(m)*fa(k))
        out.append(tot)
    return out

def theta_series(order):
    # Theta(x)=exp(sum B_{2m} x^{2m-1}/(2m(2m-1))); returns coefficients in x
    B={2:Fr(1,6),4:Fr(-1,30),6:Fr(1,42),8:Fr(-1,30),10:Fr(5,66)}
    l=[Fr(0)]*(order+1)
    for m in range(1,(order+1)//2+1):
        if 2*m-1<=order: l[2*m-1]=B[2*m]/(2*m*(2*m-1))
    r=[Fr(0)]*(order+1); r[0]=Fr(1)
    for n in range(1,order+1):
        r[n]=sum(k*l[k]*r[n-k] for k in range(1,n+1))/n
    return r

def zmul(a,b,M):
    r=[Fr(0)]*(M+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>M: break
            r[i+j]+=x*y
    return r
def zinv(a,M):
    r=[Fr(0)]*(M+1); r[0]=1/a[0]
    for n in range(1,M+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,min(n,len(a)-1)+1))/a[0]
    return r

def hyper_side(al,i,Y,nz):
    # e^{-U/z} sum_d y^d c_d(z), theta U = Y+al_i; returns [y^d][z^n], n<=nz; asserts that poles cancel
    M=nz+N+2
    a=al[i]
    others=[al[l] for l in range(3) if l!=i]
    cz=[]  # z^d c_d(z) as power series in z
    for d in range(0,N+1):
        num=[Fr(1)]+[Fr(0)]*M
        for m in range(3*d):
            num=zmul(num,[3*a,Fr(-m)],M)
        den=[Fr(fa(d))]+[Fr(0)]*M
        for b in others:
            for m in range(1,d+1):
                den=zmul(den,[b-a,Fr(m)],M)
        cz.append(zmul(num,zinv(den,M),M))
    U=[Fr(0)]+[ (Y+a).c[k]/k for k in range(1,N+1)]
    E=[{0:Fr(1)}]+[{} for _ in range(N)]
    for k in range(1,N+1):
        for j in range(1,k+1):
            for p,x in E[k-j].items():
                E[k][p+1]=E[k].get(p+1,Fr(0))+Fr(j,k)*(-U[j])*x
    res=[[Fr(0)]*(nz+1) for _ in range(N+1)]
    for tot in range(N+1):
        lau={}
        for d in range(tot+1):
            for p,x in E[tot-d].items():
                for n,c in enumerate(cz[d]):
                    if c==0: continue
                    e=n-d-p
                    if e>nz: continue
                    lau[e]=lau.get(e,Fr(0))+x*c
        for e,v in lau.items():
            if e<0: assert v==0, ("pole not cancelled",tot,e,v)
            else: res[tot][e]=v
    return res

def run(al):
    for i in range(3):
        j,k=[l for l in range(3) if l!=i]
        a=al[j]-al[i]; b=al[k]-al[i]
        Y=root(al,i)
        assert Y.c[0]==-al[i]
        W=Wf(al,Y)
        assert W.c[0]==al[i]*a*b, "W(-alpha_i)"
        # Lemma 2.5(i) on a test function g
        g=S([Fr(1),Fr(2),Fr(-3),Fr(5)])
        def T(f):
            r=[S() for _ in range(len(f)+1)]
            for n,x in enumerate(f):
                r[n]=r[n]+Y*x
                r[n+1]=r[n+1]+x.th()
            return r
        def add(f,h,cz=1):
            L=max(len(f),len(h)); r=[S() for _ in range(L)]
            for n,x in enumerate(f): r[n]=r[n]+x
            for n,x in enumerate(h): r[n]=r[n]+x*cz
            return r
        def shift(f): return [S()]+f
        def prodal(f):
            r=f
            for l in range(3):
                r=add(T(r),r,al[l])
            return r
        def ymult(f): return [Yv*x for x in f]
        def pyterm(f):
            r=add([3*x for x in T(f)],shift(f),2)   # (3T+2z) f
            r=add([3*x for x in T(r)],shift(r),1)   # (3T+z)
            r=[3*x for x in T(r)]
            return ymult(r)
        def Pop(f): return add(prodal(f),pyterm(f))
        Pg=Pop([g])
        assert all(x==0 for x in Pg[0].c), "P_0"
        wy=(W*Y.inv())
        expect=-(wy*g.th()+ (W.th()*Y.inv()*g)*Fr(1,2))
        assert all(x==0 for x in (Pg[1]-expect).c), "P_1 formula"
        # Phi
        Ph=Phi(al,Y,NZ)
        Ph0=[x.c[0] for x in Ph]
        M=NZ
        prod=[Fr(1)]+[Fr(0)]*M
        th=theta_series(M)
        for aa in (3*al[i],a,b):
            ser=[th[n]*(Fr(-1)/aa)**n for n in range(M+1)]
            prod=zmul(prod,ser,M)
        assert Ph0==prod, ("value at y=0",Ph0,prod)
        Wn=W*(1/W.c[0])
        sq=Wn.sqrt1().inv()
        Phinv=zinv(Ph0,M)
        lhs=[S() for _ in range(M+1)]
        for n in range(M+1):
            for l in range(n+1):
                lhs[n]=lhs[n]+Ph[l]*Phinv[n-l]
            lhs[n]=lhs[n]*sq
        rhs=hyper_side(al,i,Y,NZ)
        for tot in range(N+1):
            for n in range(NZ+1):
                assert lhs[n].c[tot]==rhs[tot][n], ("mismatch",i,tot,n,lhs[n].c[tot],rhs[tot][n])
        print("colour",i,"ok: W(-a_i), P_0, P_1, Phi(-a_i)=prod Theta, kernel ratio through y^%d z^%d"%(N,NZ))
        print("   sample [y^1 z^1]",lhs[1].c[1],"[y^2 z^%d]"%min(2,NZ),lhs[min(2,NZ)].c[2])
run((Fr(1),Fr(2),Fr(5)))
run((Fr(3,2),Fr(-2,7),Fr(5,3)))
print("9/1451520 == 1/161280:",Fr(9,1451520)==Fr(1,161280))
