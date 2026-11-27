# Fresh implementation of s.1 (t_c, rows of Def 1.3, edges (1.1), dilaton) over F_p, p = 2^61-1 (contains omega).
# Evidence only.
from fractions import Fraction as Fr
P=(1<<61)-1
def inv(a): return pow(a%P,P-2,P)
def fp(x):
    x=Fr(x); return x.numerator%P*inv(x.denominator)%P
def find_omega():
    e=(P-1)//3
    for g in range(2,100):
        w=pow(g,e,P)
        if w!=1: return w
OM=find_omega()

class L:
    """Laurent series in t_c: val v, coefs c (F_p); absolute precision v+len(c) (exclusive)"""
    __slots__=("v","c")
    def __init__(s,v,c): s.v=v; s.c=[x%P for x in c]
    @staticmethod
    def const(x,prec): return L(0,[x%P]+[0]*(prec-1)) if prec>0 else L(0,[])
    def prec(s): return s.v+len(s.c)
    def __add__(s,o):
        if not isinstance(o,L): o=L.const(fp(o) if isinstance(o,Fr) else o,s.prec())
        pr=min(s.prec(),o.prec()); v=min(s.v,o.v)
        c=[0]*(pr-v)
        for i,x in enumerate(s.c):
            if s.v+i<pr: c[s.v+i-v]+=x
        for i,x in enumerate(o.c):
            if o.v+i<pr: c[o.v+i-v]+=x
        return L(v,c)
    __radd__=__add__
    def __neg__(s): return L(s.v,[-x for x in s.c])
    def __sub__(s,o): return s+(-o if isinstance(o,L) else L.const(-fp(o) if isinstance(o,Fr) else -o,s.prec()))
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,L):
            k=fp(o) if isinstance(o,Fr) else o%P
            return L(s.v,[x*k for x in s.c])
        n=min(len(s.c),len(o.c))
        c=[0]*n
        for i in range(n):
            a=s.c[i]
            if a:
                for j in range(n-i):
                    c[i+j]+=a*o.c[j]
        return L(s.v+o.v,c)
    __rmul__=__mul__
    def norm(s):
        i=0
        while i<len(s.c) and s.c[i]==0: i+=1
        return L(s.v+i,s.c[i:])
    def inverse(s):
        t=s.norm()
        a0=inv(t.c[0]); n=len(t.c); b=[0]*n; b[0]=a0
        for k in range(1,n):
            acc=0
            for i in range(1,k+1): acc+=t.c[i]*b[k-i]
            b[k]=-acc%P*a0%P
        return L(-t.v,b)
    def __truediv__(s,o):
        if isinstance(o,L): return s*o.inverse()
        return s*inv(fp(o) if isinstance(o,Fr) else o)
    def __pow__(s,e):
        if e<0: return s.inverse()**(-e)
        r=L.const(1,len(s.c)+10)
        for _ in range(e): r=r*s
        return r
    def d(s):
        return L(s.v-1,[(s.v+i)*x for i,x in enumerate(s.c)])
    def coef(s,k):
        if k<s.v: return 0
        if k>=s.prec(): raise IndexError((k,s.prec()))
        return s.c[k-s.v]
    def __repr__(s): return "L(%d,%s..,prec %d)"%(s.v,s.c[:5],s.prec())

def series_setup(N):
    M=N+6
    wc=[Fr(1)]
    for n in range(0,M):
        wc.append(wc[-1]*(Fr(1,3)+n)*(Fr(2,3)+n)/((n+1)**2))
    dt=[sum(wc[:k+1]) for k in range(M+1)]
    tc=[Fr(0)]+[dt[k]/(k+1) for k in range(M)]
    def mul(a,b):
        c=[Fr(0)]*(M+1)
        for i,ai in enumerate(a):
            if ai:
                for j,bj in enumerate(b):
                    if i+j>M: break
                    c[i+j]+=ai*bj
        return c
    def compose(f,g):
        r=[Fr(0)]*(M+1)
        for co in reversed(f[:M+1]):
            r=mul(r,g); r[0]+=co
        return r
    d=[Fr(0),Fr(1)]+[Fr(0)]*(M-1)
    for k in range(2,M+1):
        e=compose(tc,d)
        d[k]-=e[k]
    def compose_t(f): return compose(f,d)
    return dict(M=M,wc=compose_t(wc),Delta=d,compose_t=compose_t)

def frL(lst,prec):
    return L(0,[fp(x) for x in lst[:prec]])

class Theory:
    def __init__(s,e1,e2,N=10,K=12,c1factor=None,rperturb=None):
        s.e1=Fr(e1); s.e2=Fr(e2); s.E=s.e1+s.e2; s.N=N; s.K=K
        E=s.E
        ss=series_setup(N+8)
        M=ss['M']; prec=N+8
        s.prec=prec
        s.Del=frL(ss['Delta'],prec)
        s.wc=frL(ss['wc'],prec)
        s.g=-s.wc
        b=[Fr(1)]
        for n in range(0,M):
            b.append(b[-1]*(-(Fr(1,3)-n))/(n+1))
        t_of_D=[-x/3 for x in b]
        tser=ss['compose_t'](t_of_D)
        s.t=frL(tser,prec)
        s.c1=-s.wc
        if c1factor is not None:
            f=[Fr(x) for x in c1factor]+[Fr(0)]*(M+1-len(c1factor))
            s.c1=s.c1*frL(ss['compose_t'](f),prec)
        s.C1i=s.c1.inverse()
        s.d2=(fp(E)*fp(Fr(1,2))%P)*((s.C1i*s.C1i-1)/s.Del)
        r=[Fr(1)]
        for m in range(0,K+2):
            nxt=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))/(27*(m+1))
            r.append(nxt)
        if rperturb:
            for l,val in rperturb.items(): r[l]=Fr(val)
        s.rFr=r
        s.r=[fp(x) for x in r]
        s._den={}
        s.rows={}; s.Delta={}; s.v={}; s.Y={}; s.dil={}
        for i in range(3):
            s.build(i)
        s.eta={(0,0):fp(s.e1*s.e2*E**3/9),(1,1):fp(-3*s.e1*s.e2*E),(1,2):fp(-3*s.e1*s.e2),(2,1):fp(-3*s.e1*s.e2)}
        s._ed={}
    def theta(s,f): return s.g*f.d()
    def den(s,Lb,n):
        key=(Lb,n)
        if key in s._den: return s._den[key]
        co=[0]*(s.K+1); co[0]=1
        Einv=inv(fp(s.E))
        for j in range(1,Lb+1):
            a=j*Einv%P
            for k in range(1,s.K+1):
                co[k]=(co[k]-a*co[k-1])%P
        for k in range(s.K+1): s._den[(Lb,k)]=co[k]
        return s._den[key]
    def build(s,i):
        E=fp(s.E); w=pow(OM,i,P)
        ti=s.t*w
        onep=ti*3+1
        Y=onep.inverse()*(-(E*inv(3)%P))
        v=onep/(ti*E)
        s.Y[i]=Y; s.v[i]=v
        s.Delta[i]=(v*v).inverse()*fp(3*s.e1*s.e2*s.E)
        vp=[L.const(1,s.prec)]
        for l in range(1,s.K+1): vp.append(vp[-1]*v)
        S=[]
        for n in range(s.K+1):
            acc=None
            for l in range(n+1):
                term=vp[l]*(s.r[l]*s.den(l+1,n-l)%P)
                acc=term if acc is None else acc+term
            S.append(acc)
        tl=s.theta(v)/v
        def Dop(F):
            out=[]
            for n in range(s.K+1):
                o=Y*F[n]
                if n>=1: o=o+tl*F[n-1]+s.theta(F[n-1])
                out.append(o)
            return out
        P1=[s.C1i*x for x in Dop(S)]
        DP1=Dop(P1)
        pref=s.Del*s.c1*s.c1
        P2=[pref*(DP1[n]-s.d2*P1[n]) for n in range(s.K+1)]
        s.rows[i]=[S,P1,P2]
        s.dil[i]=[None]+[ (S[d-1]*((-1)**d)) for d in range(1,s.K+2)]
    def pair(s,i,j,m,n):
        """[psi^m w^n] of -<P^(i)(-psi),P^(j)(-w)>"""
        sg=(-1)**(m+n)
        acc=None
        for (a,b),h in s.eta.items():
            t=s.rows[i][a][m]*s.rows[j][b][n]*h
            acc=t if acc is None else acc+t
        return acc*(-sg)
    def edge(s,i,j,m,n):
        key=(i,j,m,n)
        if key in s._ed: return s._ed[key]
        acc=None
        for k in range(n+1):
            t=s.pair(i,j,m+k+1,n-k)*((-1)**k)
            acc=t if acc is None else acc+t
        s._ed[key]=acc
        return acc
