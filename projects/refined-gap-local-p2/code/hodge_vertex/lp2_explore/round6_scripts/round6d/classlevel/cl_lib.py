# Fresh helper library for the round-6d class-level test.  Imports no project script.
from fractions import Fraction as Fr
import math

# ---------- fields ----------
class Fp:
    """Element of F_p; p set per subclass by make_fp."""
    __slots__=('v',)
    p=None
    def __init__(s,x=0):
        if isinstance(x,Fp): s.v=x.v
        elif isinstance(x,Fr): s.v=(x.numerator%s.p)*pow(x.denominator,-1,s.p)%s.p
        else: s.v=int(x)%s.p
    def _c(s,o): return o if isinstance(o,Fp) else type(s)(o)
    def __add__(s,o): return type(s)(s.v+s._c(o).v)
    __radd__=__add__
    def __sub__(s,o): return type(s)(s.v-s._c(o).v)
    def __rsub__(s,o): return type(s)(s._c(o).v-s.v)
    def __mul__(s,o): return type(s)(s.v*s._c(o).v)
    __rmul__=__mul__
    def __neg__(s): return type(s)(-s.v)
    def __truediv__(s,o):
        o=s._c(o)
        if o.v==0: raise ZeroDivisionError
        return type(s)(s.v*pow(o.v,-1,s.p))
    def __rtruediv__(s,o): return s._c(o)/s
    def __pow__(s,n): return type(s)(pow(s.v,n,s.p)) if n>=0 else (type(s)(1)/type(s)(pow(s.v,-n,s.p)))
    def __eq__(s,o): return s.v==s._c(o).v
    def __ne__(s,o): return s.v!=s._c(o).v
    def __hash__(s): return hash(s.v)
    def __repr__(s): return f"{s.v}"
    def __bool__(s): return s.v!=0

def make_fp(p):
    return type('F%d'%(p%100000),(Fp,),{'p':p,'__slots__':()})

def is_prime(n):
    if n<2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%q==0: return n==q
    d=n-1;s=0
    while d%2==0: d//=2;s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

def primes3mod4(start,count):
    n=start|3; out=[]
    while len(out)<count:
        if is_prime(n): out.append(n)
        n+=4
    return out

def fp_sqrt(x):
    p=type(x).p
    assert p%4==3
    r=type(x)(pow(x.v,(p+1)//4,p))
    return r if r*r==x else None

def rat_sqrt(x):
    x=Fr(x)
    if x<0: return None
    a=math.isqrt(x.numerator); b=math.isqrt(x.denominator)
    if a*a==x.numerator and b*b==x.denominator: return Fr(a,b)
    return None

def fsqrt(x):
    return rat_sqrt(x) if isinstance(x,Fr) else fp_sqrt(x)

# ---------- Bernoulli ----------
_B=[Fr(1)]
def bern(n):
    while len(_B)<=n:
        m=len(_B)
        _B.append(-sum(math.comb(m+1,k)*_B[k] for k in range(m))/(m+1))
    return _B[n]

# ---------- bi-series in y (or Q) and z: dict {(d,n):c}, truncated d<=D, n<=K ----------
class Ctx:
    def __init__(s,F,D,K):
        s.F=F;s.D=D;s.K=K;s.zero=F(0);s.one=F(1)

class Bi:
    __slots__=('c','t')
    def __init__(s,ctx,d=None):
        s.t=ctx; s.c=d if d is not None else {}
    def clean(s):
        z=s.t.zero
        s.c={k:v for k,v in s.c.items() if v!=z}
        return s
    def __add__(s,o):
        if not isinstance(o,Bi): o=const(s.t,o)
        r=dict(s.c)
        for k,v in o.c.items():
            r[k]=r[k]+v if k in r else v
        return Bi(s.t,r).clean()
    __radd__=__add__
    def __neg__(s): return Bi(s.t,{k:-v for k,v in s.c.items()})
    def __sub__(s,o):
        if not isinstance(o,Bi): o=const(s.t,o)
        return s+(-o)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,Bi):
            if o==s.t.zero: return Bi(s.t)
            return Bi(s.t,{k:v*o for k,v in s.c.items()}).clean()
        D,K=s.t.D,s.t.K
        r={}
        for (d1,n1),v1 in s.c.items():
            for (d2,n2),v2 in o.c.items():
                d=d1+d2;n=n1+n2
                if d>D or n>K: continue
                k=(d,n)
                r[k]=r[k]+v1*v2 if k in r else v1*v2
        return Bi(s.t,r).clean()
    __rmul__=__mul__
    def is_zero(s): return not any(v!=s.t.zero for v in s.c.values())
    def theta(s):  # y d/dy
        return Bi(s.t,{(d,n):v*d for (d,n),v in s.c.items() if d>0}).clean()
    def zshift(s,m):
        K=s.t.K
        return Bi(s.t,{(d,n+m):v for (d,n),v in s.c.items() if n+m<=K})
    def zneg(s):
        return Bi(s.t,{(d,n):(v if n%2==0 else -v) for (d,n),v in s.c.items()})
    def zcoef(s,n):
        """coefficient of z^n as Bi in y only"""
        return Bi(s.t,{(d,0):v for (d,m),v in s.c.items() if m==n})
    def get(s,d,n): return s.c.get((d,n),s.t.zero)
    def minz(s):
        return min((n for (d,n) in s.c),default=0)
    def trunc(s,K):
        return Bi(s.t,{k:v for k,v in s.c.items() if k[1]<=K})

def const(ctx,x):
    x=ctx.F(x)
    return Bi(ctx,{(0,0):x}) if x!=ctx.zero else Bi(ctx)
def ymono(ctx,d,x=1): return Bi(ctx,{(d,0):ctx.F(x)}) if d<=ctx.D else Bi(ctx)

def binv(a):
    """inverse of Bi with non-zero constant term (0,0)"""
    t=a.t
    a0=a.get(0,0)
    assert a0!=t.zero
    inv=const(t,t.one/a0)
    for _ in range(8):
        inv=inv*(const(t,2)-a*inv)
    return inv

def bsqrt1(a):
    """sqrt of a Bi with a(0,0)=1, sqrt(0,0)=1"""
    t=a.t
    assert a.get(0,0)==t.one
    r=const(t,1)
    half=t.one/t.F(2)
    for _ in range(8):
        r=(r+a*binv(r))*half
    return r

def bexp(a):
    """exp of Bi with zero (0,0) term, nilpotent (may have negative z powers but positive y powers)"""
    t=a.t
    r=const(t,1); term=const(t,1); m=1
    while True:
        term=term*a*(t.one/t.F(m))
        if term.is_zero(): break
        r=r+term; m+=1
        if m>400: raise RuntimeError
    return r

# ---------- matrices of Bi ----------
def mmul(A,B):
    n=len(A);m=len(B[0]);l=len(B)
    t=A[0][0].t
    out=[[Bi(t) for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s=Bi(t)
            for k in range(l):
                if A[i][k].c and B[k][j].c:
                    s=s+A[i][k]*B[k][j]
            out[i][j]=s
    return out
def madd(A,B): return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def msub(A,B): return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def mscal(A,x): return [[a*x for a in r] for r in A]
def mT(A): return [list(r) for r in zip(*A)]
def mid(t,n): return [[const(t,1 if i==j else 0) for j in range(n)] for i in range(n)]
def mconst(t,M): return [[const(t,x) for x in r] for r in M]
def mmap(f,A): return [[f(a) for a in r] for r in A]
def miszero(A): return all(a.is_zero() for r in A for a in r)

def inv_field(M,one):
    n=len(M)
    zero=one-one
    A=[list(r)+[one if i==j else zero for j in range(n)] for i,r in enumerate(M)]
    for c in range(n):
        p=next(i for i in range(c,n) if A[i][c]!=zero)
        A[c],A[p]=A[p],A[c]
        iv=one/A[c][c]
        A[c]=[x*iv for x in A[c]]
        for i in range(n):
            if i!=c and A[i][c]!=zero:
                f=A[i][c]
                A[i]=[x-f*y for x,y in zip(A[i],A[c])]
    return [r[n:] for r in A]

def minv(A):
    """inverse of a matrix over the local ring F[[y,z]]: Newton from the (0,0) value"""
    t=A[0][0].t
    n=len(A)
    A0=[[a.get(0,0) for a in r] for r in A]
    X=mconst(t,inv_field(A0,t.one))
    I2=mscal(mid(t,n),t.F(2))
    for _ in range(6):
        X=mmul(X,msub(I2,mmul(A,X)))
    return X

def mzero_report(A):
    """list of (i,j,d,n,value) of non-zero entries"""
    out=[]
    for i,r in enumerate(A):
        for j,a in enumerate(r):
            for (d,n),v in sorted(a.c.items()):
                if v!=a.t.zero: out.append((i,j,d,n,v))
    return out

# ---------- generic defect of A2.8 ----------
def defect_endo(EndB,Psi,eta,Hs,thetaQ):
    """Defect M(psi) of A2.8 from End B(psi) in the flat basis, given Psi (columns Psi_{aj}=eta(H^a,f_j)),
    eta, H* in flat coordinates (columns = images) and the derivation thetaQ.
    a^i(psi;H^a)= f_i^*(B(psi)H^a) = (Psi^-1 eta End B)_{ia}.
    theta a^i - f_i^*((1/psi)[H*,B]x) = M^i_k a^k.  Returns (M, commutator has no psi^0)."""
    t=EndB[0][0].t
    Psii=minv(Psi)
    A=mmul(mmul(Psii,eta),EndB)
    comm=msub(mmul(Hs,EndB),mmul(EndB,Hs))
    z0=all(a.zcoef(0).is_zero() for r in comm for a in r)
    commdiv=mmap(lambda a:Bi(t,{(d,n-1):v for (d,n),v in a.c.items() if n>=1}),comm)
    term=mmul(mmul(Psii,eta),commdiv)
    Kt=t.K-1
    diff=msub(mmap(thetaQ,A),term)
    diff=mmap(lambda a:a.trunc(Kt),diff)
    M=mmul(diff,minv(A))
    return mmap(lambda a:a.trunc(Kt),M),z0
