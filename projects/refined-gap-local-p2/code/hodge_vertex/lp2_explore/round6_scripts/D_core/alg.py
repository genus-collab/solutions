# Exact algebra: polys over Q, rational functions in psi, Laurent series in t_c, Q(omega).
from fractions import Fraction as Fr

def ptrim(p):
    p=list(p)
    while p and p[-1]==0: p.pop()
    return p
def padd(a,b):
    n=max(len(a),len(b)); return ptrim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)])
def pneg(a): return [-x for x in a]
def pmul(a,b):
    if not a or not b: return []
    r=[Fr(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b): r[i+j]+=x*y
    return ptrim(r)
def pdivmod(a,b):
    a=list(a); q=[Fr(0)]*max(len(a)-len(b)+1,0)
    lb=b[-1]
    for i in range(len(a)-len(b),-1,-1):
        c=a[i+len(b)-1]/lb; q[i]=c
        if c!=0:
            for j,y in enumerate(b): a[i+j]-=c*y
    return ptrim(q),ptrim(a[:len(b)-1] if len(b)>1 else [])
def pgcd(a,b):
    a,b=ptrim(a),ptrim(b)
    while b:
        a,b=b,pdivmod(a,b)[1]
    if a:
        l=a[-1]; a=[x/l for x in a]
    return a
def peval(p,x):
    r=Fr(0)
    for c in reversed(p): r=r*x+c
    return r

class R:
    __slots__=('n','d')
    def __init__(s,n=None,d=None,reduce=True):
        n=ptrim([Fr(c) for c in (n or [])]); d=ptrim([Fr(c) for c in (d or [1])])
        if not n: s.n=[]; s.d=[Fr(1)]; return
        if reduce and len(d)>1:
            g=pgcd(n,d)
            if len(g)>1:
                n=pdivmod(n,g)[0]; d=pdivmod(d,g)[0]
        l=d[-1]
        s.n=[c/l for c in n]; s.d=[c/l for c in d]
    @staticmethod
    def lift(x):
        return x if isinstance(x,R) else R([x])
    def __add__(s,o):
        o=R.lift(o)
        if s.d==o.d: return R(padd(s.n,o.n),s.d)
        return R(padd(pmul(s.n,o.d),pmul(o.n,s.d)),pmul(s.d,o.d))
    __radd__=__add__
    def __neg__(s): return R(pneg(s.n),s.d,reduce=False)
    def __sub__(s,o): return s+(-R.lift(o))
    def __rsub__(s,o): return R.lift(o)-s
    def __mul__(s,o):
        o=R.lift(o)
        return R(pmul(s.n,o.n),pmul(s.d,o.d))
    __rmul__=__mul__
    def __truediv__(s,o):
        o=R.lift(o)
        return R(pmul(s.n,o.d),pmul(s.d,o.n))
    def __eq__(s,o):
        o=R.lift(o)
        return s.n==o.n and s.d==o.d
    def __ne__(s,o): return not s==o
    def __hash__(s): return hash(tuple(s.n))
    def is0(s): return not s.n
    def negarg(s):
        return R([c*(-1)**i for i,c in enumerate(s.n)],[c*(-1)**i for i,c in enumerate(s.d)])
    def taylor(s,N):
        # coefficients of psi^0..psi^N, requires d[0]!=0 (else error)
        if s.d[0]==0: raise ValueError('pole at 0')
        out=[]
        for k in range(N+1):
            c=(s.n[k] if k<len(s.n) else Fr(0))
            for i in range(1,min(k,len(s.d)-1)+1):
                c-=s.d[i]*out[k-i]
            out.append(c/s.d[0])
        return out
    def __repr__(s): return 'R(%s / %s)'%(s.n,s.d)

PSI=R([0,1])

def leaves(f,E,Jmax=60):
    """partial fractions in the basis 1/(E-j psi), j>=1.
    returns (poly coefficient list, {j:c_j}, ok) ; ok False if there are other poles, multiple poles or j outside 1..Jmax."""
    f=R.lift(f)
    q,r=pdivmod(f.n,f.d) if f.n else ([],[])
    out={}; den=list(f.d); recon=R(q)
    for j in list(range(-Jmax,0))+list(range(1,Jmax+1)):
        pt=E/j
        if len(f.d)>1 and peval(f.d,pt)==0:
            dj,_=pdivmod(f.d,[-pt,Fr(1)])
            if peval(dj,pt)==0: return q,out,False
            res=peval(r,pt)/peval(dj,pt)
            cj=-j*res
            out[j]=cj
            recon=recon+R([cj],[E,-j])
    if not (recon==f): return q,out,False
    return q,out,True

def in_Lam(f,E,m):
    """True iff f lies in span{lam_1..lam_m} (no poly part)."""
    q,o,ok=leaves(f,E)
    return ok and not q and all(1<=j<=m for j in o)

class W:
    """a+b omega, omega^2=-1-omega"""
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    @staticmethod
    def lift(x): return x if isinstance(x,W) else W(x)
    def __add__(s,o):
        o=W.lift(o); return W(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return W(-s.a,-s.b)
    def __sub__(s,o):
        o=W.lift(o); return W(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return W.lift(o)-s
    def __mul__(s,o):
        o=W.lift(o)
        return W(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return W((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o): return s*W.lift(o).inv()
    def __rtruediv__(s,o): return W.lift(o)*s.inv()
    def __eq__(s,o):
        o=W.lift(o); return s.a==o.a and s.b==o.b
    def __ne__(s,o): return not s==o
    def __hash__(s): return hash((s.a,s.b))
    def __repr__(s): return '(%s+%s w)'%(s.a,s.b)
OMEGA=W(0,1)

def iszero(c):
    if isinstance(c,R): return c.is0()
    return c==0

class LS:
    """Laurent series in t_c: coefficients dict, truncated: known for orders < prec."""
    __slots__=('c','prec')
    def __init__(s,c=None,prec=10**6):
        s.prec=prec
        s.c={k:v for k,v in (c or {}).items() if k<prec and not iszero(v)}
    def lo(s): return min(s.c) if s.c else s.prec
    def __add__(s,o):
        if not isinstance(o,LS): o=LS({0:o})
        p=min(s.prec,o.prec); d=dict(s.c)
        for k,v in o.c.items():
            if k<p: d[k]=d[k]+v if k in d else v
        return LS(d,p)
    __radd__=__add__
    def __neg__(s): return LS({k:-v for k,v in s.c.items()},s.prec)
    def __sub__(s,o): return s+(-o if isinstance(o,LS) else LS({0:-o}))
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,LS): return LS({k:v*o for k,v in s.c.items()},s.prec)
        p=min(s.lo()+o.prec,o.lo()+s.prec)
        d={}
        for i,a in s.c.items():
            for j,b in o.c.items():
                if i+j<p:
                    t=a*b
                    d[i+j]=d[i+j]+t if (i+j) in d else t
        return LS(d,p)
    def __rmul__(s,o): return LS({k:o*v for k,v in s.c.items()},s.prec)
    def deriv(s): return LS({k-1:v*k for k,v in s.c.items() if k!=0},s.prec-1)
    def inv(s):
        l=s.lo(); a0=s.c[l]; n=s.prec-l
        one=Fr(1)
        out={}
        for k in range(n):
            t=(one if k==0 else 0)
            for i in range(1,k+1):
                if (l+i) in s.c and k-i in out: t=t-s.c[l+i]*out[k-i]
            out[k]=(t/a0) if (not isinstance(a0,W)) else (W.lift(t)/a0)
        return LS({k-l:v for k,v in out.items()},n-l)
    def shift(s,m): return LS({k+m:v for k,v in s.c.items()},s.prec+m)
    def get(s,k):
        if k>=s.prec: raise IndexError('order %d beyond precision %d'%(k,s.prec))
        return s.c.get(k,0)
