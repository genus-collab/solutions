from fractions import Fraction as F
class W:
    """Q(omega), omega^2=-1-omega"""
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=F(a);s.b=F(b)
    @staticmethod
    def _c(o): return o if isinstance(o,W) else W(o)
    def __add__(s,o): o=W._c(o);return W(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return W(-s.a,-s.b)
    def __sub__(s,o): o=W._c(o);return W(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return W._c(o)-s
    def __mul__(s,o):
        o=W._c(o)
        return W(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return W((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o): return s*W._c(o).inv()
    def __rtruediv__(s,o): return W._c(o)*s.inv()
    def __eq__(s,o): o=W._c(o);return s.a==o.a and s.b==o.b
    def __bool__(s): return s.a!=0 or s.b!=0
    def __repr__(s): return "(%s+%s w)"%(s.a,s.b)
class L:
    """Laurent series: val, coeffs list, absolute precision prec (exclusive)"""
    def __init__(s,val,c,prec):
        c=list(c)
        c=c[:max(prec-val,0)]
        i=0
        while i<len(c) and not c[i]: i+=1
        s.val=val+i;s.c=c[i:];s.prec=prec
        if not s.c: s.val=prec
    def coef(s,n):
        if n>=s.prec: raise ValueError("beyond prec %d>=%d"%(n,s.prec))
        i=n-s.val
        return s.c[i] if 0<=i<len(s.c) else F(0)
    def __add__(s,o):
        if not isinstance(o,L): o=L(0,[o],10**6)
        p=min(s.prec,o.prec);v=min(s.val,o.val)
        return L(v,[ s.coef(n)+o.coef(n) for n in range(v,p)],p)
    __radd__=__add__
    def __neg__(s): return L(s.val,[-x for x in s.c],s.prec)
    def __sub__(s,o): return s+(-o if isinstance(o,L) else -1*o)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,L): return L(s.val,[x*o for x in s.c],s.prec)
        p=min(s.prec+o.val,o.prec+s.val);v=s.val+o.val
        r=[F(0)]*max(p-v,0)
        for i,x in enumerate(s.c):
            if not x: continue
            for j,y in enumerate(o.c):
                if i+j<len(r): r[i+j]=r[i+j]+x*y
        return L(v,r,p)
    __rmul__=__mul__
    def inv(s):
        n=s.prec-s.val
        a0=s.c[0];i0=1/a0;r=[i0]
        for k in range(1,n):
            t=F(0)
            for i in range(1,k+1):
                if i<len(s.c): t=t+s.c[i]*r[k-i]
            r.append(-1*t*i0)
        return L(-s.val,r,-s.val+n)
    def __truediv__(s,o):
        if isinstance(o,L): return s*o.inv()
        return s*(1/o)
    def __rtruediv__(s,o): return s.inv()*o
    def __pow__(s,n):
        if n<0: return (s.inv())**(-n)
        r=L(0,[F(1)],10**6)
        for _ in range(n): r=r*s
        return r
    def d(s):
        return L(s.val-1,[(s.val+i)*x for i,x in enumerate(s.c)],s.prec-1)
    def __repr__(s): return "L(val=%d,prec=%d,%s)"%(s.val,s.prec,s.c[:6])
def compose(a,b):
    """a(x) power series list (coeffs from 0) composed with L b (val>=1)"""
    r=L(0,[F(0)],b.prec)
    for c in reversed(a):
        r=r*b+c
    return r
def rec(n):
    r=[F(1)]
    for m in range(n):
        r.append((3*(3*m*m+3*m+1)*r[m]-(m**3)*(r[m-1] if m>=1 else 0))/(27*(m+1)))
    return r
class Env:
    def __init__(s,e1,e2,N):
        s.e1=F(e1);s.e2=F(e2);s.E=s.e1+s.e2;s.N=N
        wc=[F(1)]
        for n in range(N+4):
            wc.append(wc[-1]*(F(1,3)+n)*(F(2,3)+n)/((n+1)**2))
        s.wc=wc
        pre=[F(0)]*(N+4)
        acc=F(0)
        for n in range(N+4):
            acc+=wc[n];pre[n]=acc
        tc=[F(0)]+[pre[n]/(n+1) for n in range(N+4)]
        s.tc_of_D=tc
        x=L(1,[F(1)],N+4)
        D=x
        for it in range(N+5):
            D=x-(compose(tc[2:],D)*D*D)
        s.D=D;s.x=x
        binom=[F(1)]
        for n in range(N+4):
            binom.append(binom[-1]*(F(1,3)-n)/(n+1)*(-1))
        one_third=compose(binom,D)
        s.s=1-one_third
        s.t=(s.s-1)*F(1,3)
        s.Wc=compose(wc,D)
        s.g=-1*s.Wc
        s.theta=lambda f: s.g*f.d()
    def v_of(s,col=0):
        if col==0:
            t=s.t
        else:
            om=W(0,1) if col==1 else W(-1,-1)
            t=L(s.t.val,[om*c for c in s.t.c],s.t.prec)
        return (1+3*t)/(s.E*t)
    def rows(s,z,col=0,c1=None):
        """normalised rows (Phat0,Phat1,Phat2) at numeric z, Def 1.3"""
        E=s.E
        if c1 is None: c1=s.g
        C1i=c1.inv()
        v=s.v_of(col)
        Y=1/v-E/3
        S=s.S(v,z)
        th=s.theta
        Ds=s.D
        d2=E*(C1i*C1i-1)/(2*Ds)
        P0B=v*S
        Q=(Y*P0B+z*th(P0B))*C1i
        P2=(Y*Q+z*th(Q)-d2*Q)*Ds*c1*c1
        return (S,Q/v,P2/v)
    def S(s,v,z,m=0):
        E=s.E
        rr=rec(s.N+14)
        vp=L(0,[F(1)],10**6)
        out=L(0,[F(0)],10**6)
        prec=s.N+3
        for l in range(0,len(rr)):
            den=F(1)
            for j in range(1,l+2-m):
                den*=1+j*z/E
            term=vp*(rr[l]*z**l/den)
            out=out+term
            vp=vp*v
            if v.val>=1 and vp.val>=prec+4 : break
        if v.val<1:
            pass
        return out
