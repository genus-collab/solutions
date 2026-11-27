# Evidence only. Stored F20, F30 (polynomials in XX over Q(y,eps)), F40, at XX = F_cf = theta log wc: Laurent series in xi = 1+27y,
# re-expanded in t_c; compares with (-3)^(g-1) beta_(2g-2)/t_c^(2g-2). Also t_cf numerically vs t_c, and hand-formula tests.
import os
import sys, math, re
from fractions import Fraction as Fr
from fc import F20F30, load_mma
from hodge_lib import BN

PREC=26
class L:
    """Laurent series in xi: xi^val * sum c[k] xi^k, truncated at relative length PREC"""
    __slots__=('v','c')
    def __init__(s,v,c):
        c=list(c)[:PREC]; c+= [Fr(0)]*(PREC-len(c))
        s.v=v; s.c=c
    @staticmethod
    def num(x): return L(0,[Fr(x)])
    def norm(s):
        k=0
        while k<PREC and s.c[k]==0: k+=1
        if k==PREC: return L(s.v+PREC,[Fr(0)])
        if k==0: return s
        return L(s.v+k,s.c[k:])
    def __add__(s,o):
        if not isinstance(o,L): o=L.num(o)
        v=min(s.v,o.v); n=PREC
        c=[Fr(0)]*n
        for k in range(n):
            if 0<=k-(s.v-v)<PREC: c[k]+=s.c[k-(s.v-v)]
            if 0<=k-(o.v-v)<PREC: c[k]+=o.c[k-(o.v-v)]
        return L(v,c)
    __radd__=__add__
    def __neg__(s): return L(s.v,[-x for x in s.c])
    def __sub__(s,o): return s+(-o if isinstance(o,L) else -Fr(o))
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,L): return L(s.v,[x*o for x in s.c])
        c=[Fr(0)]*PREC
        for i,x in enumerate(s.c):
            if x==0: continue
            for j,y in enumerate(o.c):
                if i+j>=PREC: break
                c[i+j]+=x*y
        return L(s.v+o.v,c)
    __rmul__=__mul__
    def inv(s):
        s=s.norm()
        assert s.c[0]!=0
        b=[Fr(0)]*PREC; b[0]=1/s.c[0]
        for n in range(1,PREC):
            b[n]=-sum(s.c[k]*b[n-k] for k in range(1,n+1))/s.c[0]
        return L(-s.v,b)
    def __truediv__(s,o):
        if not isinstance(o,L): return L(s.v,[x/o for x in s.c])
        return s*o.inv()
    def __rtruediv__(s,o): return s.inv()*o
    def __pow__(s,n):
        n=int(n)
        if n<0: return s.inv()**(-n)
        r=L.num(1)
        x=s
        while n:
            if n&1: r=r*x
            x=x*x; n>>=1
        return r
    def d(s):
        # d/dxi
        return L(s.v-1,[(s.v+k)*x for k,x in enumerate(s.c)])
    def coef(s,e):
        k=e-s.v
        if k<0: return Fr(0)
        if k>=PREC: return None
        return s.c[k]

def series_data():
    w=[Fr(1)]
    for n in range(0,PREC+2):
        w.append(w[-1]*(Fr(1,3)+n)*(Fr(2,3)+n)/(n+1)**2)
    wc=L(0,w)
    xi=L(1,[Fr(1)])
    # dt/dxi = wc/(1-xi)
    dt=wc*(1/(1-xi)) if False else wc*(L.num(1)-xi).inv()
    tc=L(1,[dt.c[k]/(k+1) for k in range(PREC)])
    XX=(xi-1)*wc.d()*wc.inv()
    y=(xi-1)*Fr(1,27)
    return wc,tc,XX,y,xi

def beta(n,e1,e2):
    # (n-1)! [t^n] 1/(4 sinh(e1 t/2) sinh(e2 t/2))
    def cosech(a,N):
        # 1/sinh(a t) = (1/(a t)) * sum_k c_k (a t)^(2k), series x/sinh x
        # x/sinh x = sum_k -(2^(2k)-2) B_2k x^(2k)/(2k)!
        out={}
        for k in range(0,N+1):
            out[2*k-1]=-(2**(2*k)-2)*BN[2*k]/math.factorial(2*k)*a**(2*k-1)
        return out
    A=cosech(e1/2,n+2); B=cosech(e2/2,n+2)
    tot=Fr(0)
    for i,x in A.items():
        if n-i in B: tot+=x*B[n-i]
    return math.factorial(n-1)*tot/4

def G_of(expr,e1,e2,XX,y):
    return eval(expr,{},dict(e1=Fr(e1),e2=Fr(e2),XX=XX,y=y,Fr=Fr))

def F40expr():
    return load_mma(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_F40.m"))

if __name__=="__main__":
    wc,tc,XX,y,xi=series_data()
    print("t_c series",[str(x) for x in tc.c[:4]],"XX(0)=",XX.c[0],"XX'(0)=",XX.c[1],"wc'(0)=",wc.c[1])
    f2,f3=F20F30(); f4=F40expr()
    print("beta_2(3,7) vs hand", beta(2,Fr(3),Fr(7)), (7*81+10*9*49+7*2401)/Fr(5760*21))
    samples=[(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7)),(Fr(1,3),Fr(1)),(Fr(4,5),Fr(1))]
    for e1,e2 in samples:
        for g,ex in ((2,f2),(3,f3),(4,f4)):
            G=G_of(ex,e1,e2,XX,y)
            target=(tc**(2-2*g))*((-3)**(g-1)*beta(2*g-2,e1,e2))
            diff=G-target
            polar=[diff.coef(-j) for j in range(2*g-2,0,-1)]
            regular=[diff.coef(j) for j in range(0,2)]
            # leading coefficient of G in t_c
            lead=(G*(tc**(2*g-2))).coef(0)
            # control: wrong propagator, XX -> XX + 1/9, and XX -> theta log A(large radius) replaced by XX + xi/5
            Gc=G_of(ex,e1,e2,XX+Fr(1,9),y); dc=Gc-target
            Gc2=G_of(ex,e1,e2,XX+xi*Fr(1,5),y); dc2=Gc2-target
            print(f"g={g} eps=({e1},{e2}): lead={lead}  (-3)^(g-1)beta={(-3)**(g-1)*beta(2*g-2,e1,e2)}  polar(G-target)={[str(p) for p in polar]}  next={regular[0]}")
            print("   controls polar nonzero:",any(dc.coef(-j)!=0 for j in range(2*g-2,0,-1)),any(dc2.coef(-j)!=0 for j in range(2*g-2,0,-1)))
    # t_cf constant: numerical comparison at y=-1/54 and y=-1/81
    for yv in (-1/54.0,-1/81.0,-1/270.0):
        Dc=1+27*yv
        # t_c by Taylor in Dc
        w=[1.0]
        for n in range(0,6000):
            w.append(w[-1]*(1/3+n)*(2/3+n)/(n+1)**2)
        # dt/dDc = wc/(1-Dc) : integrate coefficient series: partial sums
        tc_val=0.0; s=0.0
        for k in range(0,6000):
            s+=w[k]            # coefficient of Dc^k in wc/(1-Dc)
            tc_val+=s*Dc**(k+1)/(k+1)
        # BS expression (real part form) with l=log(-y)
        l=math.log(-yv)
        g=g2=0.0
        term=1.0
        fac=1.0
        H3=0.0;H1=0.0
        for d in range(1,4000):
            # (3d-1)!/(d!)^3 (-y)^d  via recurrence
            if d==1: c=2.0
            else: c=c*(3*d-1)*(3*d-2)*(3*d-3)/d**3
            pw=(-yv)**d
            # psi(3d)-psi(d+1) = H_{3d-1}-H_d
            H3=sum(1.0/i for i in range(1,3*d)) if d<60 else H3+1/(3*d-1)+1/(3*d-2)+1/(3*d-3)
            H1=sum(1.0/i for i in range(1,d+1)) if d<60 else H1+1/d
            g+=3*c*pw
            g2+=9*c*(H3-H1)*pw
            if pw*c<1e-18: break
        val=math.sqrt(3)/(2*math.pi)*(l*l/2+g*l+g2-math.pi**2/2)
        print(f"y={yv:.6f} t_c={tc_val:.12f} BS t_cf={val:.12f} diff={tc_val-val:.2e}; wrong branch would add {math.sqrt(3)/(2*math.pi)*(2*math.pi**2):.3f}")
