# Evidence only. Fresh code. Exact in Q(omega). Symplecticity, telescoped edge symmetry and S1 with c1, F free.
from fractions import Fraction as Fr
import random,sys
class W:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o):
        o=o if isinstance(o,W) else W(o); return W(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return W(-s.a,-s.b)
    def __sub__(s,o):
        o=o if isinstance(o,W) else W(o); return W(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return W(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,W) else W(o)
        return W(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return W((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o):
        o=o if isinstance(o,W) else W(o); return s*o.inv()
    def __rtruediv__(s,o): return W(o)*s.inv()
    def __pow__(s,n):
        if n<0: return (s.inv())**(-n)
        r=W(1)
        for _ in range(n): r=r*s
        return r
    def iszero(s): return s.a==0 and s.b==0
    def __repr__(s): return f"({s.a}+{s.b}w)"
om=W(0,1)
def rseq(N):
    r=[Fr(1)]
    for m in range(N):
        r.append((3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))/(27*(m+1)))
    return r
# Laurent polys in v
def ladd(a,b,c=1):
    r=dict(a)
    for k,x in b.items(): r[k]=r.get(k,Fr(0))+c*x
    return {k:x for k,x in r.items() if x!=0}
def lmul(a,b):
    r={}
    for k,x in a.items():
        for l,y in b.items(): r[k+l]=r.get(k+l,Fr(0))+x*y
    return {k:x for k,x in r.items() if x!=0}
def lth(a,E):
    r={}
    for e,x in a.items():
        if e==0: continue
        r[e-1]=r.get(e-1,Fr(0))+x*e/E
        r[e]=r.get(e,Fr(0))-x*e/3
    return {k:x for k,x in r.items() if x!=0}
def leval(a,v):
    t=W(0)
    for e,x in a.items(): t=t+(v**e)*x
    return t
def rows(E,N):
    r=rseq(N+2)
    S=[{} for _ in range(N+1)]
    for l in range(N+1):
        ser=[Fr(0)]*(N+1); ser[0]=Fr(1)
        for j in range(1,l+2):
            new=[Fr(0)]*(N+1)
            for n in range(N+1):
                new[n]=ser[n]-(Fr(j)/E*new[n-1] if n>=1 else 0)
            ser=new
        for n in range(N+1-l):
            if ser[n]!=0: S[n+l]=ladd(S[n+l],{l:r[l]*ser[n]})
    Y={-1:Fr(1),0:-E/3}
    tau={-1:-2/E,0:Fr(2,3)}
    def T(f):
        out=[{} for _ in range(N+1)]
        for n in range(N+1):
            out[n]=ladd(out[n],lmul(Y,f[n]))
            if n>=1:
                out[n]=ladd(out[n],lth(f[n-1],E))
                out[n]=ladd(out[n],lmul(tau,f[n-1]),Fr(-1,2))
        return out
    Q=T(S); TQ=T(Q)
    return S,Q,TQ
def run(e1,e2,t,c1,F,N):
    e1,e2=Fr(e1),Fr(e2); c1=c1 if isinstance(c1,W) else W(c1); F=F if isinstance(F,W) else W(F); E=e1+e2; t=W(t)
    S,Q,TQ=rows(E,N)
    y=t*t*t; Dc=1+27*y
    d2=W(E)*(c1**(-2)-1)/(2*Dc)
    eta00=e1*e2*E**3/9; eta11=-3*e1*e2*E; eta12=-3*e1*e2
    col=[]
    for i in range(3):
        ti=t*om**i; v=(1+3*ti)/(E*ti)
        s=[leval(S[n],v) for n in range(N+1)]
        q=[leval(Q[n],v) for n in range(N+1)]
        tq=[leval(TQ[n],v) for n in range(N+1)]
        p0=s; p1=[x/c1 for x in q]
        p2=[]
        for n in range(N+1):
            x=tq[n]-d2*q[n]
            if n>=1: x=x-F*q[n-1]
            p2.append(Dc*c1*x)
        col.append(dict(p=[p0,p1,p2],v=v,D=3*e1*e2*E/(v*v)))
    def pair(i,j,M,b):
        P=col[i]['p'];R=col[j]['p']
        return eta00*P[0][M]*R[0][b]+eta11*P[1][M]*R[1][b]+eta12*(P[1][M]*R[2][b]+P[2][M]*R[1][b])
    bad=0
    for i in range(3):
        for j in range(3):
            for n in range(N+1):
                s=W(0)
                for b in range(n+1):
                    s=s+((-1)**b)*pair(i,j,n-b,b)
                tgt=col[i]['D'] if (i==j and n==0) else W(0)
                if not (s-tgt).iszero(): bad+=1; print('sympl fail',i,j,n,s,tgt)
    # telescoped edge symmetry
    def n_(i,j,M,b): return ((-1)**(M+b))*pair(i,j,M,b)
    def Ed(i,j,a,b):
        r=W(0)
        for k in range(b+1): r=r-((-1)**k)*n_(i,j,a+k+1,b-k)
        return r
    M=N-1
    for i in range(3):
        for j in range(3):
            for a in range(M):
                for b in range(M-a):
                    if not (Ed(i,j,a,b)-Ed(j,i,b,a)).iszero(): bad+=1; print('edge sym fail',i,j,a,b)
    # (Pair) closed form with no c1
    for i in range(3):
        for j in range(3):
            for m in range(N+1):
                for b in range(N+1-m):
                    qi=[leval(Q[n],col[i]['v']) for n in range(N+1)]; qj=[leval(Q[n],col[j]['v']) for n in range(N+1)]
                    ti_=[leval(TQ[n],col[i]['v']) for n in range(N+1)]; tj_=[leval(TQ[n],col[j]['v']) for n in range(N+1)]
                    si=col[i]['p'][0]; sj=col[j]['p'][0]
                    rhs=eta00*si[m]*sj[b]-3*e1*e2*E*qi[m]*qj[b]-3*e1*e2*Dc*(qi[m]*tj_[b]+ti_[m]*qj[b])
                    if m>=1: rhs=rhs+3*e1*e2*Dc*F*qi[m-1]*qj[b]
                    if b>=1: rhs=rhs+3*e1*e2*Dc*F*qi[m]*qj[b-1]
                    if not (pair(i,j,m,b)-rhs).iszero(): bad+=1; print('pair fail')
        if N>8: break
    return bad
if __name__=='__main__':
    print('N=10 (3,7)',run(3,7,Fr(1,5),Fr(11,7),Fr(-2,9),10))
    print('N=10 (5/2,-3/7)',run(Fr(5,2),Fr(-3,7),Fr(3,11),W(Fr(2,3),Fr(1,5)),W(Fr(-1,7),Fr(4,3)),10))
    print('N=14 (1,2)',run(1,2,Fr(-2,7),Fr(5,3),Fr(7,5),14))
