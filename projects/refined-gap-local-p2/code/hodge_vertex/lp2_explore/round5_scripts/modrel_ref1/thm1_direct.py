# Referee script: read I-function DIRECTLY from eq:Ifun (no harmonic-number closed form), generic rational alpha.
# Truncated polys in (H,s) over Q; L=log y tracked as polynomial variable.  Then A,D,E,X from eq:ADEFXdef, test modrel.
from fractions import Fraction as Fr
import sys
Dmax=int(sys.argv[1]) if len(sys.argv)>1 else 24
def pmul(a,b,smax=3,Hmax=4):
    r={}
    for (h1,s1),x in a.items():
        for (h2,s2),y in b.items():
            if s1+s2<=smax and h1+h2<=Hmax:
                k=(h1+h2,s1+s2); r[k]=r.get(k,0)+x*y
    return {k:v for k,v in r.items() if v!=0}
def lin_inv(c,hc,sc,smax=3,Hmax=4):
    # inverse of c + (hc H + a) s form given as dict: c + u, u nilpotent in s
    pass
def inv_factor(m,alpha):
    # 1/(m+(H+alpha)s) = (1/m) sum_k (-(H+alpha)s/m)^k
    base={(0,0):Fr(1,m)}
    res=dict(base); term=dict(base)
    u={(1,1):Fr(-1,m),(0,1):Fr(-alpha,1)/m}
    for k in range(1,4):
        term=pmul(term,u)
        for kk,v in term.items(): res[kk]=res.get(kk,0)+v
    return res
def coeffs(alpha):
    terms={0:{(0,0):Fr(1)}}
    for d in range(1,Dmax+1):
        cur={(0,0):Fr(1)}
        for m in range(0,3*d):
            cur=pmul(cur,{(0,0):Fr(-m),(1,1):Fr(-3)} if m>0 else {(1,1):Fr(-3)})
        for j in range(3):
            for m in range(1,d+1):
                cur=pmul(cur,inv_factor(m,alpha[j]))
        terms[d]=cur
    return terms
def analyse(alpha,control=None):
    t=coeffs(alpha)
    e1=sum(alpha)
    # I = z exp(H s L) sum_d y^d term_d;  [H^a z^-n] I = [H^a s^(n+1)] exp(HsL) sum
    # series coefficient as function of L: exp(HsL)=sum (HsL)^k/k!
    from math import factorial
    def Icoef(a,n):
        # return dict Lpow-> list of y coeffs
        out={}
        for k in range(0,a+1):
            lp=Fr(1,factorial(k))
            ser=[Fr(0)]*(Dmax+1)
            for d in range(Dmax+1):
                ser[d]=t[d].get((a-k,n+1-k),Fr(0))*lp if n+1-k>=0 else Fr(0)
            if any(ser): out[k]=ser
        return out
    I01=Icoef(1,0); I11=Icoef(1,1); I21=Icoef(2,1)
    # check against eq:Ian shape
    def theta(c):
        out={}
        for k,ser in c.items():
            o=out.setdefault(k,[Fr(0)]*(Dmax+1))
            for d in range(Dmax+1): o[d]+=d*ser[d]
            if k>=1:
                o=out.setdefault(k-1,[Fr(0)]*(Dmax+1))
                for d in range(Dmax+1): o[d]+=k*ser[d]
        return {k:v for k,v in out.items() if any(v)}
    def smul(a,b):
        r=[Fr(0)]*(Dmax+1)
        for i,x in enumerate(a):
            if x:
                for j in range(Dmax+1-i): r[i+j]+=x*b[j]
        return r
    def sinv(a):
        r=[Fr(0)]*(Dmax+1); r[0]=1/a[0]
        for n in range(1,Dmax+1):
            r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))*r[0]
        return r
    A=theta(I01)
    assert set(A.keys())=={0},A.keys()
    A=A[0]
    Ai=sinv(A)
    tI11=theta(I11); tI21=theta(I21)
    def divA(c): return {k:smul(v,Ai) for k,v in c.items()}
    Edef=theta(divA(tI11))
    Ddef=theta(divA(tI21))
    assert set(Ddef.keys())<={0},Ddef.keys()
    assert set(Edef.keys())<={0},Edef.keys()
    D=Ddef[0]; Ee=[x/e1 for x in Edef[0]]
    X=[Fr(0)]*(Dmax+1)
    for d in range(Dmax+1): X[d]=Fr((-27)**d)
    r1=[a-b for a,b in zip(smul(smul(A,A),D),X)]
    r2=[2*a-b+c for a,b,c in zip(Ee,D,X)]
    # eq:DEXrel
    # 2/e1 theta I11 - theta I21 + A log(y/(1+27y))
    log1=[Fr(0)]*(Dmax+1)
    for n in range(1,Dmax+1): log1[n]=Fr((-1)**(n-1)*27**n,n)
    t11=tI11.get(0,[Fr(0)]*(Dmax+1))
    t21=tI21
    c0=[2*a/e1-b-c for a,b,c in zip(t11,t21.get(0,[Fr(0)]*(Dmax+1)),smul(A,log1))]
    c1=[-b+a for a,b in zip(A,t21.get(1,[Fr(0)]*(Dmax+1)))]
    return A,D,Ee,all(x==0 for x in r1),all(x==0 for x in r2),all(x==0 for x in c0),all(x==0 for x in c1),[ (i,x) for i,x in enumerate(r2) if x!=0][:2]
for alpha in [(Fr(3,7),Fr(11,5),Fr(-2,9)),(Fr(1),Fr(1),Fr(1)),(Fr(1,3),Fr(1,3),Fr(1,3)),(Fr(5),Fr(-13,2),Fr(17,3))]:
    A,D,Ee,a,b,c,d,w=analyse(alpha)
    print(alpha,'A2D=X',a,'2E=D-X',b,'DEXrel const',c,'DEXrel log',d,w,'E[:3]',[str(x) for x in Ee[:3]])
# closed-form extended test y^200 of (K)/limitconv; control perturbs one harmonic term
from math import comb
N=int(sys.argv[2]) if len(sys.argv)>2 else 200
a=[Fr(1)]
for d in range(1,N+1): a.append(a[-1]*3*(3*d-1)*(3*d-2)/d**2)
H=[Fr(0)]
for n in range(1,3*N+1): H.append(H[-1]+Fr(1,n))
ok=True
for d in range(1,N+1):
    lhs=sum(a[k]*Fr(27)**(d-k)/(d-k) for k in range(d))
    rhs=a[d]*(3*H[3*d]-H[d])
    if lhs!=rhs: ok=False; print('fail',d); break
print('(K) through y^%d'%N,ok)
d=7
print('control: replacing 3H_{3d} by 3H_{3d-1} fails:',sum(a[k]*Fr(27)**(d-k)/(d-k) for k in range(d))!=a[d]*(3*H[3*d-1]-H[d]))
