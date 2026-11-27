# Evidence only.
# Part 1 (exact): colour 0, Laurent series in t. Def 1.3 at c1=A equals eq:RHb (d2 = E*E(y), c2 = D) and the closure/row-2 ODE
#   (c1=1 gauge)  Th c = -(E/Delta_c) c - (E^2/3 + 6 y z^2) b - (E^3/27) a,  with b=Th a, c=Delta_c Th b.
# Part 2 (float): generic alpha, z^0 Frobenius route: eta_ab = Sum_i Phi_ai Phi_bi / Delta_i with Phi_2i=Y(Y-e1 E)/(A D).
from fractions import Fraction as Fr
from math import factorial as f
import cmath
TN=24     # truncation order in t
class L:  # Laurent series in t: dict exp->Fr, truncated at exp>TN
    __slots__=("d",)
    def __init__(s,d=None): s.d={k:v for k,v in (d or {}).items() if v!=0 and k<=TN}
    def __add__(s,o):
        r=dict(s.d)
        for k,v in o.d.items(): r[k]=r.get(k,Fr(0))+v
        return L(r)
    def __neg__(s): return L({k:-v for k,v in s.d.items()})
    def __sub__(s,o): return s+(-o)
    def __mul__(s,o):
        if isinstance(o,(int,Fr)): return L({k:v*o for k,v in s.d.items()})
        r={}
        for k,v in s.d.items():
            for l,u in o.d.items():
                if k+l<=TN: r[k+l]=r.get(k+l,Fr(0))+v*u
        return L(r)
    def th(s): return L({k:v*Fr(k,3) for k,v in s.d.items() if k!=0})
    def inv(s):
        m=min(s.d); c=s.d[m]
        r={-m:1/c}
        for n in range(1,TN+1+m+m if False else TN+abs(m)+m+40):
            e=-m+n
            if e>TN: break
            acc=Fr(0)
            for k,v in s.d.items():
                if k!=m and e-(k-m) in r: acc+=v*r[e-(k-m)]
            r[e]=-acc/c
        return L(r)
    def low(s): return min(s.d) if s.d else None
one=L({0:Fr(1)})
a=lambda d: Fr(f(3*d),f(d)**3)
Hn=[Fr(0)]
for n in range(1,100): Hn.append(Hn[-1]+Fr(1,n))
Dm=TN//3
A=L({3*d:a(d)*(-1)**d for d in range(Dm+1)})
G=L({3*d:-a(d)*Hn[d]*(-1)**d for d in range(Dm+1)})
h=L({3*d:3*a(d)*(Hn[3*d]-Hn[d])*(-1)**d for d in range(Dm+1)})
Ai=A.inv()
D=one+(h*Ai).th()
Ey=(G*Ai).th()
Dc=L({0:Fr(1),3:Fr(27)})
y=L({3:Fr(1)})
Xs=Dc.inv()
assert (A*A*D-Xs).d=={}
for (e1,e2) in [(Fr(5),Fr(2)),(Fr(1),Fr(1))]:
    Ew=e1+e2; alpha=Ew/3
    # dictionary
    C1i=Ai
    d2=(C1i*C1i-one)*(Ew/2)*Xs
    assert (d2-Ey*Ew).d=={}, "d2 = e1 E(y)"
    c2=Dc*(A*A)
    assert (c2*D-one).d=={}, "Delta_c c1^2 = 1/D"
    NZ=4
    r=[Fr(1)]
    for m in range(0,NZ+3):
        r.append((3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))/(27*(m+1)))
    Y=L({0:Fr(1),1:Fr(3)}).inv()*(-alpha)
    v=L({-1:1/Ew,0:3/Ew})
    def zfrac(l):
        s=[Fr(0)]*(NZ+1)
        if l<=NZ: s[l]=Fr(1)
        for j in range(1,l+2):
            for k in range(1,NZ+1):
                s[k]-=s[k-1]*Fr(j)/Ew
        return s
    S=[L() for _ in range(NZ+1)]; vp=one
    for l in range(NZ+1):
        zf=zfrac(l)
        for k in range(NZ+1): S[k]=S[k]+vp*(zf[k]*r[l])
        vp=vp*v
    def Z(a,b=None): return a
    def Th(u):   # Y + z theta on z-series
        return [Y*u[k]+(u[k-1].th() if k>=1 else L()) for k in range(NZ+1)]
    P0=[v*s for s in S]
    b=Th(P0)
    c=[Dc*x for x in Th(b)]
    lhs=Th(c)
    # Note: lhs etc truncated; the truncation order differs because of poles in t. Test only through t^(TN-2*NZ-6).
    rhs=[]
    XsE=Xs*Ew
    for k in range(NZ+1):
        t=-(XsE*c[k])-(b[k]*(Ew*Ew/3))-(P0[k]*(Ew**3/27))
        if k>=2: t=t-(y*b[k-2])*Fr(6)
        rhs.append(t)
    cut=TN-3*NZ-8
    for k in range(NZ+1):
        resid=lhs[k]-rhs[k]
        badk=[e for e in resid.d if e<=cut]
        assert not badk,(k,badk[:5])
    # eq:RHb versus Def 1.3 at c1=A
    R1=[Ai*x for x in Th(P0)]
    R2=[]
    tr=Th(R1)
    for k in range(NZ+1):
        R2.append((tr[k]-(Ey*Ew)*R1[k])*D.inv())
    D13=[]
    for k in range(NZ+1):
        D13.append((Dc*(A*A))*(tr[k]-d2*R1[k]))
    for k in range(NZ+1):
        resid=R2[k]-D13[k]
        assert not [e for e in resid.d if e<=cut]
print("Part 1 OK: dictionary c2=D, d2=e1 E(y), row-2 ODE through z^4, t^%d"%(TN-3*4-8))
# Part 2 float
alph=[0.37,1.91,-0.83]; e_1=sum(alph); e_2=alph[0]*alph[1]+alph[0]*alph[2]+alph[1]*alph[2]; e_3=alph[0]*alph[1]*alph[2]
eps=1.7
yv=0.004
def ser(cf):
    return sum(cf(d)*(-yv)**d for d in range(0,60))
aa=lambda d: f(3*d)/f(d)**3
Af=ser(lambda d: float(aa(d)))
Hf=[0.0]
for n in range(1,400): Hf.append(Hf[-1]+1.0/n)
tAf=ser(lambda d: d*float(aa(d)))
Gf=ser(lambda d: -float(aa(d))*Hf[d]); tGf=ser(lambda d: -d*float(aa(d))*Hf[d])
hf=ser(lambda d: 3*float(aa(d))*(Hf[3*d]-Hf[d])); thf=ser(lambda d: 3*d*float(aa(d))*(Hf[3*d]-Hf[d]))
Ef=(tGf*Af-Gf*tAf)/Af**2
Df=1+(thf*Af-hf*tAf)/Af**2
Xf=1/(1+27*yv)
print("float: A^2D-X =",Af**2*Df-Xf,"  2E-(D-X) =",2*Ef-(Df-Xf))
# roots of (Y+a0)(Y+a1)(Y+a2)+27 y Y^3
coef=[1+27*yv,e_1,e_2,e_3]
rt=[complex(-a) for a in alph]
for it in range(200):
    new=[]
    for i,zr in enumerate(rt):
        p=coef[0]*zr**3+coef[1]*zr**2+coef[2]*zr+coef[3]
        q=coef[0]
        for j,w in enumerate(rt):
            if j!=i: q*=(zr-w)
        new.append(zr-p/q)
    rt=new
W=[e_1*Y**2+2*e_2*Y+3*e_3 for Y in rt]
Del=[3*eps*w for w in W]
Phi=[[1,Y/Af,Y*(Y-e_1*Ef)/(Af*Df)] for Y in rt]
eta=[[sum(Phi[i][a_]*Phi[i][b_]/Del[i] for i in range(3)) for b_ in range(3)] for a_ in range(3)]
target=[[1/(3*eps*e_3),0,0],[0,0,-1/(3*eps)],[0,-1/(3*eps),e_1/(3*eps)]]
print("float: max |eta_ab - target| =",max(abs(eta[i][j]-target[i][j]) for i in range(3) for j in range(3)))
print("float: trace check: sum Y_i + e1 X =",sum(rt)+e_1*Xf)
