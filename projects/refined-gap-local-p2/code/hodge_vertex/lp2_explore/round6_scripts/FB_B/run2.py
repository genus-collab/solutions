# Evidence only. (a) continuity of Tr(x)Tr total at fixed y as alpha -> diagonal; (b) generic alpha, tiny y: total vs localisation d=0,1
# (c) polynomiality of eps1 eps2 GW_{2,1} and eps1 eps2 e3 GW_{2,0} in alpha (finite differences, exact), and exact value at T_F by interpolation
# (d) y-expansion of the T_F total against GW^{T_F}_{2,d}, d=0,1,2  (e) N(W): den at y=0 vanishes on the diagonal, reduction ~ y^2
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','FB_A'))
from trace_total import *
import g2d1
from g2d1 import loc_d0,loc_d1,refined_d1
from g2d2_limit import refined_d2

class P: pass
def mk(al,e1,e2):
    p=P(); p.al=[Fr(x) for x in al]; p.e1=Fr(e1); p.e2=Fr(e2)
    p.wt=[[3*p.al[i]]+[p.al[j]-p.al[i] for j in range(3) if j!=i]+[-p.e1,-p.e2] for i in range(3)]
    return p

print("(a) near-diagonal continuity at fixed y, F",flush=True)
for (e1,e2,y,F) in [(Fr(3),Fr(7),Fr(1,5),Fr(1)),(Fr(5,2),Fr(-3,7),Fr(3,11),Fr(-2,3))]:
    E=e1+e2
    base=Total([E/3]*3,e1,e2,y).F2(F)
    print(" eps",e1,e2,"y",y,"F",F,"T_F value",float(base),"BS F2ref equal:",base==bs_F2ref(e1,e2,y,F),flush=True)
    for k in [1,2,3,4]:
        dl=Fr(1,10**k)
        al=[E/3+dl,E/3+2*dl,E/3-3*dl]
        val=Total(al,e1,e2,y).F2(F)
        print("   delta 1e-%d: total - T_F value = %.6e"%(k,float(val-base)),flush=True)

print("(b) generic alpha, y0=1e-25: (total - GW_{2,0})/y0 vs GW_{2,1}",flush=True)
def Fser(y):
    # F = theta log A truncated at y^3 (A=1-6y+90y^2-1680y^3+...)
    a=[Fr(1),Fr(-6),Fr(90),Fr(-1680),Fr(34650)]
    th=[k*a[k] for k in range(5)]
    A=sum(a[k]*y**k for k in range(5)); dA=sum(th[k]*y**k for k in range(5))
    return dA/A
for al,e1,e2 in [((1,2,5),3,7),((2,-3,7),1,4),((1,3,8),Fr(5,2),Fr(-3,7))]:
    y0=Fr(1,10**25)
    T=Total(al,e1,e2,y0)
    val=T.F2(Fser(y0))
    p=mk(al,e1,e2); g0=loc_d0(p); g1=loc_d1(p)
    print(" alpha",al,"eps",e1,e2,": total-GW0 =%.3e ; (total-GW0)/y0 - GW1 = %.3e (GW1=%.6f)"%(float(val-g0),float((val-g0)/y0-g1),float(g1)),flush=True)
    Tc=Total(al,e1,e2,y0,hsign=-1)
    print("    control hsign=-1: total-GW0 = %.3e"%float(Tc.F2(Fser(y0))-g0),flush=True)

print("(c) polynomiality in alpha of eps1 eps2 GW^T_{2,1}, eps1 eps2 e3 GW^T_{2,0} along lines through the diagonal",flush=True)
def findiff(vals,order):
    v=list(vals)
    for _ in range(order):
        v=[v[i+1]-v[i] for i in range(len(v)-1)]
    return v
for (e1,e2) in [(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))]:
    E=e1+e2; k=e1*e2
    d=(Fr(1),Fr(2,3),Fr(-5,7))
    ss=[Fr(j) for j in range(1,14)]
    v1=[]; v0=[]
    for s in ss:
        al=[E/3+s*x for x in d]
        p=mk(al,e1,e2)
        v1.append(k*loc_d1(p))
        e3=al[0]*al[1]*al[2]
        v0.append(k*e3*loc_d0(p))
    fd1=findiff(v1,5); fd0=findiff(v0,8)
    print(" eps",e1,e2,": 5th differences of eps1eps2 GW21 all 0:",all(x==0 for x in fd1),"; 8th differences of eps1eps2 e3 GW20 all 0:",all(x==0 for x in fd0),flush=True)
    # controls: without the eps factor nothing changes; without e3 the d=0 function has a pole on the line => differences nonzero
    v0c=[]
    for s in ss:
        al=[E/3+s*x for x in d]
        v0c.append(k*loc_d0(mk(al,e1,e2)))
    print("   control, no e3: 8th differences all 0:",all(x==0 for x in findiff(v0c,8)),flush=True)
    # Lagrange at s=0
    def lag0(xs,ys):
        tot=Fr(0)
        for i,xi in enumerate(xs):
            t=ys[i]
            for j,xj in enumerate(xs):
                if j!=i: t*=(Fr(0)-xj)/(xi-xj)
            tot+=t
        return tot
    g21=lag0(ss[:5],v1[:5])/k
    g20=lag0(ss[:8],v0[:8])/(k*(E/3)**3)
    print("   exact T_F values by interpolation: GW21 =",g21,"(refined formula:",refined_d1(e1,e2),") GW20 =",g20,"(E^2/640 =",E**2/640,")",flush=True)

print("(d) y-expansion of eq:F2ref (=T_F total by run1) vs GW_{2,d}",flush=True)
N=4
def smul(a,b):
    r=[Fr(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N: r[i+j]+=x*y
    return r
def sinv(a):
    r=[Fr(0)]*(N+1); r[0]=1/a[0]
    for n in range(1,N+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
Aser=[Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d for d in range(N+1)]
Fs=smul([k*Aser[k] for k in range(N+1)],sinv(Aser))
Xs=[Fr((-27)**d) for d in range(N+1)]
iX=[Fr(1),Fr(27)]+[Fr(0)]*(N-1)
def spow(a,n):
    r=[Fr(1)]+[Fr(0)]*N
    for _ in range(n): r=smul(r,a)
    return r
def sadd(*xs):
    return [sum(x[i] for x in xs) for i in range(N+1)]
def sc(c,a): return [c*x for x in a]
f02=sadd(sc(Fr(5,8),smul(spow(Fs,3),iX)),smul(sadd(sc(Fr(5,8),iX),sc(Fr(-1,2),[Fr(1)]+[Fr(0)]*N)),spow(Fs,2)),sc(Fr(1,96),smul(sadd(sc(Fr(13),Xs),sc(Fr(20),iX),sc(Fr(-32),[Fr(1)]+[Fr(0)]*N)),Fs)),sc(Fr(-1,2160),sadd(sc(Fr(27),smul(spow(Xs,3),iX)),sc(Fr(-98),smul(spow(Xs,2),iX)),sc(Fr(121),smul(Xs,iX)),sc(Fr(-50),iX))))
f11=sadd(sc(Fr(1,16),spow(Fs,2)),sc(Fr(1,96),smul(sadd(sc(Fr(8),[Fr(1)]+[Fr(0)]*N),sc(Fr(-7),Xs)),Fs)),sc(Fr(1,5760),sadd(sc(Fr(84),spow(Xs,2)),sc(Fr(-196),Xs),sc(Fr(121),[Fr(1)]+[Fr(0)]*N))))
f20=sadd(sc(Fr(1,384),smul(Xs,Fs)),sc(Fr(-1,5760),sadd(sc(Fr(21),spow(Xs,2)),sc(Fr(-29),Xs),sc(Fr(8),[Fr(1)]+[Fr(0)]*N))))
for (e1,e2) in [(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))]:
    E=e1+e2;k=e1*e2
    ser=sadd(sc(E**4/k,f20),sc(E**2,f11),sc(k,f02))
    print(" eps",e1,e2,": y^0 == E^2/640:",ser[0]==E**2/640,"; y^1 == refined d=1:",ser[1]==refined_d1(e1,e2),"; y^2+6y^1 == refined d=2:",ser[2]+6*ser[1]==refined_d2(e1,e2),flush=True)

print("(e) Norm of W in C: den at y=0 vs reduction",flush=True)
def normW(al,y):
    R=Ring(al,y)
    W=R.add(R.add(R.scal(R.e1,R.mul(R.Y,R.Y)),R.scal(2*R.e2,R.Y)),R.scal(3*R.e3,R.one))
    M=R.matrix(W)
    det=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    return det,R
al=[Fr(1),Fr(2),Fr(5)]
d0,R=normW(al,Fr(0))
disc=((al[0]-al[1])*(al[0]-al[2])*(al[1]-al[2]))**2
print(" generic alpha: N(W)|_{y=0} =",d0,"; -e3*Disc(alpha) =",-R.e3*disc)
E=Fr(10)
for y in [Fr(1,5),Fr(-2,7)]:
    d,R=normW([E/3]*3,y)
    print(" T_F: N(W) =",d,"; E^9 y^2/Dc^2 =",E**9*y**2/(1+27*y)**2)
