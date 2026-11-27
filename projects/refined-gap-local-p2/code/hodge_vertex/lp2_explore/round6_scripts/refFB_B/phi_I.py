# Referee: Phi vs stripped I-function; truncation order K varied, delta varied. Evidence only.
import sys,os,math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from wick_ref import Phi,mul
def smul(a,b,K):
    r=[Fr(0)]*(K+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=K: r[i+j]+=x*y
    return r
def series_Y(al,j,K):
    e1=sum(al);e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2];e3=al[0]*al[1]*al[2]
    Y=[-al[j]]+[Fr(0)]*K
    Lp0=3*al[j]**2-2*e1*al[j]+e2
    for n in range(1,K+1):
        Y2=smul(Y,Y,K);Y3=smul(Y2,Y,K)
        L=[Y3[i]+(27*Y3[i-1] if i>=1 else 0)+e1*Y2[i]+e2*Y[i]+(e3 if i==0 else 0) for i in range(K+1)]
        Y[n]=-L[n]/Lp0
    return Y
def Cd(al,j,d,Nz):
    # Laurent in z, powers -d..Nz
    Lz=Nz+d
    ser=[Fr(1)]+[Fr(0)]*Lz
    for m in range(0,3*d):
        new=[Fr(0)]*(Lz+1)
        for n in range(Lz+1):
            new[n]+=3*al[j]*ser[n]
            if n+1<=Lz: new[n+1]+=-m*ser[n]
        ser=new
    for l in range(3):
        if l==j: continue
        for m in range(1,d+1):
            a=al[l]-al[j]; c=Fr(m)/a
            t=[Fr(0)]*(Lz+1);t[0]=ser[0]/a
            for n in range(1,Lz+1): t[n]=ser[n]/a-c*t[n-1]
            ser=t
    fac=Fr(1,math.factorial(d))
    return dict((n-d,ser[n]*fac) for n in range(Lz+1))
def fstripped(al,j,K,N):
    Y=series_Y(al,j,K)
    ut=[Fr(0)]*(K+1)
    for n in range(1,K+1): ut[n]=Y[n]/n
    Es=[{0:Fr(1)}]
    for n in range(1,K+1):
        cur={}
        for m in range(1,n+1):
            co=Fr(m,n)*(-ut[m])
            for p,c in Es[n-m].items():
                cur[p-1]=cur.get(p-1,0)+co*c
        Es.append(cur)
    out=[]
    for k in range(K+1):
        cur={}
        for d in range(k+1):
            cd=Cd(al,j,d,N+K+1) if d>0 else {0:Fr(1)}
            for p,c in cd.items():
                for q,x in Es[k-d].items():
                    if p+q<=N: cur[p+q]=cur.get(p+q,0)+c*x
        out.append(cur)
    return out,Y
al=[Fr(1),Fr(2),Fr(5)]
N=4
for K in (2,3,4,5):
  for j in (0,):
    fs,Yser=fstripped(al,j,K,N)
    neg=[(k,p,c) for k,x in enumerate(fs) for p,c in x.items() if p<0 and c!=0]
    print("K",K,"colour",j,"nonzero negative z powers:",neg[:3],flush=True)
    dl=Fr(1,10**7)
    Y0=-al[j]+dl
    e=[Y0+x for x in al]
    y=-(e[0]*e[1]*e[2])/(27*Y0**3)
    print("  Y series check: Y_series(y)-Y0 = %.3e ; y=%.3e"%(float(sum(Yser[n]*y**n for n in range(K+1))-Y0),float(y)))
    P=Phi(al,Y0,N); P00=Phi(al,-al[j],N)
    inv=[Fr(0)]*(N+1);inv[0]=Fr(1)
    for n in range(1,N+1):
        inv[n]=-sum(P00[k]*inv[n-k] for k in range(1,n+1))
    lhs=mul(P,inv,N)
    f=[sum(fs[k].get(n,Fr(0))*y**k for k in range(K+1)) for n in range(N+1)]
    rhs=[x/f[0] for x in f]
    print("  diffs:",["%.2e"%float(lhs[n]-rhs[n]) for n in range(N+1)],"  y^(K+1)=%.2e"%float(y**(K+1)),flush=True)
