# Spot check 1: closed form h_w(z) v S(v,z) of Definition 1.3 vs stored cPP_z9.m (exact), incl. sign of h_w. Control: h_w -> 1/h_w, E->-E
import os
import re,sys
from fractions import Fraction as F
from math import comb
path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus4", "cPP_z9.m")
s=open(path).read()
s=re.sub(r'\s+','',s)
assert s[0]=='{' and s[-1]=='}'
s=s[1:-1]
parts=[];d=0;cur=''
for ch in s:
    if ch in '[(':d+=1
    if ch in '])':d-=1
    if ch==',' and d==0: parts.append(cur);cur=''
    else: cur+=ch
parts.append(cur)
print("stored z-orders:",len(parts))
def conv(x):
    x=x.replace(r'\[Epsilon][1]','e1').replace(r'\[Epsilon][2]','e2').replace('LLL[0][y]','L').replace('Sqrt[3]','s3')
    x=re.sub(r'(?<![a-zA-Z0-9_])(\d+)(?![a-zA-Z0-9_])',lambda m:'F(%s)'%m.group(1),x)
    # exponents: F(n) exponent must be int
    x=re.sub(r'\^F\((\d+)\)',lambda m:'**%s'%m.group(1),x)
    return x
codes=[compile(conv(p),'c','eval') for p in parts]
def bern(n):
    A=[F(0)]*(n+1)
    B=[]
    for m in range(n+1):
        A[m]=F(1,m+1)
        for j in range(m,0,-1):
            A[j-1]=j*(A[j-1]-A[j])
        B.append(A[0])
    return B
B=bern(30);B[1]=F(-1,2)
K=9
def smul(a,b):
    r=[F(0)]*(K+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if i+j<=K: r[i+j]+=x*y
    return r
def sexp(a):
    # a[0]=0
    r=[F(1)]+[F(0)]*K
    for n in range(1,K+1):
        r[n]=sum(k*a[k]*r[n-k] for k in range(1,n+1))/n
    return r
def Theta_exp(c):  # log Theta(c z)
    a=[F(0)]*(K+1)
    m=1
    while 2*m-1<=K:
        a[2*m-1]=B[2*m]*c**(2*m-1)/(2*m*(2*m-1));m+=1
    return a
def hw(e1,e2,sign=1):
    E=e1+e2
    a=[x+y+w for x,y,w in zip(Theta_exp(-1/E),Theta_exp(1/e1),Theta_exp(1/e2))]
    return sexp([sign*x for x in a])
r=[F(1)]
rm1=F(0)
for m in range(0,K+2):
    rp=(3*(3*m*m+3*m+1)*r[m]-(m**3)*(r[m-1] if m>=1 else F(0)))/(27*(m+1))
    r.append(rp)
print("r:",r[:5])
def Sser(v,E):
    # S(v,z) to z^K
    tot=[F(0)]*(K+1)
    for l in range(0,K+1):
        den=[F(1)]+[F(0)]*K
        for j in range(1,l+2):
            # 1/(1+j z/E)
            g=[(-F(j)/E)**n for n in range(K+1)]
            den=smul(den,g)
        for n in range(0,K+1-l):
            tot[n+l]+=r[l]*v**l*den[n]
    return tot
res=True
for (e1v,e2v) in [(F(7),F(3)),(F(5),F(2)),(F(4),F(9))]:
  for Lv in [F(1,7),F(-3,11)]:
    env=dict(F=F,e1=e1v,e2=e2v,L=Lv,sq=F(1))
    vals=[]
    for c in codes:
        env['s3']=F(1);a=eval(c,env)
        env['s3']=F(2);b=eval(c,env)
        p={F(2):1,F(1,2):-1}[b/a]
        # sqrt3 * stored (sq=1) : power (p+1): p=1 ->3*a ; p=-1 -> a
        vals.append(a*(3 if p==1 else 1))
    E=e1v+e2v
    v=3/(E+3*Lv)
    # expected: P0 * (sqrt3 sq E^2) = h_w v S   since P0 = v h S/sqrt(3 e1 e2 E), sq = sqrt(e1e2/E^3)
    h=hw(e1v,e2v);S=Sser(v,E)
    rhs=[v*x for x in smul(h,S)]
    lhs=[x*E**2 for x in vals]
    ok=lhs==rhs
    ctrl=[v*x for x in smul(hw(e1v,e2v,-1),S)]
    print((e1v,e2v),Lv,"match z^0..z^9:",ok," control h->1/h match:",lhs==ctrl, "first",lhs[:2],rhs[:2])
    res=res and ok
print("ALL",res)
