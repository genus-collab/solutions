# misc controls: genus-1 BS data vs localisation; Pair at free c1,F on T_F; symplecticity at free c1,F; modular claim u->1-u by q-series.
import math
from fractions import Fraction as Fr
from tf import Closed,TF
from loc import TF_limit
# genus one: Q dQ F_{1,0}=-X/(24A), Q dQ F_{0,1}=-F/(2A)+(X-2)/(12A)
N=6
a=[Fr(1)]+[Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d for d in range(1,N)]
def mul(u,v):
    o=[Fr(0)]*N
    for i in range(N):
        for j in range(N-i): o[i+j]+=u[i]*v[j]
    return o
def inv(u):
    o=[Fr(0)]*N;o[0]=1/u[0]
    for n in range(1,N): o[n]=-sum(u[k]*o[n-k] for k in range(1,n+1))/u[0]
    return o
Ai=inv(a); th=[k*a[k] for k in range(N)]; F=mul(th,Ai)
X=[Fr((-27)**k) for k in range(N)]
one=[Fr(1)]+[Fr(0)]*(N-1)
f10=[-x/24 for x in mul(X,Ai)]
f01=[ -x/2+y/12 for x,y in zip(mul(F,Ai),mul([p-2*q for p,q in zip(X,one)],Ai))]
# Q=y exp(g): convert: coefficient extraction via y(Q) not needed: d=1,2 by Q-series
g=[Fr(0)]+[3*Fr(math.factorial(3*d-1),math.factorial(d)**3)*(-1)**d for d in range(1,N)]
ex=[Fr(0)]*N; ex[0]=Fr(1)
for n in range(1,N): ex[n]=sum(k*g[k]*ex[n-k] for k in range(1,n+1))/n
Qy=[Fr(0)]+ex[:N-1]
# series in y of a function = Sum c_d Q^d => solve triangular
def toQ(f):
    c=[Fr(0)]*N; r=list(f); pw=one[:]
    for d in range(N):
        c[d]=r[d]/ (pw[d]) if pw[d]!=0 else Fr(0)
        r=[x-c[d]*y for x,y in zip(r,pw)]
        pw=mul(pw,Qy)
    return c
c10=toQ(f10); c01=toQ(f01)
for (e1,e2) in ((3,7),(1,2)):
    E=Fr(e1+e2); k=Fr(e1*e2)
    for d in (1,2,3):
        pred=(E**2/k*c10[d]+c01[d])/d
        v,ok=TF_limit(1,d,e1,e2)
        print("g=1 eps",(e1,e2),"d",d,"BS genus-1 data:",pred,"localisation:",v,"match",pred==v)
# symplecticity and Pair on T_F at free c1,F
E=Fr(10); cl=Closed(E,8)
T=TF(cl,3,7,Fr(1,5),Fr(-7,3),Fr(11,5),'hv',0,[E,-3,-7])
bad=0
for i in range(3):
    for j in range(3):
        for K in range(0,9):
            x=None
            from tf import Om
            x=Om(0)
            for m in range(K+1):
                x=x+T.pair(i,j,m,K-m)*Fr((-1)**(K-m))
            if K==0 and i==j: x=x-T.Del[i]
            if not x.zero(): bad+=1
print("T_F symplecticity violations at free c1,F:",bad)
# modular: u(tau)=27 r/(1+27 r), r=eta(3tau)^12/eta(tau)^12 on tau=i s; u(-1/(3 i s))
def eta_ratio(s):
    q=math.exp(-2*math.pi*s)
    # r=q*prod (1-q^{3n})^12/(1-q^n)^12
    p=1.0
    for n in range(1,400):
        qn=q**n
        if qn<1e-30: break
        p*=((1-q**(3*n))/(1-qn))**12
    return q*p
for s in (0.3,0.5,0.8):
    r=eta_ratio(s); u=27*r/(1+27*r)
    r2=eta_ratio(1/(3*s)); u2=27*r2/(1+27*r2)
    print("tau=i*%.1f: u=%.12f u(-1/3tau)=%.12f sum=%.12f"%(s,u,u2,u+u2))
