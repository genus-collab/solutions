from fractions import Fraction as Fr
from math import comb
N=10
def mul(a,b): return [sum(a[k]*b[n-k] for k in range(n+1)) for n in range(N)]
def inv(a):
    h=[Fr(0)]*N; h[0]=1/a[0]
    for n in range(1,N): h[n]=-sum(a[k]*h[n-k] for k in range(1,n+1))/a[0]
    return h
def binpow(a,r):
    # a in 1+uQ[[u]], (a)^r via exp/log-free binomial: (1+x)^r
    x=a[:]; x[0]=Fr(0)
    res=[Fr(0)]*N; res[0]=Fr(1)
    xp=[Fr(1)]+[Fr(0)]*(N-1)
    coef=Fr(1)
    for k in range(1,N):
        xp=mul(xp,x); coef=coef*(r-k+1)/k
        res=[res[i]+coef*xp[i] for i in range(N)]
    return res
S=[Fr(4,5*n+4) for n in range(N)]
tau=inv(S)
v=[Fr(0)]+binpow(tau,Fr(-5,4))[:N-1]   # u*tau^{-5/4}
vp=[(n+1)*v[n+1] for n in range(N-1)]+[Fr(0)]
lhs=mul(mul(mul(mul(mul(vp,vp),vp),vp),vp),[Fr(0),Fr(1)]+[Fr(0)]*(N-2))
om=[Fr(comb(5,k))*(-1)**k for k in range(6)]+[Fr(0)]*(N-6)
lhs=mul(lhs,om)
print("v'^5 u (1-u)^5 - v:",[lhs[n]-v[n] for n in range(N-1)])
print("v:",v[:4])
# c_g
from fractions import Fraction
B={4:Fr(-1,30),6:Fr(1,42),8:Fr(-1,30),10:Fr(5,66)}
for g in range(2,6):
    cg=(-1)**(g-1)*B[2*g]*Fr(5)**(1-g)/(2*g*(2*g-2))
    print(g,cg)
