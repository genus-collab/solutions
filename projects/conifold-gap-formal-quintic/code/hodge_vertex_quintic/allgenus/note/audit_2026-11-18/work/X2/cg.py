from fractions import Fraction as F
from math import comb
# Bernoulli numbers B_n (B_1=-1/2)
def bern(N):
    B=[F(0)]*(N+1); B[0]=F(1)
    for m in range(1,N+1):
        B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/F(m+1)
    return B
B=bern(20)
for h in range(2,9):
    ftr=B[2*h]/(2*h*(2*h-2))
    c1=F((-5))**(1-h)*ftr
    c2=F((-1)**(h-1))*B[2*h]/(2*h*(2*h-2)*5**(h-1))
    c3=abs(B[2*h])/(2*h*(2*h-2)*5**(h-1))
    print(h,ftr,c1,c1==c2==c3)
# tau, v, p00, d, delta to order u^2 for generic alpha: check d*t_c - v in u^2 Q[[u]] at order u^2 (needs t_c to u^2 only)
# series in u truncated at u^3
N=4
def mul(a,b): 
    c=[F(0)]*N
    for i in range(N):
        for j in range(N-i): c[i+j]+=a[i]*b[j]
    return c
def inv(a):
    b=[F(0)]*N; b[0]=1/a[0]
    for n in range(1,N): b[n]=-sum(a[k]*b[n-k] for k in range(1,n+1))/a[0]
    return b
def powr(a,r): # a[0]=1, binomial
    # exp(r log a)
    # log a
    x=[F(0)]+a[1:]
    lg=[F(0)]*N; xp=[F(1)]+[F(0)]*(N-1)
    for k in range(1,N):
        xp=mul(xp,x); lg=[lg[i]+F((-1)**(k+1),k)*xp[i] for i in range(N)]
    y=[r*t for t in lg]; e=[F(1)]+[F(0)]*(N-1); yp=[F(1)]+[F(0)]*(N-1); fact=1
    for k in range(1,N):
        yp=mul(yp,y); fact*=k; e=[e[i]+yp[i]/fact for i in range(N)]
    return e
tauinv=[F(4,5*n+4) for n in range(N)]
tau=inv(tauinv)
print('tau',tau)
v=[F(0)]+powr(tau,F(-5,4))[:N-1]
print('v',v)
