from fractions import Fraction as F
from math import comb, floor
# Bernoulli numbers
def bern(n):
    B=[F(0)]*(n+1); B[0]=F(1)
    for m in range(1,n+1):
        B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/F(m+1)
    return B
B=bern(12)
for g in range(2,6):
    cg=F((-1)**(g-1))*B[2*g]/(2*g*(2*g-2)*5**(g-1))
    print("c_%d ="%g, cg)
nu=[0,-1,3,2,1]
print("mu table (b,k): exponent 4k+nu_b-nu_{b+k}, must be div by 5")
zero=[]
for k in range(0,12):
    row=[]
    for b in range(5):
        e=4*k+nu[b]-nu[(b+k)%5]
        assert e%5==0
        mu=e//5
        assert mu>=max(0,(4*k-4)/5)
        if mu==0: zero.append((b,k))
        row.append(mu)
    print(k,row)
print("mu=0 at",zero)
# degree bound: at u^n, k with mu(b,k)<=n ; max k vs floor((5n+4)/4)
for n in range(0,8):
    ks=[k for k in range(0,60) for b in range(5) if (4*k+nu[b]-nu[(b+k)%5])//5<=n]
    print(n,max(ks),(5*n+4)//4)
# lem:mai-Psi: m_a m_{3-a}=1 symbolic check with random rationals
import random
C0=F(random.randint(2,50),random.randint(2,50)); C1=F(random.randint(2,50),random.randint(2,50)); L=F(random.randint(2,50),random.randint(2,50))
m=[L/C0, L**2/(C0*C1), C0*C1/L**2, C0/L, F(1)]
print([m[a]*m[(3-a)%5] for a in range(5)])
# forward recursion at z^0 from X0=L/C0
ck=[C1, L**5/(C0**2*C1**2), C1, C0, C0]
X=[L/C0]
for a in range(4): X.append(L*X[a]/ck[a])  # colour 0
print([X[a]==m[a] for a in range(5)], "closing:", L*X[4]==C0*X[0])
# B6(1/5)/6
x=F(1,5); B6=x**6-3*x**5+F(5,2)*x**4-F(1,2)*x**2+F(1,42)
print(B6/6, F(4537,3937500))
# c(z)c(-z) coefficient
m5,m10=F(3,7),F(5,11)
print("z^10 coeff of c(z)c(-z):", 2*m10-m5**2, "= -(m5^2-2m10):", -(m5**2-2*m10))
