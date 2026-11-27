from fractions import Fraction as Fr
from math import comb
# Bernoulli numbers, convention x/(e^x-1)
def bern(n):
    B=[Fr(1)]
    for m in range(1,n+1):
        B.append(-sum(comb(m+1,k)*B[k] for k in range(m))/Fr(m+1))
    return B
B=bern(12)
def Bpoly(n,x): return sum(comb(n,k)*B[k]*x**(n-k) for k in range(n+1))
print("B1",B[1],"B6(1/5)/6 =",Bpoly(6,Fr(1,5))/6, "claim", Fr(4537,3937500))
# c_g
for g in range(2,6):
    print(g, Fr((-1)**(g-1))*B[2*g]/(2*g*(2*g-2)*5**(g-1)))
# h_1 via weights
import cmath
z=cmath.exp(2j*cmath.pi/5)
s=sum(1/(1-z**j) for j in range(1,5))+1/(-5)
print("sum a^-1 over w0:",s, " h1=",-s/12)
