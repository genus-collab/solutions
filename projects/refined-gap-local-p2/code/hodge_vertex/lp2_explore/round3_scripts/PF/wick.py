# r_l from BS Wick definition (lines 5257-5262), s3=1:
# r_l = sum_k (l+k)!/(3^(l+k) k!) [s1^(l-2k) s2^(l+k)] (1-X/27)^(-1/2),
# X = 18 s1 s2 - 4 s1^3 - 4 s2^3 + s1^2 s2^2   (D = -27 + X)
from fractions import Fraction as F
from math import factorial
import sys
L=int(sys.argv[1])
P1=L; P2=L+L//2
def mul(a,b):
    c={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if i+k<=P1 and j+l<=P2:
                c[(i+k,j+l)]=c.get((i+k,j+l),0)+x*y
    return {k:v for k,v in c.items() if v}
X={(1,1):F(18),(3,0):F(-4),(0,3):F(-4),(2,2):F(1)}
Xs={k:v/27 for k,v in X.items()}
G={(0,0):F(1)}; term={(0,0):F(1)}; n=0
while True:
    n+=1
    term=mul(term,Xs)
    if not term: break
    # coefficient binom(-1/2,n)(-1)^n = (2n)!/(4^n n!^2)
    c=F(factorial(2*n),4**n*factorial(n)**2)
    for k,v in term.items(): G[k]=G.get(k,0)+c*v
r=[]
for l in range(L+1):
    s=F(0)
    for k in range(l//2+1):
        s+=F(factorial(l+k),3**(l+k)*factorial(k))*G.get((l-2*k,l+k),0)
    r.append(s)
for l,x in enumerate(r): print(l,x)
