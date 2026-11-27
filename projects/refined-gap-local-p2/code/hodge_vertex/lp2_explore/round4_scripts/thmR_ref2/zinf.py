from fractions import Fraction as F
from math import factorial
N=18
r=[F(1)]; rm1=F(0)
for m in range(N+2):
    prev=r[m-1] if m>=1 else F(0)
    r.append((3*(3*m*m+3*m+1)*r[m]-m**3*prev)/(27*(m+1)))
# series in V
def mul(a,b):
    c=[F(0)]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]): c[i+j]+=x*y
    return c
def der(a): return [ (i+1)*a[i+1] for i in range(N)]+[F(0)]
def theta(a): return mul([F(1),F(-1,3)]+[F(0)]*(N-1), der(a))
S=[F(0)]*(N+1)
for l in range(N): S[l+1]=r[l]/factorial(l+1)
# y = 1/(V-3)^3 = -1/27 (1-V/3)^-3
inv=[F(0)]*(N+1)
for k in range(N+1):
    from math import comb
    inv[k]=F(-1,27)*comb(k+2,2)*F(1,3)**k
Dc=[F(1)+27*inv[0]]+[27*x for x in inv[1:]]
tS=theta(S)
L=[a+b+c for a,b,c in zip(mul(Dc,theta(tS)),mul([27*x for x in inv],tS),mul([6*x for x in inv],[F(1)]+[F(0)]*N))]
L2=[a+b for a,b in zip(mul(Dc,theta(tS)),mul([27*x for x in inv],tS))]
L2=[a+b for a,b in zip(L2,mul([6*x for x in inv],S))]
print("L2b(theta S)=",[str(x) for x in L[:N-2]])
print("alt (6yS instead of 6y thetaS... ) =",[str(x) for x in L2[:6]])
print("Dc series",[str(x) for x in Dc[:4]])
f=tS
y=inv
L=[a+b+c for a,b,c in zip(mul(Dc,theta(theta(f))),mul([27*x for x in y],theta(f)),mul([6*x for x in y],f))]
print("fixed L2b(theta S)=",[str(x) for x in L[:N-3]])
# control: perturb S by r_2 change
S2=S[:]; S2[3]+=F(1,100)
f=theta(S2)
L=[a+b+c for a,b,c in zip(mul(Dc,theta(theta(f))),mul([27*x for x in y],theta(f)),mul([6*x for x in y],f))]
print("control",[str(x) for x in L[:5]])
