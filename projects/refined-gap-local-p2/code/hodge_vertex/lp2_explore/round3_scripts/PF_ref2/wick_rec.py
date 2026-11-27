from fractions import Fraction as F
from math import factorial
L=18
P1,P2=L,L+L//2+1
def mul(a,b):
    c={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if i+k<=P1 and j+l<=P2:
                c[(i+k,j+l)]=c.get((i+k,j+l),0)+x*y
    return c
X={(2,2):F(1),(3,0):F(-4),(1,1):F(18),(0,3):F(-4)}
mX={k:-v/27 for k,v in X.items()}
G={(0,0):F(1)}; pw={(0,0):F(1)}; bc=F(1)
for m in range(1,(P1+P2)//2+2):
    pw=mul(pw,mX); bc=bc*(F(-1,2)-m+1)/m
    for k,v in pw.items(): G[k]=G.get(k,0)+bc*v
def wick(l):
    s=F(0)
    for k in range(0,l//2+1):
        s+=F(factorial(l+k),3**(l+k)*factorial(k))*G.get((l-2*k,l+k),0)
    return s
r=[F(1)]; rm=F(0)
rec=[F(1)]
prev=F(0)
for m in range(0,L):
    nxt=(3*(3*m*m+3*m+1)*rec[m]-m**3*(rec[m-1] if m>=1 else 0))/(27*(m+1))
    rec.append(nxt)
W=[wick(l) for l in range(L+1)]
for l in range(L+1): print(l,W[l],rec[l],W[l]==rec[l])
# control: different recurrence variant (e.g. sign of m^3 term)
alt=[F(1)]
for m in range(0,6):
    alt.append((3*(3*m*m+3*m+1)*alt[m]+m**3*(alt[m-1] if m>=1 else 0))/(27*(m+1)))
print("control alt",alt[:6],[W[i]==alt[i] for i in range(6)])
