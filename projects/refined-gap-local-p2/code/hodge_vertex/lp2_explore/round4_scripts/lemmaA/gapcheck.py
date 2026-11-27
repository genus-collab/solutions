# exact: gap_g(x) = (2g-3)! [t^{2g}] t^2/(4 sinh(x t/2) sinh(t/2)), compare (-3)^(g-1) gap_g with stored local-P2 vertex poles
from fractions import Fraction as F
from math import factorial
N=12  # t-order
def sinh_over(a):  # series of sinh(a t)/(a t) in t, coefficients as dict power->poly in x (poly = list)
    pass
# represent polys in x as dict exp->Fraction; series in t as list of polys
def pmul(p,q):
    r={}
    for i,a in p.items():
        for j,b in q.items():
            r[i+j]=r.get(i+j,0)+a*b
    return {k:v for k,v in r.items() if v!=0}
def padd(p,q,s=1):
    r=dict(p)
    for k,v in q.items(): r[k]=r.get(k,0)+s*v
    return {k:v for k,v in r.items() if v!=0}
def smul(A,B):
    C=[{} for _ in range(N+1)]
    for i in range(N+1):
        for j in range(N+1-i):
            C[i+j]=padd(C[i+j],pmul(A[i],B[j]))
    return C
# S1 = sinh(x t/2)/(x t/2) = sum (x/2)^{2k} t^{2k}/(2k+1)!
S1=[{} for _ in range(N+1)]; S2=[{} for _ in range(N+1)]
for k in range(0,N//2+1):
    S1[2*k]={2*k:F(1,2**(2*k)*factorial(2*k+1))}
    S2[2*k]={0:F(1,2**(2*k)*factorial(2*k+1))}
P=smul(S1,S2)  # = sinh sinh /(x t^2/4)
# invert P (P[0]=1)
Inv=[{} for _ in range(N+1)]; Inv[0]={0:F(1)}
for n in range(1,N+1):
    acc={}
    for i in range(1,n+1): acc=padd(acc,pmul(P[i],Inv[n-i]))
    Inv[n]={k:-v for k,v in acc.items()}
# t^2/(4 sinh sinh) = (1/x) * Inv ; [t^{2g}] of t^2/(4 sinh sinh) -> Inv[2g-2]/x
for g in (2,3,4):
    gap={k-1:v*factorial(2*g-3) for k,v in Inv[2*g].items()}  # exponents shifted by 1/x
    lp2={k:v*(-3)**(g-1) for k,v in gap.items()}
    print("g",g,"conifold gap_g:",sorted(gap.items()),"\n   (-3)^(g-1) gap_g:",sorted(lp2.items()))
# stored local-P2 single-vertex poles (hodge_vertex summaries), as x^-1 * poly
stored={2:("-(7+10x^2+7x^4)/(1920x)",{-1:F(-7,1920),1:F(-10,1920),3:F(-7,1920)}),
        3:("-(1+x^2)(31+18x^2+31x^4)/(17920x)",None),
        4:("-(381+620x^2+686x^4+620x^6+381x^8)/(143360x)",{-1:F(-381,143360),1:F(-620,143360),3:F(-686,143360),5:F(-620,143360),7:F(-381,143360)})}
# expand g=3: (1+x^2)(31+18x^2+31x^4) = 31+49x^2+49x^4+31x^6
stored[3]=(stored[3][0],{-1:F(-31,17920),1:F(-49,17920),3:F(-49,17920),5:F(-31,17920)})
for g in (2,3,4):
    gap={k-1:v*factorial(2*g-3) for k,v in Inv[2*g].items()}
    lp2={k:v*(-3)**(g-1) for k,v in gap.items()}
    print("g",g,"match (-3)^(g-1) gap_g == stored local-P2 pole:", lp2==stored[g][1])
