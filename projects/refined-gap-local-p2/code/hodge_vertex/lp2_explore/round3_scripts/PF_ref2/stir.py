from fractions import Fraction as F
from math import factorial
N=6; J=6*N+2
def run(al,be):
    # phi_ng(xi) = -al*sum_{j>=3} (-1)^{j+1} (xi/al)^j / j
    png=[F(0)]*(J+1)
    for j in range(3,J+1): png[j]=-al*F((-1)**(j+1),j)/al**j
    # exp(-png/z): series in xi with coefficients dict zpow->c
    ex=[dict() for _ in range(J+1)]; ex[0]={0:F(1)}
    term=[dict() for _ in range(J+1)]; term[0]={0:F(1)}
    for k in range(1,2*N+3):
        new=[dict() for _ in range(J+1)]
        for i in range(J+1):
            for zp,c in term[i].items():
                for j in range(3,J+1-i):
                    if png[j]: new[i+j][zp-1]=new[i+j].get(zp-1,0)-c*png[j]/k
        term=new
        for i in range(J+1):
            for zp,c in term[i].items(): ex[i][zp]=ex[i].get(zp,0)+c
    # gamma (1+xi/al)^{be-1}
    g=[F(1)]
    for j in range(1,J+1): g.append(g[-1]*(be-1-(j-1))/j/al)
    prod=[dict() for _ in range(J+1)]
    for i in range(J+1):
        for zp,c in ex[i].items():
            for j in range(J+1-i):
                prod[i+j][zp]=prod[i+j].get(zp,0)+c*g[j]
    # <xi^{2m}> = (z al)^m (2m-1)!!
    res={}
    for i in range(0,J+1,2):
        m=i//2; df=1
        for t in range(1,2*m,2): df*=t
        for zp,c in prod[i].items():
            res[zp+m]=res.get(zp+m,0)+c*df*al**m
    return [res.get(k,0) for k in range(N+1)]
def bern(n,x):
    B=[F(1)]
    for m in range(1,n+1):
        B.append(-sum(F(factorial(m+1),factorial(k)*factorial(m+1-k))*B[k] for k in range(m))/(m+1))
    return sum(F(factorial(n),factorial(k)*factorial(n-k))*B[k]*x**(n-k) for k in range(n+1))
def theta(w,x): # coefficients in z where argument = w*z
    lg=[F(0)]*(N+1)
    for k in range(1,N+1): lg[k]=-bern(k+1,x)*(-w)**k/(k*(k+1))
    out=[F(1)]+[F(0)]*N
    # exp of series
    for k in range(1,N+1):
        out[k]=sum(j*lg[j]*out[k-j] for j in range(1,k+1))/k
    return out
for al,be in [(F(2),F(-2)),(F(-3,2),F(-3)),(F(5),F(2,7))]:
    a=run(al,be); t=theta(1/al,be)
    print(al,be,a==t)
a=run(F(2),F(-2)); t=theta(F(1,2),F(-1))
print("control wrong beta:",a==t)
