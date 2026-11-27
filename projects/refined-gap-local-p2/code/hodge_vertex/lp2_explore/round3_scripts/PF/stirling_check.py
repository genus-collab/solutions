# F5 check: formal Asym of e^{-(q-log q)/z} q^{beta-1} dq at q=1 (alpha=1) vs Theta(z,beta)
from fractions import Fraction as Fr
from math import comb, factorial
K=8; beta=Fr(2,7); D=6*K+2   # xi-degree needed
def mul(a,b):  # series in xi, coefficients dict zpow->Fr
    c=[{} for _ in range(D+1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j>D: break
            for p,u in x.items():
                for q,v in y.items():
                    c[i+j][p+q]=c[i+j].get(p+q,0)+u*v
    return c
png=[{} for _ in range(D+1)]
for k in range(3,D+1): png[k]={-1:Fr(-(-1)**k,k)}   # -phi_ng/z
E=[{} for _ in range(D+1)]; E[0]={0:Fr(1)}; term=[dict(x) for x in E]
for n in range(1,D//3+1):
    term=mul(term,png); term=[{p:v/n for p,v in x.items()} for x in term]
    for i in range(D+1):
        for p,v in term[i].items(): E[i][p]=E[i].get(p,0)+v
g=[{0:Fr(1)}]; c=Fr(1)
for i in range(1,D+1):
    c=c*(beta-1-(i-1))/i; g.append({0:c})
f=mul(E,g)
res={}
for i in range(0,D+1,2):
    m=i//2; dfac=1
    for j in range(1,2*m,2): dfac*=j
    for p,v in f[i].items():
        if p+m<=K: res[p+m]=res.get(p+m,0)+v*dfac
# Theta(z,beta)=exp(-sum_k B_{k+1}(beta)(-z)^k/(k(k+1)))
B=[Fr(1)]
for n in range(1,K+3): B.append(-sum(comb(n+1,k)*B[k] for k in range(n))/(n+1))
def Bpoly(n,x): return sum(comb(n,k)*B[k]*x**(n-k) for k in range(n+1))
L=[Fr(0)]+[-Bpoly(k+1,beta)*(-1)**k/(k*(k+1)) for k in range(1,K+1)]
th=[Fr(1)]+[Fr(0)]*K
# exp of series L
for n in range(1,K+1):
    th[n]=sum(k*L[k]*th[n-k] for k in range(1,n+1))/n
print([res.get(k,0)==th[k] for k in range(K+1)])
