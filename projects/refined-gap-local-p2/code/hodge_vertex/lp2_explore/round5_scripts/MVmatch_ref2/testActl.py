# Controls for (A): each wrong variant must fail at g=1..5 on the same Schur data.
from fractions import Fraction as F
from hodge import *
from schur import *
from collections import Counter
tau=F(3,2);sig=tau*(1+tau)
mus=[(3,),(2,1),(2,2),(1,1)]
res={}
def rec(name,ok):
    a,b=res.get(name,(0,0));res[name]=(a+(0 if ok else 1),b+1)
for mu in mus:
    l=len(mu);n=sum(mu);G=5
    c=hatc(mu,tau,2*G-2+l+n+1)
    pr=F(1);prbar=F(1)
    for m in mu: pr*=Pi(m,tau);prbar*=Pi(m,-1-tau)
    aut=1
    for v in Counter(mu).values(): aut*=factorial(v)
    for g in range(1,G+1):
        lhs=c[2*g-2+l+n]
        Mt=M(g,tau,mu)
        rec("true",lhs==sig**(g+l-1)*pr*Mt)
        rec("extra (-1)^g",lhs==(-1)**g*sig**(g+l-1)*pr*Mt)
        rec("drop sigma^g",lhs==sig**(l-1)*pr*Mt)
        rec("Pi(-1-tau)",lhs==sig**(g+l-1)*prbar*Mt)
        rec("labelled 1/Aut",lhs==sig**(g+l-1)*pr*Mt/aut)
        rec("weights flipped",lhs==sig**(g+l-1)*pr*M(g,tau,mu,[F(1),1/tau,-1/(1+tau)]))
        rec("leaf 1/(1+mu psi)",lhs==sig**(g+l-1)*pr*M(g,tau,tuple(-m for m in mu)))
for k,v in res.items(): print(k,"mismatches/tested",v[0],"/",v[1])
