# Test (A) at genera/taus beyond the author: g up to 6, special tau where Pi vanishes; controls.
from fractions import Fraction as F
from hodge import *
from schur import *
from collections import Counter
import time,sys

def M0(mu):
    n=sum(mu);l=len(mu)
    return F(n)**(l-3)

def Mg(g,mu,tau,ts=None):
    if g==0: return M0(mu)
    return M(g,tau,mu,ts)

def run(gmax,mus,taus):
    nt=0;bad=0
    for tau in taus:
        tau=F(tau);sig=tau*(1+tau)
        for mu in mus:
            l=len(mu);n=sum(mu)
            gm=gmax if isinstance(gmax,int) else gmax(mu)
            P=2*gm-2+l+n+1
            c=hatc(mu,tau,P)
            # forbidden powers
            for i in range(P):
                e=-n+i
                # allowed: e=2g-2+l for g>=0
                ok=(e-l+2)>=0 and (e-l+2)%2==0
                if not ok and c[i]!=0:
                    print("FORBIDDEN",mu,tau,e,c[i]);bad+=1
            prpi=F(1)
            for m in mu: prpi*=Pi(m,tau)
            for g in range(gm+1):
                e=2*g-2+l
                lhs=c[e+n]
                rhs=sig**(g+l-1)*prpi*Mg(g,mu,tau)
                nt+=1
                if lhs!=rhs:
                    bad+=1
                    print("MISMATCH",g,mu,tau,lhs,rhs)
        print("tau",tau,"done",nt,"bad",bad,time.strftime("%X"),flush=True)
    return nt,bad

if __name__=="__main__":
    taus=[F(3,2),F(-7,5),F(-1,2),F(-1,3),F(-2,3),F(1),F(5,11)]
    mus5=[(1,),(2,),(3,),(4,),(5,),(1,1),(2,1),(2,2),(3,2),(3,1),(1,1,1),(2,1,1),(2,2,1),(2,2,2),(3,2,1)]
    print(run(lambda mu:5 if len(mu)<=2 else 4,mus5,taus))
    print(run(6,[(1,),(2,),(3,)],[F(3,2),F(-1,3)]))
