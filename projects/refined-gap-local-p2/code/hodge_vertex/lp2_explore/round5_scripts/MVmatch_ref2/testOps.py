# Operators for R+1, R+0 of (B) as polynomial identities in s (coefficientwise), genus 2,3 incl. non-integers.
from fractions import Fraction as F
from hodge import Qk
def Qlist(g,tau,R,extra=2):
    N=3*g-3+len(R)
    return [Qk(g,tau,R,k) for k in range(N+extra+1)]
def padd(a,b):
    n=max(len(a),len(b));return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
def pmul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            r[i+j]+=x*y
    return r
def der(a): return [i*a[i] for i in range(1,len(a))] or [F(0)]
def strip(a):
    a=list(a)
    while a and a[-1]==0:a.pop()
    return a
tot=0;bad=0
for g,R,tau in [(2,(3,),F(5,3)),(2,(F(5,2),2),F(-9,7)),(3,(4,),F(2,9)),(3,(),F(5,3)),(2,(),F(-9,7)),(1,(5,),F(5,3)),(1,(F(1,2),3),F(7,3)),(0,(2,3,4),F(5,3))]:
    r=len(R);J=sum(R) if R else F(0);n=2*g-2+r
    Q=Qlist(g,tau,R)
    o1=padd(pmul([F(1),F(2),F(1)],der(Q)),pmul([J+n,F(n)],Q))
    o0=padd(pmul([J,F(n)],Q),pmul([F(0),F(1),F(1)],der(Q)))
    q1=Qlist(g,tau,tuple(R)+(1,))
    q0=Qlist(g,tau,tuple(R)+(0,))
    a=strip(o1)==strip(q1);b=strip(o0)==strip(q0)
    # control
    c=strip(padd(pmul([F(1),F(2),F(1)],der(Q)),pmul([J,F(n)],Q)))==strip(q1)
    print(g,R,tau,"R+1",a,"R+0",b,"control(J+n->J) passes:",c,flush=True)
