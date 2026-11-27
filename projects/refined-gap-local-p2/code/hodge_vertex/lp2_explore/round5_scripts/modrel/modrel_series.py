# Evidence only (exact, Fractions). Checks, from the DEFINITION eq:Ifun:
#  (a) eq:Ian coefficients, (b) L2b computations, (c) A^2 D = X, 2E = D - X, to y^N, (d) convolution eq:limitconv.
from fractions import Fraction as Fr
from math import factorial as f
N=60
Hn=[Fr(0)]
for n in range(1,3*N+2): Hn.append(Hn[-1]+Fr(1,n))
a=[Fr(f(3*d),f(d)**3) for d in range(N+1)]
sg=[(-1)**d for d in range(N+1)]
def mul(p,q):
    r=[Fr(0)]*(N+1)
    for i,x in enumerate(p):
        if x==0: continue
        for j,y in enumerate(q):
            if i+j>N: break
            r[i+j]+=x*y
    return r
def inv(p):
    r=[Fr(0)]*(N+1); r[0]=1/p[0]
    for n in range(1,N+1):
        r[n]=-sum(p[k]*r[n-k] for k in range(1,n+1))/p[0]
    return r
def th(p): return [k*x for k,x in enumerate(p)]
A=[a[d]*sg[d] for d in range(N+1)]
G=[-a[d]*Hn[d]*sg[d] for d in range(N+1)]                # theta I_1^[1]/e_1
h=[3*a[d]*(Hn[3*d]-Hn[d])*sg[d] for d in range(N+1)]     # P = A log y + h
Delta=[Fr(0)]*(N+1); Delta[0]=Fr(1); Delta[1]=Fr(27)
X=inv(Delta)
def L2b(p):
    tp=th(p); ttp=th(tp)
    r=mul(Delta,ttp)
    for k in range(N):
        r[k+1]+=27*tp[k]+6*p[k]
    return r
# (b)
assert L2b(A)==[0]*(N+1)
logD=[Fr(0)]+[Fr((-1)**(k+1)*27**k,k) for k in range(1,N+1)]   # log(1+27y)
f_=[a[d]*(3*Hn[3*d]-Hn[d])*sg[d] for d in range(N+1)]
assert f_==[-x for x in mul(A,logD)], "f = -A log Delta"
Rhs=[Fr(0)]*(N+1)
for d in range(1,N+1): Rhs[d]=27*(2*d-1)*a[d-1]*sg[d]
assert L2b(f_)==Rhs
assert L2b([-x for x in mul(A,logD)])==Rhs
tA=th(A); t=mul(Delta,tA); t=[2*x for x in t]
for k in range(N): t[k+1]+=27*A[k]
lh=L2b(h)
assert [x+y for x,y in zip(t,lh)]==[0]*(N+1), "L2b P=0"
assert L2b(G)==[-x for x in tA], "L2b G=-theta A"
# (c)
Ai=inv(A)
hA=mul(h,Ai); tt=th(hA)
D=[x for x in tt]; D[0]+=1        # D = theta(log y + h/A)
E=th(mul(G,Ai))
A2=mul(A,A)
assert mul(A2,D)==X, "A^2 D = X"
assert [2*x for x in E]==[x-y for x,y in zip(D,X)], "2E=D-X"
# (d) eq:limitconv
for d in range(1,N+1):
    lhs=sum(Fr(1,d-k)*a[k]/27**k for k in range(d))
    rhs=a[d]/27**d*(2*Hn[d]+Fr(1,d)+3*Hn[3*d-1]-3*Hn[d])
    assert lhs==rhs
print("series identities OK to y^%d"%N, "E[1..3]=",E[1:4],"D[0..3]=",D[:4])
# (a) eq:Ian from eq:Ifun: d-th summand of I/(z y^{H/z}) in s=1/z to s^2, H-degree<=2, generic alpha
alpha=[Fr(2,7),Fr(5,3),Fr(-11,4)]
e1=sum(alpha)
def pm(p,q):
    r={}
    for (a1,b1),x in p.items():
        for (a2,b2),y in q.items():
            if a1+a2<=2 and b1+b2<=2:
                r[(a1+a2,b1+b2)]=r.get((a1+a2,b1+b2),Fr(0))+x*y
    return r
def pinv_lin(c,al):
    u={(1,1):Fr(-1)/c,(1,0):-al/c}
    t={(0,0):Fr(1)}; acc={(0,0):Fr(1)}
    for k in range(2):
        t=pm(t,u)
        for kk,v in t.items(): acc[kk]=acc.get(kk,Fr(0))+v
    return {kk:v/c for kk,v in acc.items()}
for d in range(1,7):
    term={(0,0):Fr(1)}
    for m in range(0,3*d):
        fac={(1,1):Fr(-3)}
        if m>0: fac[(0,0)]=Fr(-m)
        term=pm(term,fac)
    for j in range(3):
        for m in range(1,d+1):
            term=pm(term,pinv_lin(Fr(m),alpha[j]))
    c=Fr(f(3*d-1),f(d)**3)*(-1)**d
    assert term.get((1,1),0)==3*c
    assert term.get((1,0),0)==0
    assert term.get((2,0),0)==0
    assert term.get((2,1),0)==-3*e1*c*Hn[d]
    assert term.get((2,2),0)==9*c*(Hn[3*d-1]-Hn[d])
print("eq:Ian coefficients reproduced from eq:Ifun for d<=6 at generic alpha")
