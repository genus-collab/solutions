# Direct check: (Ev-3)^3 * e^{-u/z} L_pf e^{u/z} applied to F = v*S(v,z) vanishes through z^K.
# theta = (1/E - v/3) d/dv,  Y+alpha = 1/v, alpha=E/3, 27y = 27/(Ev-3)^3.  Uses recurrence r_l.
from fractions import Fraction as Fr
import sys
K=int(sys.argv[1]); E=Fr(sys.argv[2])
r=[Fr(1)]
for n in range(K+2):
    r.append((3*(3*n*n+3*n+1)*r[n]-(n**3*r[n-1] if n else 0))/(27*(n+1)))
def ser_inv_prod(l):  # coefficients of prod_{j=1}^{l+1} (1+j z/E)^{-1} up to z^K
    c=[Fr(0)]*(K+1); c[0]=Fr(1)
    for j in range(1,l+2):
        a=Fr(j)/E  # divide by (1+a z)
        for k in range(1,K+1): c[k]-=a*c[k-1]
    return c
F={k:{} for k in range(K+1)}
for l in range(K+1):
    c=ser_inv_prod(l)
    for k in range(l,K+1):
        F[k][l+1]=F[k].get(l+1,0)+r[l]*c[k-l]
def ladd(p,q,s=1):
    o=dict(p)
    for e,x in q.items(): o[e]=o.get(e,0)+s*x
    return {e:x for e,x in o.items() if x}
def lmul(p,q):
    o={}
    for e,x in p.items():
        for f,y in q.items(): o[e+f]=o.get(e+f,0)+x*y
    return {e:x for e,x in o.items() if x}
def theta(p):
    o={}
    for e,x in p.items():
        if e: 
            o[e-1]=o.get(e-1,0)+x*e/E
            o[e]=o.get(e,0)-x*e/3
    return {e:x for e,x in o.items() if x}
def sadd(A,B,s=1): return {k:ladd(A.get(k,{}),B.get(k,{}),s) for k in range(K+1)}
def T(A,shift=Fr(0)):   # (z theta + 1/v + shift(z)) ; shift given as (c0,c1): c0 + c1 z
    pass
def Top(A,c0,c1):  # (z theta + 1/v + c0 + c1 z) A
    o={k:{} for k in range(K+1)}
    for k in range(K+1):
        p=A[k]
        o[k]=ladd(o[k],{e-1:x for e,x in p.items()})
        o[k]=ladd(o[k],{e:c0*x for e,x in p.items()})
        if k+1<=K:
            o[k+1]=ladd(o[k+1],theta(p))
            o[k+1]=ladd(o[k+1],{e:c1*x for e,x in p.items()})
    return o
t1=Top(Top(Top(F,0,0),0,0),0,0)
cub={0:Fr(-27),1:27*E,2:-9*E*E,3:E**3}  # (Ev-3)^3
t1={k:lmul(cub,p) for k,p in t1.items()}
t2=Top(Top(Top(F,-E/3,Fr(2,3)),-E/3,Fr(1,3)),-E/3,Fr(0))
t2={k:{e:27*x for e,x in p.items()} for k,p in t2.items()}
res=sadd(t1,t2)
bad=[k for k in range(K+1) if res[k]]
print("E=%s K=%d nonzero orders:"%(E,K), bad)
