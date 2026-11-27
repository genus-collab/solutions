# Symbolic shift-operator computation.
# Basis e_n = V^n/(a)_n, V=E v, a=1+E/z.  An operator O acts on G=sum_n b_n e_n by
# (O G)_m = sum_j c_j(m) b_{m+j};  O is stored as {j: poly in (m,z,E)}.
from fractions import Fraction as F
def padd(p,q):
    r=dict(p)
    for k,v in q.items():
        r[k]=r.get(k,0)+v
    return {k:v for k,v in r.items() if v}
def pmul(p,q):
    r={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(i+j for i,j in zip(a,b)); r[k]=r.get(k,0)+x*y
    return {k:v for k,v in r.items() if v}
def pscal(c,p): return {k:c*v for k,v in p.items() if c*v}
from math import comb
def pshift(p,j):  # m -> m+j
    r={}
    for (a,b,c),x in p.items():
        for i in range(a+1):
            k=(i,b,c); r[k]=r.get(k,0)+x*comb(a,i)*F(j)**(a-i)
    return {k:v for k,v in r.items() if v}
ONE={(0,0,0):F(1)}; M={(1,0,0):F(1)}; Z={(0,1,0):F(1)}; EE={(0,0,1):F(1)}
def oadd(A,B):
    r=dict(A)
    for j,p in B.items(): r[j]=padd(r.get(j,{}),p)
    return {j:p for j,p in r.items() if p}
def omul(A,B):  # A after B
    r={}
    for j,a in A.items():
        for k,b in B.items():
            r[j+k]=padd(r.get(j+k,{}),pmul(a,pshift(b,j)))
    return {j:p for j,p in r.items() if p}
def oscal(p,A): return {j:pmul(p,a) for j,a in A.items()}
def opow(A,n):
    R={0:ONE}
    for _ in range(n): R=omul(A,R)
    return R
sigma={1:ONE}                  # sigma e_n = e_{n-1}  => (sigma G)_m = b_{m+1}
N={0:M}
# mu = zV : mu e_n = (E+(n+1)z) e_{n+1} => (mu G)_m = (E+m z) b_{m-1}
mu={-1:padd(EE,pmul(M,Z))}
tau=oadd(oscal({(0,0,0):F(3)},sigma),oscal({(0,0,0):F(-1)},N))
def c(x): return {(0,0,0):F(x)} if x else {}
def lin(pz,pE,op=None):  # returns z*pz*op + ... helper
    pass
I={0:ONE}
A=oadd(mu,{0:pscal(F(-3),Z)})          # mu-3z
term1=omul(opow(A,3),opow(tau,3))
def fac(j):  # z tau + j z - E
    return oadd(oscal(Z,tau),{0:padd(pscal(F(j),Z),pscal(F(-1),EE))})
term2=oscal(c(27),omul(fac(0),omul(fac(1),fac(2))))
Q=oadd(term1,term2)
# recurrence operator Rec: (Rec b)_m = 27 m b_{m+1} - 3(3m^2-3m+1) b_m + (m-1)^3 b_{m-1}
m2=pmul(M,M)
Rec={1:pscal(F(27),M), 0:pscal(F(-3),padd(padd(pscal(F(3),m2),pscal(F(-3),M)),c(1))),
     -1:pmul(padd(M,c(-1)),pmul(padd(M,c(-1)),padd(M,c(-1))))}
import pickle
pickle.dump((Q,Rec),open('ops.pkl','wb'))
for j in sorted(Q): print(j, len(Q[j]))
