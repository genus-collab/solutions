# Sanity check (not a proof): Conc(C_0, I_1) = 0 in Q[[q]][l] to order q^N,
# with the jet matrix of lem:gau-conc (iii), u = 1 - 5^5 q, D = q d/dq.
from fractions import Fraction as Fr
from math import factorial
N = 10
def D(a): return [k*a[k] for k in range(len(a))]
def mul(a,b):
    c=[Fr(0)]*N
    for i in range(N):
        if a[i]==0: continue
        for j in range(N-i): c[i+j]+=a[i]*b[j]
    return c
def add(*xs): return [sum(x[k] for x in xs) for k in range(N)]
def sc(s,a): return [s*x for x in a]
C0=[Fr(factorial(5*d),factorial(d)**5) for d in range(N)]
It=[Fr(0)]+[5*Fr(factorial(5*d),factorial(d)**5)*sum(Fr(1,m) for m in range(d+1,5*d+1)) for d in range(1,N)]
u=[Fr(1),Fr(-5**5)]+[Fr(0)]*(N-2); one_u=[Fr(0),Fr(5**5)]+[Fr(0)]*(N-2)
c=[C0]; i=[It]
for k in range(3): c.append(D(c[-1])); i.append(D(i[-1]))
zero=[Fr(0)]*N
J=[[zero,sc(Fr(2,5),one_u),one_u,sc(-1,u)],
   [sc(Fr(-2,5),one_u),zero,u,zero],
   [sc(-1,one_u),sc(-1,u),zero,zero],
   [u,zero,zero,zero]]
# antisymmetry of J (l-part of Conc(C0,I1) is l*Conc(C0,C0))
print("J antisymmetric:", all(J[a][b]==sc(-1,J[b][a]) for a in range(4) for b in range(4)))
# l-free part: jets of I1 without l: b*c_{b-1} + i_b
jet=[add(sc(b,c[b-1]) if b>0 else zero, i[b]) for b in range(4)]
tot=zero
for a in range(4):
    for b in range(4):
        tot=add(tot,mul(mul(c[a],J[a][b]),jet[b]))
print("Conc(C0,I1) l-free part mod q^%d:"%N, tot)
