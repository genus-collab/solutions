# Reconcile BS eq:betagiv (5116-5124, graph sum with pure psi vertices, Upsilon edges/dilatons) with Z[g,0] of the draft.
# genus 2 (7 graphs) and genus 3 by Mumford equivalence: direct sum at g=2. Exact rationals. Evidence only.
from fractions import Fraction as Fr
from math import factorial
from hodge import wk,bern
import sys
e1,e2=Fr(sys.argv[1]),Fr(sys.argv[2]); E=e1+e2
Nz=10
def mul(a,b):
    c=[Fr(0)]*(Nz+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=Nz: c[i+j]+=x*y
    return c
def exp_ser(a):
    # a[0]=0
    r=[Fr(0)]*(Nz+1); r[0]=Fr(1)
    for n in range(1,Nz+1):
        r[n]=sum(k*a[k]*r[n-k] for k in range(1,n+1))/n
    return r
# Theta(x)=exp sum_{m>=1} B_{2m} x^{2m-1}/(2m(2m-1)); Upsilon(z/e1,z/e2)=h_w(z)=Theta(-z/E)Theta(z/e1)Theta(z/e2)
def log_hw(sign):
    # coefficient series of log h_w(sign*z)
    a=[Fr(0)]*(Nz+1)
    for m in range(1,Nz):
        k=2*m-1
        if k>Nz: break
        co=bern(2*m)/(2*m*(2*m-1))
        a[k]=co*((-1/E)**k+(1/e1)**k+(1/e2)**k)*(sign**k)
    return a
Ups_minus=exp_ser(log_hw(-1))    # Upsilon(-z/e1,-z/e2)=h_w(-z)
geo=[(1/E)**k for k in range(Nz+1)]   # (1-z/E)^{-1}
A=mul(Ups_minus,geo)
Dil={d:( -A[d-1] if d>=2 else None) for d in range(2,Nz+1)}  # Cont(d)=[z^{a-1}](1-Ups (1-z/E)^-1); a>=2 only: a=1 gives 0
Dil[1]=Fr(0)
# edge: [z^m w^n](1-Ups(z)Ups(w))/(z+w)
num={}
for m in range(Nz+1):
    for n in range(Nz+1):
        num[(m,n)]=(1 if m==n==0 else 0)-Ups_minus[m]*Ups_minus[n]
def edge(m,n):
    s=Fr(0)
    for k in range(n+1):
        if m+k+1<=Nz: s+=(-1)**k*num[(m+k+1,n-k)]
    return s
def vertex(g,exps):
    # Sum_m 1/m! Sum_d <psi^exps psi^d> prod Dil(d), pure psi
    n=len(exps); total=3*g-3+n-sum(exps)
    res=Fr(0)
    m=0
    while True:
        # sum (d_k) = total+m, d>=2
        if 2*m>total+m: break
        def rec(m_,bud):
            if m_==0:
                if bud==0: yield ()
                return
            for d in range(2,bud-2*(m_-1)+1):
                for r in rec(m_-1,bud-d): yield (d,)+r
        for t in rec(m,total+m):
            w=Fr(1,factorial(m))
            for d in t: w*=Dil[d]
            res+=w*wk(g,tuple(sorted(tuple(exps)+t)))
        m+=1
    return res
def graph(genera,edges):
    nh=[0]*len(genera)
    for u,v in edges: nh[u]+=1; nh[v]+=1
    dims=[3*genera[v]-3+nh[v] for v in range(len(genera))]
    he=[]
    for u,v in edges: he+= [u,v]
    tot=Fr(0)
    def rec(i,cur,used):
        nonlocal tot
        if i==len(he):
            term=Fr(1)
            for e in range(len(edges)):
                term*=edge(cur[2*e],cur[2*e+1])
                if term==0: return
            for v in range(len(genera)):
                term*=vertex(genera[v],tuple(cur[k] for k in range(len(he)) if he[k]==v))
                if term==0: return
            tot+=term
            return
        v=he[i]
        for a in range(0,dims[v]-used[v]+1):
            used[v]+=a; cur.append(a); rec(i+1,cur,used); cur.pop(); used[v]-=a
    rec(0,[],[0]*len(genera))
    return tot
G2=[([2],[],1),([1],[(0,0)],2),([0],[(0,0),(0,0)],8),([1,1],[(0,1)],2),([1,0],[(0,1),(1,1)],2),([0,0],[(0,1),(0,1),(0,1)],12),([0,0],[(0,0),(0,1),(1,1)],8)]
s=sum(graph(g,ed)/a for g,ed,a in G2)
# beta_2
beta2=(7*e1**4+10*e1**2*e2**2+7*e2**4)/(5760*e1*e2)
print("BS graph sum g=2:",s," *( -E^3 e1 e2)=",s*(-E**3*e1*e2)," beta_2=",beta2)
# draft: Z[2,0] = E^{-3} Z_2(tau), Z_2=-(24+28s+7s^2)/5760, s=1/sigma
sig=-e1*e2/E**2; sinv=1/sig
Z=-(24+28*sinv+7*sinv**2)/5760/E**3
print("Z[2,0] draft:",Z," equal to BS graph sum:",Z==s)
