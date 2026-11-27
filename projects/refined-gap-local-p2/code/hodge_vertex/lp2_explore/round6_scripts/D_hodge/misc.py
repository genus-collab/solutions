# Extra checks: E00 loop support and t_c=0 value, chi_j, Dil(d;0), specialisations incl. diagonal and sign choices. Evidence only.
import sys
from fractions import Fraction as Fr
from itertools import product
from rows import *
from graphs_check import Vtx,graph_value,recon,GR3
from leaf_check_lib import solve_mod
def run(e1,e2,verbose=True):
    th=Theory(e1,e2,N=8,K=14)
    E=fp(th.E); tc=L(1,[1]+[0]*(th.prec-1))
    rho=(th.v[0]*(-E))/tc
    lam=lambda j,a: pow(j,a,P)*inv(pow(E,a+1,P))%P
    M=9
    ok_all=True
    for k in range(0,5):
        idx=[(i,l) for i in range(1,k+2) for l in range(1,k+2) if (i-1)+(l-1)<=k]
        A=[];b=[]
        for m in range(M):
            for n in range(M):
                if m+n>12: continue
                A.append([lam(i,m)*lam(l,n)%P for i,l in idx])
                b.append((rho*rho*th.edge(0,0,m,n)).coef(k))
        x=solve_mod(A,b)
        if x is None: ok_all=False
        if k==0:
            print(" E00(t_c=0) coefficient:",recon(x[0]),"expected",-th.e1*th.e2*th.E**4/27)
    print(" rho^2 E00 support (i-1)+(l-1)<=k, k<=4:",ok_all)
    # control: support k-1 fails
    bad=True
    k=3
    idx=[(i,l) for i in range(1,k+2) for l in range(1,k+2) if (i-1)+(l-1)<=k-1]
    A=[];b=[]
    for m in range(M):
        for n in range(M):
            if m+n>12: continue
            A.append([lam(i,m)*lam(l,n)%P for i,l in idx])
            b.append((rho*rho*th.edge(0,0,m,n)).coef(k))
    print(" control: smaller support solvable?",solve_mod(A,b) is not None)
    V=Vtx(th)
    worst=10**9
    for name,gen,ed in GR3:
        for col in product(range(3),repeat=len(gen)):
            worst=min(worst,graph_value(th,V,gen,ed,col).norm().v)
    print(" genus-3 sample graphs: min valuation all colourings:",worst)
    sv=V.vertex(3,0,())
    print(" single vertex g=3 polar:",[recon(sv.coef(k)) for k in range(-4,0)])
    a1,a2=th.e1,th.e2
    return
for e in [("4/5","1"),("1","1"),("5","-2"),("-7","3"),("1/8","1")]:
    print("eps",e); run(Fr(e[0]),Fr(e[1]))
