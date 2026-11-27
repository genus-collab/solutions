# Extensions of rstar.py, exact Fractions.
# (a) genus-3 divisor identity for R and R^* through y^1.
# (b) arbitrary diagonal constants: R~ X with X = diag(exp(a_1 z + a_3 z^3 + a_5 z^5)) still satisfies the divisor identity
#     (this is the ambiguity of A2.9, and shows that flatness alone cannot fix a_k); R~ X with a t-dependent X, not tested.
# (c) left multiplication control: rows with the wrong sign in the theta term of the rows (R(z) -> R(-z)) fail.
# (d) genus 2, y^1 with X arbitrary: value differs from localisation, so the constants are really seen.
import rstar as rs
from fractions import Fraction as Fr
from rstar import *

def theory_X(Z,aa):
    "rows of the true theory, multiplied on the right by diag(x_j(z)/h...) with extra odd exponent aa[j][k]"
    C=[rows(j,Z) for j in range(3)]
    Ds=[c['Delta'] for c in C]; P=[c['P'] for c in C]
    Pn=[]
    for j in range(3):
        e=[Fr(0)]*(Z+1)
        for k,v in aa[j].items():
            if k<=Z: e[k]=v
        x=[Fr(0)]*(Z+1); x[0]=Fr(1)
        for n in range(1,Z+1): x[n]=sum(k*e[k]*x[n-k] for k in range(1,n+1))/n
        row=[]
        for a in range(3):
            r=[zero[:] for _ in range(Z+1)]
            for n in range(Z+1):
                for m in range(n+1):
                    r[n]=add(r[n],scal(x[m],P[j][a][n-m]))
            row.append(r)
        Pn.append(row)
    return wrap(Pn,Ds,Z)
def wrap(P,Ds,Z):
    t=Th(); t.P=P; t.D=Ds; t.Z=Z
    t.E={}
    for i in range(3):
        for j in range(3):
            num={}
            for m in range(Z+1):
                for b in range(Z+1-m):
                    s=zero[:]
                    for (c,d),et in ETA.items():
                        s=add(s,scal(et*(-1)**(m+b),mul(P[i][c][m],P[j][d][b])))
                    s=scal(Fr(-1),s)
                    if m==0 and b==0 and i==j: s=add(s,Ds[i])
                    num[(m,b)]=s
            assert all(x==0 for x in num[(0,0)])
            e={}
            for tot in range(Z):
                for a in range(tot,-1,-1):
                    b=tot-a
                    e[(a,b)]=add(num[(a+1,b)],scal(Fr(-1),e[(a+1,b-1)])) if b>0 else num[(a+1,0)]
            t.E[(i,j)]=e
    t.vcache={}
    return t
def flipped(Z):
    "R(-z): rows with z -> -z. Symplectic, but not flat for the true divisor."
    C=[rows(j,Z) for j in range(3)]
    Ds=[c['Delta'] for c in C]; P=[c['P'] for c in C]
    Pn=[[[scal(Fr((-1)**n),P[j][a][n]) for n in range(Z+1)] for a in range(3)] for j in range(3)]
    return wrap(Pn,Ds,Z)
def report(name,t,gmax2=True):
    print("=====",name)
    C4=gsum(t,0,[1,1,1,1]); C5=gsum(t,0,[1,1,1,1,1])
    print(" g=0 defect:",[str(x) for x in add(mul(A,C5),scal(Fr(-1),th(C4)))],flush=True)
    F11=gsum(t,1,[1]); F12=gsum(t,1,[1,1])
    print(" g=1 defect:",[str(x) for x in add(mul(A,F12),scal(Fr(-1),th(F11)))],flush=True)
    F2=gsum(t,2,[]); F21=gsum(t,2,[1])
    print(" g=2 F_2:",[str(x) for x in F2])
    print(" g=2 defect:",[str(x) for x in add(mul(A,F21),scal(Fr(-1),th(F2)))],flush=True)
if __name__=="__main__":
    Z=ZM
    aa=[{1:Fr(2,7),3:Fr(-5,11),5:Fr(3,13)},{1:Fr(-1,3),3:Fr(1,2),5:Fr(7,5)},{1:Fr(4,9),3:Fr(2,5),5:Fr(-3,8)}]
    report("R X, X arbitrary odd diagonal constants",theory_X(Z,aa))
    report("R(-z)",flipped(Z))
    R=theory(Z+2,False) if False else None
