# Critic: L1 to t^10 and L3 (bi-degree form) to t^6 at a fresh sample directly from the Wick rows; no closed form used.
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
from cylib import *

a0,a1,a2,e1 = Fr(2,9),Fr(11,4),Fr(-7,6),Fr(3,5)
e2 = a0+a1+a2-e1
w = weights(a0,a1,a2,e1,e2)
N,L = 18,13
t0=time.time()
R = build_rows(w,N,L)
print('rows',time.time()-t0,flush=True)
cls=[w.Ee,(-w.e1)%P,(-w.e2)%P]
rows,_ = framed_rows(R,cls)
E=w.Ee
# rho: Delta_1 t^2/(-e1 e2 E^3) = rho^-2
vD,sD=R.D1
assert vD==-2
c=inv((-w.e1*w.e2%P*pow(E,3,P))%P)
rm2=scal(c,sD)
assert rm2[0]==1
rho=sinv(ssqrt(rm2,1,L),L)
print('rho_1 == u/(2E^2):', rho[1]==((E-2*w.a1)*(E-2*w.a2)%P*inv(2*E*E%P))%P)
# rho q_t: psi^d coefficient, d = 1..N+1
Q=[]
for d in range(1,N+2):
    s=mul(rho,rows[1][0][d-1][1],L)
    Q.append([((-1)**d*x)%P for x in s])
def annih(k):
    pol=[1]
    for j in range(1,k+2):
        new=[0]*(len(pol)+1)
        for i,c_ in enumerate(pol):
            new[i]=(new[i]+E*c_)%P
            new[i+1]=(new[i+1]-j*c_)%P
        pol=new
    return pol
for k in range(0,11):
    ser=[0]+[Q[d-1][k] for d in range(1,N+2)]  # psi^0..psi^(N+1)
    pol=annih(k)
    out=[sum(pol[i]*ser[m-i] for i in range(min(m,len(pol)-1)+1))%P for m in range(N+2)]
    top=max([m for m,x in enumerate(out) if x],default=-1)
    print('k',k,'[t^k](rho q_t) x prod(E-j psi): degree',top,'(bound k =',k,') tested through psi^%d'%(N+1),' limit const part ok' if (top<=k or k==0) else 'FAIL',flush=True)
# L3 bi-degree
M=14
Ed=edge_table(R,rows,1,1,M)
r2=mul(rho,rho,L)
for k in range(0,7):
    pol=annih(k)
    def e(m,n):
        v=lmul((0,r2),Ed[(m,n)],L)
        return lcoef(v,k)
    # multiply in psi and psi', test vanishing outside [0,k]x[0,k]
    bad=0;tested=0
    for m in range(0,M+1):
        for n in range(0,M+1-m):
            s=0
            for i,a in enumerate(pol):
                for j,b in enumerate(pol):
                    if m-i>=0 and n-j>=0:
                        s+=a*b*e(m-i,n-j)
            s%=P
            if m>k or n>k:
                tested+=1
                if s: bad+=1
    print('k',k,'E_11 bi-degree test: violations',bad,'of',tested,flush=True)
print('time',time.time()-t0)
