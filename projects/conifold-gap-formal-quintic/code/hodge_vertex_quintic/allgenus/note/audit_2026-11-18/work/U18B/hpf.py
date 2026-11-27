# Exact HPF recursion (lem:gau-HPF) with ph_0=1 and free coefficient f_m = [z^{5m}]G_inf; compare Phi(z,1) with h(z).
from fractions import Fraction as F
import math
NMAX=10
def padd(a,b):
    n=max(len(a),len(b)); return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
def pscale(a,c): return [c*x for x in a]
def pmulu(a,k=1): return [F(0)]*k+a   # times u^k
def theta(a): return [i*x for i,x in enumerate(a)]
def one_minus_u(a): return padd(a,pmulu(pscale(a,-1)))
# series in zh: list of polys
def O(e,Phi):
    # O_e = 1 + zh(1-u)(4/5 theta_zh - e/5 - theta_u)
    out=[p[:] for p in Phi]+[[F(0)]]
    for m,p in enumerate(Phi):
        q=padd(pscale(p,F(4*m,5)-F(e,5)),pscale(theta(p),-1))
        q=one_minus_u(q)
        out[m+1]=padd(out[m+1],q)
    return out
def shift(Phi,k,ucoef,upow):  # ucoef*u^upow*zh^k*Phi
    out=[[F(0)] for _ in range(k)]+[pscale(pmulu(p,upow),ucoef) for p in Phi]
    return out
def lhs(Phi):
    P=[Phi]
    for j in range(1,6): P.append(O(-j,P[-1]))
    terms=[P[5],[pscale(p,-1) for p in Phi]]
    br=[shift(P[4],1,F(3),0),shift(P[3],2,F(17,5),1),shift(P[2],3,F(9,5),2),shift(P[1],4,F(274,625),3),shift(Phi,5,F(24,625),4)]
    tot=[]
    for t in terms:
        for i,p in enumerate(t):
            while len(tot)<=i: tot.append([F(0)])
            tot[i]=padd(tot[i],p)
    for b in br:
        for i,p in enumerate(b):
            while len(tot)<=i: tot.append([F(0)])
            tot[i]=padd(tot[i],pscale(one_minus_u(p),-1))
    return tot
# Bernoulli and G_inf free coefficients
def bern(n):
    B=[F(0)]*(n+1); B[0]=F(1)
    for m in range(1,n+1):
        B[m]=-sum(math.comb(m+1,k)*B[k] for k in range(m))/F(m+1)
    return B
B=bern(60)
def bernpoly(n,x): return sum(math.comb(n,k)*B[k]*x**(n-k) for k in range(n+1))
# log G_inf = sum_l (-1)^{l+1} B_{5l+1}(1/5)/(l(5l+1)) z^{5l}
lg=[F(0)]*(NMAX+1)
for l in range(1,NMAX//5+1): lg[5*l]=F((-1)**(l+1))*bernpoly(5*l+1,F(1,5))/(l*(5*l+1))
G=[F(0)]*(NMAX+1); G[0]=F(1)
for n in range(1,NMAX+1): G[n]=sum(m*lg[m]*G[n-m] for m in range(1,n+1))/n
free={m:G[5*m] for m in range(1,NMAX//5+1)}  # free coeff f_m = [z^{5m}] G_inf (cor:gam-free)
Phi=[[F(1)]]
for n in range(1,NMAX+1):
    Phi.append([F(0)])
    L=lhs(Phi)
    r=L[n] if n<len(L) else [F(0)]
    # coefficient zh^n of lhs = (1-u)[(4(n-1)-5theta)ph_{n-1}+r_{n-1}] ... we instead solve for ph_n from coefficient zh^{n+1}? use direct:
    # coefficient of zh^{n+1} involves ph_n linearly: (1-u)(4n-5theta)ph_n + (1-u) r_n
    deg=(4*n)//5+2
    # determine r_n by setting ph_n=0
    r=L[n+1] if n+1<len(L) else [F(0)]
    # divide by (1-u): r = (1-u)*rn
    rn=[]; rem=r[:]+[F(0)]*5
    acc=F(0)
    for i in range(len(rem)):
        acc+=rem[i]; rn.append(acc)   # (1-u)^{-1} series; should terminate
    # check termination
    while rn and rn[-1]==0: rn.pop()
    ph=[]
    for k in range(len(rn)+2):
        rk=rn[k] if k<len(rn) else F(0)
        if 4*n==5*k:
            assert rk==0, ("resonance fails",n,k,rk)
            ph.append(free[n//5])
        else:
            ph.append(-rk/(4*n-5*k))
    while ph and ph[-1]==0: ph.pop()
    Phi[n]=ph if ph else [F(0)]
# h(z)
# T_k = 5^{1-2k} + sum_j (zeta^j-1)^{1-2k}; power sums of y=1/(zeta^j-1), roots of 5y^4+10y^3+10y^2+5y+1
e=[F(1),F(-10,5),F(10,5),F(-5,5),F(1,5)]  # monic y^4+2y^3+2y^2+y+1/5 : coefficients c1..c4 = 2,2,1,1/5
c=[F(1),F(2),F(2),F(1),F(1,5)]
ps={0:F(4)}
for m in range(1,2*NMAX+2):
    s=F(0)
    for i in range(1,5):
        if m-i>0: s-=c[i]*ps[m-i]
        elif m-i==0: s-=i*c[i]
    ps[m]=s
lh=[F(0)]*(NMAX+1)
for k in range(1,NMAX//2+2):
    if 2*k-1<=NMAX:
        Tk=F(5)**(1-2*k)+ps[2*k-1]
        lh[2*k-1]=B[2*k]/(2*k*(2*k-1))*Tk
h=[F(0)]*(NMAX+1); h[0]=F(1)
for n in range(1,NMAX+1): h[n]=sum(m*lh[m]*h[n-m] for m in range(1,n+1))/n
for n in range(0,NMAX+1):
    v=sum(Phi[n]); print(n,"Phi_n(u)=",[str(x) for x in Phi[n]],"Phi_n(1)=",v,"h_n=",h[n],"equal" if v==h[n] else "DIFF")
