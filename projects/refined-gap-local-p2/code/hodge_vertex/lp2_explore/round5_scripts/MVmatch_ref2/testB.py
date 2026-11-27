# (B) at genus 4 incl. non-integer leaves and parts 0,1; (E) Laurent degree + special tau + E-scaling (literal E-class).
from fractions import Fraction as F
from hodge import *
from testD import lag
from math import comb,factorial
import time

from collections import Counter
def M_units(g,tau,R,a):
    ts=omega_ts(tau);r=len(R);dim=3*g-3+r+a
    tot=F(0)
    for S,c in omega_monos(g,ts,dim).items():
        if c==0: continue
        rem=dim-sum(S)
        for tdeg in range(rem+1):
            for pt in parts_min(tdeg,1):
                if len(pt)>a: continue
                pad=pt+(0,)*(a-len(pt))
                mult=factorial(a)
                for v in Counter(pad).values(): mult//=factorial(v)
                for d in compositions(rem-tdeg,r):
                    w=F(1)
                    for l,e in zip(R,d): w*=F(l)**e
                    if w==0: continue
                    tot+=c*w*mult*ch(g,S,tuple(sorted(d+pad)))
    return tot
def Mseq(g,tau,R,amax):
    return [M_units(g,tau,R,a) for a in range(amax+1)]

def lagF(n,al,x):
    s=F(0)
    for i in range(n+1):
        s+=(-1)**i*comb(n+al,n-i)*F(x)**i/factorial(i)
    return s

for g,R,tau in [(4,(F(5,2),),F(5,3)),(4,(F(1,2),F(1,3)),F(-9,7)),(4,(0,),F(5,3)),(4,(1,),F(5,3)),(4,(0,3),F(2,9)),(5,(F(5,2),),F(5,3)),(4,(),F(5,3)),(3,(F(-7,3),2,F(1,5)),F(1,3))]:
    t0=time.time()
    r=len(R);J=sum(R) if R else F(0);N=3*g-3+r
    Ms=Mseq(g,tau,R,N+2)
    ok=True;cnt=0
    Xi_lag=sum(F((-1)**a,factorial(a))*lagF(N-a,2*g-2+r+a,-J-a)*Ms[a] for a in range(N+1))
    for m in range(N+3):
        lhs=Qk(g,tau,R,m) if g<=4 or m<=N else None
        if lhs is None: continue
        rhs=sum(Ms[a]/factorial(a)*(-1)**(m-a)*lagF(m-a,2*g-3+r+a,-J-a) for a in range(m+1))
        cnt+=1
        if lhs!=rhs: ok=False;print("  diff",m,lhs,rhs)
    xd=sum((-1)**k*Qk(g,tau,R,k) for k in range(N+1))
    # control : wrong Laguerre index
    ctl=sum(F((-1)**a,factorial(a))*lagF(N-a,2*g-3+r+a,-J-a)*Ms[a] for a in range(N+1))
    print("g",g,"R",R,"tau",tau,"coeffs ok",ok,cnt,"Xi eq",Xi_lag==xd,"Xi",xd,"ctl eq",ctl==xd,round(time.time()-t0,1),flush=True)

# string/lam0 at genus 3,4 literal (J noninteger)
for g,R,tau in [(3,(F(5,2),3),F(5,3)),(4,(F(3,2),),F(-9,7)),(3,(),F(5,3)),(4,(),F(2,5))]:
    D=(sum(R) if R else 0)-len(R); J=sum(R) if R else 0
    x=Xi_def(g,tau,R)
    print("g",g,"R",R,"Xi(R+1)=J Xi:",Xi_def(g,tau,tuple(R)+(1,))==J*x,"Xi(R+0)=(D-2g+2)Xi:",Xi_def(g,tau,tuple(R)+(0,))==(D-2*g+2)*x,flush=True)

# (E): Laurent polynomial; degree bound <= h in 1/sigma; tau <-> -1-tau; specialisation
def lagr_fit(xs,ys):
    # return coefficients of interpolating polynomial
    n=len(xs)
    A=[[x**k for k in range(n)]+[y] for x,y in zip(xs,ys)]
    for i in range(n):
        p=next(j for j in range(i,n) if A[j][i]!=0)
        A[i],A[p]=A[p],A[i]
        A[i]=[v/A[i][i] for v in A[i]]
        for j in range(n):
            if j!=i and A[j][i]!=0:
                f=A[j][i]
                A[j]=[u-f*v for u,v in zip(A[j],A[i])]
    return [A[i][n] for i in range(n)]

for h,R in [(2,(3,)),(3,(5,)),(3,(4,2)),(4,(7,))]:
    taus=[F(k,3) for k in range(1,h+5)]
    # g(w)=Xi as polynomial in w=-1/sigma, degree<=h expected
    w=[-1/(t*(1+t)) for t in taus]
    y=[Xi_def(h,t,R) for t in taus]
    co=lagr_fit(w,y)
    print("E Laurent h",h,"R",R,"deg<=h:",all(c==0 for c in co[h+1:]),"top",co[h],"coefs0",co[0],flush=True)
    print("  sym tau->-1-tau",Xi_def(h,F(2,3),R)==Xi_def(h,F(-5,3),R))

# Zero of Pi: in-range Xi vanish at tau=-1/2 (Pi_2=0), -1/3 (Pi_3=0), -1/4 (Pi_4=0) from the definition
for h,R,t in [(2,(2,),F(-1,2)),(3,(2,2),F(-1,2)),(3,(4,),F(-1,2)),(3,(3,2),F(-1,3)),(3,(3,),F(-1,3)),(4,(4,3),F(-1,4)),(3,(3,),F(-2,3))]:
    print("zero of Pi: h",h,"R",R,"tau",t,"Xi_def",Xi_def(h,t,R),flush=True)

# E-scaling literally: class with w=E(1,tau,-1-tau), leaves 1/(E-j psi), T=-psi^2/(E-psi)
def Xi_E(h,tau,R,E):
    # Omega_w= prod Lambda(-1/w_a)
    ts=[-1/(F(E)*w) for w in (F(1),tau,-1-tau)]
    r=len(R);N=3*h-3+r
    tot=F(0)
    from collections import Counter
    for k in range(N+1):
        dim=N+k
        for S,c in omega_monos(h,ts,dim).items():
            if c==0: continue
            rem=dim-sum(S)
            for tdeg in range(2*k,rem+1):
                ptl=[()] if (k==0 and tdeg==0) else ([p for p in parts_min(tdeg,2) if len(p)==k] if k>0 else [])
                for pt in ptl:
                    mult=factorial(k)
                    for v in Counter(pt).values(): mult//=factorial(v)
                    wT=F(1)
                    for b in pt: wT*=-F(E)**(1-b)
                    for d in compositions(rem-tdeg,r):
                        wl=F(1)
                        for j,e in zip(R,d):
                            wl*=F(j)**e*F(E)**(-1-e)
                        if wl==0: continue
                        tot+=c*wT*wl*mult*ch(h,S,tuple(sorted(d+pt)))/factorial(k)
    return tot
for h,R,tau,E in [(2,(3,),F(5,3),F(7,2)),(3,(5,),F(2,9),F(-4,3)),(3,(4,2),F(5,3),F(3))]:
    print("E-scale",h,R,Xi_E(h,tau,R,E)==E**(-(3*h-3+2*len(R)))*Xi_def(h,tau,R),
      "ctl3-2r",Xi_E(h,tau,R,E)==E**(3-2*len(R))*Xi_def(h,tau,R),flush=True)
    J=sum(R);D=J-len(R)
    print("  S_P2",Xi_E(h,tau,R+(1,),E)==F(J)/E**2*Xi_E(h,tau,R,E),"L0_P2",Xi_E(h,tau,R+(0,),E)==F(D-2*h+2)/E**2*Xi_E(h,tau,R,E),"ctl J/E",Xi_E(h,tau,R+(1,),E)==F(J)/E*Xi_E(h,tau,R,E),flush=True)
