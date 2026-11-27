from fractions import Fraction as Fr
import random
random.seed(7)
def rnd(): return Fr(random.randint(-9,9) or 1, random.randint(1,7))
ok=True
for trial in range(6):
    a0,a1,a2,z=rnd(),rnd(),rnd(),rnd()
    E=a0+a1+a2; A=a0+a2; B=a0+a1; u=(E-2*a1)*(E-2*a2)
    if 0 in (a1,a2,A,B,E,E+z,E+2*z,z): continue
    # M1
    H1=(E*z+E*(a2+A)-2*a2*A)/(E*(E+2*z))
    G1=(E/(E+z))*(H1-Fr(1,2))
    assert G1==u/(2*(E+z)*(E+2*z))
    r1=(E/(E+z))*(H1-(E**2+u)/(2*E**2))     # [t^1] r_{0,1}
    assert r1==-u*z/(E*(E+z)*(E+2*z))
    assert G1==r1+(u/(2*E**2))*(E/(E+z))
    # M2, M4: need a rational point with rational delta: choose s* rational and solve tau from N(s*)=0
    s=rnd()
    # N(s)=B tau s^2-(E+a0 tau)s+A=0 -> tau (B s^2 - a0 s) = E s - A
    if B*s*s-a0*s==0: continue
    tau=(E*s-A)/(B*s*s-a0*s)
    if tau==0 or s==1 or tau*s==1: continue
    delta=E+a0*tau-2*B*tau*s
    assert delta**2==(E+a0*tau)**2-4*A*B*tau
    L2=B*s-a0; U1=1-tau*s
    sp=s*L2/delta           # s*_tau
    U1p=-s-tau*sp
    deltap=a0-2*B*s-2*B*tau*sp
    Rtr=a0-2*B*s-B*tau*sp
    # 2 (log alg)_tau * delta
    two_logalg_delta = deltap + delta*U1p/U1 - delta*sp/s + delta*sp/(1-s)
    assert 2*Rtr==two_logalg_delta, "R_tr identity"
    # corrected line
    lhs=-s*(delta+tau*L2)/U1 - L2 + L2*s/(1-s) - deltap
    assert lhs==2*B*tau*s*L2/delta, "M2 corrected display"
    # printed (wrong) line would assert lhs + 2Btau s L2/delta == 2 B tau s L2/delta i.e. lhs==0
    assert lhs!=0 or True
    # final printed line
    assert -s*(delta+tau*L2)/U1 - L2 + L2*s/(1-s) - a0 + 2*B*s == 0
    # M4
    alg2=a1*A*delta*U1/(E**3*s*(1-s))
    assert U1/(1-s)==(E-B*tau*s)/a1
    assert alg2==A*delta*(E-B*tau*s)/(E**3*s)
    assert E-B*tau*s==delta+tau*L2
    # N at s=E/(B tau)
    s0=E/(B*tau)
    N=lambda x: B*tau*x*x-(E+a0*tau)*x+A
    assert N(s0)==a1*a2/B
    # M5: LofU vs L(U)
    Q=1/(1-tau)
    U=rnd()
    LofU=-(Q*B*U+Q*a2)/(1-Q)
    assert LofU==(B*U+a2)/tau
    LU=-a1*Q*U/(1-Q*U)
    LeqU=B*Q*U*U-(a0+Q*(a1-a2))*U-a2
    assert LofU-LU==-LeqU/(tau*(1-Q*U))
    # LeqU at 1/Q
    assert B*Q/Q**2-(a0+Q*(a1-a2))/Q-a2 == -a1*tau
    # M3
    assert E**2+u==2*a0*E+4*a1*a2
    # sharpness
    def P(k,E,a2,A): 
        p=Fr(1)
        for i in range(1,k+1): p*=((k+1)*a2-i*E)*((k+1)*A-i*E)
        return p
    print("trial ok", trial)
# sharpness checks: u=0 (E=2a2), odd k
a1_,a2_=Fr(3),Fr(5); E_=2*a2_; A_=E_-a1_
for k in range(1,12):
    p=Fr(1)
    for i in range(1,k+1): p*=((k+1)*a2_-i*E_)*((k+1)*A_-i*E_)
    print("u=0,E=2a2 k=",k,"P_k==0:",p==0)
a2_=Fr(7); E_=Fr(9); a1_=Fr(2); A_=E_-a1_
p=Fr(1)
for i in range(1,9): p*=(9*a2_-i*E_)*(9*A_-i*E_)
print("E=9a2/7, P_8==0:",p==0)
