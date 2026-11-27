from fractions import Fraction as Fr
# weights (generic, not CY-restricted; e1,e2 irrelevant here)
a0,a1,a2 = Fr(5,3),Fr(-4,7),Fr(13,6)
E=a0+a1+a2; A=a0+a2; B=a0+a1; u=(E-2*a1)*(E-2*a2)
z=Fr(3,11)
H1=(E*z+E*(a2+A)-2*a2*A)/(E*(E+2*z))
G1_correct=(E/(E+z))*(H1-Fr(1,2))
G1_printed=(E/(E+z))*(H1-(E*E+u)/(2*E*E))
G1_final=u/(2*(E+z)*(E+2*z))
t1row=-u*z/(E*(E+z)*(E+2*z))
print("G1 (H1-1/2) == u/(2(E+z)(E+2z)):", G1_correct==G1_final)
print("printed (H1-(E^2+u)/2E^2) == u/(2(E+z)(E+2z)):", G1_printed==G1_final)
print("printed == [t^1] r01 = -uz/(E(E+z)(E+2z)):", G1_printed==t1row)
# R_tr intermediate display: pick tau with rational delta
# delta^2 = (E+a0 tau)^2 - 4AB tau ; search small rationals
import itertools
def sqrt_fr(q):
    from math import isqrt
    if q<0: return None
    n,d=q.numerator,q.denominator
    rn,rd=isqrt(n),isqrt(d)
    return Fr(rn,rd) if rn*rn==n and rd*rd==d else None
found=None
for p in range(-40,41):
    for q in range(1,41):
        tau=Fr(p,q)
        if tau==0: continue
        d2=(E+a0*tau)**2-4*A*B*tau
        r=sqrt_fr(d2)
        if r is not None and r!=0: found=(tau,r); break
    if found: break
tau,delta=found
s=(E+a0*tau-delta)/(2*B*tau)
L2=B*s-a0; U1=1-tau*s
s_tau=s*L2/delta
Rtr=a0-2*B*s-B*tau*s*L2/delta
delta_tau=a0-2*B*s-2*B*tau*s_tau
lhs_core=-s*(delta+tau*L2)/U1-L2+L2*s/(1-s)
print("tau,delta=",tau,delta)
print("printed display: core - delta_tau + X == X  (i.e. core - delta_tau == 0):", lhs_core-delta_tau==0)
print("correct display: core - delta_tau == 2B tau s L2/delta:", lhs_core-delta_tau==2*B*tau*s*L2/delta)
print("final line core - a0 + 2Bs == 0:", lhs_core-a0+2*B*s==0)
# identity R_tr = delta (log alg)_tau check via 2 R_tr/delta = (log alg^2)_tau using finite formula
twologalg_tau = delta_tau/delta + (-s-tau*s_tau)/U1 - s_tau/s + s_tau/(1-s)
print("2R_tr/delta == (log alg^2)_tau:", 2*Rtr/delta==twologalg_tau)
