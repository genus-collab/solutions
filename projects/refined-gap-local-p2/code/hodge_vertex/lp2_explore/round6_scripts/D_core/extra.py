# Extra exact tests: E-scaling of the rows (so "non-special eps" is one test), Lemma G gauge statement, Prop 3.1, sharpness, E00 telescoping.
import sys
from fractions import Fraction as Fr
from alg import *
from core import *

K=6; prec=K+7
r=rec_r(prec+6)

def rescale(f,lam):
    # f(psi) -> f(lam psi)
    f=R.lift(f)
    return R([c*lam**i for i,c in enumerate(f.n)],[c*lam**i for i,c in enumerate(f.d)])

# 1. scaling: sigma_a(psi;E)= E^a sigma_a^(E=1)(psi/E)
B1=Base(Fr(1),prec); _,_,sig1=rows0(B1,r)
ok=True
for E in (Fr(2,3),Fr(3,5),Fr(-2),Fr(-4,3),Fr(11,14)):
    B=Base(E,prec); _,_,sig=rows0(B,r)
    for a_ in range(3):
        for k in range(-2,K+2):
            if not (R.lift(sig[a_].get(k))==rescale(sig1[a_].get(k),1/E)*(E**a_)): ok=False; print('scale fail',E,a_,k)
print('1. sigma_a(psi;E) = E^a sigma_a(psi/E;1) through t_c^%d at 5 values of E:'%(K+1),ok)

# 2. Prop 3.1 row 2 (general c1, off and on shell)
for label,fs in (('on-shell',None),('off-shell f=1+Delta/10',[Fr(1),Fr(1,10)])):
    E=Fr(3,5)
    B=Base(E,prec,fsec=fs); S,a,sig=rows0(B,r)
    z=R([0,Fr(-1)])
    Th=lambda F:(B.Y0*F)+(B.theta(F)*z)
    XX=B.theta(B.c1)*B.C1i
    Rd=(B.Delta*B.theta(XX))+(B.Delta*XX*XX)+((B.Delta-LS({0:Fr(1)}))*XX)+((B.Delta-LS({0:Fr(1)}))*Fr(2,9))
    Di=B.Delta.inv()
    r0=-E**3/27
    s=-1*(B.c1*((B.d2*E)+(B.Delta*B.d2*B.d2)+LS({0:E*E/3})))
    t=-1*(B.C1i*((Di*E)+B.d2))
    lhs=B.C1i*Th(a[2])
    rhs=(a[0]*r0)+(s*a[1])+(t*a[2])-((B.c1*Rd*a[1])*(z*z))
    d=lhs-rhs
    print('2. Prop 3.1 row 2,',label,': difference zero through order',d.prec-1,':',all(iszero(d.get(k)) for k in range(-3,d.prec)),' Rd lowest order',Rd.lo() if Rd.c else 'zero to prec %d'%Rd.prec)

# 3. Lemma G: hybrid gauge rows; pairing equal, but H2(b) "no constant" fails
e1,e2=Fr(-2,5),Fr(1); E=e1+e2
B=Base(E,prec); S,a,sig=rows0(B,r); eta=eta_of(e1,e2)
h1=B.c1*sig[1]
h2=(B.C1i*sig[2])+(((B.C1i*B.C1i-LS({0:Fr(1)}))*(E/2))*B.c1*sig[1])
sh=[sig[0],h1,h2]
print('3. hybrid gauge: sigma^hy^(-1)=',[R.lift(sh[a_].get(-1)).n for a_ in range(3)])
cons=[(a_,k,leaves(R.lift(sh[a_].get(k)),E)[0]) for a_ in (1,2) for k in range(0,3)]
print('   polynomial parts of hybrid rows k=0,1,2:',[(a_,k,[str(x) for x in q]) for a_,k,q in cons])
# edge pairing invariance at two independent psi, psi' : compare coefficient by coefficient
okp=True
for n in range(-2,3):
    tot=R([0])
    for k1 in range(-1,n+2):
        k2=n-k1
        for (p,qq),m in eta.items():
            A=R.lift(sig[p].get(k1)); Bb=R.lift(sig[qq].get(k2)).taylor(3)
            A2=R.lift(sh[p].get(k1)); Bb2=R.lift(sh[qq].get(k2)).taylor(3)
            for b in range(4):
                pass
    for b in range(4):
        x=R([0]); y=R([0])
        for k1 in range(-1,n+2):
            k2=n-k1
            for (p,qq),m in eta.items():
                x=x+R.lift(sig[p].get(k1))*(m*R.lift(sig[qq].get(k2)).taylor(b)[b])
                y=y+R.lift(sh[p].get(k1))*(m*R.lift(sh[qq].get(k2)).taylor(b)[b])
        if not (x==y): okp=False
print('   <sigma,sigma\'> equal in both gauges (orders -2..2, psi\'^0..3):',okp)

# 4. Sharpness: top leaves of rho q_t and sigma at k<=7, and nonvanishing
E=Fr(3,5)
B=Base(E,prec); S,a,sig=rows0(B,r)
print('4. r_k >0 for k<=40:',all(x>0 for x in rec_r(40)))
for k in range(1,6):
    q,o,ok=leaves(R([0,-1])*sig[0].get(k),E)
    print('   rho q_t k=%d support %s'%(k,sorted(o)))
