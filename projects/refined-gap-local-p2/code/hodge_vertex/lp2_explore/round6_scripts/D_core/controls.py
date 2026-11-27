# Controls: every membership tester must be able to fail.  Exact.
import sys
from fractions import Fraction as Fr
from alg import *
from core import *

def scan(sig,E,K,label):
    out=[]
    for k in range(0,K+1):
        for comp,bound in ((0,k+1),(1,k+1),(2,k)):
            q,o,ok=leaves(R.lift(sig[comp].get(k)),E)
            bad=(not ok) or bool(q) or (o and (min(o)<1 or max(o)>bound))
            if bad: out.append((comp,k,'poly',q[:3] if q else [],'ok',ok,'top',max(o) if o else 0))
    print(' ',label,': failures of H2(b) =',len(out),'first:',out[:3])
    return out

def K1(sig,eta,e1,e2,K):
    E=e1+e2
    bad=[]
    for n in range(-3,K+1):
        tot=R([0])
        for k1 in range(-2,n+3):
            k2=n-k1
            if k2<-2: continue
            for (p,qq),m in eta.items():
                tot=tot+R.lift(sig[p].get(k1))*R.lift(sig[qq].get(k2)).negarg()*m
        want=3*e1*e2*E**3 if n==-2 else Fr(0)
        if not (tot==want): bad.append(n)
    return bad

def E00_low(sig,eta,E,nrange,BW=5):
    # rho^2 E_00 at order n, [psi'^0..BW]: returns dict n -> list of Rs (from (1.1))
    tayl={(a,k):R.lift(sig[a].get(k)).taylor(BW) for a in range(3) for k in range(-3,max(nrange)+3)}
    out={}
    for n in nrange:
        Ns=[]
        for b in range(BW+1):
            x=R([0])
            for k1 in range(-3,n+4):
                k2=n-k1
                if k2<-3 or k2>max(nrange)+2: continue
                for (p,qq),m in eta.items():
                    f=sig[p].get(k1)
                    if iszero(f): continue
                    x=x+R.lift(f)*(-m*tayl[qq,k2][b])
            Ns.append(x)
        Hs=[];Hp=R([0])
        for b in range(BW+1):
            val=Fr(0)
            for m_ in range(b+1): val+=(-1)**(b-m_)*Ns[m_].taylor(b-m_)[b-m_]
            H=(Ns[b]-val-Hp)/R([0,1]); Hs.append(H); Hp=H
        out[n]=Hs
    return out

e1,e2=Fr(-2,5),Fr(1); E=e1+e2; K=4; prec=K+8
eta=eta_of(e1,e2)
print('eps',e1,e2,'E',E)
# baseline
B=Base(E,prec); r=rec_r(prec+6); S,a,sig=rows0(B,r)
print('C0 baseline:'); scan(sig,E,K,'on-shell'); print('  K=1 failures',K1(sig,eta,e1,e2,K))
# C1 off-shell consistent section c1=C1cf(1+Q Delta)
Q=Fr(1,10)
B1=Base(E,prec,fsec=[Fr(1),Q]); S1,a1,sig1=rows0(B1,r)
print('C1 off-shell section f=1+Q Delta, Q=%s:'%Q); bad=scan(sig1,E,K,'off-shell'); print('  K=1 failures',K1(sig1,eta,e1,e2,K))
# formula E Delta C1cf^2 theta(f)/t_c
thf=B1.theta(B1.f)
form=(((B1.Delta*B1.C1cf*B1.C1cf*thf)*E)).shift(-1)
okf=True
for k in range(-1,K+1):
    q,o,ok=leaves(R.lift(sig1[2].get(k)),E)
    lin=q[1] if len(q)>1 else Fr(0)
    if lin!=form.get(k): okf=False; print('  mismatch at',k,lin,form.get(k))
print('  H2(d) psi-linear coefficient = E Delta C1cf^2 theta f/t_c for k=-1..%d:'%K,okf,'values',[form.get(k) for k in range(0,3)])
# C2 perturbed r_2
r2=rec_r(prec+6,perturb={2:Fr(1,7)})
S2,a2,sig2=rows0(B,r2)
print('C2 r_2 perturbed:'); scan(sig2,E,K,'r2+1/7'); print('  K=1 failures',K1(sig2,eta,e1,e2,K))
qz=[]
for k in range(0,K+1):
    q,o,ok=leaves(R.lift(sig2[0].get(k)),E); qz.append(sum((c/j for j,c in o.items()),Fr(0)))
print('  sum beta_j/j (want E,0,0,..):',qz)
# C3 XX shift
d=Fr(1,10)
S3,a3,sig3=rows0(B,r,xxshift=d)
print('C3 XX shifted by %s (no change of c1):'%d); scan(sig3,E,K,'XXshift'); print('  K=1 failures',K1(sig3,eta,e1,e2,K))
tb=E00_low(sig3,eta,E,[-2,-1,0])
print('  rho^2 E00 [t_c^-1] psi^0 psi\'^0..2:',[tb[-1][b].taylor(1)[0] for b in range(3)],' expected pole 3 e1 e2 E^2 delta =',3*e1*e2*E*E*d)
# C4 z->+psi
S4,a4,sig4=rows0(B,r,sgn=+1)
print('C4 rows at z=+psi instead of -psi:'); scan(sig4,E,K,'z=+psi')
# C5 wrong normalisation of section, f=2
B5=Base(E,prec,fsec=[Fr(2)])
S5,a5,sig5=rows0(B5,r)
print('C5 c1 = 2 C1cf (d2 not analytic): d2 lowest order',B5.d2.lo(),' sigma_2 lowest order',sig5[2].lo())
scan(sig5,E,K,'f=2')
# C6 wrong leaf scale
print('C6 leaf scale E+1 used in tester:')
scan(sig,E+1,K,'scale E+1')
# C7 one level below the bound
fails=0;tot=0
for k in range(0,K+1):
    for comp,bound in ((0,k+1),(1,k+1),(2,k)):
        if bound==0: continue
        tot+=1
        if not in_Lam(R.lift(sig[comp].get(k)),E,bound-1): fails+=1
print('C7 membership in the next lower space fails in %d of %d cases'%(fails,tot))
# C8 sign c1 -> -c1: f=-1. H2 up to sign of rows 1,2 ; edges unchanged
B8=Base(E,prec,fsec=[Fr(-1)]); S8,a8,sig8=rows0(B8,r)
print('C8 c1=-C1cf: sigma_1 -> -sigma_1, sigma_2 -> -sigma_2:',all(R.lift(sig8[1].get(k))==-R.lift(sig[1].get(k)) and R.lift(sig8[2].get(k))==-R.lift(sig[2].get(k)) for k in range(-1,K+1)))
# C9 sharper Z-infty control: replace C1cf-coordinate t_c by Delta itself in rho: rho'=-Ev/Delta => const appears
print('C9 rho built with Delta_c instead of t_c: q=rho\'q_t constants')
rhop=(B.v0*(-E))*B.Delta.inv()
Sp=rhop*S
for k in range(0,4):
    f=R([0,-1])*Sp.get(k)
    q,o,ok=leaves(f,E)
    print('   k',k,'const part',q[:1])
