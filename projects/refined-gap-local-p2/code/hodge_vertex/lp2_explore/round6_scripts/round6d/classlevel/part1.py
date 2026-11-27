# Part 1: defect M of A2.8 computed directly, local P2, generic alpha.
# Field: Q (samples with all W_i^0 perfect squares, so sqrt W_i^0 is rational) or F_p with sqrt adjoined as elements of F_p.
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from rows_p2 import *
from fractions import Fraction as Fr
from itertools import product

def run(F,alpha,eps,D,K,label,verbose=True):
    t0=time.time()
    P=P2Data(F,alpha,eps,D,K)
    t=P.tt
    roots=[fsqrt(P.W0(i)) for i in range(3)]
    if any(r is None for r in roots):
        return None
    rr=P.rows_rational()
    R=[[rr[a][j]*(F(1)/roots[j]) for j in range(3)] for a in range(3)]
    Psi=[[R[a][j].zcoef(0) for j in range(3)] for a in range(3)]
    res={}
    Lam=[[ (P.Ainv*P.yb(P.Y[i],t) if i==j else const(t,0)) for j in range(3)] for i in range(3)]
    res['C Psi=Psi Lambda']=miszero(msub(mmul(P.Cm,Psi),mmul(Psi,Lam)))
    etai=mconst(t,P.etai_f)
    res['Psi^T etai Psi=1']=miszero(msub(mmul(mmul(mT(Psi),etai),Psi),mid(t,3)))
    Rm=mmap(lambda a:a.zneg(),R)
    res['R(-z)^T etai R(z)=1']=miszero(msub(mmul(mmul(mT(Rm),etai),R),mid(t,3)))
    ok=True
    for j in range(3):
        hh=P.h(j)
        for n in range(K+1):
            for a in range(3):
                if R[a][j].get(0,n)!=Psi[a][j].get(0,0)*hh.get(0,n): ok=False
    res['R|Q=0=Psi|Q=0 h']=ok
    res['Psi0=Delta^-1/2']=all((Psi[0][j]*Psi[0][j]*P.Delta[j]-const(t,1)).is_zero() for j in range(3))
    De=P.A*P.yb(P.C11[2],t); eE=P.A*P.yb(P.C11[1],t)
    X=binv(const(t,1)+ymono(t,1,27))
    res['A^2 D=X (comparison)']=(P.A*P.A*De-X).is_zero()
    res['2 e1 E=e1(D-X) (comparison)']=(eE*F(2)-(De-X)*P.e[0]).is_zero()
    res['C^0_11=0']=all(x==0 for x in P.C11[0])
    Hs=P.Hs
    eta=P.eta
    Psii=minv(Psi)
    Rend=mmul(mmul(etai,R),mmul(Psii,eta))   # End R(z) in the flat basis
    B=minv(Rend)
    M,z0=defect_endo(B,Psi,eta,Hs,P.thetaQ)
    res['[H*,B] has no psi^0']=z0
    res['M=0']=miszero(M)
    # a^i(psi;H^a) = R_{ai}(-psi): A2.7
    A=mmul(mmul(Psii,eta),B)
    res['a^i(psi;H^a)=R_ai(-psi)']=miszero(msub(A,mT(Rm)))
    # equation in the form closing A2.8: z thQ(R Psi)=(H*R-RH*)Psi, Psi as vectors = etai Psi
    Pv=mmul(etai,Psi)
    lhs=mmap(lambda a:a.zshift(1),mmap(P.thetaQ,mmul(Rend,Pv)))
    rhs=mmul(msub(mmul(Hs,Rend),mmul(Rend,Hs)),Pv)
    res['z thQ(R Psi)=[H*,R]Psi']=miszero(mmap(lambda a:a.trunc(K),msub(lhs,rhs)))
    ctr={}
    def report(name,Mx):
        nz=mzero_report(Mx)
        if not nz:
            ctr[name]='M=0'
        else:
            ctr[name]='M!=0 (%d nonzero coeffs; psi-orders %s; min y-order %d)'%(len(nz),sorted(set(e[3] for e in nz)),min(e[2] for e in nz))
    # R^*(z)=R^-1(-z) in place of R: legs B*(psi)=R^*(psi)^-1=End R(-psi)
    Bstar=mmap(lambda a:a.zneg(),Rend)
    Mx,z=defect_endo(Bstar,Psi,eta,Hs,P.thetaQ); report('R^*(z)=R^-1(-z)',Mx)
    # order of first failure in the z^2 row
    Bneg=mmap(lambda a:a.zneg(),B)
    Mx,z=defect_endo(Bneg,Psi,eta,Hs,P.thetaQ); report('R(-z)',Mx)
    Mx,z=defect_endo(B,Psi,eta,mscal(Hs,F(-1)),P.thetaQ); report('H*->-H*',Mx)
    Mx,z=defect_endo(B,Psi,eta,Hs,lambda a:-P.thetaQ(a)); report('theta->-theta',Mx)
    Psi00=mconst(t,[[Psi[a][j].get(0,0) for j in range(3)] for a in range(3)])
    Mx,z=defect_endo(B,Psi00,eta,Hs,P.thetaQ); report('constant frame Psi(0) in M',Mx)
    for hs,nm in ((-1,'h->h^-1 (expected M=0, A2.9)'),(0,'no h (expected M=0, A2.9)')):
        rr2=P.rows_rational(hsign=hs)
        R2=[[rr2[a][j]*(F(1)/roots[j]) for j in range(3)] for a in range(3)]
        R2end=mmul(mmul(etai,R2),mmul(Psii,eta))
        Mx,z=defect_endo(minv(R2end),Psi,eta,Hs,P.thetaQ); report(nm,Mx)
    if verbose:
        print('---',label,'alpha',alpha,'eps',eps,'D',D,'K',K,'field',F.__name__,'time %.1f'%(time.time()-t0),flush=True)
        print('  checks:',res)
        for k,v in ctr.items(): print('  ctl',k,':',v)
    return res,ctr

if __name__=='__main__':
    D=int(sys.argv[1]) if len(sys.argv)>1 else 3
    K=int(sys.argv[2]) if len(sys.argv)>2 else 5
    found=[]
    rng=[Fr(n) for n in range(-14,15) if n!=0]
    es=[(Fr(a),Fr(b)) for a in range(1,8) for b in list(range(-7,0))+list(range(1,8)) if a+b!=0 and (b<0 or a<=b)]
    for al in product(rng,repeat=3):
        if not (al[0]<al[1]<al[2]): continue
        base=[3*al[i]*(al[(i+1)%3]-al[i])*(al[(i+2)%3]-al[i]) for i in range(3)]
        for ep in es:
            c=ep[0]*ep[1]
            if all(rat_sqrt(v*c) is not None for v in base):
                found.append((al,ep))
    print('rational samples with all W_i^0 squares:',len(found),[(tuple(map(str,a)),tuple(map(str,e))) for a,e in found[:8]],flush=True)
    for al,ep in found[:3]:
        run(Fr,al,ep,D,K,'Q (square sample)')
    samples=[((1,2,5),(3,7)),((2,-3,8),(1,4)),((1,3,8),(Fr(5,2),Fr(-3,7))),((Fr(3),Fr(7),Fr(-11,2)),(Fr(5,3),Fr(-2,7))),((2,-3,7),(1,4))]
    plist=primes3mod4(2**89,400)
    for al,ep in samples:
        cnt=0
        for p in plist:
            F=make_fp(p)
            P=P2Data(F,al,ep,1,1)
            if any(fsqrt(P.W0(i)) is None for i in range(3)): continue
            cnt+=1
            if cnt==1:
                run(F,al,ep,D,K,'F_p, p=%d'%p)
            else:
                out=run(F,al,ep,D,K,'F_p',verbose=False)
                print('  second prime p=%d: M=0 is %s; controls %s'%(p,out[0]['M=0'],[v[:4] for v in out[1].values()]),flush=True)
            if cnt==2: break
