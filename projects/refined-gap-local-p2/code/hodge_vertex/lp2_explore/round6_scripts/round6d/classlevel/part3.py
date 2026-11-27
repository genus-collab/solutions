# Part 3: rank 2. Resolved conifold x C^2, rows of leaf_structure_proof.md s.13 Prop 13.2(b), generic weights; defect M of A2.8 in rank 2.
# Field: Q (samples with e(wt_+), e(wt_-) squares) or F_p with the sqrt's as elements of F_p. Fresh code.
import sys,os,time,math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from fractions import Fraction as Fr
from itertools import product

def p_polys(N):
    """p_k(x), k=1..N, with p_1=x, p_{k+1}=x(1+x)p_k'; x=X/(1-X) under theta_X. coefficient lists."""
    out={1:[0,1]}
    for k in range(1,N):
        p=out[k]
        dp=[i*p[i] for i in range(1,len(p))]
        q=[0]*(len(dp)+2)
        for i,c in enumerate(dp):
            q[i+1]+=c;q[i+2]+=c
        out[k+1]=q
    return out

class Coni:
    def __init__(s,F,a0,a1,a2,e1,e2,D,K):
        s.F=F;a0,a1,a2,e1,e2=[F(x) for x in (a0,a1,a2,e1,e2)]
        s.a=(a0,a1,a2);s.e=(e1,e2);s.D=D;s.K=K
        s.t=Ctx(F,D,K);t=s.t
        Q=ymono(t,1,1)
        s.wtp=[a0,a1,a2,-e1,-e2];s.wtm=[-a0,a0+a1,a0+a2,-e1,-e2]
        s.ew=[]
        for w in (s.wtp,s.wtm):
            r=F(1)
            for x in w: r=r*x
            s.ew.append(r)
        # U_+ in F[[Q]], V=Q U_- in F[[Q]]
        U=const(t,-a2/a0)
        for _ in range(D+2):
            # Newton for Q(a0+a1)U^2-(a0+Q(a1-a2))U-a2=0
            f=Q*(a0+a1)*U*U-(const(t,a0)+Q*(a1-a2))*U-const(t,a2)
            df=Q*U*(2*(a0+a1))-(const(t,a0)+Q*(a1-a2))
            U=U-f*binv(df)
        V=const(t,a0/(a0+a1))
        for _ in range(D+2):
            f=V*V*(a0+a1)-(const(t,a0)+Q*(a1-a2))*V-Q*a2
            df=V*(2*(a0+a1))-(const(t,a0)+Q*(a1-a2))
            V=V-f*binv(df)
        s.U=U;s.V=V
        one=const(t,1)
        # P=QU/(1-QU), T=U/(1-U)
        Pp=(Q*U)*binv(one-Q*U); Tp=U*binv(one-U)
        Pm=V*binv(one-V)
        # T at U_-=V/Q: U/(1-U)=V/(Q-V)=-1/(1-Q/V)
        Tm=-binv(one-Q*binv(V))
        s.PT=[(Pp,Tp),(Pm,Tm)]
        # L_i
        Lp=-((Q*U)*(a0+a1)+Q*a2)*binv(one-Q)
        Lm=-(V*(a0+a1)+Q*a2)*binv(one-Q)
        s.L=[Lp,Lm]
        # Frobenius data, (i) of Prop 13.2(b)
        s.Lsum=-(const(t,a0)+Q*(a1+a2))*binv(one-Q)
        s.Lprod=-(Q*(a1*a2))*binv(one-Q)
        sig=a0+a1+a2
        s.Delta=[const(t,e1*e2)*(const(t,a0*a1*a2)+L*(2*a1*a2)-L*L*sig) for L in s.L]
        pref=F(1)/((a0+a1)*(a0+a2)*e1*e2)
        s.eta_f=[[sig/(a1*a2)*pref,pref],[pref,-a0*pref]]
        s.eta=mconst(t,s.eta_f)
        s.etai_f=inv_field(s.eta_f,F(1))
        s.etai=mconst(t,s.etai_f)
        s.Hs=[[const(t,0),-s.Lprod],[const(t,1),s.Lsum]]   # columns: images of 1, H
        s.Cm=mT(s.Hs)
        s.thetaQ=lambda x:x.theta()
        s.rows()
    def hv(s,v,sign=1):
        F=s.F;t=s.t
        ex=Bi(t)
        k=1
        while 2*k-1<=t.K:
            b=F(bern(2*k)/(2*k*(2*k-1)))
            p=F(0)
            for a in v: p=p+F(1)/a**(2*k-1)
            ex=ex+Bi(t,{(0,2*k-1):-b*p*sign})
            k+=1
        return bexp(ex)
    def wick(s,i):
        """r^pt_{0,i}=Sum_n ((2n)!/(n!2^n)) (-z/W_2)^n [v^{2n}] Exp[Sum_{k>=3} W_k v^k/(k! z)], (C7)"""
        F=s.F;t=s.t;K=s.K
        a0,a1,a2=s.a
        P,T=s.PT[i]
        Nk=3*K+6
        polys=p_polys(Nk)
        Ppow=[const(t,1)];Tpow=[const(t,1)]
        for m in range(Nk+2):
            Ppow.append(Ppow[-1]*P);Tpow.append(Tpow[-1]*T)
        W={}
        for k in range(2,Nk+1):
            # W_k=theta^k W=-a1 p_k(P)-(a0+a2) p_k(T), k>=2 (W_1 has the extra constant -a2)
            pk=polys[k]
            val=Bi(t)
            for m,c in enumerate(pk):
                if c:
                    val=val+(Ppow[m]*F(c))*(-a1)+(Tpow[m]*F(c))*(-(a0+a2))
            W[k]=val
        W2inv=binv(W[2])
        s.W2=getattr(s,'W2',{});s.W2[i]=W[2]
        # j factors: S_j(v)=(Sum_{k>=3} W_k v^k/k!)^j/j!: dict v-power->Bi
        res=const(t,1)   # z^0 term
        out=[const(t,1)]+[Bi(t) for _ in range(K)]
        # coefficient of z^m: Sum_{n} (2n)!/(n! 2^n) (-1/W_2)^n [v^{2n}] S_{n-m}
        Sj={0:{0:const(t,1)}}
        Jmax=2*K
        base={k:W[k]*(F(1)/F(math.factorial(k))) for k in range(3,Nk+1)}
        for j in range(1,Jmax+1):
            cur={}
            for vp,x in Sj[j-1].items():
                for k,y in base.items():
                    if vp+k>2*(K+j): continue   # 2n with n=m+j, m<=K
                    if vp+k in cur: cur[vp+k]=cur[vp+k]+x*y
                    else: cur[vp+k]=x*y
            cur={vp:x*(F(1)/F(j)) for vp,x in cur.items()}
            Sj[j]=cur
        pw=[const(t,1)]
        for n in range(3*K+2):
            pw.append(pw[-1]*(W2inv*F(-1)))
        for m in range(0,K+1):
            tot=Bi(t)
            for j in range(0,2*m+1):
                n=m+j
                if 2*n in Sj.get(j,{}):
                    co=F(math.factorial(2*n))/F(math.factorial(n)*2**n)
                    tot=tot+Sj[j][2*n]*pw[n]*co
            out[m]=tot
        r=Bi(t)
        for m,x in enumerate(out):
            r=r+x.zshift(m)
        return r
    def rows(s):
        F=s.F;t=s.t
        a0,a1,a2=s.a;e1,e2=s.e
        nw=[a1,a0+a2,-e1,-e2]
        hnw=s.hv(nw)
        s.v=[]
        s.rpt=[]
        for i in range(2):
            r0pt=s.wick(i)
            s.rpt.append(r0pt)
            r0=hnw*r0pt
            Dl=s.Delta[i]
            logd=Dl.theta()*binv(Dl)*(F(1)/F(2))
            r1=s.L[i]*r0+(r0.theta()-logd*r0).zshift(1)
            s.v.append([r0,r1])
def run(F,w,D,K,label,verbose=True):
    t0=time.time()
    C=Coni(F,*w,D,K)
    t=C.t
    roots=[fsqrt(C.ew[0]),fsqrt(C.ew[1])]
    if any(r is None for r in roots): return None
    res={}
    # Delta^{1/2}=sqrt(e(wt)) (Delta/e(wt))^{1/2}
    Dh=[bsqrt1(C.Delta[i]*(F(1)/C.ew[i]))*roots[i] for i in range(2)]
    Dhi=[binv(x) for x in Dh]
    R=[[C.v[j][a]*Dhi[j] for j in range(2)] for a in range(2)]    # R_{aj}, v_j=Delta_j^(1/2) R_{.j}
    Psi=[[R[a][j].zcoef(0) for j in range(2)] for a in range(2)]
    Lam=[[C.L[i] if i==j else const(t,0) for j in range(2)] for i in range(2)]
    res['C Psi=Psi Lambda']=miszero(msub(mmul(C.Cm,Psi),mmul(Psi,Lam)))
    res['Psi^T etai Psi=1']=miszero(msub(mmul(mmul(mT(Psi),C.etai),Psi),mid(t,2)))
    Rm=mmap(lambda a:a.zneg(),R)
    res['R(-z)^T etai R(z)=1']=miszero(msub(mmul(mmul(mT(Rm),C.etai),R),mid(t,2)))
    res['L_i is a root (L^2-Lsum L+Lprod=0)']=all((L*L-C.Lsum*L+C.Lprod).is_zero() for L in C.L)
    res['Delta_i(0)=e(wt_i)']=all(C.Delta[i].get(0,0)==C.ew[i] for i in range(2))
    # L_i=-a1 Q U/(1-QU) (note), Delta_i=1/eta(e_i,1) (13.2)
    a0,a1,a2=C.a
    ok=True
    one=const(t,1)
    P0=C.PT[0][0];P1=C.PT[1][0]
    res['L_i=-a1 QU/(1-QU)']=(C.L[0]+P0*a1).is_zero() and (C.L[1]+P1*a1).is_zero()
    # row 0: at Q=0 v_i=h_{wt_i}(1,L_i(0))
    ok=True
    wts=[C.wtp,C.wtm]
    for i in range(2):
        hh=C.hv(wts[i])
        for n in range(K+1):
            if C.v[i][0].get(0,n)!=hh.get(0,n): ok=False
            if C.v[i][1].get(0,n)!=hh.get(0,n)*C.L[i].get(0,0): ok=False
    res['v_i|Q=0=h_wt_i (1,L_i(0))']=ok
    # QDE row in the form (b)(ii): second-order equation for r_0 at each root; follows from the first order system
    Hs=C.Hs;eta=C.eta;etai=C.etai
    Psii=minv(Psi)
    Rend=mmul(mmul(etai,R),mmul(Psii,eta))
    B=minv(Rend)
    M,z0=defect_endo(B,Psi,eta,Hs,C.thetaQ)
    res['[H*,B] has no psi^0']=z0
    res['M=0']=miszero(M)
    A=mmul(mmul(Psii,eta),B)
    res['a^i(psi;H^a)=R_ai(-psi)']=miszero(msub(A,mT(Rm)))
    Pv=mmul(etai,Psi)
    lhs=mmap(lambda a:a.zshift(1),mmap(C.thetaQ,mmul(Rend,Pv)))
    rhs=mmul(msub(mmul(Hs,Rend),mmul(Rend,Hs)),Pv)
    res['z thQ(R Psi)=[H*,R]Psi']=miszero(mmap(lambda a:a.trunc(K),msub(lhs,rhs)))
    # (F8)
    res['F8 ratio']=(C.ew[0]/C.ew[1])==-(a1*a2)/((a0+a1)*(a0+a2))
    ctr={}
    def report(name,Mx):
        nz=mzero_report(Mx)
        ctr[name]='M=0' if not nz else 'M!=0 (%d nonzero coeffs; psi-orders %s; min Q-order %d)'%(len(nz),sorted(set(e[3] for e in nz)),min(e[2] for e in nz))
    Bstar=mmap(lambda a:a.zneg(),Rend)
    Mx,z=defect_endo(Bstar,Psi,eta,Hs,C.thetaQ); report('R^*(z)=R^-1(-z)',Mx)
    Bneg=mmap(lambda a:a.zneg(),B)
    Mx,z=defect_endo(Bneg,Psi,eta,Hs,C.thetaQ); report('R(-z)',Mx)
    Mx,z=defect_endo(B,Psi,eta,mscal(Hs,F(-1)),C.thetaQ); report('H*->-H*',Mx)
    Psi00=mconst(t,[[Psi[a][j].get(0,0) for j in range(2)] for a in range(2)])
    Mx,z=defect_endo(B,Psi00,eta,Hs,C.thetaQ); report('constant frame Psi(0) in M',Mx)
    # hnw dropped: constants, expected M=0
    hnwi=C.hv([a1,a0+a2,-C.e[0],-C.e[1]],-1)
    Rc=[[R[a][j]*hnwi for j in range(2)] for a in range(2)]
    Rcend=mmul(mmul(etai,Rc),mmul(Psii,eta))
    Mx,z=defect_endo(minv(Rcend),Psi,eta,Hs,C.thetaQ); report('no h_nw (expected M=0, ambiguity)',Mx)
    if verbose:
        print('---',label,'(a0,a1,a2,e1,e2)',tuple(map(str,w)),'D',D,'K',K,'field',F.__name__,'time %.1f'%(time.time()-t0),flush=True)
        print('  checks:',res)
        for k,v in ctr.items(): print('  ctl',k,':',v)
    return res,ctr

if __name__=='__main__':
    D=int(sys.argv[1]) if len(sys.argv)>1 else 4
    K=int(sys.argv[2]) if len(sys.argv)>2 else 5
    found=[]
    rng=[Fr(n) for n in range(-12,13) if n!=0]
    for a0,a1,a2 in product(rng,repeat=3):
        if a0+a1==0 or a0+a2==0 or a1>a2: continue
        r=-a1*a2/((a0+a1)*(a0+a2))
        if rat_sqrt(r) is None: continue
        # e(wt_+)=a0 a1 a2 e1 e2: take e1=1, e2=a0 a1 a2
        c=a0*a1*a2
        found.append((a0,a1,a2,Fr(1),c))
    print('rational samples with both e(wt) squares: at least',len(found),[tuple(map(str,f)) for f in found[:6]],flush=True)
    for w in found[:3]:
        run(Fr,w,D,K,'Q (square sample)')
    samples=[(1,2,3,5,7),(Fr(2),Fr(-3),Fr(5,2),Fr(1,3),Fr(4)),(Fr(3,7),Fr(11),Fr(-2,5),Fr(5,3),Fr(-2,7)),(1,1,1,1,1),(2,3,5,1,1)]
    plist=primes3mod4(2**89,400)
    for w in samples:
        cnt=0
        for p in plist:
            F=make_fp(p)
            C=Coni(F,*w,0,0)
            if any(fsqrt(x) is None for x in C.ew): continue
            cnt+=1
            if cnt==1:
                run(F,w,D,K,'F_p, p=%d'%p)
            else:
                out=run(F,w,D,K,'F_p',verbose=False)
                print('  second prime p=%d: M=0 is %s; controls %s'%(p,out[0]['M=0'],[v[:4] for v in out[1].values()]),flush=True)
            if cnt==2: break
