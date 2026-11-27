# Colour 0: H2, D0, Z-infinity, K=1, Lemma 4.1, (I1),(I2), E00.  Exact.  Evidence only.
import sys,time
from fractions import Fraction as Fr
from alg import *
from core import *

def gauss(A):
    n=len(A); A=[row[:] for row in A]
    for i in range(n):
        p=next(j for j in range(i,n) if A[j][i]!=0)
        A[i],A[p]=A[p],A[i]
        A[i]=[x/A[i][i] for x in A[i]]
        for j in range(n):
            if j!=i and A[j][i]!=0:
                fac=A[j][i]
                A[j]=[x-fac*y for x,y in zip(A[j],A[i])]
    return [A[i][n] for i in range(n)]

def run(e1,e2,K=6,BW=14,verbose=True):
    E=e1+e2; tau=-e1/E
    print('--- eps=(%s,%s) E=%s tau=%s K=%d'%(e1,e2,E,tau,K)); sys.stdout.flush()
    prec=K+6
    B=Base(E,prec)
    r=rec_r(prec+4)
    eta=eta_of(e1,e2)
    S,a,sig=rows0(B,r)
    res={}
    z=R([0,Fr(-1)])
    pf=stripped_PF(B,a[0],z)
    res['PF_stripped(prec %d)'%pf.prec]=all(iszero(pf.get(k)) for k in range(-3,pf.prec))
    # S_infty=-t_c
    Sinf=LS({},prec); vp=B.v0*E; fact=Fr(1)
    for l in range(prec):
        fact*=(l+1)
        Sinf=Sinf+(vp*(r[l]/fact))
        vp=vp*(B.v0*E)
    res['Zinf_period']=all((Sinf.get(k)==(-1 if k==1 else 0)) for k in range(0,prec))
    # (a)
    ok=True
    for a_ in range(3):
        for k in range(-4,-1):
            if not iszero(sig[a_].get(k)): ok=False
    res['H2a_vanish_below_-1']=ok
    res['H2a_sigma_minus1']= (iszero(sig[0].get(-1)) and R.lift(sig[1].get(-1))==E and R.lift(sig[2].get(-1))==-E*E)
    # (b)
    tops=[]
    okb=True
    for k in range(0,K+2):
        for comp,bound in ((0,k+1),(1,k+1),(2,k)):
            f=R.lift(sig[comp].get(k))
            q,o,lk=leaves(f,E)
            if not lk or q or (o and (min(o)<1 or max(o)>bound)):
                okb=False; print('H2b FAIL',comp,k,q,o,lk)
            tops.append((comp,k,max(o) if o else 0,bound))
    res['H2b']=okb
    if verbose: print('  top leaf < bound (comp,k,top,bound):',[t for t in tops if t[2]<t[3]] or 'all bounds attained')
    # (c)
    res['H2c']=(R.lift(sig[0].get(0))==R([E],[E,-1]) and R.lift(sig[1].get(0))==R([E*E/9],[E,-1]) and iszero(sig[2].get(0)))
    res['sigma0_1']=(R.lift(sig[0].get(1))==R([-E/18],[E,-1])+R([E/9],[E,-2]))
    # Z-infinity as sum beta_j/j
    okz=True
    for k in range(0,K+2):
        q,o,lk=leaves(R.lift(sig[0].get(k)),E)
        s=sum((c/j for j,c in o.items()),Fr(0))
        if s!=(E if k==0 else 0): okz=False
    res['Zinf_beta']=okz
    # K=1
    okK=True
    for n in range(-4,K+1):
        tot=R([0])
        for k1 in range(-2,n+3):
            k2=n-k1
            if k2<-2: continue
            for (p,qq),m in eta.items():
                x=R.lift(sig[p].get(k1)); y=R.lift(sig[qq].get(k2)).negarg()
                tot=tot+x*y*m
        want=3*e1*e2*E**3 if n==-2 else Fr(0)
        if not (tot==want): okK=False; print('K=1 FAIL',n,tot)
    res['K=1']=okK
    # Lemma 4.1 (direct comparison with the rows)
    v=B.v0; thv=B.theta(v)
    res['thetav']=all((thv-(LS({0:1/E})-(v*Fr(1,3)))).get(k)==0 for k in range(0,prec-2))
    def Sm(m):
        out=LS({},B.prec); vp=LS({0:Fr(1)},B.prec)
        for l in range(B.prec+1):
            den=R([1])
            for j in range(1,l+1-m+1):
                den=den*R([1,Fr(-j)/E])
            zp=R([1])
            for _ in range(l): zp=zp*z
            out=out+(vp*(zp/den*r[l])); vp=vp*v
        return out
    S1=Sm(1); S2=Sm(2)
    Phat1=((B.C1i*thv*B.vinv0)*E)*S1
    CT=B.C1i*thv
    logder=B.theta(CT)*CT.inv()
    Phat2_=((thv*E)*(B.vinv0*B.vinv0))*S2 + ((logder*z)-(B.d2*R([1])))*(B.vinv0*S1)
    Phat2=((B.c1*B.Delta*thv)*E)*Phat2_
    d1=Phat1-(a[1]*B.vinv0); d2_=Phat2-(a[2]*B.vinv0)
    res['Lemma4.1']=all(iszero(d1.get(k)) for k in range(-2,prec-4)) and all(iszero(d2_.get(k)) for k in range(-2,prec-5))
    # (I1) (I2)
    Yhat=B.C1i*B.Y0
    Yreg=Yhat-LS({-1:E})
    res['Yhat(0)=E/18']=(Yreg.get(0)==E/18 and Yhat.get(-1)==E)
    def Dc(F):
        return (Yhat*F)+(F*LS({-1:z}))+(F.deriv()*z)
    i1=sig[1]-Dc(sig[0])
    res['I1']=all(iszero(i1.get(k)) for k in range(-3,K+1))
    kappa=B.c1*B.c1*B.c1*B.Delta; mu=(B.c12-LS({0:Fr(1)}))*(E/2)
    res['kappa,mu']=(kappa.get(1)==-1 and kappa.get(0)==0 and mu.get(0)==0)
    i2=sig[2]-((kappa*Dc(sig[1]))+(mu*sig[1]))
    res['I2']=all(iszero(i2.get(k)) for k in range(-3,K+1))
    # Lemma 4.3 coefficient formula
    ok43=True
    DC=Dc(sig[0])
    for n in range(-2,K):
        val=(R([E,-(n+2)])*sig[0].get(n+1) if not iszero(sig[0].get(n+1)) else R([0]))
        for m in range(0,n+3):
            if not iszero(sig[0].get(n-m)): val=val+sig[0].get(n-m)*Yreg.get(m)
        if not (val==R.lift(DC.get(n))): ok43=False
    res['Lemma4.3']=ok43
    # Dil, L1
    okl=True
    for k in range(0,K+2):
        f=R([0,-1])*sig[0].get(k)
        if k==0:
            if not (f==R([E])-R([E*E],[E,-1])): okl=False
        else:
            if not in_Lam(f,E,k+1): okl=False; print('L1 FAIL',k,f)
    res['L1 (rho q_t)']=okl
    # Dil^(d;0)
    tay=R.lift(S.get(0)).taylor(8)
    res['D0(iv)']=all(-tay[d-1]==-E**(1-d) for d in range(1,9))
    # E00 from definition (1.1), leaves in psi at each [psi'^b]; rho^2 E_00 ; then fit in psi'
    t0=time.time()
    tayl={}
    for a_ in range(3):
        for k in range(-1,K+2):
            tayl[a_,k]=R.lift(sig[a_].get(k)).taylor(BW)
    def N(n,b):
        out=R([0])
        for k1 in range(-1,n+2):
            k2=n-k1
            if k2<-1 or k2>K+1: continue
            for (p,qq),m in eta.items():
                f=sig[p].get(k1)
                if iszero(f): continue
                out=out+R.lift(f)*(-m*tayl[qq,k2][b])
        return out
    okneg=True; okv=True; okiii=True; oksym=True; okformula=True
    table={}
    for n in range(-2,K+1):
        Ns=[N(n,b) for b in range(BW+1)]
        Htab=[]
        Hprev=R([0])
        for b in range(BW+1):
            val=Fr(0)
            for m_ in range(b+1):
                val+=(-1)**(b-m_)*Ns[m_].taylor(b-m_)[b-m_]
            G=Ns[b]-val
            H=(G-Hprev)/R([0,1])
            Htab.append(H); Hprev=H
        table[n]=Htab
        if n<0:
            if any(not H.is0() for H in Htab): okneg=False; print('E00 negative order FAIL',n)
            continue
        cmap={}
        for b,H in enumerate(Htab):
            q_,o,lk=leaves(H,E)
            if not lk or q_ or (o and (min(o)<1 or max(o)>n+1)):
                okv=False; print('E00(v) FAIL',n,b,q_,o,lk)
            for j,c in o.items(): cmap.setdefault(j,{})[b]=c
        for i,cs in cmap.items():
            lmax=n+2-i
            if lmax<1: okiii=False; print('E00(iii) support FAIL',n,i); continue
            A=[[Fr(l)**b/E**(b+1) for l in range(1,lmax+1)]+[cs.get(b,Fr(0))] for b in range(lmax)]
            sol=gauss(A)
            for b in range(BW+1):
                if sum(sol[l-1]*Fr(l)**b/E**(b+1) for l in range(1,lmax+1))!=cs.get(b,Fr(0)):
                    okiii=False; print('E00(iii) FAIL',n,i,b)
            table[n,'C',i]=sol
        nn=n+2
        for i in range(1,nn):
            for l in range(1,nn-i+1):
                C=Fr(0)
                for k1 in range(0,nn-1):
                    k2=nn-2-k1
                    for (p,qq),m in eta.items():
                        _,o1,_=leaves(R.lift(sig[p].get(k1)),E)
                        _,o2,_=leaves(R.lift(sig[qq].get(k2)),E)
                        C-=m*o1.get(i,Fr(0))*o2.get(l,Fr(0))
                val=Fr(i*l)/(E*(i+l))*C
                got=table[n,'C',i][l-1] if (n,'C',i) in table and l-1<len(table[n,'C',i]) else Fr(0)
                if val!=got: okformula=False; print('E00(i) mismatch',n,i,l,val,got)
                got2=table[n,'C',l][i-1] if (n,'C',l) in table and i-1<len(table[n,'C',l]) else Fr(0)
                if got!=got2: oksym=False
        if n==0:
            res['E00(iv)']=(Htab[0]==R([-e1*e2*E**4/27/E],[E,-1]))
    res['E00 no pole (n<0)']=okneg
    res['E00(v)']=okv; res['E00(iii)']=okiii; res['E00(i)']=okformula; res['E00 sym']=oksym
    print('  E00 time %.1fs'%(time.time()-t0))
    for k,v in res.items(): print('  %-24s %s'%(k,v))
    sys.stdout.flush()
    return B,r,sig,table

if __name__=='__main__':
    K=int(sys.argv[1]) if len(sys.argv)>1 else 5
    pts=[(Fr(-1,3),Fr(1)),(Fr(-2,5),Fr(1)),(Fr(-3),Fr(1)),(Fr(-7,3),Fr(1)),(Fr(-5,7),Fr(3,2)),(Fr(1,3),Fr(1))]
    for e1,e2 in pts:
        run(e1,e2,K=K)
