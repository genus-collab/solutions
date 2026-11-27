# Colours 1,2: rows from Def 1.3 in Q(omega), Theorem M, (L2), orthogonality. Exact. Evidence only.
import sys,time
from fractions import Fraction as Fr
from alg import *
from core import *

def colour_rows(B,r,j,Zmax):
    E=B.E
    om_inv=OMEGA*OMEGA if j==1 else OMEGA   # omega^{-j}
    n=B.prec+2
    v=LS({k:(W(c)*om_inv*(-3/E)) for k,c in B.u.c.items()},B.u.prec)+LS({0:W(3/E)},B.u.prec)
    vinv=v.inv()
    Y=vinv-LS({0:W(E/3)},B.u.prec)
    # powers of v
    vp=[LS({0:W(1)},B.prec)]
    for l in range(Zmax+2): vp.append(vp[-1]*v)
    # coef(l,n)=[z^n] prod_{i<=l+1}(1+i z/E)^{-1}
    def coef(l):
        ser=[Fr(1)]+[Fr(0)]*Zmax
        for i in range(1,l+2):
            new=[Fr(0)]*(Zmax+1)
            for k in range(Zmax+1):
                new[k]=ser[k]-(Fr(i)/E*new[k-1] if k>=1 else 0)
            ser=new
        return ser
    cf=[coef(l) for l in range(Zmax+1)]
    a0=[]
    for m in range(Zmax+1):
        t=LS({},B.prec)
        for l in range(m+1):
            t=t+(vp[l+1]*(r[l]*cf[l][m-l]))
        a0.append(t)
    a1=[]
    for m in range(Zmax+1):
        t=Y*a0[m]
        if m>=1: t=t+B.theta(a0[m-1])
        a1.append(B.C1i*t)
    a2=[]
    for m in range(Zmax+1):
        t=(Y-B.d2)*a1[m]
        if m>=1: t=t+B.theta(a1[m-1])
        a2.append((B.Delta*B.c12)*t)
    P=[[vinv*x for x in row] for row in (a0,a1,a2)]
    return v,vinv,Y,P

def run(e1,e2,K=4,BW=6):
    E=e1+e2
    print('--- eps=(%s,%s) E=%s tau=%s K=%d BW=%d'%(e1,e2,E,-e1/E,K,BW)); sys.stdout.flush()
    prec=K+BW+8
    B=Base(E,prec)
    r=rec_r(prec+BW+6)
    eta=eta_of(e1,e2)
    S,a,sig=rows0(B,r)
    res={}
    Zmax=BW+4
    # leaf decomposition of sigma
    Sl={}
    for k in range(0,K+1):
        for l in range(1,k+2):
            Sl[k,l]=[Fr(0)]*3
        for p in range(3):
            _,o,ok=leaves(R.lift(sig[p].get(k)),E)
            for l,c in o.items(): Sl[k,l][p]=c
    tayl={(p,k):R.lift(sig[p].get(k)).taylor(BW+3) for p in range(3) for k in range(-1,K+2)}
    rows={}
    for j in (1,2):
        v,vinv,Y,P=colour_rows(B,r,j,Zmax)
        rows[j]=(v,vinv,Y,P)
    # (R): no negative powers, component 2 zero at t_c=0
    res['R_no_neg']=all(P_lo>=0 for j in (1,2) for row in rows[j][3] for x in row for P_lo in [x.lo()])
    res['v_j0 comp2=0']=all(rows[j][3][2][m].get(0)==0 for j in (1,2) for m in range(Zmax+1))
    # Delta_j(0)
    okd=True
    for j in (1,2):
        v=rows[j][0]
        Dj=(v*v).inv()*(3*e1*e2*E)
        om=OMEGA if j==1 else OMEGA*OMEGA
        want=(W(e1*e2*E**3/3)*om*om)/((W(1)-om)*(W(1)-om))
        if Dj.get(0)!=want: okd=False
    res['Delta_j(0)']=okd
    # orthogonality j,j' : <P^(j)(z),P^(j')(-z)>=delta Delta_j
    oko=True
    Z=BW+3
    for j in (1,2):
        for jj in (1,2):
            Pj=rows[j][3]; Pjj=rows[jj][3]
            v=rows[j][0]; Dj=(v*v).inv()*(3*e1*e2*E)
            for m in range(Z+1):
                tot=LS({},prec)
                for mm in range(m+1):
                    for (p,q),c in eta.items():
                        tot=tot+((Pj[p][mm]*Pjj[q][m-mm])*(c*(-1)**(m-mm)))
                want=Dj if (j==jj and m==0) else LS({},prec)
                d=tot-want
                for k in range(0,K+1):
                    if not iszero(d.get(k)): oko=False; print('orth fail',j,jj,m,k)
    res['orth jj']=oko
    # orthogonality 0,j : <sigma(psi),P^(j)(z=psi)>=0 (Taylor in psi)
    oko=True
    for j in (1,2):
        Pj=rows[j][3]
        for n in range(-1,K+1):
            for m in range(Z+1):
                tot=W(0)
                for k1 in range(-1,n+1):
                    k2=n-k1
                    for (p,q),c in eta.items():
                        for mm in range(m+1):
                            tot=tot+W(tayl[p,k1][mm]*c)*Pj[q][m-mm].get(k2)
                if tot!=0: oko=False; print('orth0j fail',j,n,m,tot)
    res['orth 0j']=oko
    # Ehat_0j from (1.1)
    okmc=True; okneg=True; okmb=True; okrank=True; okmd=True; okme=True
    H={}
    for j in (1,2):
        Pj=rows[j][3]
        def vjk(c,k,b):
            return W.lift(Pj[c][b].get(k))*((-1)**b)
        for k in range(-1,K+1):
            Ns=[]
            for b in range(BW+1):
                re=R([0]); im=R([0])
                for k1 in range(-1,k+1):
                    k2=k-k1
                    if k2<0: continue
                    for (p,q),c in eta.items():
                        f=sig[p].get(k1)
                        if iszero(f): continue
                        w=vjk(q,k2,b)*(-c)
                        re=re+R.lift(f)*w.a; im=im+R.lift(f)*w.b
                Ns.append((re,im))
            Hprev=(R([0]),R([0]))
            for b in range(BW+1):
                out=[]
                for part in (0,1):
                    val=Fr(0)
                    for m_ in range(b+1):
                        val+=(-1)**(b-m_)*Ns[m_][part].taylor(b-m_)[b-m_]
                    G=Ns[b][part]-val
                    out.append((G-Hprev[part])/R([0,1]))
                Hprev=tuple(out)
                H[j,k,b]=Hprev
                if k<0:
                    if not (out[0].is0() and out[1].is0()): okneg=False; print('E0j neg',j,k,b)
                    continue
                for part in (0,1):
                    if not in_Lam(out[part],E,k+1):
                        okmc=False; print('M-c FAIL',j,k,b,part,leaves(out[part],E))
                # M-b
                exp=[R([0]),R([0])]
                for k1 in range(0,k+1):
                    for l in range(1,k1+2):
                        tot=W(0)
                        for n_ in range(b+1):
                            pair=W(0)
                            for (p,q),c in eta.items():
                                pair=pair+W(Sl[k1,l][p]*c)*vjk(q,k-k1,b-n_)
                            tot=tot+pair*(Fr(-l)**n_/E**(n_+1))
                        tot=tot*(-l)
                        exp[0]=exp[0]+R([tot.a],[E,-l]); exp[1]=exp[1]+R([tot.b],[E,-l])
                if not (exp[0]==out[0] and exp[1]==out[1]):
                    okmb=False; print('M-b FAIL',j,k,b)
        # rank one at t_c=0
        for b in range(BW+1):
            # chi_j(w) coefficient
            chi=W(0)
            for n_ in range(b+1):
                pair=W(0)
                for (p,q),c in eta.items():
                    bb=(E,E*E/9,Fr(0))[p]
                    pair=pair+W(bb*c)*vjk(q,0,b-n_)
                chi=chi+pair*(Fr(-1)**n_/E**(n_+1))
            chi=chi*(-1)
            if not (H[j,0,b][0]==R([chi.a],[E,-1]) and H[j,0,b][1]==R([chi.b],[E,-1])):
                okrank=False
        # (M-d): rho E_j0(m,n) via code telescoping against rho E_0j(n,m)
        M_=3; N_=3
        def Pjpsi(c,k,p): return W.lift(Pj[c][p].get(k))*((-1)**p)
        def pre(p_,q_,k):
            tot=W(0)
            for k1 in range(0,k+2):
                k2=k-k1
                if k2<-1: continue
                for (p,q),c in eta.items():
                    tot=tot-Pjpsi(p,k1,p_)*W(c*tayl[q,k2][q_])
            return tot
        for k in range(0,K+1):
            for m_ in range(M_+1):
                for n_ in range(N_+1):
                    ed=W(0)
                    for kk in range(n_+1):
                        ed=ed+pre(m_+kk+1,n_-kk,k)*((-1)**kk)
                    hre=H[j,k,m_][0].taylor(n_)[n_]; him=H[j,k,m_][1].taylor(n_)[n_]
                    if ed!=W(hre,him): okmd=False; print('M-d FAIL',j,k,m_,n_,ed,hre,him)
    # M-e
    for j in (1,2):
        for jj in (1,2):
            Pj=rows[j][3]; Pjj=rows[jj][3]
            def pre(i1,i2,p_,q_):
                A=rows[i1][3]; Bb=rows[i2][3]
                tot=LS({},prec)
                for (p,q),c in eta.items():
                    tot=tot-((A[p][p_]*Bb[q][q_])*(c*(-1)**(p_+q_)))
                return tot
            for m_ in range(3):
                for n_ in range(3):
                    e1_=LS({},prec); e2_=LS({},prec)
                    for kk in range(n_+1):
                        e1_=e1_+(pre(j,jj,m_+kk+1,n_-kk)*((-1)**kk))
                    for kk in range(m_+1):
                        e2_=e2_+(pre(jj,j,n_+kk+1,m_-kk)*((-1)**kk))
                    d=e1_-e2_
                    if any(not iszero(d.get(k)) for k in range(0,K+1)) or e1_.lo()<0:
                        okme=False; print('M-e FAIL',j,jj,m_,n_)
    res['M-a no neg']=okneg; res['M-c']=okmc; res['M-b']=okmb; res['rank one chi_j']=okrank; res['M-d']=okmd; res['M-e sym']=okme
    for k,v in res.items(): print('  %-20s %s'%(k,v))
    sys.stdout.flush()

if __name__=='__main__':
    K=int(sys.argv[1]) if len(sys.argv)>1 else 4
    BW=int(sys.argv[2]) if len(sys.argv)>2 else 6
    pts=[(Fr(-1,3),Fr(1)),(Fr(-3),Fr(1)),(Fr(-5,7),Fr(3,2))]
    for e1,e2 in pts:
        t=time.time(); run(e1,e2,K=K,BW=BW); print('  time %.1f'%(time.time()-t))
