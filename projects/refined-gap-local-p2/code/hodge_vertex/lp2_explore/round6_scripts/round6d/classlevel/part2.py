# Part 2: class-level Omega=R.omega with psi insertions at d>0: GP localisation against the R.omega graph sum with legs.
import sys,os,time,json
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cl_lib import *
from rows_p2 import P2Data
from gsum import Rows,ymul,yadd,yscal,yinv
from gp import GP,compositions
import math

def to_rows(P,Pm,D,K):
    """Pm[a][j] Bi matrix"""
    t=P.tt
    PL=[[[[Pm[a][j].get(d,n) for d in range(D+1)] for n in range(K+1)] for j in range(3)] for a in range(3)]
    Dl=[[P.Delta[j].get(d,0) for d in range(D+1)] for j in range(3)]
    return Rows(3,PL,Dl,P.etai_f,D,K)

class Setup:
    def __init__(s,alpha,eps,D,K):
        s.alpha=alpha;s.eps=eps;s.D=D;s.K=K
        s.P=P2Data(Fr,alpha,eps,D,K)
        s.G=GP(alpha,eps)
        P=s.P;t=P.tt
        Pm=P.Prows()
        s.R=to_rows(P,Pm,D,K)
        # Q(y)
        L=[Fr(0)]*(D+1)
        for d in range(1,D+1): L[d]=Fr(P.Ay[d].numerator,P.Ay[d].denominator)/d
        ex=[Fr(0)]*(D+1);ex[0]=Fr(1)
        for n in range(1,D+1):
            ex[n]=sum(m*L[m]*ex[n-m] for m in range(1,n+1))/n
        s.Q=[Fr(0)]+ex[:D]
        s.Qp=[[Fr(1)]+[Fr(0)]*D]
        for d in range(1,D+1): s.Qp.append(ymul(s.Qp[-1],s.Q,D))
        s.Ainv=yinv(list(P.Ay),D)
        s.gcache={}
        s.Cm=[[[Fr(0)]*(D+1) for _ in range(3)] for _ in range(3)]
        s.Cm[0][1]=[Fr(1)]+[Fr(0)]*D
        s.Cm[1]=[list(x) for x in P.C11]
        s.Cm[2]=[list(x) for x in P.C12]
        s.Pm=Pm
    def rows_variant(s,name):
        P=s.P;t=P.tt;D=s.D;K=s.K
        if name=='R^*':
            Pm0=mmap(lambda a:a.zcoef(0),s.Pm)
            Pneg=mmap(lambda a:a.zneg(),s.Pm)
            return to_rows(P,mmul(mmul(Pm0,minv(Pneg)),Pm0),D,K)
        if name=='R(-z)':
            return to_rows(P,mmap(lambda a:a.zneg(),s.Pm),D,K)
        if name=='h^-1':
            return to_rows(P,P.Prows(hsign=-1),D,K)
        if name=='no h':
            return to_rows(P,P.Prows(hsign=0),D,K)
        raise ValueError
    def gw(s,g,ins,mode='anc'):
        """list N_d for d=0..D"""
        key=(g,tuple(ins),mode)
        if key in s.gcache: return s.gcache[key]
        out=[s.G.value(g,list(ins),d,mode) for d in range(s.D+1)]
        s.gcache[key]=out
        return out
    def gw_y(s,g,ins,mode='anc',wt=None):
        N=s.gw(g,ins,mode)
        r=[Fr(0)]*(s.D+1)
        for d,x in enumerate(N):
            f=x*(wt(d) if wt else 1)
            r=yadd(r,yscal(s.Qp[d],f))
        return r

def thetaQ_y(S,v):
    return ymul(S.Ainv,[d*x for d,x in enumerate(v)],S.D)

def main_compare(S,cases,log):
    """cases: list of (g,n). Compare all (k,a) insertions. Returns stats."""
    stats={}
    for g,n in cases:
        dim=3*g-3+n
        tests=0;bad=0;nonzero=0;t0=time.time();ks=0
        for tot in range(dim+1):
            for k in compositions(tot,n):
                for aa in __import__('itertools').product(range(3),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    gsum=S.R.integral(g,list(ins))
                    gw=S.gw_y(g,ins)
                    tests+=1
                    if any(x!=0 for x in gw): nonzero+=1
                    if max(k)>=2: ks+=1
                    if gsum!=gw:
                        bad+=1
                        print('MISMATCH',S.alpha,S.eps,g,ins,'gsum',gsum,'gw',gw,flush=True)
        stats[(g,n)]=dict(tests=tests,mismatch=bad,nonzero=nonzero,with_k_ge2=ks,secs=round(time.time()-t0,1))
        print('weights',S.alpha,S.eps,'(g,n)=',(g,n),'D=',S.D,stats[(g,n)],flush=True)
    return stats

def controls(S,cases):
    out={}
    for g,n in cases:
        dim=3*g-3+n
        for tot in sorted(set([0,1,2,dim])):
            if tot>dim: continue
            for k in compositions(tot,n):
                for aa in __import__('itertools').product(range(3),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    gw=S.gw_y(g,ins)
                    # legs +psi
                    v=S.R.integral(g,list(ins),legsign=1)
                    out.setdefault(('legs +psi',g,n),[0,0])
                    out[('legs +psi',g,n)][0]+=1
                    if v!=gw: out[('legs +psi',g,n)][1]+=1
                    for nm in ('R^*','R(-z)','h^-1','no h'):
                        if nm not in S.__dict__.setdefault('vars',{}):
                            S.vars[nm]=S.rows_variant(nm)
                        v=S.vars[nm].integral(g,list(ins))
                        out.setdefault((nm,g,n),[0,0,None])
                        out[(nm,g,n)][0]+=1
                        if v!=gw:
                            out[(nm,g,n)][1]+=1
                            fd=min(i for i in range(S.D+1) if v[i]!=gw[i])
                            out[(nm,g,n)][2]=fd if out[(nm,g,n)][2] is None else min(fd,out[(nm,g,n)][2])
    for k,v in out.items():
        print('control',k,'tests,fails'+(',first failing y-order' if len(v)>2 else ''),v,flush=True)
    return out

def divisor_checks(S,g):
    """ancestor divisor and unit equations:
    <tb_k(H^a) tb_0(H)>(Q) = theta_Q <tb_k(H^a)> + Sum_c C^c_{1a} <tb_{k-1}(H^c)>,   <tb_k(H^a) tb_0(1)> = <tb_{k-1}(H^a)>.
    Tested on the GP side and on the graph-sum side, separately."""
    dim=3*g-2
    res={'GP div':[0,0],'GS div':[0,0],'GP unit':[0,0],'GS unit':[0,0],'GP div without bubble (control)':[0,0]}
    z=[Fr(0)]*(S.D+1)
    for k in range(0,dim+2):
        for a in range(3):
            # graph sum side
            I2=S.R.integral(g,[(k,a),(0,1)])
            I1=S.R.integral(g,[(k,a)]) if k<=dim else z
            rhs=thetaQ_y(S,I1)
            if k>0:
                for c in range(3):
                    rhs=yadd(rhs,ymul(S.Cm[a][c],S.R.integral(g,[(k-1,c)]),S.D))
            res['GS div'][0]+=1
            if I2!=rhs: res['GS div'][1]+=1; print('GS div fail',g,k,a)
            # GP side
            G2=S.gw_y(g,((k,a),(0,1)))
            G1=S.gw_y(g,((k,a),),wt=lambda d:d) if k<=dim else z
            rhsg=G1
            if k>0:
                for c in range(3):
                    rhsg=yadd(rhsg,ymul(S.Cm[a][c],S.gw_y(g,((k-1,c),)),S.D))
            res['GP div'][0]+=1
            if G2!=rhsg: res['GP div'][1]+=1; print('GP div fail',g,k,a,G2,rhsg)
            res['GP div without bubble (control)'][0]+=1
            if G2!=G1: res['GP div without bubble (control)'][1]+=1
            # unit
            if k>0:
                U2=S.R.integral(g,[(k,a),(0,0)]);U1=S.R.integral(g,[(k-1,a)])
                res['GS unit'][0]+=1
                if U2!=U1: res['GS unit'][1]+=1
                V2=S.gw_y(g,((k,a),(0,0)));V1=S.gw_y(g,((k-1,a),))
                res['GP unit'][0]+=1
                if V2!=V1: res['GP unit'][1]+=1
    print('divisor/unit checks',S.alpha,S.eps,'g=',g,'D=',S.D,res,flush=True)
    return res

def descendant_test(S,g,n,kmax):
    """Descendants at all k (not limited by dim Mbar): GP against ancestors from the GRAPH SUM, through S_z from genus-0 GP.
    <tau_k(H^a)...> = Sum <tb_m(S_r H^a)...> ."""
    D=S.D
    eta=S.P.etai_f
    # S_r[a][c][d]: coefficient of Q^d H^c in S_r H^a
    Scache={}
    def Sr(r,a):
        if (r,a) in Scache: return Scache[(r,a)]
        out=[[Fr(0)]*(D+1) for _ in range(3)]
        if r==0:
            out[a][0]=Fr(1)
        else:
            for d in range(1,D+1):
                for b in range(3):
                    v=S.G.value(0,[(r-1,a),(0,b)],d,'desc')
                    for c in range(3):
                        out[c][d]+=v*eta[b][c]
        Scache[(r,a)]=out
        return out
    def toy(lst):
        r=[Fr(0)]*(D+1)
        for d,x in enumerate(lst): r=yadd(r,yscal(S.Qp[d],x))
        return r
    tests=0;bad=0;t0=time.time();nz=0
    dim=3*g-3+n
    if n==1:
        for k in range(0,kmax+1):
            for a in range(3):
                lhs=S.gw_y(g,((k,a),),mode='desc')
                rhs=[Fr(0)]*(D+1)
                for m in range(0,min(k,dim)+1):
                    r=k-m
                    SS=Sr(r,a)
                    for c in range(3):
                        rhs=yadd(rhs,ymul(toy(SS[c]),S.R.integral(g,[(m,c)]),D))
                tests+=1
                if any(lhs): nz+=1
                if lhs!=rhs:
                    bad+=1;print('desc mismatch',g,k,a,lhs,rhs,flush=True)
    else:
        for k in range(0,kmax+1):
          for l in range(0,kmax+1):
            for a in range(3):
              for b in range(3):
                if (k,a)>(l,b): continue
                lhs=S.gw_y(g,((k,a),(l,b)),mode='desc')
                rhs=[Fr(0)]*(D+1)
                for m in range(0,k+1):
                    for mm in range(0,l+1):
                        if m+mm>dim: continue
                        SS1=Sr(k-m,a);SS2=Sr(l-mm,b)
                        for c in range(3):
                            f1=toy(SS1[c])
                            if not any(f1): continue
                            for e in range(3):
                                f2=toy(SS2[e])
                                if not any(f2): continue
                                val=S.R.integral(g,[(m,c),(mm,e)]) if (m,c)<=(mm,e) else S.R.integral(g,[(mm,e),(m,c)])
                                rhs=yadd(rhs,ymul(ymul(f1,f2,D),val,D))
                tests+=1
                if any(lhs): nz+=1
                if lhs!=rhs:
                    bad+=1;print('desc mismatch',g,k,l,a,b,lhs,rhs,flush=True)
    print('descendant test weights',S.alpha,S.eps,'(g,n)',(g,n),'k<=',kmax,'D',D,'tests',tests,'nonzero',nz,'mismatch',bad,'secs',round(time.time()-t0,1),flush=True)

if __name__=='__main__':
    which=sys.argv[1]
    if which=='quick':
        S=Setup((1,2,5),(3,7),1,6)
        main_compare(S,[(1,1),(1,2),(2,1)],None)
    if which=='main':
        idx=int(sys.argv[2]);D=int(sys.argv[3]);cases=eval(sys.argv[4])
        W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))][idx]
        S=Setup(W[0],W[1],D,7)
        main_compare(S,cases,None)
    if which=='ctl':
        idx=int(sys.argv[2]);D=int(sys.argv[3]);cases=eval(sys.argv[4])
        W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))][idx]
        S=Setup(W[0],W[1],D,7)
        controls(S,cases)
    if which=='div':
        idx=int(sys.argv[2]);D=int(sys.argv[3])
        W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))][idx]
        S=Setup(W[0],W[1],D,7)
        for g in (1,2): divisor_checks(S,g)
    if which=='desc':
        idx=int(sys.argv[2]);D=int(sys.argv[3])
        W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))][idx]
        S=Setup(W[0],W[1],D,7)
        descendant_test(S,1,1,5)
        descendant_test(S,2,1,6)
        descendant_test(S,1,2,3)
        if len(sys.argv)>4: descendant_test(S,2,2,int(sys.argv[4]))
