# Extensions of part 2: lambda_1 insertions, divisor equation at n=2->3, A1.7 leg coefficient, rank-2 class-level test.
import sys,os,time,itertools
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cl_lib import *
from gsum import Rows,ymul,yadd,yscal,yinv,graphs
from gp import GP,compositions,fixed_graphs
from part2 import Setup
from part3 import Coni
import math

W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))]

def th(S,v):
    return ymul(S.Ainv,[d*x for d,x in enumerate(v)],S.D)

def lam_test(S,cases,lams=(1,2,3)):
    for g,n in cases:
        dim=3*g-3+n
        tests=bad=nz=0;t0=time.time()
        for lam in lams:
            if lam>dim: continue
            for tot in range(dim-lam+1):
                for k in compositions(tot,n):
                    for aa in itertools.product(range(3),repeat=n):
                        ins=tuple((k[i],aa[i]) for i in range(n))
                        if tuple(sorted(ins))!=ins: continue
                        gs=S.R.integral(g,list(ins),lam=lam)
                        N=[S.G.value(g,list(ins),d,'anc',lam=lam) for d in range(S.D+1)]
                        gw=[Fr(0)]*(S.D+1)
                        for d,x in enumerate(N): gw=yadd(gw,yscal(S.Qp[d],x))
                        tests+=1
                        if any(gw): nz+=1
                        if gs!=gw:
                            bad+=1;print('lam MISMATCH',g,ins,lam,gs,gw,flush=True)
        print('lambda_1 insertions',S.alpha,S.eps,(g,n),'D',S.D,'tests',tests,'nonzero',nz,'mismatch',bad,'secs',round(time.time()-t0,1),flush=True)

def div3(S,g):
    """<tb_k(a) tb_l(b) tb_0(H)> = theta_Q <tb_k(a) tb_l(b)> + Sum C^c_{1a}<tb_{k-1}(c) tb_l(b)> + Sum C^c_{1b}<tb_k(a) tb_{l-1}(c)>."""
    dim=3*g-1
    res={'GS':[0,0],'GP':[0,0]}
    for tot in range(0,dim+2):
        for k in range(0,tot+1):
            l=tot-k
            for a in range(3):
                for b in range(a,3):
                    def gs2(k1,a1,k2,a2):
                        return S.R.integral(g,[(k1,a1),(k2,a2)])
                    def gw2(k1,a1,k2,a2,w=None):
                        return S.gw_y(g,((k1,a1),(k2,a2)),wt=w)
                    L=S.R.integral(g,[(k,a),(l,b),(0,1)])
                    rhs=th(S,gs2(k,a,l,b))
                    L2=S.gw_y(g,((k,a),(l,b),(0,1)))
                    rhs2=gw2(k,a,l,b,lambda d:d)
                    for c in range(3):
                        if k>0:
                            rhs=yadd(rhs,ymul(S.Cm[a][c],gs2(k-1,c,l,b),S.D))
                            rhs2=yadd(rhs2,ymul(S.Cm[a][c],gw2(k-1,c,l,b),S.D))
                        if l>0:
                            rhs=yadd(rhs,ymul(S.Cm[b][c],gs2(k,a,l-1,c),S.D))
                            rhs2=yadd(rhs2,ymul(S.Cm[b][c],gw2(k,a,l-1,c),S.D))
                    res['GS'][0]+=1;res['GP'][0]+=1
                    if L!=rhs: res['GS'][1]+=1;print('div3 GS fail',g,k,l,a,b,flush=True)
                    if L2!=rhs2: res['GP'][1]+=1;print('div3 GP fail',g,k,l,a,b,flush=True)
    print('divisor n=2->3 (A1.5 at n=2 with psi on both legs)',S.alpha,S.eps,'g',g,'D',S.D,'[tests,fails]',res,flush=True)

def a17():
    """A1.7: with omega_{0,n}(f_i^n)=Delta_i^{(n-2)/2}. Rank 2, Delta_i=delta_i^2, R=exp(r z^m). Value of Int (R.omega)_{0,m+3}(f_i,f_j,...,f_j)
    should be -(R_m)_{ji} Delta_j^{(m+1)/2}."""
    print('A1.7 leg coefficient:')
    for m,r,dl in [(1,[[Fr(2),Fr(3)],[Fr(3),Fr(-1)]],(Fr(2),Fr(3))),(2,[[Fr(0),Fr(5,2)],[Fr(-5,2),Fr(0)]],(Fr(2),Fr(3))),(3,[[Fr(1),Fr(-7)],[Fr(-7),Fr(2)]],(Fr(3,2),Fr(1,2))),(4,[[Fr(0),Fr(1,3)],[Fr(-1,3),Fr(0)]],(Fr(2),Fr(5)))]:
        # R_{aj}: in the f frame P_a^{(j)}=Delta_j^{1/2} R_{aj}(z), flat basis := f basis, eta=1
        K=m+3
        R=[[[Fr(0)]*(K+1) for _ in range(2)] for _ in range(2)]
        for a in range(2):
            R[a][a][0]=Fr(1)
            for j in range(2):
                R[a][j][m]+=r[a][j]
                # second order term of the exponential, so that R is symplectic through 2m
                if 2*m<=K:
                    R[a][j][2*m]+=sum(r[a][c]*r[c][j] for c in range(2))/2
                if 3*m<=K:
                    R[a][j][3*m]+=sum(r[a][c]*r[c][e]*r[e][j] for c in range(2) for e in range(2))/6
                if 4*m<=K:
                    R[a][j][4*m]+=sum(r[a][c]*r[c][e]*r[e][f]*r[f][j] for c in range(2) for e in range(2) for f in range(2))/24
        P=[[[[R[a][j][n]*dl[j]] for n in range(K+1)] for j in range(2)] for a in range(2)]
        Ro=Rows(2,P,[[dl[0]**2],[dl[1]**2]],[[Fr(1),Fr(0)],[Fr(0),Fr(1)]],0,K)
        for i,j in ((0,1),(1,0)):
            ins=[(0,i)]+[(0,j)]*(m+2)
            v=Ro.integral(0,ins)[0]
            exp=-r[j][i]*dl[j]**(m+1)
            print('  m',m,'(i,j)',(i,j),'value',v,'predicted',exp,'OK' if v==exp else 'FAIL',flush=True)

class GPconi(GP):
    ncol=2
    def __init__(s,a0,a1,a2,e1,e2):
        s.a=[Fr(x) for x in (a0,a1,a2)];s.ep=[Fr(e1),Fr(e2)];s.al=None
        s.vcache={};s.ecache={};s.scache={}
    def wt(s,i):
        a0,a1,a2=s.a
        if i==0: return [a0,a1,a2,-s.ep[0],-s.ep[1]]
        return [-a0,a0+a1,a0+a2,-s.ep[0],-s.ep[1]]
    def hval(s,c): return Fr(0) if c==0 else -s.a[0]
    def tF(s,i,j,d): return (s.a[0] if i==0 else -s.a[0])/d
    def edge(s,i,j,d):
        key=(i,j,d)
        if key in s.ecache: return s.ecache[key]
        a0,a1,a2=s.a
        t=s.tF(i,j,d)
        wi=s.wt(i)
        r=Fr(1)
        # TP^1=O(2): H^0 weights w_i-l t, l=0..2d, drop the zero
        for l in range(0,2*d+1):
            w=wi[0]-l*t
            if w!=0: r/=w
        # O(-1)'s: H^1 weights w+l t, l=1..d-1
        for w in (wi[1],wi[2]):
            for l in range(1,d):
                r*=w+l*t
        r/=s.ep[0]*s.ep[1]
        s.ecache[key]=r
        return r

class ConiSetup:
    def __init__(s,w,D,K):
        s.w=w;s.D=D;s.K=K
        s.C=Coni(Fr,*w,D,K)
        C=s.C
        s.G=GPconi(*w)
        PL=[[[[C.v[j][a].get(d,n) for d in range(D+1)] for n in range(K+1)] for j in range(2)] for a in range(2)]
        Dl=[[C.Delta[j].get(d,0) for d in range(D+1)] for j in range(2)]
        s.R=Rows(2,PL,Dl,C.etai_f,D,K)
    def gw(s,g,ins,lam=0,mode='anc'):
        return [s.G.value(g,list(ins),d,mode,lam=lam) for d in range(s.D+1)]
    def variant(s,name):
        C=s.C;t=C.t;D=s.D;K=s.K
        Pm=[[C.v[j][a] for j in range(2)] for a in range(2)]
        if name=='R^*':
            P0=mmap(lambda a:a.zcoef(0),Pm)
            Pn=mmap(lambda a:a.zneg(),Pm)
            Pm=mmul(mmul(P0,minv(Pn)),P0)
        elif name=='R(-z)':
            Pm=mmap(lambda a:a.zneg(),Pm)
        PL=[[[[Pm[a][j].get(d,n) for d in range(D+1)] for n in range(K+1)] for j in range(2)] for a in range(2)]
        Dl=[[C.Delta[j].get(d,0) for d in range(D+1)] for j in range(2)]
        return Rows(2,PL,Dl,C.etai_f,D,K)

def coni_test(S,cases):
    # genus-0 check
    e1e2=S.G.ep[0]*S.G.ep[1]
    print('conifold <H,H>_{0,d} from GP (expected 1/(d e1 e2)):',[S.G.value(0,[(0,1),(0,1)],d,'desc')*d*e1e2 for d in range(1,S.D+1)],flush=True)
    for g,n in cases:
        dim=3*g-3+n
        tests=bad=nz=0;t0=time.time()
        for tot in range(dim+1):
            for k in (compositions(tot,n) if n>0 else [()]):
                for aa in itertools.product(range(2),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    if n==0 and tot>0: continue
                    gs=S.R.integral(g,list(ins))
                    gw=S.gw(g,ins)
                    tests+=1
                    if any(gw): nz+=1
                    if gs!=gw:
                        bad+=1;print('conifold MISMATCH',g,ins,gs,gw,flush=True)
        print('conifold class-level',S.w,(g,n),'D',S.D,'tests',tests,'nonzero',nz,'mismatch',bad,'secs',round(time.time()-t0,1),flush=True)
    # controls
    Rs={'R^*':S.variant('R^*'),'R(-z)':S.variant('R(-z)')}
    for g,n in cases:
        if n==0: continue
        dim=3*g-3+n
        cnt={}
        for tot in range(dim+1):
            for k in compositions(tot,n):
                for aa in itertools.product(range(2),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    gw=S.gw(g,ins)
                    v=S.R.integral(g,list(ins),legsign=1)
                    cnt.setdefault('legs +psi',[0,0]);cnt['legs +psi'][0]+=1
                    if v!=gw: cnt['legs +psi'][1]+=1
                    for nm,Rr in Rs.items():
                        v=Rr.integral(g,list(ins))
                        cnt.setdefault(nm,[0,0]);cnt[nm][0]+=1
                        if v!=gw: cnt[nm][1]+=1
        print('conifold controls',S.w,(g,n),'[tests,fails]',cnt,flush=True)

if __name__=='__main__':
    which=sys.argv[1]
    if which=='a17': a17()
    if which=='lam':
        idx=int(sys.argv[2]);D=int(sys.argv[3])
        S=Setup(W[idx][0],W[idx][1],D,7)
        lam_test(S,[(1,1),(1,2),(2,0),(2,1),(2,2)])
    if which=='div3':
        idx=int(sys.argv[2]);D=int(sys.argv[3])
        S=Setup(W[idx][0],W[idx][1],D,7)
        div3(S,1);div3(S,2)
    if which=='coni':
        D=int(sys.argv[2])
        for w in [(1,2,3,5,7),(Fr(2),Fr(-3),Fr(5,2),Fr(1,3),Fr(4))]:
            S=ConiSetup(w,D,7)
            coni_test(S,[(1,1),(1,2),(2,0),(2,1),(2,2)])
