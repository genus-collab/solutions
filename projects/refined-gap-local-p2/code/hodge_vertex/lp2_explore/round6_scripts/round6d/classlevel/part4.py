# Part 4: general Hodge classes through ch_{2k-1} (genus 3 and re-check of genus 1,2); Hodge-vertex form with legs; sensitivity tables; recorded-value comparisons.
import sys,os,time,itertools,math,pickle
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cl_lib import *
from gsum import Rows,ymul,yadd,yscal,yinv,graphs,partitions,stable_graphs,ypow
import gsum
from gp import GP,compositions,fixed_graphs
from part2 import Setup,to_rows
from intersect import anc,chint,wk,ch1int

W=[((1,2,5),(3,7)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7)))]

def hodge_ch_poly(weights,g,normalise=False):
    """prod_w Lambda^vee(w) (or /w^g) as a dict {rs:coef}, rs exponents of ch_1,ch_3,...; uses Lambda^vee(w)/w^g=exp(-Sum_k (2k-2)! ch_{2k-1} w^(1-2k))."""
    if g==0: return {():Fr(1)}
    top=1 if g==1 else 3*g-3
    nch=g
    # exponent e=Sum_k c_k ch_{2k-1}, c_k=-(2k-2)! p_{2k-1}(1/w)
    c=[-math.factorial(2*k-2)*sum(Fr(1)/w**(2*k-1) for w in weights) for k in range(1,nch+1)]
    pref=Fr(1)
    if not normalise:
        for w in weights: pref*=w**g
    out={}
    ranges=[range(0,top//(2*k+1)+1) for k in range(nch)]
    for rs in itertools.product(*ranges):
        deg=sum(r*(2*i+1) for i,r in enumerate(rs))
        if deg>top: continue
        co=pref
        for i,r in enumerate(rs):
            co*=c[i]**r/math.factorial(r)
        t=list(rs)
        while t and t[-1]==0: t.pop()
        out[tuple(t)]=out.get(tuple(t),Fr(0))+co
    return out

class GPch(GP):
    def vint(s,g,col,F_S,F_X,Mk,rx=0):
        assert rx==0
        key=(g,col,tuple(sorted(F_S)),tuple(sorted(F_X)),tuple(sorted(Mk)),rx)
        if key in s.vcache: return s.vcache[key]
        F_S=key[2];F_X=key[3];Mk=key[4]
        nS=len(F_S);nX=len(F_X)
        n=nS+nX+len(Mk)
        dim=3*g-3+n
        hod=s.__dict__.setdefault('hc',{})
        if (col,g) not in hod: hod[(col,g)]=hodge_ch_poly(s.wt(col),g)
        base=sum(b for t,b in F_S)+sum(a+b for a,b in Mk)
        total=Fr(0)
        for rs,co in hod[(col,g)].items():
            deg=sum(r*(2*i+1) for i,r in enumerate(rs))
            rem=dim-base-deg
            if rem<0: continue
            for exps in compositions(rem,nS+nX):
                ex=exps[:nS];ey=exps[nS:]
                w=Fr(1)
                for e,(t,b) in zip(ex,F_S): w/=t**(e+1)
                for e,t in zip(ey,F_X): w/=t**(e+1)
                SS=tuple(sorted([(e,b) for e,(t,b) in zip(ex,F_S)]+list(Mk)))
                if nX==0:
                    val=chint(g,tuple(sorted(p[0]+p[1] for p in SS)),rs)
                else:
                    val=anc(g,SS,tuple(sorted(ey)),rs)
                total+=co*w*val
        s.vcache[key]=total
        return total

class RowsHodge(Rows):
    """vertex classes Delta_i^{g-1} Hodge_{wt_i}: Hodge-vertex form of A3.6. Rows are R~ (no h)."""
    def set_weights(s,wts):
        s.wts=wts;s.hc={}
        return s
    def V(s,i,g,exps,rv=0):
        assert rv==0
        exps=tuple(sorted(exps))
        key=(i,g,exps,0)
        if key in s.vc: return s.vc[key]
        n=len(exps)
        dim=3*g-3+n
        if (i,g) not in s.hc: s.hc[(i,g)]=hodge_ch_poly(s.wts[i],g,normalise=True)
        r=[Fr(0)]*(s.D+1)
        for rs,co in s.hc[(i,g)].items():
            deg=sum(c*(2*j+1) for j,c in enumerate(rs))
            ex=dim-sum(exps)-deg
            if ex<0: continue
            for part in partitions(ex):
                ds=tuple(sorted(p+1 for p in part))
                w=chint(g,tuple(sorted(exps+ds)),rs)
                if w==0: continue
                mult={}
                for x in ds: mult[x]=mult.get(x,0)+1
                coef=Fr(1)
                for x,m in mult.items(): coef/=math.factorial(m)
                term=[Fr(1)]+[Fr(0)]*s.D
                for x in ds: term=ymul(term,s.dil(i,x),s.D)
                r=yadd(r,yscal(term,w*coef*co))
        r=ymul(r,ypow(s.Delta[i],g-1,s.D),s.D)
        s.vc[key]=r
        return r

def run_g3(idx,D,n):
    S=Setup(W[idx][0],W[idx][1],D,3*3-3+n+1)
    S.G=GPch(W[idx][0],W[idx][1])
    g=3
    dim=3*g-3+n
    t0=time.time()
    print('graphs (3,%d):'%n,len(graphs(3,n)),'secs',round(time.time()-t0,1),flush=True)
    tests=bad=nz=0
    for tot in range(dim+1):
        if n==0 and tot>0: break
        for k in (compositions(tot,n) if n>0 else [()]):
            for aa in itertools.product(range(3),repeat=n):
                ins=tuple((k[i],aa[i]) for i in range(n))
                if tuple(sorted(ins))!=ins: continue
                t1=time.time()
                gs=S.R.integral(g,list(ins))
                t2=time.time()
                gw=S.gw_y(g,ins)
                tests+=1
                if any(gw): nz+=1
                ok=(gs==gw)
                if not ok: bad+=1
                print('g=3',ins,'D',D,'ok' if ok else 'MISMATCH %s %s'%(gs,gw),'gs %.1fs gw %.1fs'%(t2-t1,time.time()-t2),flush=True)
    print('GENUS 3 weights',W[idx],'n',n,'D',D,'tests',tests,'nonzero',nz,'mismatch',bad,'secs',round(time.time()-t0,1),flush=True)

def recheck_g12(idx,D):
    """GP through ch-monomials against GP through lambda_1 polys (lambda_2=lambda_1^2/2); and against the graph sum"""
    S=Setup(W[idx][0],W[idx][1],D,7)
    G2=GPch(W[idx][0],W[idx][1])
    bad=tests=0
    for g,n in [(1,1),(1,2),(2,0),(2,1),(2,2)]:
        dim=3*g-3+n
        for tot in range(dim+1):
            if n==0 and tot>0: break
            for k in (compositions(tot,n) if n>0 else [()]):
                for aa in itertools.product(range(3),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    for d in range(D+1):
                        tests+=1
                        if G2.value(g,list(ins),d)!=S.G.value(g,list(ins),d): bad+=1
    print('GP(ch monomials) vs GP(lambda_1 poly): tests',tests,'mismatch',bad,flush=True)

def hodge_vertex(idx,D):
    """A3.6 with legs: R~-sum with vertex classes Delta^{g_v-1} Hodge_{wt_i} = GP, psi at d>0. Control: vertex point class with R~ rows."""
    S=Setup(W[idx][0],W[idx][1],D,7)
    P=S.P
    Pt=P.Prows(hsign=0)
    Rt=to_rows(P,Pt,D,7)
    RH=RowsHodge(3,Rt.P,Rt.Delta,Rt.etai,D,7).set_weights([P.wt(i) for i in range(3)])
    for g,n in [(1,1),(1,2),(2,0),(2,1),(2,2)]:
        dim=3*g-3+n
        tests=bad=nz=0;t0=time.time()
        for tot in range(dim+1):
            if n==0 and tot>0: break
            for k in (compositions(tot,n) if n>0 else [()]):
                for aa in itertools.product(range(3),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    a=RH.integral(g,list(ins))
                    gw=S.gw_y(g,ins)
                    tests+=1
                    if any(gw): nz+=1
                    if a!=gw:
                        bad+=1;print('Hodge-vertex MISMATCH',g,ins,flush=True)
        print('Hodge-vertex form with legs',W[idx],(g,n),'D',D,'tests',tests,'nonzero',nz,'mismatch',bad,'secs',round(time.time()-t0,1),flush=True)

def sensitivity(idx,D):
    """which total psi powers separate R from R^*, and descendants from ancestors."""
    S=Setup(W[idx][0],W[idx][1],D,7)
    Rs=S.rows_variant('R^*')
    for g,n in [(1,2),(2,1),(2,2)]:
        dim=3*g-3+n
        tab={}
        tabd={}
        for tot in range(dim+1):
            for k in compositions(tot,n):
                for aa in itertools.product(range(3),repeat=n):
                    ins=tuple((k[i],aa[i]) for i in range(n))
                    if tuple(sorted(ins))!=ins: continue
                    gw=S.gw_y(g,ins)
                    v=Rs.integral(g,list(ins))
                    tab.setdefault(tot,[0,0]);tab[tot][0]+=1
                    if v!=gw: tab[tot][1]+=1
                    # descendants in place of ancestors on the GP side, against the graph sum
                    if tot>0:
                        de=S.gw(g,ins,'desc');an=S.gw(g,ins,'anc')
                        tabd.setdefault(tot,[0,0,0]);tabd[tot][0]+=1
                        if de!=an:
                            tabd[tot][1]+=1
                            if de[0]!=an[0]: tabd[tot][2]+=1
        print('R^* control by total psi power k: {k:[tests,fails]}',(g,n),'dim',dim,tab,flush=True)
        print('descendant vs ancestor on the GP side, {k:[tests,differ,differ at d=0]}',(g,n),tabd,flush=True)

def recorded():
    """comparison with values recorded in the round-6c report (external, comparison only)."""
    S=Setup((1,2,5),(3,7),1,5)
    Rs=S.rows_variant('R^*')
    print('R^* g=2 d=1 graph sum:',Rs.integral(2,[]),'recorded 356413156283/928972800',flush=True)
    print('R g=2:',S.R.integral(2,[]),flush=True)
    G=GPch((1,2,5),(3,7))
    print('GP g=3 d=0,1:',G.value(3,[],0),G.value(3,[],1),'recorded R-sum 3032411/25401600, 392389/10584',flush=True)
    # singular sample of s.A.1
    Gs=GP((2,-3,7),(1,4))
    print('sample (2,-3,7;1,4), d=1:',Gs.value(2,[],1),flush=True)
    try:
        print('d=2:',Gs.value(2,[],2))
    except ZeroDivisionError:
        print('sample (2,-3,7;1,4), d=2: term-by-term localisation raises ZeroDivisionError, as s.A.1 warns',flush=True)

if __name__=='__main__':
    w=sys.argv[1]
    if w=='g3': run_g3(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
    if w=='recheck': recheck_g12(int(sys.argv[2]),int(sys.argv[3]))
    if w=='hv': hodge_vertex(int(sys.argv[2]),int(sys.argv[3]))
    if w=='sens': sensitivity(int(sys.argv[2]),int(sys.argv[3]))
    if w=='rec': recorded()
