# Evidence only. Frame change on T_F: point-frame sum with rows h_w P^ versus Hodge-vertex sum with rows P^,
# exact in Q(omega), with F an arbitrary rational number (the indeterminate of the polynomial in F). No y-series.
import os
import sys, math, re, time, itertools
from fractions import Fraction as Fr
from collections import Counter
from hodge_lib import *

def rseq(n):
    r=[Fr(1)]
    prev=Fr(0)
    for m in range(0,n):
        nxt=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else Fr(0)))/(27*(m+1))
        r.append(nxt)
    return r

# Laurent polys in v: dict exp->Fr
def ladd(a,b,c=Fr(1)):
    r=dict(a)
    for k,x in b.items():
        r[k]=r.get(k,Fr(0))+c*x
    return r
def lmulv(a,k,c=Fr(1)):
    return {e+k:c*x for e,x in a.items()}
def ltheta(a,E):
    # theta v = 1/E - v/3
    r={}
    for e,x in a.items():
        if e==0: continue
        r[e-1]=r.get(e-1,Fr(0))+x*e/E
        r[e]=r.get(e,Fr(0))-x*e/3
    return r
def leval(a,v):
    tot=W(0,0)
    for e,x in a.items():
        tot=tot+(v**e)*x
    return tot

def theta_series(c,N):
    l=[Fr(0)]*(N+1)
    m=1
    while 2*m-1<=N:
        l[2*m-1]=BN[2*m]/(2*m*(2*m-1))*c**(2*m-1); m+=1
    r=[Fr(0)]*(N+1); r[0]=Fr(1)
    for n in range(1,N+1):
        r[n]=sum(k*l[k]*r[n-k] for k in range(1,n+1))/n
    return r
def zmul(a,b,N):
    r=[Fr(0)]*(N+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j<=N: r[i+j]+=x*y
    return r

class Rows:
    def __init__(s,e1,e2,N):
        s.e1=Fr(e1); s.e2=Fr(e2); s.E=s.e1+s.e2; s.N=N; E=s.E
        r=rseq(N+1)
        # 1/prod_{j<=L}(1+jz/E) series
        S=[{} for _ in range(N+1)]
        for l in range(0,N+1):
            ser=[Fr(0)]*(N+1); ser[0]=Fr(1)
            for j in range(1,l+2):
                # divide by (1+jz/E)
                c=Fr(j)/E
                for n in range(1,N+1):
                    ser[n]-=c*ser[n-1]
            for n in range(l,N+1):
                S[n][l]=S[n].get(l,Fr(0))+r[l]*ser[n-l]
        s.S=S
        a=[lmulv(x,1) for x in S]
        def Th(f):
            # Y = -E/3 + 1/v
            out=[]
            for n in range(N+1):
                g=ladd(lmulv(f[n],0,-E/3),lmulv(f[n],-1))
                if n>=1: g=ladd(g,ltheta(f[n-1],E))
                out.append(g)
            return out
        b=Th(a); c=Th(b)
        s.Q=[lmulv(x,-1) for x in b]
        s.TQ=[lmulv(x,-1) for x in c]
        # h_w(z)=Theta(-z/E)Theta(z/e1)Theta(z/e2)
        h=zmul(zmul(theta_series(-1/E,N),theta_series(1/s.e1,N),N),theta_series(1/s.e2,N),N)
        hi=zmul(zmul(theta_series(1/E,N),theta_series(-1/s.e1,N),N),theta_series(-1/s.e2,N),N)
        s.h=h; s.hi=hi

class Theory:
    """mode: 'pt' rows h^s P^ with point vertices; 'hodge' rows h^s P^ with Omega_w vertices."""
    def __init__(s,rows,t,F,mode,hs,wts=None,cols=(0,1,2)):
        s.rows=rows; E=rows.E; N=rows.N
        s.t=Fr(t); s.F=Fr(F); s.mode=mode; s.kap=rows.e1*rows.e2; s.E=E
        y=s.t**3; s.Dc=1+27*y; s.y=y
        s.vs=[];s.S=[];s.Qr=[];s.TQr=[];s.Del=[]
        hh={0:[Fr(1)]+[Fr(0)]*N,1:rows.h,-1:rows.hi}[hs]
        for i in range(3):
            ti=(OM**i)*s.t
            v=(ti*3+1)/(ti*E)
            s.vs.append(v)
            Sv=[leval(x,v) for x in rows.S]
            Qv=[leval(x,v) for x in rows.Q]
            TQv=[leval(x,v) for x in rows.TQ]
            def zm(ser):
                out=[]
                for n in range(N+1):
                    tot=W(0,0)
                    for k in range(n+1):
                        if hh[k]!=0: tot=tot+ser[n-k]*hh[k]
                    out.append(tot)
                return out
            s.S.append(zm(Sv)); s.Qr.append(zm(Qv)); s.TQr.append(zm(TQv))
            s.Del.append((v**-2)*(3*s.kap*E))
        s.eta00=s.kap*E**3/9
        s.wts=wts
        s.hodge=hodge_cache(wts) if mode=='hodge' else None
        s._e={}; s._v={}; s._n={}
    def n(s,i,j,M,b):
        key=(i,j,M,b)
        if key in s._n: return s._n[key]
        kap=s.kap; E=s.E; Dc=s.Dc
        Qi,Qj,Ti,Tj=s.Qr[i],s.Qr[j],s.TQr[i],s.TQr[j]
        x=s.S[i][M]*s.S[j][b]*s.eta00
        x=x-(Qi[M]*Qj[b])*(3*kap*E)
        x=x-(Qi[M]*Tj[b]+Ti[M]*Qj[b])*(3*kap*Dc)
        if s.F!=0:
            if M>=1: x=x+(Qi[M-1]*Qj[b])*(3*kap*Dc*s.F)
            if b>=1: x=x+(Qi[M]*Qj[b-1])*(3*kap*Dc*s.F)
        s._n[key]=x
        return x
    def edge(s,i,j,a,b):
        key=(i,j,a,b)
        if key in s._e: return s._e[key]
        tot=W(0,0)
        for k in range(b+1):
            tot=tot+s.n(i,j,a+k+1,b-k)*Fr(-1)*Fr((-1)**k)*Fr((-1)**(a+b+1))
        s._e[key]=tot
        return tot
    def dil(s,i,d):
        return s.S[i][d-1]*Fr((-1)**d)
    def vert(s,i,g,hs):
        hs=tuple(sorted(hs)); key=(i,g,hs)
        if key in s._v: return s._v[key]
        n=len(hs); dim=3*g-3+n-sum(hs)
        tot=W(0,0)
        if dim>=0:
            def parts(rem,maxp):
                # multisets of dilaton degrees d = p+1 >= 2, consuming p = d-1 of the budget
                yield ()
                for p in range(min(rem,maxp),0,-1):
                    for r in parts(rem-p,p): yield (p+1,)+r
            for ds in parts(dim,dim):
                if s.mode=='pt':
                    val=wk(g,tuple(sorted(hs+ds)))
                else:
                    val=s.hodge(g,tuple(sorted(hs+ds)))
                if val==0: continue
                mult=1
                for c in Counter(ds).values(): mult*=math.factorial(c)
                term=W(val/mult,0)
                for d in ds: term=term*s.dil(i,d)
                tot=tot+term
        tot=tot*(s.Del[i]**(g-1))
        s._v[key]=tot
        return tot
    def graph(s,gen,edges):
        nv=len(gen)
        val=[0]*nv
        for u,v in edges:
            val[u]+=1; val[v]+=1
        bud=[3*gen[i]-3+val[i] for i in range(nv)]
        total=W(0,0)
        for col in itertools.product(range(3),repeat=nv):
            hs=[[] for _ in range(nv)]
            def rec(k,pref):
                nonlocal total
                if k==len(edges):
                    x=pref
                    for i in range(nv):
                        x=x*s.vert(col[i],gen[i],tuple(hs[i]))
                        if x.iszero(): return
                    total=total+x
                    return
                u,v=edges[k]
                bu=bud[u]-sum(hs[u])
                for a in range(0,bu+1):
                    bv=(bud[v]-sum(hs[v])) if u!=v else (bu-a)
                    for b in range(0,bv+1):
                        hs[u].append(a); hs[v].append(b)
                        rec(k+1,pref*s.edge(col[u],col[v],a,b))
                        hs[u].pop(); hs[v].pop()
            rec(0,W(1,0))
        return total
    def total(s,g):
        tot=W(0,0); per=[]
        for gen,edges,aut in stable_graphs(g):
            x=s.graph(gen,edges)*Fr(1,aut)
            per.append(x); tot=tot+x
        return tot,per

def bs_F2ref(e1,e2,y,F):
    e1=Fr(e1); e2=Fr(e2)
    X=1/(1+27*Fr(y)); F=Fr(F); E=e1+e2; k=e1*e2
    f02=5*F**3/(8*X)+(5/(8*X)-Fr(1,2))*F**2+Fr(1,96)*(13*X+20/X-32)*F-(27*X**3-98*X**2+121*X-50)/(2160*X)
    f11=F**2/16+Fr(1,96)*(8-7*X)*F+(84*X**2-196*X+121)/5760
    f20=X/384*F-(21*X**2-29*X+8)/5760
    return E**4/k*f20+E**2*f11+k*f02

class Sym:
    pass
def load_mma(path_or_str,isfile=True):
    s=open(path_or_str).read() if isfile else path_or_str
    s=s.replace('\\[Epsilon][1]','e1').replace('\\[Epsilon][2]','e2').replace('^','**').replace('\n',' ')
    s=re.sub(r'(?<![\w.])(\d+)(?![\w.])',lambda m:'Fr(%s)'%m.group(1),s)
    return s
def F20F30():
    s=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "ref", "F20_F30_y.m")).read()
    s=s.strip()
    assert s[0]=='{' and s[-1]=='}'
    body=s[1:-1]
    # split at top level comma
    d=0
    for i,ch in enumerate(body):
        if ch in '([{': d+=1
        elif ch in ')]}': d-=1
        elif ch==',' and d==0:
            cut=i; break
    return load_mma(body[:cut],False),load_mma(body[cut+1:],False)

if __name__=="__main__":
    G=int(sys.argv[1]) if len(sys.argv)>1 else 2
    samples=[(Fr(3),Fr(7),Fr(1,5)),(Fr(1),Fr(1),Fr(-2,7)),(Fr(5,2),Fr(-3,7),Fr(3,11))]
    if len(sys.argv)>2: samples=samples[:int(sys.argv[2])]
    Fs=[Fr(0),Fr(1),Fr(-2,3),Fr(5,3),Fr(11,7),Fr(-7,2),Fr(2,9)]
    if G==3: Fs=Fs[:int(sys.argv[3])] if len(sys.argv)>3 else Fs[:2]
    f2s,f3s=F20F30()
    N=3*G-3+4 if G==2 else 3*G+3
    for (e1,e2,t) in samples:
        rows=Rows(e1,e2,N)
        wts=[e1+e2,-e1,-e2]
        for F in Fs:
            t0=time.time()
            pt=Theory(rows,t,F,'pt',1); P,pp=pt.total(G)
            hv=Theory(rows,t,F,'hodge',0,wts); H,hh=hv.total(G)
            bad=Theory(rows,t,F,'hodge',1,wts); B,_=bad.total(G)
            ctrl2=Theory(rows,t,F,'pt',-1); C2,_=ctrl2.total(G)
            y=t**3
            env=dict(e1=e1,e2=e2,y=y,XX=F,Fr=Fr)
            ref=eval(f2s if G==2 else f3s,{},env)
            extra=""
            if G==2:
                extra=" bsF2ref-pt=%s"%(P.a-bs_F2ref(e1,e2,y,F))
            print(f"g={G} eps=({e1},{e2}) t={t} F={F}: pt-hodge={(P-H).a},{(P-H).b} pt-stored={(P-ref).a},{P.b}{extra} | ctrl hodge(s=+1 rows)-pt nonzero={not (B-P).iszero()}, pt(h^-1)-pt nonzero={not (C2-P).iszero()} [{time.time()-t0:.1f}s]",flush=True)
            if F==Fs[0]:
                print("   value",P.a, " graphs differ between frames:",sum(1 for a,b in zip(pp,hh) if not (a-b).iszero()),"of",len(pp),flush=True)
