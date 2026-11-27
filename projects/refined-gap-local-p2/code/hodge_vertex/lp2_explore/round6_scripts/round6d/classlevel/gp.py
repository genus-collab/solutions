# Graber-Pandharipande localisation for the e_T(R pi_* f^* N)^{-1} twisted theory of P2, N=O(-3)_{3 alpha}+O_{-eps1}+O_{-eps2}.
# Fresh code. Correlators with insertions (k,a): psi^k ev^*H^a, psi = ancestor (pulled back from Mbar_{g,n}) or descendant.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from itertools import product,combinations,combinations_with_replacement,permutations
from functools import lru_cache
import math
from intersect import anc,ch1int,wk,g0mult

def compositions(n,k):
    if k==0:
        if n==0: yield ()
        return
    if k==1:
        yield (n,);return
    for i in range(n+1):
        for r in compositions(n-i,k-1):
            yield (i,)+r
def poscomps(n,k):
    for c in compositions(n-k,k):
        yield tuple(x+1 for x in c)

def connected(V,edges):
    seen={0};st=[0]
    while st:
        u=st.pop()
        for e in edges:
            for x,y in ((e[0],e[1]),(e[1],e[0])):
                if x==u and y not in seen:
                    seen.add(y);st.append(y)
    return len(seen)==V

@lru_cache(maxsize=None)
def fixed_graphs(g,n,d,ncol=3):
    """list of (structure, weight) with weight=1/|Aut| (no 1/d_e). structure=(gs,cols,edges,legs). Canonical under vertex permutations."""
    out={}
    for E in range(0,d+1):
        if (d==0)!=(E==0): continue
        for degs in poscomps(d,E) if E>0 else [()]:
            for V in range(1,E+2):
                b1=E-V+1
                if b1>g: continue
                pairs=list(combinations(range(V),2))
                for ch in product(pairs,repeat=E):
                    if V>1 and not connected(V,ch): continue
                    if V==1 and E>0: continue
                    for gs in compositions(g-b1,V):
                        for cols in product(range(ncol),repeat=V):
                            if any(cols[u]==cols[v] for u,v in ch): continue
                            for legs in product(range(V),repeat=n):
                                # canonical form
                                best=None
                                for p in permutations(range(V)):
                                    gs2=[0]*V;cols2=[0]*V
                                    for i in range(V): gs2[p[i]]=gs[i];cols2[p[i]]=cols[i]
                                    ee=tuple(sorted((min(p[u],p[v]),max(p[u],p[v]),dd) for (u,v),dd in zip(ch,degs)))
                                    ll=tuple(p[x] for x in legs)
                                    key=(tuple(gs2),tuple(cols2),ee,ll)
                                    if best is None or key<best: best=key
                                out[best]=out.get(best,0)+1
                if False: pass
    # weight: number of labelled data / (V! E!) ; edges labelled through 'ch' with degree composition
    res=[]
    for key,cnt in out.items():
        V=len(key[0]);E=len(key[2])
        res.append((key,Fr(cnt,math.factorial(V)*math.factorial(E))))
    return res

class GP:
    def __init__(s,alpha,eps):
        s.al=[Fr(x) for x in alpha];s.ep=[Fr(x) for x in eps]
        s.vcache={}
        s.ecache={}
        s.scache={}
    def wt(s,i):
        a=s.al
        j,k=[l for l in range(3) if l!=i]
        return [3*a[i],a[j]-a[i],a[k]-a[i],-s.ep[0],-s.ep[1]]
    def W(s,i):
        r=Fr(1)
        for w in s.wt(i): r*=w
        return r
    def tF(s,i,j,d): return (s.al[j]-s.al[i])/d
    def edge(s,i,j,d):
        key=(i,j,d)
        if key in s.ecache: return s.ecache[key]
        a=s.al;t=s.tF(i,j,d)
        r=Fr(1)
        for l in range(1,3*d):
            r*=3*a[i]+l*t
        for m in range(3):
            for l in range(0,d+1):
                w=(a[m]-a[i])-l*t
                if w!=0: r/=w
        r/=s.ep[0]*s.ep[1]
        s.ecache[key]=r
        return r
    def hodge(s,i,g):
        """coefficients in lambda_1^r of prod_{w in wt_i} Lambda^vee_g(w). g<=2."""
        p=[Fr(1)]
        if g==0: return p
        assert g<=2
        for w in s.wt(i):
            if g==1: f=[w,Fr(-1)]
            else: f=[w*w,-w,Fr(1,2)]     # lambda_2=lambda_1^2/2
            q=[Fr(0)]*(len(p)+len(f)-1)
            for x,u in enumerate(p):
                for y,v in enumerate(f):
                    q[x+y]+=u*v
            p=q
        # truncation: lambda_1^2=0 for g=1, lambda_1^4=0 for g=2
        lim=1 if g==1 else 3
        return p[:lim+1]
    def stabilise(s,struct):
        """Returns per mark m: (vertex v, port) of the surviving component it ends on, and per vertex the surviving ports.
        Ports of a stable vertex v: ('F',e_index) flags, ('M',m) marks."""
        if struct in s.scache: return s.scache[struct]
        gs,cols,edges,legs=struct
        V=len(gs);n=len(legs)
        val=[0]*V
        inc=[[] for _ in range(V)]
        for ei,(u,v,dd) in enumerate(edges):
            inc[u].append((ei,0));inc[v].append((ei,1));val[u]+=1;val[v]+=1
        marks=[[m for m in range(n) if legs[m]==v] for v in range(V)]
        stable=[gs[v]>0 or val[v]+len(marks[v])>=3 for v in range(V)]
        comps={}   # comp id -> dict(g, ports=set)
        link={}
        def lk(p,q):
            link[p]=q;link[q]=p
        for v in range(V):
            if stable[v]:
                ports=set()
                for (ei,end) in inc[v]:
                    ports.add(('F',v,ei,end))
                for m in marks[v]:
                    ports.add(('M',v,m))
                comps[('v',v)]={'g':gs[v],'ports':ports}
        for ei,(u,v,dd) in enumerate(edges):
            comps[('e',ei)]={'g':0,'ports':set()}
        for ei,(u,v,dd) in enumerate(edges):
            for end,x in ((0,u),(1,v)):
                port=('E',ei,end)
                if stable[x]:
                    comps[('e',ei)]['ports'].add(port)
                    lk(port,('F',x,ei,end))
                else:
                    if val[x]==2 and not marks[x]:
                        comps[('e',ei)]['ports'].add(port)
                    elif val[x]==1 and len(marks[x])==1:
                        comps[('e',ei)]['ports'].add(port)
                        lk(port,('mark',marks[x][0]))
                    elif val[x]==1 and not marks[x]:
                        pass
                    else:
                        raise RuntimeError
        for x in range(V):
            if not stable[x] and val[x]==2:
                (e1,en1),(e2,en2)=inc[x]
                lk(('E',e1,en1),('E',e2,en2))
        for v in range(V):
            if stable[v]:
                for m in marks[v]:
                    lk(('M',v,m),('mark',m))
        owner={}
        for cid,c in comps.items():
            for p in c['ports']: owner[p]=cid
        changed=True
        while changed:
            changed=False
            for cid in list(comps):
                c=comps[cid]
                if c['g']==0 and len(c['ports'])<=2:
                    ports=list(c['ports'])
                    if len(ports)==1:
                        q=link.pop(ports[0])
                        link.pop(q,None)
                        if q[0]!='mark':
                            comps[owner[q]]['ports'].discard(q)
                        else:
                            raise RuntimeError('mark on a one-point tail')
                    elif len(ports)==2:
                        q1=link.pop(ports[0]);q2=link.pop(ports[1])
                        if q1==ports[1]:
                            raise RuntimeError('closed loop, no stable component')
                        link[q1]=q2;link[q2]=q1
                    del comps[cid]
                    changed=True
                    break
        res={'markport':{},'surv':{}}
        for cid,c in comps.items():
            if cid[0]!='v':
                raise RuntimeError('edge component survived')
            res['surv'][cid[1]]=set(c['ports'])
        for p,q in list(link.items()):
            if p[0]=='mark':
                res['markport'][p[1]]=(owner[q][1],q)
        s.scache[struct]=res
        return res
    def vertex_int(s,g,col,flags,Sinfo):
        """Int over Mbar_{g,val+marks}: prod_F 1/(t_F-psi_F) Hodge(col,g) with extras.
        flags: list of (t_F, tag) ; Sinfo: dict tag->(b or plain exponent) ...
        Parameters: flags=[(t,b_or_None_flag)]: for each flag: t, b ancestor power, forgotten bool; marks: list of (a,b) fixed."""
        raise NotImplementedError
    def vint(s,g,col,F_S,F_X,Mk,rx=0):
        """F_S: tuple of (t,b): flags kept in the target, carrying ancestor power b.
        F_X: tuple of t: forgotten flags.
        Mk: tuple of (a,b): mark ports with plain exponent a, ancestor power b."""
        key=(g,col,tuple(sorted(F_S)),tuple(sorted(F_X)),tuple(sorted(Mk)),rx)
        if key in s.vcache: return s.vcache[key]
        F_S=key[2];F_X=key[3];Mk=key[4]
        nS=len(F_S);nX=len(F_X)
        n=nS+nX+len(Mk)
        dim=3*g-3+n
        hod=s.hodge(col,g)
        base=sum(b for t,b in F_S)+sum(a+b for a,b in Mk)
        total=Fr(0)
        for r,co in enumerate(hod):
            if co==0: continue
            rem=dim-base-r-rx
            if rem<0: continue
            for exps in compositions(rem,nS+nX):
                ex=exps[:nS];ey=exps[nS:]
                w=Fr(1)
                for e,(t,b) in zip(ex,F_S): w/=t**(e+1)
                for e,t in zip(ey,F_X): w/=t**(e+1)
                SS=tuple(sorted([(e,b) for e,(t,b) in zip(ex,F_S)]+list(Mk)))
                if nX==0:
                    # all ancestors degenerate to plain powers since nothing is forgotten
                    val=ch1int(g,tuple(sorted(p[0]+p[1] for p in SS)),r+rx)
                else:
                    val=anc(g,SS,tuple(sorted(ey)),r+rx)
                total+=co*w*val
        s.vcache[key]=total
        return total
    ncol=3
    def hval(s,c): return -s.al[c]
    def value(s,g,ins,d,mode='anc',lam=0):
        """ins: list of (k,a). Returns the correlator at degree d, with lambda_1^lam of the stabilised curve."""
        n=len(ins)
        total=Fr(0)
        for struct,weight in fixed_graphs(g,n,d,s.ncol):
            gs=struct[0]
            pos=[v for v in range(len(gs)) if gs[v]>0]
            if lam==0:
                total+=weight*s.contrib(struct,g,ins,mode,[0]*len(gs))
            else:
                if not pos: continue
                for dist in compositions(lam,len(pos)):
                    rv=[0]*len(gs);mult=math.factorial(lam)
                    for v,r in zip(pos,dist):
                        rv[v]=r;mult//=math.factorial(r)
                    total+=weight*mult*s.contrib(struct,g,ins,mode,rv)
        return total
    def contrib(s,struct,g,ins,mode,rv):
        gs,cols,edges,legs=struct
        V=len(gs);n=len(legs)
        val=Fr(1)
        for (u,v,dd) in edges:
            val*=s.edge(cols[u],cols[v],dd)/dd
            if val==0: return Fr(0)
        for m,(k,a) in enumerate(ins):
            val*=s.hval(cols[legs[m]])**a
        if val==0: return val
        inc=[[] for _ in range(V)]
        for ei,(u,v,dd) in enumerate(edges):
            inc[u].append((ei,0,s.tF(cols[u],cols[v],dd)))
            inc[v].append((ei,1,s.tF(cols[v],cols[u],dd)))
        need_anc=(mode=='anc') and any(k>0 for k,a in ins)
        if need_anc:
            st=s.stabilise(struct)
        for v in range(V):
            c=cols[v];g_v=gs[v];vl=len(inc[v])
            mk=[m for m in range(n) if legs[m]==v]
            W=s.W(c)
            if g_v==0 and vl+len(mk)<3:
                # unstable: contracted point
                if vl==1 and not mk:
                    val*=inc[v][0][2]
                elif vl==1 and len(mk)==1:
                    k,a=ins[mk[0]]
                    if mode=='desc' and k>0:
                        val*=(-inc[v][0][2])**k
                    # ancestor powers are carried by the surviving component
                elif vl==2 and not mk:
                    val*=W/(inc[v][0][2]+inc[v][1][2])
                else:
                    raise RuntimeError
                continue
            # stable vertex
            pref=W**(vl-1)
            if mode=='desc' or not need_anc:
                Mk=tuple((ins[m][0] if mode=='desc' else 0,0) for m in mk)
                F_S=tuple((t,0) for (ei,end,t) in inc[v])
                val*=pref*s.vint(g_v,c,F_S,(),Mk,rv[v])
            else:
                # ancestor powers ending on this vertex
                pw={}
                for m,(k,a) in enumerate(ins):
                    if k>0:
                        vw,port=st['markport'][m]
                        if vw==v: pw[port]=k
                if not pw:
                    F_S=tuple((t,0) for (ei,end,t) in inc[v])
                    Mk=tuple((0,0) for m in mk)
                    val*=pref*s.vint(g_v,c,F_S,(),Mk,rv[v])
                else:
                    surv=st['surv'][v]
                    F_S=[];F_X=[]
                    for (ei,end,t) in inc[v]:
                        port=('F',v,ei,end)
                        if port in surv:
                            F_S.append((t,pw.get(port,0)))
                        else:
                            F_X.append(t)
                    Mk=tuple((0,pw.get(('M',v,m),0)) for m in mk)
                    val*=pref*s.vint(g_v,c,tuple(F_S),tuple(F_X),Mk,rv[v])
            if val==0: return val
        return val
