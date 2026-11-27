# Graph sum of Theorem FB_A (with legs allowed) and independent virtual localisation.
from fractions import Fraction as Fr
from itertools import product
from math import factorial
from collections import Counter
from engine import wk, graphs, hodge_norm, stable
from rows import *

def partitions_min2(n, mn=1):
    """multisets of d_k-1>=1 summing to n"""
    if n==0:
        yield (); return
    for p in range(mn,n+1):
        for q in partitions_min2(n-p,p):
            yield (p,)+q

class GS:
    def __init__(s,c):
        s.c=c; s.D=c.D
        D=c.D
        s.Dil=[{} for _ in range(3)]
        for i in range(3):
            for d in range(2,c.Zmax+2):
                s.Dil[i][d]=[(1 if d%2==0 else -1)*c.P[i][0][y].get(d-1,Fr(0)) for y in range(D+1)]
        # Delta powers
        s.Dpow={}
        for i in range(3):
            for k in range(-1,6):
                s.Dpow[(i,k)]=ypow(c.Delta[i],k,D)
        # numerators N_ij(k,l) for psi^k w^l
        s.E={}
        Zm=c.Zmax
        for i in range(3):
            for j in range(3):
                N={}
                for a in range(3):
                    for b in range(3):
                        ee=c.etainv[a][b]
                        if ee==0: continue
                        for y1 in range(D+1):
                            for p,x in c.P[i][a][y1].items():
                                if x==0 or p>Zm: continue
                                for y2 in range(D+1-y1):
                                    for q,w in c.P[j][b][y2].items():
                                        if w==0 or q>Zm: continue
                                        sg=1 if (p+q)%2==0 else -1
                                        key=(p,q)
                                        if key not in N: N[key]=[Fr(0)]*(D+1)
                                        N[key][y1+y2]-=ee*x*w*sg
                if i==j:
                    key=(0,0)
                    if key not in N: N[key]=[Fr(0)]*(D+1)
                    N[key]=yadd(N[key],c.Delta[i])
                # (psi+w)E=N
                Eij={}
                z0=[Fr(0)]*(D+1)
                assert all(x==0 for x in N.get((0,0),z0)), "N(0,0)!=0"
                for tot in range(1,Zm+1):
                    prev=z0
                    for k in range(0,tot):
                        l=tot-1-k
                        # N_{k,l+1}=E_{k-1,l+1}+E_{k,l}
                        e=[x-w for x,w in zip(N.get((k,l+1),z0),prev)]
                        if any(x!=0 for x in e): Eij[(k,l)]=e
                        prev=e
                    # consistency: N_{tot,0}=E_{tot-1,0}
                    assert N.get((tot,0),z0)==prev, "division failure"
                s.E[(i,j)]=Eij
        s.Lcache={}
    def L(s,i,g,h):
        """Delta^{g-1} sum_m 1/m! sum_{d_k>=2} <psi^h psi^d> prod Dil_i(d_k)"""
        key=(i,g,h)
        if key in s.Lcache: return s.Lcache[key]
        D=s.D; n=len(h)
        res=[Fr(0)]*(D+1)
        if stable(g,n):
            rem=3*g-3+n-sum(h)
            if rem>=0:
                for part in partitions_min2(rem):
                    ds=tuple(p+1 for p in part)
                    v=wk(g,tuple(sorted(h+ds)))
                    if v==0: continue
                    t=[v]+[Fr(0)]*D
                    for dd,mu in Counter(ds).items():
                        for _ in range(mu): t=ymul(t,s.Dil[i][dd],D)
                        t=[x/factorial(mu) for x in t]
                    res=yadd(res,t)
        res=ymul(res,s.Dpow[(i,g-1)],D)
        s.Lcache[key]=res
        return res
    def legvec(s,i,a):
        """dict k-> y-series for P_a^{(i)}(-psi) coefficient of psi^k"""
        c=s.c; out={}
        for y in range(s.D+1):
            for p,x in c.P[i][a][y].items():
                if x==0: continue
                out.setdefault(p,[Fr(0)]*(s.D+1))[y]=x*(1 if p%2==0 else -1)
        return out
    def total(s,g,leg=None):
        """graph sum genus g; leg=None or a in {0,1,2}: one insertion H^a."""
        D=s.D
        n=0 if leg is None else 1
        res=[Fr(0)]*(D+1)
        byg=[]
        for G in graphs(g,n):
            gv=G['gv']; legs=G['legs']; edges=G['edges']; V=len(gv)
            val=[legs[v]+sum((u==v)+(w==v) for (u,w) in edges) for v in range(V)]
            dims=[3*gv[v]-3+val[v] for v in range(V)]
            gt=[Fr(0)]*(D+1)
            for col in product(range(3),repeat=V):
                # recursive over edges
                def rec(ei,hs,cur):
                    nonlocal gt
                    if ei==len(edges):
                        # legs
                        if n==1:
                            lv=legs.index(1)
                            lv_vec=s.legvec(col[lv],leg)
                            acc=[Fr(0)]*(D+1)
                            for k,ser in lv_vec.items():
                                if sum(hs[lv])+k>dims[lv]: continue
                                hs2=[list(x) for x in hs]; hs2[lv].append(k)
                                t=ymul(cur,ser,D)
                                for v in range(V):
                                    t=ymul(t,s.L(col[v],gv[v],tuple(sorted(hs2[v]))),D)
                                acc=yadd(acc,t)
                            gt=yadd(gt,acc)
                        else:
                            t=cur
                            for v in range(V):
                                t=ymul(t,s.L(col[v],gv[v],tuple(sorted(hs[v]))),D)
                            gt=yadd(gt,t)
                        return
                    u,w=edges[ei]
                    Eij=s.E[(col[u],col[w])]
                    for (a,b),ser in Eij.items():
                        if u==w:
                            if sum(hs[u])+a+b>dims[u]: continue
                        else:
                            if sum(hs[u])+a>dims[u] or sum(hs[w])+b>dims[w]: continue
                        hs[u].append(a); hs[w].append(b)
                        rec(ei+1,hs,ymul(cur,ser,D))
                        hs[w].pop(); hs[u].pop()
                rec(0,[[] for _ in range(V)],[Fr(1)]+[Fr(0)]*D)
            gt=[x/G['aut'] for x in gt]
            res=yadd(res,gt)
        return res

# ---------- localisation ------------
def edge_factor(alpha,eps,i,j,d):
    """1/e(H^0-H^1)^mov of f^*TZ on an edge of degree d between P_i,P_j."""
    ai,aj=alpha[i],alpha[j]
    num=Fr(1)
    for l in range(1,3*d):
        num*=3*ai+Fr(l,d)*(aj-ai)
    den=eps[0]*eps[1]*Fr(1)  # (-e1)(-e2)
    for l in range(3):
        al=alpha[l]
        for m in range(0,d+1):
            w=(al-ai)+Fr(m,d)*(ai-aj)
            if w!=0: den*=w
    return num/den

def loc(alpha,eps,g,d,ins=()):
    """GW^T_{g,d}(H^{a}...) via GP. ins = tuple of exponents a (n markings). labelled-vertex sum."""
    wts=[]
    for i in range(3):
        wts.append([3*alpha[i]]+[alpha[l]-alpha[i] for l in range(3) if l!=i]+[-eps[0],-eps[1]])
    W=[]
    for i in range(3):
        p=Fr(1)
        for w in wts[i]: p*=w
        W.append(p)
    n=len(ins)
    total=Fr(0)
    Hv=[-a for a in alpha]
    if d==0:
        for i in range(3):
            v=Fr(1)
            for a in ins: v*=Hv[i]**a
            total+=v*W[i]**(g-1)*hodge_norm(g,(0,)*n,wts[i])/1
        return total
    # V vertices
    for V in range(2,d+2):
        pairs=[(u,v) for u in range(V) for v in range(u+1,V)]
        # edges: multiset of (pair,deg)
        kinds=[(p,dd) for p in pairs for dd in range(1,d+1)]
        def multisets(idx,rem):
            if rem==0:
                yield {}; return
            if idx==len(kinds): return
            p,dd=kinds[idx]
            for m in range(0,rem//dd+1):
                for ms in multisets(idx+1,rem-m*dd):
                    if m>0:
                        x=dict(ms); x[kinds[idx]]=m; yield x
                    else:
                        yield ms
        for ms in multisets(0,d):
            E=sum(ms.values())
            b1=E-V+1
            if b1<0 or b1>g: continue
            # connected
            adj={v:set() for v in range(V)}
            for ((u,v),dd),m in ms.items():
                adj[u].add(v); adj[v].add(u)
            seen={0}; st=[0]
            while st:
                x=st.pop()
                for w in adj[x]:
                    if w not in seen: seen.add(w); st.append(w)
            if len(seen)<V: continue
            sym=Fr(1,factorial(V))
            for m in ms.values(): sym/=factorial(m)
            for ((u,v),dd),m in ms.items():
                sym/=dd**m
            for col in product(range(3),repeat=V):
                if any(col[u]==col[v] for ((u,v),dd) in ms): continue
                edge_w=Fr(1)
                flags={v:[] for v in range(V)}
                for ((u,v),dd),m in ms.items():
                    ee=edge_factor(alpha,eps,col[u],col[v],dd)
                    edge_w*=ee**m
                    wu=(alpha[col[v]]-alpha[col[u]])/dd
                    wv=(alpha[col[u]]-alpha[col[v]])/dd
                    flags[u]+= [wu]*m
                    flags[v]+= [wv]*m
                # genus and marking distributions
                for gs in product(range(g-b1+1),repeat=V):
                    if sum(gs)!=g-b1: continue
                    for marks in product(range(V),repeat=n):
                        t=edge_w*sym
                        for k,mk in enumerate(marks):
                            t*=Hv[col[mk]]**ins[k]
                        for v in range(V):
                            nm=marks.count(v)
                            t*=vertex(wts[col[v]],W[col[v]],gs[v],flags[v],nm)
                            if t==0: break
                        total+=t
    return total

def vertex(wt,Wv,gv,fl,nm):
    val=len(fl)+nm
    if gv==0 and val==1:
        return fl[0]
    if gv==0 and val==2:
        if nm==0:
            return Wv/(fl[0]+fl[1])
        return Fr(1)  # one flag + one marking
    # stable
    N=len(fl)
    pref=Wv**(gv+val-1-nm)
    # W^{flags -1} * prod Lambda^vee(w)=W^{-1}W^{g}Hodge; flags count only
    dim=3*gv-3+val
    tot=Fr(0)
    def rec(k,exps,c):
        nonlocal tot
        if k==N:
            tot+=c*hodge_norm(gv,tuple(sorted(exps+[0]*nm)),tuple(wt))
            return
        for e in range(0,dim-sum(exps)+1):
            rec(k+1,exps+[e],c/fl[k]**(e+1))
    rec(0,[],Fr(1))
    return pref*tot

def q_to_series(c,gw):
    """sum_d Q^d gw[d] as y series"""
    D=c.D
    res=[Fr(0)]*(D+1)
    Qp=[Fr(1)]+[Fr(0)]*D  # Q^0
    Qy=[Fr(0)]+c.Qhat[:-1]
    for d in range(D+1):
        res=yadd(res,yscal(gw[d],Qp))
        Qp=ymul(Qp,Qy,D)
    return res
