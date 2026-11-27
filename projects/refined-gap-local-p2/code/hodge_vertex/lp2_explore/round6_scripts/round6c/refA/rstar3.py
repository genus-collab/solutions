# Discriminates R from R^* := R^{-1}(-z) at genus >= 2, exact Fractions, no GW input in test (i).
# (i)  divisor identity  theta_Q Int_{Mbar_g,n} Omega(H^n) = Int_{Mbar_g,n+1} Omega(H^{n+1}) for Omega = R.omega.
# (ii) comparison of the genus-2 sum with a hand-formula localisation at d = 0, 1.
import sys
from fractions import Fraction as Fr
from math import factorial, comb
from itertools import product, permutations, combinations
from functools import lru_cache
sys.setrecursionlimit(10000)

ALPHA=(Fr(1),Fr(2),Fr(5)); EPS=(Fr(3),Fr(7)); DY=1; ZM=8
# ---- y-series
def add(a,b): return [x+y for x,y in zip(a,b)]
def scal(c,a): return [c*x for x in a]
def mul(a,b):
    r=[Fr(0)]*(DY+1)
    for i,x in enumerate(a):
        if x:
            for j in range(DY+1-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def inv(a):
    r=[Fr(0)]*(DY+1); r[0]=1/a[0]
    for n in range(1,DY+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def sqrt1(a):
    r=[Fr(0)]*(DY+1); r[0]=Fr(1)
    for n in range(1,DY+1):
        r[n]=(a[n]-sum(r[k]*r[n-k] for k in range(1,n)))/2
    return r
def th(a): return [i*x for i,x in enumerate(a)]
one=[Fr(1)]+[Fr(0)]*DY; zero=[Fr(0)]*(DY+1); y1=[Fr(0),Fr(1)]+[Fr(0)]*(DY-1)
def const(c): return [Fr(c)]+[Fr(0)]*DY
def spow(a,n):
    if n<0: return spow(inv(a),-n)
    r=one[:]
    for _ in range(n): r=mul(r,a)
    return r
# ---- WK by DVV
def df(n):
    r=1
    while n>1: r*=n; n-=2
    return r
@lru_cache(maxsize=None)
def wk(g,a):
    a=tuple(sorted(a)); n=len(a)
    if g<0 or 2*g-2+n<=0 or sum(a)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and n==1: return Fr(1,24)
    k=a[-1]-1; rest=a[:-1]; res=Fr(0)
    for j,d in enumerate(rest):
        if d+k<0: continue
        b=list(rest); b[j]=d+k
        res+=Fr(df(2*k+2*d+1),df(2*d-1))*wk(g,tuple(sorted(b)))
    for r in range(k):
        s=k-1-r; c=Fr(df(2*r+1)*df(2*s+1),2)
        res+=c*wk(g-1,tuple(sorted(rest+(r,s))))
        idx=range(len(rest))
        for m in range(len(rest)+1):
            for I in combinations(idx,m):
                aI=tuple(rest[i] for i in I); aJ=tuple(rest[i] for i in idx if i not in I)
                for g1 in range(g+1):
                    res+=c*wk(g1,tuple(sorted(aI+(r,))))*wk(g-g1,tuple(sorted(aJ+(s,))))
    return res/df(2*k+3)
assert wk(2,(4,))==Fr(1,1152) and wk(2,(2,2,2))==Fr(7,240) and wk(3,(7,))==Fr(1,82944)
# ---- Bernoulli, Theta
def bern(n):
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1): B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/(m+1)
    return B[n]
def h_series(wts,Z):
    e=[Fr(0)]*(Z+1)
    for k in range(1,Z+2):
        if 2*k-1<=Z: e[2*k-1]=-bern(2*k)/(2*k*(2*k-1))*sum(1/w**(2*k-1) for w in wts)
    r=[Fr(0)]*(Z+1); r[0]=Fr(1)
    for n in range(1,Z+1): r[n]=sum(k*e[k]*r[n-k] for k in range(1,n+1))/n
    return r
# ---- rows
a0,a1,a2=ALPHA; e1=a0+a1+a2; e2=a0*a1+a0*a2+a1*a2; e3=a0*a1*a2; ee=EPS[0]*EPS[1]
A=[Fr(factorial(3*d),factorial(d)**3)*(-1)**d for d in range(DY+1)]; Ai=inv(A)
X=inv(add(one,scal(Fr(27),y1)))
Dfn=mul(X,mul(Ai,Ai)); e1E=scal(e1/2,add(Dfn,scal(Fr(-1),X)))
def Wser(Y): return add(add(scal(e1,mul(Y,Y)),scal(2*e2,Y)),const(3*e3))
def Lhat(Y):
    p=one
    for al in ALPHA: p=mul(p,add(Y,const(al)))
    return add(p,scal(Fr(27),mul(y1,mul(Y,mul(Y,Y)))))
def Lhatp(Y):
    s=zero[:]
    for l in range(3):
        p=one
        for m in range(3):
            if m!=l: p=mul(p,add(Y,const(ALPHA[m])))
        s=add(s,p)
    return add(s,scal(Fr(81),mul(y1,mul(Y,Y))))
def rootY(j):
    Y=const(-ALPHA[j])
    for _ in range(DY+2):
        Y=add(Y,scal(Fr(-1),mul(Lhat(Y),inv(Lhatp(Y)))))
    assert all(x==0 for x in Lhat(Y))
    return Y
# z-series of y-series: list n->series
def Dop(f,Y,c,Z):
    "(z theta + Y + c z) f, truncated at z^Z"
    out=[zero[:] for _ in range(Z+1)]
    for n in range(Z+1):
        t=mul(Y,f[n]) if n<len(f) else zero[:]
        if n>=1 and n-1<len(f): t=add(t,add(th(f[n-1]),scal(c,f[n-1])))
        out[n]=t
    return out
def Pop(f,Y,Z):
    p=[x[:] for x in f]
    for al in ALPHA: p=Dop(p,add(Y,const(al)),Fr(0),Z)
    q=[x[:] for x in f]
    for c in (Fr(2,3),Fr(1,3),Fr(0)): q=Dop(q,Y,c,Z)
    q=[scal(Fr(27),mul(y1,x)) for x in q]
    return [add(u,v) for u,v in zip(p,q)]
def row0(j,Z):
    Y=rootY(j); W=Wser(Y); w=scal(1/W[0],W); sw=sqrt1(w); swi=inv(sw)
    f=[swi]+[zero[:] for _ in range(Z)]
    YoW=mul(Y,inv(W))
    for n in range(1,Z+1):
        r=Pop(f[:n]+[zero[:]],Y,n+1)[n+1]
        rhs=mul(sw,mul(YoW,r))
        assert rhs[0]==0,("constant term",n,rhs)
        g=[Fr(0)]+[rhs[d]/d for d in range(1,DY+1)]
        f[n]=mul(swi,g)
    # check PF annihilation through z^(Z+1)
    res=Pop(f,Y,Z+1)
    assert all(all(x==0 for x in res[n]) for n in range(Z+2))
    wt=[3*ALPHA[j]]+[ALPHA[l]-ALPHA[j] for l in range(3) if l!=j]+[-EPS[0],-EPS[1]]
    h=h_series(wt,Z)
    P0=[zero[:] for _ in range(Z+1)]
    for n in range(Z+1):
        for m in range(n+1):
            P0[n]=add(P0[n],scal(h[m],mul(sw,f[n-m])))
    W5=Fr(1)
    for x in wt: W5*=x
    Delta=scal(W5,w)
    assert Delta==scal(3*ee,W)
    return Y,Delta,mul(th(w),inv(w)),P0,wt,W5
def rows(j,Z):
    Y,Delta,tl,P0,wt,W5=row0(j,Z)
    def T(f,extra):
        out=[zero[:] for _ in range(Z+1)]
        for n in range(Z+1):
            t=mul(add(Y,extra),f[n])
            if n>=1: t=add(t,add(th(f[n-1]),scal(Fr(-1,2),mul(tl,f[n-1]))))
            out[n]=t
        return out
    Q=T(P0,zero)
    P1=[mul(Ai,x) for x in Q]
    P2=[mul(inv(Dfn),x) for x in T(P1,scal(Fr(-1),e1E))]
    return dict(Delta=Delta,P=[P0,P1,P2],Y=Y,wt=wt,W5=W5)
ETA={ (0,0):3*e3*ee,(1,1):-3*e1*ee,(1,2):-3*ee,(2,1):-3*ee }
class Th: pass
def theory(Z,star=False):
    C=[rows(j,Z) for j in range(3)]
    Ds=[c['Delta'] for c in C]
    P=[c['P'] for c in C]
    if star:
        # R^*(z) = R^{-1}(-z) = R^T(z) in the f-frame:
        # P*_a^{(i)} = Sum_j eps_{aj} (Sum_{c,b} eps_{ci} eta^{cb} P_b^{(j)}) / Delta_j, eps_{aj} = [z^0] P_a^{(j)}
        eps_=[[P[j][a][0] for a in range(3)] for j in range(3)]
        Pn=[[[zero[:] for _ in range(Z+1)] for a in range(3)] for i in range(3)]
        for i in range(3):
            for a in range(3):
                for n in range(Z+1):
                    tot=zero[:]
                    for j in range(3):
                        s=zero[:]
                        for (c,b),et in ETA.items():
                            s=add(s,scal(et,mul(eps_[i][c],P[j][b][n])))
                        tot=add(tot,mul(mul(eps_[j][a],s),inv(Ds[j])))
                    Pn[i][a][n]=tot
        P=Pn
    t=Th(); t.P=P; t.D=Ds; t.Z=Z
    # symplectic defect
    for i in range(3):
        for j in range(3):
            for n in range(0,Z+1):
                s=zero[:]
                for k in range(n+1):
                    for (a,b),et in ETA.items():
                        s=add(s,scal(et*(-1)**(n-k),mul(P[i][a][k],P[j][b][n-k])))
                tgt=(Ds[i] if (i==j and n==0) else zero)
                assert s==tgt,("symp",i,j,n,s)
    # edges: E[(i,j)][(a,b)]
    N=Z
    t.E={}
    for i in range(3):
        for j in range(3):
            num={}
            for m in range(N+1):
                for b in range(N+1-m):
                    s=zero[:]
                    for (c,d),et in ETA.items():
                        s=add(s,scal(et*(-1)**(m+b),mul(P[i][c][m],P[j][d][b])))
                    s=scal(Fr(-1),s)
                    if m==0 and b==0 and i==j: s=add(s,Ds[i])
                    num[(m,b)]=s
            assert all(x==0 for x in num[(0,0)])
            e={}
            for tot in range(N):
                for a in range(tot,-1,-1):
                    b=tot-a
                    e[(a,b)]=add(num[(a+1,b)],scal(Fr(-1),e[(a+1,b-1)])) if b>0 else num[(a+1,0)]
            t.E[(i,j)]=e
    t.vcache={}
    return t
def partitions(n,mx=None):
    if mx is None: mx=n
    if n==0:
        yield ()
        return
    for p in range(min(n,mx),0,-1):
        for q in partitions(n-p,p): yield (p,)+q
def vval(t,g,pw,c):
    pw=tuple(sorted(pw)); key=(g,pw,c)
    if key in t.vcache: return t.vcache[key]
    n=len(pw); Dv=3*g-3+n-sum(pw)
    res=zero[:]
    if Dv>=0 and 2*g-2+n>0:
        for part in partitions(Dv):
            ds=tuple(p+1 for p in part)
            w=wk(g,tuple(sorted(pw+ds)))
            if w==0: continue
            coef=Fr(1)
            term=one[:]
            cnt={}
            for d in ds: cnt[d]=cnt.get(d,0)+1
            for d,m in cnt.items(): coef/=factorial(m)
            for d in ds:
                term=mul(term,scal(Fr((-1)**d),t.P[c][0][d-1]))
            res=add(res,scal(w*coef,term))
    res=mul(res,spow(t.D[c],g-1))
    t.vcache[key]=res
    return res
# ---- labelled-leg stable graphs
def graphs(g,n):
    seen={}; out=[]
    for V in range(1,max(1,2*g-2+n)+1):
        pairs=[(u,v) for u in range(V) for v in range(u+1,V)]
        def comps(k,tot):
            if k==0:
                if tot==0: yield ()
                return
            for i in range(tot+1):
                for c in comps(k-1,tot-i): yield (i,)+c
        for gv in product(range(g+1),repeat=V):
            need=g-sum(gv)+V-1
            if need<0: continue
            for c in comps(V+len(pairs),need):
                loops=c[:V]; adj=dict(zip(pairs,c[V:]))
                sn={0}; st=[0]
                while st:
                    x=st.pop()
                    for (u,v),m in adj.items():
                        if m and x in (u,v):
                            o=u+v-x
                            if o not in sn: sn.add(o); st.append(o)
                if len(sn)<V: continue
                for lg in product(range(V),repeat=n):
                    nl=[lg.count(v) for v in range(V)]
                    val=[2*loops[v]+nl[v]+sum(m for (a,b),m in adj.items() if v in (a,b)) for v in range(V)]
                    if any(2*gv[v]-2+val[v]<=0 for v in range(V)): continue
                    best=None; k0=None; keys=[]
                    for perm in permutations(range(V)):
                        pos={o:nw for nw,o in enumerate(perm)}
                        kk=(tuple((gv[perm[i]],loops[perm[i]]) for i in range(V)),
                            tuple(adj[(min(perm[i],perm[j]),max(perm[i],perm[j]))] for i in range(V) for j in range(i+1,V)),
                            tuple(pos[x] for x in lg))
                        if perm==tuple(range(V)): k0=kk
                        keys.append(kk)
                        if best is None or kk<best: best=kk
                    if best in seen: continue
                    seen[best]=1
                    aut=sum(1 for kk in keys if kk==k0)
                    for v in range(V): aut*=2**loops[v]*factorial(loops[v])
                    for p in pairs: aut*=factorial(adj[p])
                    edges=[]
                    it=iter(best[1])
                    for i in range(V):
                        for j in range(i+1,V): edges+=[(i,j)]*next(it)
                    for i in range(V): edges+=[(i,i)]*best[0][i][1]
                    out.append(dict(gv=[x[0] for x in best[0]],edges=edges,legs=list(best[2]),aut=aut))
    return out
def gsum(t,g,legrow):
    "Int_{Mbar_{g,n}} (R.omega)(H^{a_1},...,H^{a_n}), legrow = list of row indices a_l"
    n=len(legrow); tot=zero[:]
    for G in graphs(g,n):
        V=len(G['gv']); edges=G['edges']; legs=G['legs']
        sub=[zero[:]]
        dims=[]
        for v in range(V):
            nv=sum(1 for (a,b) in edges for x in (a,b) if x==v)+legs.count(v)
            dims.append(3*G['gv'][v]-3+nv)
        for col in product(range(3),repeat=V):
            def recl(il,pw,cur):
                if il==n:
                    val=cur
                    for v in range(V):
                        val=mul(val,vval(t,G['gv'][v],tuple(pw[v]),col[v]))
                        if not any(val): return
                    sub[0]=add(sub[0],val)
                    return
                v=legs[il]
                rem=dims[v]-sum(pw[v])
                for k in range(rem+1):
                    f=scal(Fr((-1)**k),t.P[col[v]][legrow[il]][k])
                    if not any(f): continue
                    pw2=[list(x) for x in pw]; pw2[v].append(k)
                    recl(il+1,pw2,mul(cur,f))
            def rec(ie,pw,cur):
                if ie==len(edges):
                    recl(0,pw,cur)
                    return
                u,v=edges[ie]
                for a in range(dims[u]-sum(pw[u])+1):
                    if u==v:
                        rb=dims[u]-sum(pw[u])-a
                    else:
                        rb=dims[v]-sum(pw[v])
                    for b in range(rb+1):
                        f=t.E[(col[u],col[v])][(a,b)]
                        if not any(f): continue
                        pw2=[list(x) for x in pw]; pw2[u].append(a); pw2[v].append(b)
                        rec(ie+1,pw2,mul(cur,f))
            rec(0,[[] for _ in range(V)],one[:])
        tot=add(tot,scal(Fr(1,G['aut']),sub[0]))
    return tot
# ---- localisation d = 0, 1 (hand formula, no R-matrix)
def wts_at(i): return [3*ALPHA[i]]+[ALPHA[l]-ALPHA[i] for l in range(3) if l!=i]+[-EPS[0],-EPS[1]]
def Vg(i,om,g):
    w=wts_at(i); pr=Fr(1)
    for x in w: pr*=x
    p1=sum(1/x for x in w); p2=sum(1/x**2 for x in w)
    e2_=sum(1/(w[a]*w[b]) for a in range(5) for b in range(a+1,5))
    e3_=sum(1/(w[a]*w[b]*w[c]) for a in range(5) for b in range(a+1,5) for c in range(b+1,5))
    S=sum(1/(w[a]*w[b]**2) for a in range(5) for b in range(5) if a!=b)
    if g==0: return Fr(1)
    if g==1:
        return pr*(Fr(1,24)/om**3-p1*Fr(1,24)/om**2)
    if g==2:
        val=Fr(1,1152)/om**4-p1*Fr(1,480)/om**3+(p2*Fr(7,5760)+e2_*Fr(7,2880))/om**2+(-e3_*Fr(1,1440)-S*Fr(1,2880))/om
        return pr**2*val/om**2
def loc_d1(g):
    tot=Fr(0)
    for i,j,k in ((0,1,2),(0,2,1),(1,2,0)):
        om_i=ALPHA[j]-ALPHA[i]; om_j=-om_i
        Ed=(3*ALPHA[i]+om_i)*(3*ALPHA[i]+2*om_i)/((ALPHA[k]-ALPHA[i])*(ALPHA[k]-ALPHA[j])*ee)
        for g1 in range(g+1):
            tot+=Ed*Vg(i,om_i,g1)*Vg(j,om_j,g-g1)
    return tot
def loc_d0_g2():
    tot=Fr(0)
    for i in range(3):
        w=wts_at(i); pr=Fr(1)
        for x in w: pr*=x
        e3_=sum(1/(w[a]*w[b]*w[c]) for a in range(5) for b in range(a+1,5) for c in range(b+1,5))
        S=sum(1/(w[a]*w[b]**2) for a in range(5) for b in range(5) if a!=b)
        tot+=pr*(-e3_*Fr(1,2880)-S*Fr(1,5760))
    return tot
if __name__=="__main__":
    print("loc: g=0 d=1 *eps1eps2 =",loc_d1(0)*ee," g=2 d=0",loc_d0_g2()," g=1 d=1",loc_d1(1)," g=2 d=1",loc_d1(2),flush=True)
    Z=ZM
    R=theory(Z,False); S=theory(Z,True)
    print("graph counts (2,0),(2,1),(1,1),(1,2),(0,4),(0,5):",[len(graphs(*x)) for x in ((2,0),(2,1),(1,1),(1,2),(0,4),(0,5))],flush=True)
    print("z^1 rows equal for R, R^*:",all(R.P[i][a][1]==S.P[i][a][1] for i in range(3) for a in range(3)))
    print("z^2 rows equal for R, R^*:",all(R.P[i][a][2]==S.P[i][a][2] for i in range(3) for a in range(3)))
    print("y=0 rows equal for R, R^*:",all(R.P[i][a][n][0]==S.P[i][a][n][0] for i in range(3) for a in range(3) for n in range(Z+1)),flush=True)
    for name,t in (("R",R),("R^*",S)):
        print("=====",name)
        C4=gsum(t,0,[1,1,1,1]); C5=gsum(t,0,[1,1,1,1,1])
        print(" g=0: A*Int(0,5)-theta Int(0,4):",[str(x) for x in add(mul(A,C5),scal(Fr(-1),th(C4)))],flush=True)
        F11=gsum(t,1,[1]); F12=gsum(t,1,[1,1])
        print(" g=1: A*Int(1,2)-theta Int(1,1):",[str(x) for x in add(mul(A,F12),scal(Fr(-1),th(F11)))],flush=True)
        F2=gsum(t,2,[]); F21=gsum(t,2,[1])
        print(" g=2: F_2 y-series:",[str(x) for x in F2])
        print(" g=2: A*Int(2,1)-theta F_2:",[str(x) for x in add(mul(A,F21),scal(Fr(-1),th(F2)))],flush=True)
    print("===== genus 3 divisor identity, y^0, y^1", flush=True)
    for name,t in (("R",R),("R^*",S)):
        F3=gsum(t,3,[]); F31=gsum(t,3,[1])
        print(name," F_3:",[str(x) for x in F3])
        print(name," defect:",[str(x) for x in add(mul(A,F31),scal(Fr(-1),th(F3)))],flush=True)
