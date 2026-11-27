# Check s.7.5 (i),(ii),(iii): leaf decompositions, the identity
#   Delta_0^{h-1} d^nFpot_h(q_t)[f...] = (3e1e2E^3/t_c^2)^{h-1} Sum_r 1/r! dd^{n+r}(rho f.., delta^r)
# with everything expanded in Xi_h(R), against the direct finite decoration sum, for h=2,3, graphs (h)--(1) and (h)+loop.
# Evidence only.
import sys,time
from fractions import Fraction as Fr
from math import factorial,comb
from itertools import product
from collections import Counter
from rows import *
from graphs_check import Vtx,graph_value,recon

def solve_mod(A,b):
    # least squares free: A m x n, exact solve with extra equations check; returns x or None
    m=len(A); n=len(A[0]) if m else 0
    M=[row[:]+[b[i]] for i,row in enumerate(A)]
    piv=[];r=0
    for c in range(n):
        p=None
        for i in range(r,m):
            if M[i][c]%P: p=i;break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        iv=inv(M[r][c])
        M[r]=[x*iv%P for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%P:
                f=M[i][c]
                M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c);r+=1
    for i in range(r,m):
        if M[i][n]%P: return None
    x=[0]*n
    for i,c in enumerate(piv): x[c]=M[i][n]
    return x

class Leaf:
    def __init__(s,th):
        s.th=th; s.E=fp(th.E)
    def lam_coef(s,j,a):
        # [psi^a] 1/(E-j psi)= j^a/E^{a+1}
        return pow(j,a,P)*inv(pow(s.E,a+1,P))%P
    def decompose(s,series,nmaxleaf):
        """series: list of F_p psi-coefs (len>=nmaxleaf+3). Fit in span{1, lam_1..lam_nmaxleaf}. returns (const,{j:c}) or None"""
        L_=len(series)
        A=[]
        for a in range(L_):
            A.append([1 if a==0 else 0]+[s.lam_coef(j,a) for j in range(1,nmaxleaf+1)])
        x=solve_mod(A,series)
        if x is None: return None
        return x[0],{j:x[j] for j in range(1,nmaxleaf+1) if x[j]}

class Xi:
    def __init__(s,th,V):
        s.th=th;s.V=V;s.E=fp(th.E);s.memo={}
        s.Tc={d:(-inv(pow(s.E,d-1,P)))%P for d in range(2,60)}
    def xi(s,h,R):
        """Xi_h(R) = Sum_m 1/m! Int Omega prod lam_j(psi_i) prod T  -- F_p; direct definition s.1.5"""
        R=tuple(sorted(R)); key=(h,R)
        if key in s.memo: return s.memo[key]
        r=len(R); tot=0
        m=0
        V=s.V
        while True:
            dim=3*h-3+r+m
            if 2*m>dim: break
            # enumerate T degrees as non-increasing partitions (with multiplicity) and leaf exponents
            def parts(m_,bud,mx):
                if m_==0:
                    yield ()
                    return
                for d in range(min(mx,bud-2*(m_-1)),1,-1):
                    for q in parts(m_-1,bud-d,d):
                        yield (d,)+q
            for ds in parts(m,dim,dim):
                cnt=Counter(ds); mult=factorial(m)
                for v in cnt.values(): mult//=factorial(v)
                # ordered count = mult; weight 1/m!
                wT=inv(factorial(m))*mult%P
                for d in ds: wT=wT*s.Tc[d]%P
                bud=dim-sum(ds)
                def leafexp(i,b):
                    if i==r:
                        yield ()
                        return
                    for a in range(0,b+1):
                        for q in leafexp(i+1,b-a):
                            yield (a,)+q
                for ex in leafexp(0,bud):
                    w=wT
                    for j,a in zip(R,ex):
                        w=w*pow(j,a,P)%P*inv(pow(s.E,a+1,P))%P
                    I=V.omega_int(h,ex+ds)
                    if I!=0:
                        tot=(tot+w*fp(I))%P
            m+=1
        s.memo[key]=tot
        return tot

def main():
    e1,e2=Fr(sys.argv[1]),Fr(sys.argv[2])
    h=int(sys.argv[3]); kmax=int(sys.argv[4])   # kmax: top Taylor order k
    th=Theory(e1,e2,N=kmax+6,K=kmax+3*h+8)
    V=Vtx(th); X=Xi(th,V); LF=Leaf(th)
    E=fp(th.E)
    # rho = -E v_0/t_c as Laurent
    tc=L(1,[1]+[0]*(th.prec-1))
    rho=(th.v[0]*(-E))/tc
    assert rho.coef(0)==1
    # q^(k): [t_c^k](rho q_t) leaves; rho q_t=-psi rho S(v,-psi): psi^a coefficient: -(-1)^{a-1} rho [z^{a-1}]S
    Sz=th.rows[0][0]
    leaves_q={}
    for k in range(0,kmax+1):
        ser=[0]+[ (-(-1)**(a-1))*(rho*Sz[a-1]).coef(k)%P for a in range(1,kmax+5)]
        dec=LF.decompose(ser,k+1)
        if dec is None: print("q^(%d): NOT in C+Lambda_%d"%(k,k+1)); return
        leaves_q[k]=dec
    print("L1: constants of [t_c^k](rho q_t), k=0..:",[recon(leaves_q[k][0]) for k in leaves_q],"(k=0 expected E=%s)"%th.E)
    print("    k=0 leaf coefficient of lam_1:",recon(leaves_q[0][1].get(1,0)),"expected -E^2=",-th.E**2)
    # half-edge: colour 0 half-edge against regular colour j, far exponent b (rho Ehat_0j)[w^b], and against colour 0 far exponent
    def edge_leaves(j,b):
        out={}
        for k in range(0,kmax+1):
            ser=[(rho*th.edge(0,j,a,b)).coef(k) for a in range(0,kmax+5)]
            dec=LF.decompose(ser,k+1)
            if dec is None: return None
            out[k]=dec
        return out
    for j in (1,2,0):
        for b in (0,1,2):
            d=edge_leaves(j,b)
            print("rho Ehat_0%d[w^%d] in Lambda_{k+1}, k<=%d, const zero: %s"%(j,b,kmax,d is not None and all(d[k][0]==0 for k in d)))
    # Rescaled vertex factor with one half-edge f_b := Ehat_01(.,w)[w^b]
    K0=fp(3*th.e1*th.e2*th.E**3)
    def taylor_vertex(f_leaves,extra=None):
        """coefficients of d^1Fpot_h(rho q_t)[rho f] up to t_c^kmax, as a list by k, using Xi; f_leaves[k]=(const,{j:c})"""
        res=[0]*(kmax+1)
        # delta order pieces
        dl={k:leaves_q[k][1] for k in range(1,kmax+1)}
        # enumerate multisets of dilaton orders (k1<=k2<=...) total + k0 <= kmax
        def comps(total,mn):
            yield ()
            for k in range(mn,total+1):
                for c in comps(total-k,k):
                    yield (k,)+c
        for k0 in f_leaves:
            for j0,c0 in f_leaves[k0][1].items():
                for dls in comps(kmax-k0,1):
                    tot=k0+sum(dls)
                    # weight 1/r! with ordered sum: multiset multiplicity
                    r=len(dls)
                    cnt=Counter(dls); mult=factorial(r)
                    for v in cnt.values(): mult//=factorial(v)
                    w=inv(factorial(r))*mult%P*c0%P
                    # expand product of leaves of the delta's
                    def rec(i,cur,wt):
                        nonlocal res
                        if i==r:
                            res[tot]=(res[tot]+wt*X.xi(h,(j0,)+tuple(cur)))%P
                            return
                        for j,c in dl[dls[i]].items():
                            rec(i+1,cur+[j],wt*c%P)
                    rec(0,[],w)
        return res
    for j,b in [(1,0),(1,1),(2,0)]:
        fl=edge_leaves(j,b)
        tv=taylor_vertex(fl)
        # direct: Delta_0^{h-1} Sum_a f_a * d/dpsi... the (0.1) finite sum
        direct=None
        dim=3*h-3+1
        for a in range(0,dim+1):
            t=V.vertex(h,0,(a,))*th.edge(0,j,a,b)
            direct=t if direct is None else direct+t
        ok=True; report=[]
        pref=pow(K0,h-1,P)
        for k in range(0,kmax+1):
            lhs=direct.coef(2-2*h+k); rhs=pref*tv[k]%P
            report.append(lhs==rhs)
        print("(i)+(ii), h=%d, half-edge to colour %d, far exponent %d: direct = Taylor in Xi at t_c^{2-2h+k}, k=0..%d:"%(h,j,b,kmax),report, "direct polar part zero:",all(direct.coef(2-2*h+k)==0 for k in range(0,2*h-2)))
    # Xi vanish table
    def parts_D(D,mn=2):
        # multisets of parts>=2 with D(R)=D
        if D==0:
            yield ()
            return
        for p in range(mn,D+2):
            for q in parts_D(D-(p-1),p):
                yield (p,)+q
    for D in range(1,2*h):
        vals=[(R,X.xi(h,R)!=0) for R in parts_D(D)]
        print("Xi_%d(R) non-zero? D=%d:"%(h,D),vals)
    print("string check Xi(1,2)=2/E^2 Xi(2) with control R=(3):", X.xi(h,(1,3))==3*inv(E*E)%P*X.xi(h,(3,))%P)
    print("Z h0 times (-e1e2E^3)^(h-1):",recon(X.xi(h,())*pow(fp(-th.e1*th.e2*th.E**3),h-1,P)%P))
main()
