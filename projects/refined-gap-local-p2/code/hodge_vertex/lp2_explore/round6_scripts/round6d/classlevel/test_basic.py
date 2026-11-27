# Unit checks of the graph-sum and GP engines against independent values.
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cl_lib import bern
from intersect import ch1int,wk
from gsum import *
from gp import *
from rows_p2 import P2Data
import math

# stable graph counts and Aut check: Sum over graphs of (-1)^E/|Aut| style is not used; print counts
print('graph counts (g,n):',{(g,n):len(stable_graphs(g,n)) for (g,n) in [(1,1),(1,2),(2,0),(2,1),(2,2),(3,0)]},flush=True)

# scalar check: h_v acting on 1 equals Hodge_v, with psi powers: Int psi^k (h_v.1)_{g,n} = Int psi^k prod Lambda^vee(v)/v^g
def scalar_rows(v,D,K):
    # h_v(z)=exp(-Sum_k B_2k/(2k(2k-1)) p_{2k-1}(1/v) z^{2k-1})
    ex=[Fr(0)]*(K+1)
    k=1
    while 2*k-1<=K:
        ex[2*k-1]=-bern(2*k)/(2*k*(2*k-1))*sum(Fr(1)/a**(2*k-1) for a in v)
        k+=1
    h=[Fr(0)]*(K+1);h[0]=Fr(1)
    # exp via recursion n h_n=Sum m ex_m h_{n-m}
    for n in range(1,K+1):
        h[n]=sum(m*ex[m]*h[n-m] for m in range(1,n+1))/n
    P=[[[[h[n]]+[Fr(0)]*D for n in range(K+1)]]]
    return Rows(1,P,[[Fr(1)]+[Fr(0)]*D],[[Fr(1)]],D,K)

def hodge_int(v,g,exps):
    # prod Lambda^vee(a)/a^g as poly in lambda_1, g<=2
    p=[Fr(1)]
    if g>0:
        for a in v:
            f=[Fr(1),Fr(-1)/a] if g==1 else [Fr(1),-Fr(1)/a,Fr(1,2)/a**2]
            q=[Fr(0)]*(len(p)+len(f)-1)
            for x,u in enumerate(p):
                for y,w in enumerate(f): q[x+y]+=u*w
            p=q
    tot=Fr(0)
    for r,c in enumerate(p):
        tot+=c*ch1int(g,tuple(sorted(exps)),r)
    return tot

v=[Fr(3),Fr(-5,2),Fr(7)]
R=scalar_rows(v,0,8)
bad=0;cnt=0
for g,n in [(1,1),(1,2),(2,0),(2,1),(2,2),(2,3)]:
    dim=3*g-3+n
    for k in compositions(0,0) if False else [c for tot in range(dim+1) for c in compositions(tot,n)] if n>0 else [()]:
        ins=[(x,0) for x in k]
        a=R.integral(g,ins)[0]
        b=hodge_int(v,g,list(k))
        cnt+=1
        if a!=b:
            bad+=1;print('scalar mismatch',g,n,k,a,b)
print('scalar h.1=Hodge with psi powers: tests',cnt,'mismatches',bad,flush=True)

# Control of the same: legs +psi
bad2=0;c2=0
for g,n in [(2,1),(2,2)]:
    dim=3*g-3+n
    for k in [c for tot in range(dim+1) for c in compositions(tot,n)]:
        ins=[(x,0) for x in k]
        a=R.integral(g,ins,legsign=1)[0]
        b=hodge_int(v,g,list(k))
        c2+=1
        if a!=b: bad2+=1
print('scalar control legs +psi: tests',c2,'mismatches',bad2,flush=True)

# GP: values recorded in round 6 reports (external, used only as comparison): g=2,n=0 at (1,2,5;3,7)
G=GP((1,2,5),(3,7))
t0=time.time()
print('GP g2 d0',G.value(2,[],0),'expected 19189/100800')
print('GP g2 d1',G.value(2,[],1),'expected 34231/1680','time',time.time()-t0,flush=True)
# genus 0 check: <H^2,H^2,H^2?> three-point d=1 of the twisted theory against the QDE data
P=P2Data(Fr,(1,2,5),(3,7),3,1)
# C^c_{1a}: eta(H*H^a,H^b)=<H,H^a,H^b>; at Q^d by divisor <H,H^a,H^b>_{0,d}=d<H^a,H^b>_{0,d}
Qy=[Fr(0)]*5
eta=P.eta_f
for a in range(3):
    for b in range(3):
        # GP two-point series
        vals=[G.value(0,[(0,a),(0,b)],d) if d>0 else None for d in range(0,4)]
        # H*H^a matrix in Q series: Sum_c C^c_{1a} eta_{cb}
        Cm=[None,P.C11,P.C12][a]
        if Cm is None:
            lhs=[eta[1][b]]+[Fr(0)]*3
        else:
            lhs=[sum(Cm[c][dd]*eta[c][b] for c in range(3)) for dd in range(4)]
        # convert to y: Q(y)
        Qs=Q_of_y=None
        A=P.Ay
        # log Q=log y+Sum A_d/d y^d
        import itertools
        L=[Fr(0)]*4
        for d in range(1,4): L[d]=Fr(A[d])/d
        ex=[Fr(0)]*4;ex[0]=Fr(1)
        for n in range(1,4):
            ex[n]=sum(m*L[m]*ex[n-m] for m in range(1,n+1))/n
        Q1=[Fr(0)]+ex[:3]
        rhs=[Fr(0)]*4
        cl=sum(eta[1][c]*0 for c in range(3))
        # classical part eta(H H^a,H^b)=eta[a+1][b] (with reduction): use interpolation over nodes
        nodes=[-x for x in P.al]
        rhs[0]=sum(nodes[i]**(1+a+b)/P.W0(i) for i in range(3))
        Qp=[Fr(1),Fr(0),Fr(0),Fr(0)]
        for d in range(1,4):
            Qp=ymul(Qp,Q1,3)
            for e in range(4): rhs[e]+=d*vals[d]*Qp[e]
        if a>0 and any(x!=y for x,y in zip(lhs,rhs)):
            print('genus-0 QDE vs GP mismatch',a,b,lhs,rhs)
        elif a>0:
            print('genus-0 two-point GP = QDE matrix, a,b=',a,b,'through y^3')
