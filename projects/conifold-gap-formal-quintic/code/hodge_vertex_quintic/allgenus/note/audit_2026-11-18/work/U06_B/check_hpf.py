# Sanity checks for unit U06 (lens B): HPF recursion, T-system, WKB operator.
from fractions import Fraction as Fr
import sys

# ---------- polynomials in u: dict {k: coeff}
def padd(a,b):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,Fr(0))+v
        if r[k]==0: del r[k]
    return r
def pscal(c,a):
    c=Fr(c)
    return {k:c*v for k,v in a.items() if c*v!=0}
def pmul(a,b):
    r={}
    for k1,v1 in a.items():
        for k2,v2 in b.items():
            r[k1+k2]=r.get(k1+k2,Fr(0))+v1*v2
    return {k:v for k,v in r.items() if v!=0}
def ptheta(a):
    return {k:k*v for k,v in a.items() if k!=0}
ONE={0:Fr(1)}
ONEMU={0:Fr(1),1:Fr(-1)}   # 1-u
def umon(j): return {j:Fr(1)}

# series in zh: list of polys, truncated at N
def O(e,F,N):
    # O_e F = F + zh (1-u) (4/5 theta_zh - e/5 - theta_u) F
    G=[dict(F[m]) if m<len(F) else {} for m in range(N+1)]
    for m in range(1,N+1):
        if m-1<len(F):
            f=F[m-1]
            t=padd(pscal(Fr(4*(m-1)-e,5),f), pscal(-1,ptheta(f)))
            G[m]=padd(G[m],pmul(ONEMU,t))
    return G
def shift(F,j,N):  # multiply by zh^j
    return [ {} ]*j + F[:N+1-j]
def smulpoly(p,F):
    return [pmul(p,f) for f in F]
def sadd(F,G):
    n=max(len(F),len(G))
    return [padd(F[i] if i<len(F) else {}, G[i] if i<len(G) else {}) for i in range(n)]

def HPF(Phi,N):
    P=[Phi]
    for j in range(1,6):
        P.append(O(-j,P[-1],N))
    res=sadd(P[5],smulpoly({0:Fr(-1)},P[0]))
    sig=[Fr(3),Fr(17,5),Fr(9,5),Fr(274,625),Fr(24,625)]
    for k in range(1,6):
        term=shift(P[5-k],k,N)
        term=smulpoly(pmul(pscal(-sig[k-1],ONEMU),umon(k-1)),term)
        res=sadd(res,term)
    return res

def pdiv_1mu(a):
    # divide polynomial a by (1-u) exactly
    if not a: return {}
    d=max(a.keys())
    # a = (1-u) b ; b_k = sum_{j<=k} a_j
    b={}
    s=Fr(0)
    for k in range(0,d):
        s+=a.get(k,Fr(0))
        if s!=0: b[k]=s
    s+=a.get(d,Fr(0))
    assert s==0, "not divisible by 1-u"
    return b

def solve_hpf(N,free):
    ph=[ONE]
    rs=[{}]
    for n in range(1,N+1):
        Phi=ph+[{}]            # ph_n=0
        res=HPF(Phi,n+1)
        rn=pdiv_1mu(res[n+1])
        rs.append(rn)
        p={}
        for k,v in rn.items():
            if 5*k==4*n:
                raise Exception("resonance violated at n=%d: %s"%(n,v))
            p[k]=-v/(4*n-5*k)
        if n%5==0:
            m=n//5
            f=free.get(m,Fr(0))
            if f!=0: p[4*m]=f
        ph.append(p)
    return ph,rs

N=11
ph,rs=solve_hpf(N,{})
print("r_1 =",rs[1])
for n in range(0,8):
    print("ph_%d ="%n, dict(sorted(ph[n].items())))
print("degrees:",[max(p.keys()) if p else 0 for p in ph])
# verify the whole HPF vanishes through zh^{N+1}? (coefficient n+1 uses ph_n) check up to N
res=HPF(ph,N)
print("HPF residual through zh^N zero:",all(not r for r in res[:N+1]))

# ---------- T-system, colour 0
def Tk(k,F,N):
    # T_k F = -F + (t/5)(k - 4 t d/dt) F ; on t^n: -t^n + (k-4n)/5 t^{n+1}
    G=[-F[n] for n in range(N+1)]
    for n in range(N):
        G[n+1]+=Fr(k-4*n,5)*F[n]
    return G
def tsystem(N):
    F0=[Fr(1)]+[Fr(0)]*N
    # closing: F0 = -T0 T1 T2 T_{-2} T_{-1} F0 ; solve order by order
    for n in range(1,N+1):
        # coefficient of t^{n+1} of closing relation determines F0[n]; do linear solve by two evaluations
        def resid(x):
            G=list(F0); G[n]=x
            G=G+[Fr(0)]
            H=Tk(-1,G,N+1); H=Tk(-2,H,N+1); H=Tk(2,H,N+1); H=Tk(1,H,N+1); H=Tk(0,H,N+1)
            return G[n+1]+H[n+1] if n+1<=N else None
        # extend arrays to N+1
        def resid2(x):
            G=list(F0)+[Fr(0)]
            G[n]=x
            H=Tk(-1,G,N+1); H=Tk(-2,H,N+1); H=Tk(2,H,N+1); H=Tk(1,H,N+1); H=Tk(0,H,N+1)
            return G[n+1]+H[n+1]
        r0=resid2(Fr(0)); r1=resid2(Fr(1))
        F0[n]=-r0/(r1-r0)
    G=list(F0)
    F1=Tk(-1,G,N); F2=[-x for x in Tk(-2,F1,N)]; F3=Tk(2,F2,N); F4=[-x for x in Tk(1,F3,N)]
    back=[-x for x in Tk(0,F4,N)]
    return [F0,F1,F2,F3,F4],back
Fs,back=tsystem(N)
print("Fs_0 =",Fs[0][:6])
print("closing check (first N coeffs):",all(back[n]==Fs[0][n] for n in range(N)))
print("Fs_b(0) =",[F[0] for F in Fs])
print("ph_n(0) == [t^n]Fs_0 :",[ph[n].get(0,Fr(0))==Fs[0][n] for n in range(N+1)])

# ---------- h = h_{w_0}
def bern(n):
    B=[Fr(0)]*(n+1)
    A=[Fr(0)]*(n+1)
    for m in range(n+1):
        A[m]=Fr(1,m+1)
        for j in range(m,0,-1):
            A[j-1]=j*(A[j-1]-A[j])
        B[m]=A[0]
    return B
BB=bern(2*N+2)
BB[1]=Fr(-1,2)
# power sums of x=1/(1-zeta^j): roots of 5x^4-10x^3+10x^2-5x+1
e=[Fr(1),Fr(2),Fr(2),Fr(1),Fr(1,5)]
p=[Fr(4)]
for m in range(1,2*N+2):
    s=Fr(0)
    for i in range(1,min(m,4)+1):
        s+=(-1)**(i-1)*e[i]*(p[m-i] if m-i>0 else Fr(0))
    if m<=4: s+=(-1)**(m-1)*m*e[m]
    p.append(s)
S=lambda m: p[m]+Fr(1,(-5)**m)
logh=[Fr(0)]*(N+1)
k=1
while 2*k-1<=N:
    logh[2*k-1]=-BB[2*k]/(2*k*(2*k-1))*S(2*k-1)
    k+=1
# exp
h=[Fr(1)]+[Fr(0)]*N
for n in range(1,N+1):
    h[n]=sum(Fr(j)*logh[j]*h[n-j] for j in range(1,n+1))/n
print("h =",h[:6])
val1=lambda P: sum(P.values(),Fr(0))
print("ph_n(1) (free coeffs 0):",[val1(ph[n]) for n in range(N+1)])
print("ph_n(1)==h_n:",[val1(ph[n])==h[n] for n in range(N+1)])
print("free coefficient needed at n=5 for ph_5(1)=h_5:",h[5]-val1(ph[5]))
