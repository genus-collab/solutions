#!/usr/bin/env python
"""Exact finite-order evidence for the MV/mKP proof of (VI). Evidence only.
 Part 1: (*) on the Schur side, p_1 slice, series in p and x.
 Part 2: Z_h extracted from the MV Schur side by the w-polynomial fit, compared with beta_{2h-2}.
 Part 3: order-eps^1 equation: L[B_claim] = RHS as series in 1/N, tau-derivative exact via poly in tau.
"""
from fractions import Fraction as Fr
import sys, math

# --- truncated Laurent series in x: (val, list) with fixed absolute precision
class LS:
    __slots__=('v','c')
    def __init__(s,v,c): s.v=v; s.c=c   # sum c[i] x^(v+i); truncated at len
    def prec(s): return s.v+len(s.c)
def ls_add(a,b):
    v=min(a.v,b.v); p=min(a.prec(),b.prec()); n=max(p-v,0)
    c=[Fr(0)]*n
    for i in range(n):
        k=v+i
        if a.v<=k<a.prec(): c[i]+=a.c[k-a.v]
        if b.v<=k<b.prec(): c[i]+=b.c[k-b.v]
    return LS(v,c)
def ls_scale(a,f): return LS(a.v,[f*t for t in a.c])
def ls_mul(a,b):
    n=min(len(a.c),len(b.c)); c=[Fr(0)]*n
    for i,t in enumerate(a.c[:n]):
        if t:
            for j,u in enumerate(b.c[:n-i]): c[i+j]+=t*u
    return LS(a.v+b.v,c)
def ls_inv(a):
    n=len(a.c); c=[Fr(0)]*n; c[0]=1/a.c[0]
    for i in range(1,n):
        s=Fr(0)
        for j in range(1,i+1): s+=a.c[j]*c[i-j]
        c[i]=-s/a.c[0]
    return LS(-a.v,c)
def ls_exp_lin(lam,n):  # exp(lam x) to n terms
    c=[Fr(1)]
    for i in range(1,n): c.append(c[-1]*lam/i)
    return LS(0,c)
def two_sinh(h,n):  # 2 sinh(h x/2)= h x (1+...)  n terms from x^1
    c=[]
    for i in range(n):
        c.append(Fr(h)**(2*i+1)/(Fr(2)**(2*i)*math.factorial(2*i+1)) if True else 0)
    # stored in x^(1+2i); expand
    d=[Fr(0)]*n
    for i in range(n):
        if 2*i<n: d[2*i]=c[i]
    return LS(1,d)

def partitions(n,m=None):
    if m is None: m=n
    if n==0: yield (); return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k): yield (k,)+r
def hooks(nu):
    conj=[sum(1 for t in nu if t>j) for j in range(nu[0])] if nu else []
    return [nu[i]-j-1+conj[j]-i-1+1 for i in range(len(nu)) for j in range(nu[i])]
def kappa(nu): return sum(t*(t-2*(i+1)+1) for i,t in enumerate(nu))

def schur_coeffs(tau,amax,prec):
    """z[n] = [p^n] Z as LS with val -n and `prec` terms; also dz[n]=[p^n] dZ/dtau. Z= sum_nu e^{(tau+1/2)kappa x/2} prod 1/(h 2sinh(hx/2)) p^|nu|."""
    z=[None]*(amax+1); dz=[None]*(amax+1)
    for n in range(0,amax+1):
        tot=LS(-n,[Fr(0)]*prec); dtot=LS(-n,[Fr(0)]*prec)
        for nu in partitions(n):
            ser=ls_exp_lin((tau+Fr(1,2))*Fr(kappa(nu),2),prec)
            for h in hooks(nu):
                ser=ls_mul(ser,ls_scale(ls_inv(two_sinh(h,prec)),Fr(1,h)))
            tot=ls_add(tot,ser)
            # d/dtau multiplies by x kappa/2
            dser=LS(ser.v+1,[t*Fr(kappa(nu),2) for t in ser.c[:-1]])
            dser=LS(dser.v,dser.c)
            dtot=ls_add(dtot,dser) if n>0 else dtot
        z[n]=tot; dz[n]=dtot
    return z,dz

def log_series(z,amax,prec):
    """f[n]=[p^n] log Z (n>=1) via n f_n = n z_n - sum_{k=1}^{n-1} k f_k z_{n-k}"""
    f=[None]*(amax+1)
    for n in range(1,amax+1):
        s=ls_scale(z[n],Fr(n))
        for k in range(1,n):
            s=ls_add(s,ls_scale(ls_mul(f[k],z[n-k]),Fr(-k)))
        f[n]=ls_scale(s,Fr(1,n))
    return f

def dlog_series(z,dz,f,amax):
    """df[n]=[p^n] d/dtau log Z: Z' = Z F' => dz_n = sum_{k=1}^n df_k z_{n-k}"""
    df=[None]*(amax+1)
    for n in range(1,amax+1):
        s=dz[n]
        for k in range(1,n):
            s=ls_add(s,ls_scale(ls_mul(df[k],z[n-k]),Fr(-1)))
        df[n]=s
    return df

def part1(tau,amax=7,prec=26):
    z,dz=schur_coeffs(tau,amax,prec)
    f=log_series(z,amax,prec); df=dlog_series(z,dz,f,amax)
    # (*): (2/x)(F_t(T)-F_t(T+beta))+(D^2-D)(F(T)+F(T+beta))+(F'(T)-F'(T+beta))^2 = 0
    # p^n coefficient: e^{n beta} shift. D -> n.
    bad=0; report=[]
    for n in range(1,amax+1):
        eb=ls_exp_lin(tau*n,prec)
        # T1 = (2/x)(1-e^{n beta}) df_n
        one_minus=LS(0,[Fr(1)-eb.c[0]]+[-t for t in eb.c[1:]])
        T1=ls_scale(ls_mul(one_minus,df[n]),Fr(2)); T1=LS(T1.v-1,T1.c)
        onep=LS(0,[Fr(1)+eb.c[0]]+eb.c[1:])
        T2=ls_scale(ls_mul(onep,f[n]),Fr(n*n-n))
        T3=None
        for k in range(1,n):
            ek=ls_exp_lin(tau*k,prec); el=ls_exp_lin(tau*(n-k),prec)
            a=LS(0,[Fr(1)-ek.c[0]]+[-t for t in ek.c[1:]]); b=LS(0,[Fr(1)-el.c[0]]+[-t for t in el.c[1:]])
            t=ls_scale(ls_mul(ls_mul(ls_mul(a,b),f[k]),f[n-k]),Fr(k*(n-k)))
            T3=t if T3 is None else ls_add(T3,t)
        tot=ls_add(T1,T2)
        if T3 is not None: tot=ls_add(tot,T3)
        nz=[i for i,t in enumerate(tot.c) if t!=0]
        report.append((n,tot.v,tot.prec(),len(nz)))
        if nz: bad+=1
    print("Part1 tau=%s (p^n, val, prec, #nonzero):"%tau,report," fails:",bad)
    # control: flip sign of T1
    return bad

def S_coef(n):
    """[y^{2n}] (y/2)/sinh(y/2)"""
    c=[Fr(0)]*(n+1); sh=[Fr(1,(4**k)*math.factorial(2*k+1)) for k in range(n+1)]
    c[0]=Fr(1)
    for i in range(1,n+1):
        c[i]=-sum(sh[j]*c[i-j] for j in range(1,i+1))
    return c
def beta(h,e1,e2):
    """beta_{2h-2}(e1,e2) = (2h-3)! [t^{2h-2}] 1/(4 sinh sinh) = (2h-3)! /(e1 e2) [t^{2h}] S(e1 t) S(e2 t)"""
    s=S_coef(h); tot=sum(s[i]*s[h-i]*e1**(2*i)*e2**(2*(h-i)) for i in range(h+1))
    return Fr(math.factorial(2*h-3))*tot/(e1*e2)

def tree_series(N):
    # u = T(y) = sum a^(a-1) y^a/a!
    return [Fr(0)]+[Fr(a**(a-1),math.factorial(a)) for a in range(1,N+1)]
def pser_mul(a,b,N):
    c=[Fr(0)]*(N+1)
    for i,t in enumerate(a):
        if t:
            for j,u in enumerate(b):
                if i+j>N: break
                c[i+j]+=t*u
    return c
def solve(M,r):
    n=len(M); A=[row[:]+[r[i]] for i,row in enumerate(M)]
    for i in range(n):
        p=next(k for k in range(i,n) if A[k][i]!=0)
        A[i],A[p]=A[p],A[i]
        A[i]=[t/A[i][i] for t in A[i]]
        for k in range(n):
            if k!=i and A[k][i]!=0:
                A[k]=[s-A[k][i]*t for s,t in zip(A[k],A[i])]
    return [A[i][n] for i in range(n)]

def part2(tau,hmax):
    sig=tau*(1+tau)
    amax=3*hmax-2
    prec=2*hmax+2*amax+6
    z,dz=schur_coeffs(tau,amax,prec)
    f=log_series(z,amax,prec)
    # M_g(1^a)= a! sigma^{-(g+a-1)} [x^{2g-2+a}] f[a]
    res=[]
    # w series
    u=tree_series(amax)
    w=[Fr(1)]+[Fr(0)]*amax
    # 1/(1-u)
    w=[Fr(1)]+[Fr(0)]*amax
    for i in range(1,amax+1):
        w[i]=sum(u[j]*w[i-j] for j in range(1,i+1))
    for h in range(2,hmax+1):
        na=3*h-2
        G=[Fr(0)]*(na+1)
        for a in range(1,na+1):
            k=2*h-2+a
            cf=f[a].c[k-f[a].v]
            G[a]=cf/sig**(h+a-1)   # y^a coefficient of G_h = M/a!
        # G_h - G_h(0) = sum_{k=0}^{3h-3} c_k (w^{2h-2+k} - 1)
        pw=[Fr(1)]+[Fr(0)]*na
        for _ in range(2*h-2): pw=pser_mul(pw,w[:na+1],na)
        cols=[]
        cur=pw
        wn=w[:na+1]+[Fr(0)]*(max(0,na-amax))
        for k in range(3*h-2):
            cols.append(cur[:]); cur=pser_mul(cur,wn,na)
        Mat=[[cols[k][a] for k in range(3*h-2)] for a in range(1,na+1)]
        c=solve(Mat,G[1:])
        Zh=c[0]
        lhs=sig**(h-1)*Zh
        rhs=beta(h,-tau,1+tau)
        res.append((h,str(lhs),str(rhs),lhs==rhs))
        # Also: <Omega>_h = sum c_k = classical value check is not used
    print("Part2 tau=%s:"%tau,res)
    return res

# --- Part 3: polynomials in tau as dict-free lists; rational functions handled by exact eval + derivative via dual numbers
class Dual:
    __slots__=('a','b')
    def __init__(s,a,b=Fr(0)): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o):
        o=o if isinstance(o,Dual) else Dual(o); return Dual(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return Dual(-s.a,-s.b)
    def __sub__(s,o):
        o=o if isinstance(o,Dual) else Dual(o); return Dual(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return Dual(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,Dual) else Dual(o); return Dual(s.a*o.a,s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def inv(s): return Dual(1/s.a,-s.b/(s.a*s.a))
    def __truediv__(s,o):
        o=o if isinstance(o,Dual) else Dual(o); return s*o.inv()
    def __rtruediv__(s,o): return Dual(o)*s.inv()
    def __pow__(s,n):
        r=Dual(1)
        for _ in range(n): r=r*s
        return r

def part3(tau0,hmax=14):
    """Check L[B_claim]=RHS through N^{-(2hmax-1)}, with Dual(tau,1) giving d/dtau.
       Series in 1/N: dict n-> coefficient of N^{-n}."""
    t=Fr(tau0); tau=Dual(t,1); e=1+tau; sig=tau*e
    b={2*h-2:beta_dual(h,-tau,e) for h in range(2,hmax+1)}
    nmax=2*hmax-2
    # B' , B_tau
    def shift(ser):  # E f: f(N+tau), f=sum c_n N^{-n}, n>=1
        out={}
        for n,c in ser.items():
            for k in range(0,nmax-n+1):
                bin_=Fr(comb_neg(n,k))
                out[n+k]=out.get(n+k,Fr(0))+c*bin_*t**k
        return out
    Bt={n:c.b for n,c in b.items()}
    Bp={n+1:-n*c.a for n,c in b.items()}
    EBt=shift(Bt); EBp=shift(Bp)
    L={}
    def add(n,v): L[n]=L.get(n,Fr(0))+v
    for n,c in EBt.items(): add(n,-2*c)
    for n,c in Bt.items(): add(n,2*c)
    for n,c in Bp.items(): add(n,-c)
    for n,c in EBp.items(): add(n,-c)
    # -((2N+tau)/e)(EB'-B')
    ea=e.a
    D={}
    for n,c in EBp.items(): D[n]=D.get(n,Fr(0))+c
    for n,c in Bp.items(): D[n]=D.get(n,Fr(0))-c
    for n,c in D.items():
        add(n-1,-2*c/ea); add(n,-t*c/ea)
    # RHS: (1/sig^2){ -l (MN - sig'/12) + tau (M+N)/2 - tau(1+2sig)/24 (1/M+1/N)}
    sg=sig.a; sp=sig.b
    ell={k:(-1)**(k-1)*t**k/k for k in range(1,nmax+4)}
    R={}
    def addR(n,v): R[n]=R.get(n,Fr(0))+v
    for k,c in ell.items():
        addR(k-2,-c); addR(k-1,-t*c); addR(k,sp/12*c)
    addR(-1,t); addR(0,t*t/2)  # tau(M+N)/2 = tau N + tau^2/2
    onM={k+1:(-t)**k for k in range(0,nmax+2)}
    for n,c in onM.items(): addR(n,-t*(1+2*sg)/24*c)
    addR(1,-t*(1+2*sg)/24)
    for n in list(R): R[n]=R[n]/sg**2
    bad=[n for n in range(-1,nmax+1) if L.get(n,Fr(0))!=R.get(n,Fr(0))]
    print("Part3 tau=%s: orders N^{-n}, n=-1..%d: mismatches %s; first RHS coefficient(s):"%(tau0,nmax,bad),[ (n,str(R.get(n,0))) for n in range(-1,5)])
    return bad

def comb_neg(n,k):
    # coefficient of (tau/N)^k in (1+tau/N)^{-n} is (-1)^k C(n+k-1,k)
    return (-1)**k*math.comb(n+k-1,k)
def beta_dual(h,e1,e2):
    s=S_coef(h); tot=Dual(0)
    for i in range(h+1):
        tot=tot+s[i]*s[h-i]*(e1**(2*i))*(e2**(2*(h-i)))
    return Dual(math.factorial(2*h-3))*tot/(e1*e2)

if __name__=="__main__":
    part1(Fr(2,7)); part1(Fr(-5,3))
    part3(Fr(2,7)); part3(Fr(-5,3)); part3(Fr(3,2))
    for tau in (Fr(2,7),Fr(-5,3)):
        part2(tau,int(sys.argv[1]) if len(sys.argv)>1 else 4)
