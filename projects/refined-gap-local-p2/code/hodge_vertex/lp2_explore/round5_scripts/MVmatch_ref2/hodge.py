# Referee own Hodge engine: WK by DVV; ch classes by Mumford; independent of author code
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
from math import factorial
import sys
sys.setrecursionlimit(100000)

def dfact(n):
    r=1
    while n>1:
        r*=n;n-=2
    return r

@lru_cache(maxsize=None)
def wk(g,d):
    n=len(d)
    if g<0: return F(0)
    if sum(d)!=3*g-3+n: return F(0)
    if 2*g-2+n<=0: return F(0)
    if n==0: return F(0)
    if g==0 and n==3: return F(1)
    if g==1 and d==(1,): return F(1,24)
    if d[0]==0:
        rest=d[1:]
        s=F(0)
        for i in range(len(rest)):
            if rest[i]>0:
                t=list(rest);t[i]-=1
                s+=wk(g,tuple(sorted(t)))
        return s
    if 1 in d:
        i=d.index(1); rest=d[:i]+d[i+1:]
        return (2*g-2+len(rest))*wk(g,rest)
    k=d[-1]; rest=list(d[:-1])
    s=F(0)
    # DVV: <tau_k prod tau_d>: (2k+1)!! <...> = sum_j (2k+2d_j-1)!!/(2d_j-1)!! <tau_{k+d_j-1}...>
    #  + 1/2 sum_{a+b=k-2} (2a+1)!!(2b+1)!! [ <tau_a tau_b ...>_{g-1} + sum splits ]
    for j in range(len(rest)):
        t=rest[:j]+rest[j+1:]+[k+rest[j]-1]
        s+=F(dfact(2*k+2*rest[j]-1),dfact(2*rest[j]-1))*wk(g,tuple(sorted(t)))
    tot=F(0)
    for a in range(k-1):
        b=k-2-a
        c=dfact(2*a+1)*dfact(2*b+1)
        tot+=c*wk(g-1,tuple(sorted(rest+[a,b])))
        m=len(rest)
        for mask in range(1<<m):
            I=[rest[i] for i in range(m) if mask>>i&1]
            J=[rest[i] for i in range(m) if not mask>>i&1]
            num=sum(I)+a-(len(I)+1)+3
            if num%3: continue
            g1=num//3
            g2=g-g1
            if g1<0 or g2<0: continue
            v1=wk(g1,tuple(sorted(I+[a])))
            if v1==0: continue
            tot+=c*v1*wk(g2,tuple(sorted(J+[b])))
    return (s+tot/2)/dfact(2*k+1)

def bern(N):
    B=[F(0)]*(N+1);A=[F(0)]*(N+1)
    for m in range(N+1):
        A[m]=F(1,m+1)
        for j in range(m,0,-1):
            A[j-1]=j*(A[j-1]-A[j])
        B[m]=A[0]
    return B
BB=bern(60)

@lru_cache(maxsize=None)
def ch(g,S,d):
    """<prod_{k in S} ch_k prod tau_{d_i}>_{g,n}; S sorted tuple of odd ints, d sorted tuple"""
    n=len(d)
    if g<0 or 2*g-2+n<=0: return F(0)
    if sum(d)+sum(S)!=3*g-3+n: return F(0)
    if not S:
        return wk(g,d)
    if g==0: return F(0)
    if S[-1]>2*g-1: return F(0)
    k=S[-1]; R=S[:-1]; l=(k+1)//2
    pref=BB[2*l]/factorial(2*l)
    tot=ch(g,R,tuple(sorted(d+(2*l,))))
    for i in range(n):
        t=list(d);t[i]+=2*l-1
        tot-=ch(g,R,tuple(sorted(t)))
    bd=F(0)
    mR=len(R)
    for a in range(2*l-1):
        b=2*l-2-a
        sg=(-1)**a
        bd+=sg*ch(g-1,R,tuple(sorted(d+(a,b))))
        for mr in range(1<<mR):
            R1=tuple(R[i] for i in range(mR) if mr>>i&1)
            R2=tuple(R[i] for i in range(mR) if not mr>>i&1)
            for mask in range(1<<n):
                I=tuple(d[i] for i in range(n) if mask>>i&1)
                J=tuple(d[i] for i in range(n) if not mask>>i&1)
                num=sum(I)+a+sum(R1)-(len(I)+1)+3
                if num%3: continue
                g1=num//3;g2=g-g1
                if g1<0 or g2<0: continue
                if 2*g1-2+len(I)+1<=0 or 2*g2-2+len(J)+1<=0: continue
                v1=ch(g1,R1,tuple(sorted(I+(a,))))
                if v1==0: continue
                bd+=sg*v1*ch(g2,R2,tuple(sorted(J+(b,))))
    return pref*(tot+bd/2)

def omega_monos(g,ts,maxdeg):
    """prod_{t in ts} Lambda(t) as dict S->coef; log Lambda(t)=sum_k (-1)^(k-1)(k-1)! t^k ch_k, even ch = 0"""
    res={():F(1)}
    for k in range(1,2*g,2):
        ck=factorial(k-1)*sum(t**k for t in ts)
        new={}
        for S,c in res.items():
            m=0;deg=sum(S);term=c
            while deg+m*k<=maxdeg:
                key=tuple(sorted(S+(k,)*m))
                new[key]=new.get(key,F(0))+term
                m+=1
                term=term*ck/m
        res=new
    return res

def compositions(total,n):
    if n==0:
        if total==0: yield ()
        return
    if n==1:
        yield (total,);return
    for i in range(total+1):
        for r in compositions(total-i,n-1):
            yield (i,)+r

def omega_ts(tau):
    return [F(-1),-1/F(tau),1/(1+F(tau))]

def M(g,tau,leaves,ts=None):
    """Int Omega prod 1/(1-l psi) over Mbar_{g,n}; genus 0 unstable convention NOT here"""
    if ts is None: ts=omega_ts(tau)
    n=len(leaves)
    dim=3*g-3+n
    mon=omega_monos(g,ts,dim)
    tot=F(0)
    for S,c in mon.items():
        if c==0: continue
        rem=dim-sum(S)
        for d in compositions(rem,n):
            w=F(1)
            for l,e in zip(leaves,d):
                w*=F(l)**e
            if w==0: continue
            tot+=c*w*ch(g,S,tuple(sorted(d)))
    return tot

def parts_min(n,minp):
    if n==0:
        yield ();return
    for p in range(minp,n+1):
        for r in parts_min(n-p,p):
            yield (p,)+r

def Qk(g,tau,R,k,sgn=1):
    """(1/k!) Int Omega lam_R Tt^k over Mbar_{g,r+k}, literal T-legs Tt=psi^2/(1-psi)"""
    ts=omega_ts(tau)
    r=len(R)
    N=3*g-3+r
    dim=N+k
    mon=omega_monos(g,ts,dim)
    tot=F(0)
    for S,c in mon.items():
        if c==0: continue
        rem=dim-sum(S)
        if rem<2*k: continue
        for tdeg in range(2*k,rem+1):
            ptl=[()] if (k==0 and tdeg==0) else ([p for p in parts_min(tdeg,2) if len(p)==k] if k>0 else [])
            for pt in ptl:
                mult=factorial(k)
                for v in Counter(pt).values(): mult//=factorial(v)
                for d in compositions(rem-tdeg,r):
                    w=F(1)
                    for l,e in zip(R,d):
                        w*=F(l)**e
                    if w==0: continue
                    tot+=c*w*mult*ch(g,S,tuple(sorted(d+pt)))
    return tot/factorial(k)

def Xi_def(g,tau,R):
    N=3*g-3+len(R)
    return sum((-1)**k*Qk(g,tau,R,k) for k in range(N+1))

if __name__=="__main__":
    print("wk",wk(2,(4,)),wk(3,(7,)),wk(2,(2,3)), "expect 1/1152 1/82944 29/5760")
    print("lam1 M11",ch(1,(1,),(0,)),"1/24")
    print("lam1^3 M2",ch(2,(1,1,1),()),"1/2880")
    print("lam_g psi^(2g-2): expect b_g")
    for g,exp in [(1,F(1,24)),(2,F(7,5760)),(3,F(31,967680)),(4,F(127,154828800)),(5,F(73,3503554560))]:
        mon=omega_monos(g,[F(1)],3*g-2)
        v=F(0)
        for S,c in mon.items():
            if sum(S)==g:
                v+=c*ch(g,S,(2*g-2,))
        print(g,v,exp,v==exp)
    # A_g: lambda_{g-2}lambda_{g-1}lambda_g
    def A(g):
        B=BB
        return F(1,2*factorial(2*g-2))*abs(B[2*g-2])/(2*g-2)*abs(B[2*g])/(2*g)
    for g in range(2,6):
        for tau in [F(2,7),F(-5,3)]:
            ts=omega_ts(tau); dim=3*g-3
            mon=omega_monos(g,ts,dim)
            v=sum(c*ch(g,S,()) for S,c in mon.items() if sum(S)==dim)
            sig=tau*(1+tau)
            print(g,tau,v==-A(g)*sig**(1-g))
