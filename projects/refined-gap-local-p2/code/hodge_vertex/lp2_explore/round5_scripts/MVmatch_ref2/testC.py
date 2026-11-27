# (C): FF computed literally from the Schur sum in K((x))[[s,P]] and compared with the table of QQ.
# Schur-only tests: termination of s-series, A_g constant (g up to 6), genus 0/1 unstable entries.
from fractions import Fraction as F
from schur import *
from hodge import BB,Qk,Xi_def
from math import comb,factorial
import sys,time

def smul(a,b,M):
    r=[F(0)]*(M+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>M: break
            r[i+j]+=x*y
    return r

def binom_series(k,M):
    # (1+s)^k any integer k
    out=[F(1)]
    for m in range(1,M+1):
        out.append(out[-1]*F(k-m+1,m))
    return out

def e_minus_u(M):
    # e^{-u}, u=s/(1+s) : exp series
    u=[F(0)]+[F((-1)**(m-1)) for m in range(1,M+1)]
    res=[F(0)]*(M+1);res[0]=F(1)
    # exp(-u): n e_n = sum k (-u_k) e_{n-k}
    for n in range(1,M+1):
        res[n]=sum(k*(-u[k])*res[n-k] for k in range(1,n+1))/n
    return res

class FF:
    def __init__(s,R,tau,M,G,xshift=True,eu=True,sgn=1,usePi=True):
        s.R=tuple(R);l=len(R);s.l=l;tau=F(tau);s.tau=tau;sig=tau*(1+tau);s.sig=sig;s.M=M;s.J=sum(R);s.r=l
        P=2*G-2+2*M+s.J+l+1;s.P=P
        em=e_minus_u(M)
        epow=[None]*(M+s.J+2)
        epow[0]=[F(1)]+[F(0)]*M
        for i in range(1,len(epow)):
            epow[i]=smul(epow[i-1],em,M)
        s.Z={}
        for m in range(1<<l):
            sub=tuple(R[i] for i in range(l) if m>>i&1)
            cnt=len(sub);JS=sum(sub)
            pr=1
            piS=F(1)
            for q in sub:
                pr*=q
                if usePi: piS*=Pi(q,tau)
            sigc=sig**cnt*piS
            arr=[[F(0)]*P for _ in range(M+1)]  # arr[mm][idx], exponent = -JS-cnt-2mm+idx
            for a in range(M+1):
                n=JS+a
                A=[F(0)]*P
                for nu in partitions(n):
                    c=chi(nu,tuple(sorted(sub,reverse=True))+(1,)*a) if n>0 else 1
                    if c==0: continue
                    y=Y(nu,tau,P)
                    for i in range(P): A[i]+=c*y[i]
                f=F(sgn**a,pr*factorial(a))/(sig**a*sigc)
                for i in range(P):
                    if A[i]==0: continue
                    if xshift:
                        ser=binom_series(n-i,M-a)
                    else:
                        ser=[F(1)]+[F(0)]*(M-a)
                    if eu:
                        ser=smul(ser,epow[n][:M-a+1],M-a)
                    for dm,v in enumerate(ser):
                        idx=i+2*dm
                        if idx<P:
                            arr[a+dm][idx]+=A[i]*f*v
            s.Z[m]=arr
        # log
        Z0=s.Z[0]
        c0=[None]*(M+1)
        c0[0]=[F(0)]*P
        for mm in range(1,M+1):
            r=[mm*v for v in Z0[mm]]
            for b in range(1,mm):
                pp=mul(c0[b],Z0[mm-b],P)
                for i in range(P): r[i]-=b*pp[i]
            c0[mm]=[v/mm for v in r]
        s.cum={0:c0}
        for m in range(1,1<<l):
            low=m&-m;rest=m^low
            rows=[list(v) for v in s.Z[m]]
            for sp in submasks(rest):
                B=low|sp
                if B==m: continue
                o=m^B
                for a in range(M+1):
                    for b in range(a+1):
                        pp=mul(s.cum[B][b],s.Z[o][a-b],P)
                        for i in range(P): rows[a][i]-=pp[i]
            out=[None]*(M+1)
            for a in range(M+1):
                r=list(rows[a])
                for b in range(a):
                    pp=mul(out[b],Z0[a-b],P)
                    for i in range(P): r[i]-=pp[i]
                out[a]=r
            s.cum[m]=out
    def QQ(s,g):
        """list of [s^m] sigma^{1-g}[x^{2g-2}] d_P_R FF|_{P=0}, m=0..M"""
        m=(1<<s.l)-1
        res=[]
        for mm in range(s.M+1):
            idx=2*g-2+s.J+s.r+2*mm
            res.append(s.cum[m][mm][idx]*s.sig**(1-g))
        return res
    def forbidden(s):
        # x exponents that are not even >= -2 : should vanish
        m=(1<<s.l)-1
        bad=0
        for mm in range(s.M+1):
            for idx,v in enumerate(s.cum[m][mm]):
                e=-s.J-s.r-2*mm+idx
                if (e<-2 or e%2) and v!=0:
                    bad+=1
        return bad

def A(g):
    return F(1,2*factorial(2*g-2))*abs(BB[2*g-2])/(2*g-2)*abs(BB[2*g])/(2*g)

def series_rat_invpow(k,M):
    return binom_series(-k,M)

if __name__=="__main__":
    M=int(sys.argv[1]);G=int(sys.argv[2])
    for tau in [F(5,3),F(-1,2)]:
        sig=tau*(1+tau)
        t0=time.time()
        ff=FF((),tau,M,G)
        print("tau",tau,"forbidden",ff.forbidden(),round(time.time()-t0,1),flush=True)
        # genus 0
        q0=ff.QQ(0)
        ex=[F(0)]*(M+1)
        ex[1]+=1
        if M>=2: ex[2]+=F(1,4)
        inv=binom_series(-1,M)
        for m in range(M+1):
            if m>=3: ex[m]+=inv[m-3]*F(1,6)
        print(" (0,empty) ok",q0==ex)
        q1=ff.QQ(1)
        ex=[F(0)]*(M+1)
        for m in range(1,M+1):
            ex[m]=F((-1)**(m-1),m)*F(1,24)+F(1,24)*(1+1/sig)*inv[m]
        print(" (1,empty) ok",q1==ex)
        # control: wrong (1+1/sigma) -> 1/sigma
        ex2=[F(0)]*(M+1)
        for m in range(1,M+1):
            ex2[m]=F((-1)**(m-1),m)*F(1,24)+F(1,24)*(1/sig)*inv[m]
        print(" control(1,empty) must fail:",q1==ex2)
        for g in range(2,G+1):
            q=ff.QQ(g)
            N=3*g-3
            om=-A(g)*sig**(1-g)
            ser=series_rat_invpow(2*g-2,M)
            full=[q[m]+om*ser[m] for m in range(M+1)]
            tail=full[N+1:]
            ctl=[q[m]-om*ser[m] for m in range(M+1)][N+1:]
            ctl0=q[N+1:]
            print(" g",g,"N",N,"tail(with <Omega>) zero:",all(v==0 for v in tail),"len",len(tail),"| control -<Omega> zero:",all(v==0 for v in ctl),"| control none zero:",all(v==0 for v in ctl0),flush=True)
            if g<=4 and tau==F(5,3):
                agree=all(full[m]==Qk(g,tau,(),m) for m in range(N+1))
                print("   coefficients vs literal Hodge def:",agree,flush=True)
            # value at s=-1 of polynomial part and the pole
            print("   Xi_g(empty)=",sum((-1)**m*full[m] for m in range(N+1)),"; <Omega>=",om)
