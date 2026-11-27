# Schur side of MV mod a 61-bit prime; own implementation (log S exponents, hook power sums).
# (i) (*) coefficientwise in p_1, x, with controls. (ii) w-polynomial form of G_h + Z_h vs beta, h up to hmax.
import sys
P=(1<<61)-1
def inv(a): return pow(a%P,P-2,P)
def fr(n,d=1): return n%P*inv(d)%P
from math import comb,factorial
def bern(n):
    B=[0]*(n+1);B[0]=1
    for m in range(1,n+1):
        B[m]=(-sum(comb(m+1,k)*B[k] for k in range(m))*inv(m+1))%P
    return B
BN=bern(260)
def partitions(n,m=None):
    if m is None:m=n
    if n==0: yield ();return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k): yield (k,)+r
def hooks(nu):
    conj=[sum(1 for t in nu if t>j) for j in range(nu[0])] if nu else []
    return [nu[i]-j+conj[j]-i-1 for i in range(len(nu)) for j in range(nu[i])]
def kappa(nu): return sum(t*(t-2*(i+1)+1) for i,t in enumerate(nu))
def logS(L):
    # log S(y)=log((y/2)/sinh(y/2)) = -sum_{k>=1} B_{2k} y^{2k}/(2k (2k)!)
    c=[0]*(L+1)
    fac=1
    for n in range(1,L+1):
        fac=fac*n%P
        if n%2==0: c[n]=(-BN[n]*inv(n)%P*inv(fac))%P
    return c
def series_exp(a,L):
    # exp of a(x) with a[0]=0, length L+1
    o=[0]*(L+1);o[0]=1
    for n in range(1,L+1):
        s=0
        for k in range(1,n+1):
            if a[k]: s+=k*a[k]%P*o[n-k]
        o[n]=s%P*inv(n)%P
    return o
def mul(a,b,L):
    o=[0]*(L+1)
    for i,x in enumerate(a):
        if x and i<=L:
            for j in range(0,L+1-i):
                y=b[j]
                if y:o[i+j]+=x*y
    return [t%P for t in o]
def build(tau,nmax,hmaxx):
    # F[a] stored relative: Zt[n]=x^n Z_n series in x length Ln
    Lmax=2*nmax+2*hmaxx+4
    ls=logS(Lmax)
    Z=[None]*(nmax+1); Zd=[None]*(nmax+1)
    for n in range(nmax+1):
        tot=[0]*(Lmax+1);dt=[0]*(Lmax+1)
        for nu in partitions(n):
            H=hooks(nu);ka=kappa(nu)
            a=[0]*(Lmax+1)
            a[1]=(tau+fr(1,2))*fr(ka,2)%P
            for k in range(2,Lmax+1,2):
                if ls[k]:
                    ps=sum(pow(h,k,P) for h in H)%P
                    a[k]=(a[k]+ls[k]*ps)%P
            ex=series_exp(a,Lmax)
            pre=1
            for h in H: pre=pre*inv(h*h%P)%P
            for i in range(Lmax+1): tot[i]=(tot[i]+pre*ex[i])%P
            cc=pre*fr(ka,2)%P
            for i in range(Lmax): dt[i+1]=(dt[i+1]+cc*ex[i])%P
        Z[n]=tot;Zd[n]=dt
    return Z,Zd,Lmax
def logZ(Z,Zd,nmax,L):
    # F_n, Fd_n with relative shift: Z_n=x^-n Z[n]; F_n=x^-n f[n]
    f=[None]*(nmax+1);fd=[None]*(nmax+1)
    for n in range(1,nmax+1):
        s=[n*t%P for t in Z[n]]
        for k in range(1,n):
            m=mul(f[k],Z[n-k],L)
            for i in range(L+1): s[i]=(s[i]-k*m[i])%P
        f[n]=[t*inv(n)%P for t in s]
    # Fd_n = Zd_n - sum_{k=1}^{n-1} (k/n?) ... use d log Z = dZ/Z: n-grading: Fd_n = Zd_n - sum_{k=1}^{n-1} Fd_k Z_{n-k}
    for n in range(1,nmax+1):
        s=list(Zd[n])
        for k in range(1,n):
            m=mul(fd[k],Z[n-k],L)
            for i in range(L+1): s[i]=(s[i]-m[i])%P
        fd[n]=s
    return f,fd
def expx(c,L):
    o=[1]*(L+1)
    for i in range(1,L+1): o[i]=o[i-1]*c%P*inv(i)%P
    return o
def star_check(tau,nmax,L,sign_quad=1,sign_exp=1):
    Z,Zd,Lm=build(tau,nmax,L)
    f,fd=logZ(Z,Zd,nmax,Lm)
    bad=0;tot=0
    # coefficient of p^n: all multiplied by x^{n}: shift bookkeeping
    # E_n=e^{n tau x} (or sign_exp); terms: (2/x)(1-E_n)Fd_n + (n^2-n)(1+E_n)F_n + sum_{i+j=n} (1-E_i) i F_i (1-E_j) j F_j
    Es={i:expx(sign_exp*i*tau%P,Lm) for i in range(1,nmax+1)}
    for n in range(1,nmax+1):
        E=Es[n]
        oneminus=[(-t)%P for t in E];oneminus[0]=(1-E[0])%P
        # (1-E_n)/x series
        om1=oneminus[1:]+[0]
        t1=mul(om1,fd[n],Lm)
        t1=[2*t%P for t in t1]
        onep=[t for t in E];onep[0]=(1+E[0])%P
        t2=mul(onep,f[n],Lm);t2=[(n*n-n)*t%P for t in t2]
        t3=[0]*(Lm+1)
        for i in range(1,n):
            j=n-i
            ai=[(-t)%P for t in Es[i]];ai[0]=(1-Es[i][0])%P
            aj=[(-t)%P for t in Es[j]];aj[0]=(1-Es[j][0])%P
            m=mul(mul(ai,f[i],Lm),mul(aj,f[j],Lm),Lm)
            for k in range(Lm+1): t3[k]=(t3[k]+i*j*m[k])%P
        for k in range(L+1):
            v=(t1[k]+t2[k]+sign_quad*t3[k])%P
            tot+=1
            if v: bad+=1
    return bad,tot
def beta_modp(n,a,b):
    N=n+4
    def S(c):
        ii=[0]*(N+1)
        for i in range(0,N+1,2): ii[i]=pow(c*fr(1,2)%P,i,P)*inv(factorial(i+1)%P)%P
        o=[0]*(N+1);o[0]=1
        for i in range(1,N+1):
            o[i]=(-sum(ii[j]*o[i-j] for j in range(1,i+1)))%P
        return o
    A=S(a);B=S(b)
    co=sum(A[i]*B[n+2-i] for i in range(n+3))%P*inv(a*b%P)%P
    return factorial(n-1)%P*co%P
def yser_w(b,Y):
    # [y^a] w^b, w=1/(1-u), u=T(y): Lagrange: [y^a] H(u)= (1/a)[u^{a-1}] H'(u) e^{a u}; H=(1-u)^-b
    o=[0]*(Y+1);o[0]=1
    for a in range(1,Y+1):
        s=0
        for i in range(a):
            # H' coefficient u^i: b*(1-u)^(-b-1) -> b*C(b+i,i)
            c=b*comb(b+i,i)%P
            s+=c*pow(a,a-1-i,P)%P*inv(factorial(a-1-i)%P)
        o[a]=s%P*inv(a)%P
    return o
def extract(tau,hmax):
    sig=tau*(1+tau)%P
    ymax=3*hmax+1
    Z,Zd,Lm=build(tau,ymax,hmax)
    f,fd=logZ(Z,Zd,ymax,Lm)
    res={}
    for h in range(2,hmax+1):
        # [y^a] sigma^{h-1} x^{2h-2} G_h = [x^{2h-2}] F_a/(sigma x)^a ; F_a = x^-a f[a]
        data={}
        for a in range(1,3*h+2):
            if a>ymax: break
            co=f[a][2*h-2+2*a]
            data[a]=co*inv(pow(sig,a,P))%P*inv(pow(sig,h-1,P))%P
        bs=list(range(2*h-2,5*h-4))
        ys={b:yser_w(b,max(data)) for b in bs}
        rows=[[ys[b][a] for b in bs]+[data[a]] for a in sorted(data)]
        # solve
        m=len(bs);piv=[];r=0
        for c in range(m):
            p=None
            for i in range(r,len(rows)):
                if rows[i][c]: p=i;break
            if p is None: continue
            rows[r],rows[p]=rows[p],rows[r]
            iv=inv(rows[r][c]);rows[r]=[t*iv%P for t in rows[r]]
            for i in range(len(rows)):
                if i!=r and rows[i][c]:
                    fct=rows[i][c];rows[i]=[(x-fct*y)%P for x,y in zip(rows[i],rows[r])]
            piv.append(c);r+=1
        over=[rows[i][-1] for i in range(r,len(rows))]
        consistent=all(v==0 for v in over)
        Zh=rows[0][-1] if piv and piv[0]==0 else None
        beta=beta_modp(2*h-2,(-tau)%P,(1+tau)%P)
        lhs=pow(sig,h-1,P)*Zh%P
        res[h]=(r==m,consistent,len(over),lhs==beta,lhs==(beta*sig%P))
        print(f" h={h}: full rank {r==m}, extra eqs {len(over)} consistent {consistent}, sigma^(h-1)Z_h==beta: {lhs==beta}; control beta*sigma: {lhs==(beta*sig%P)}");sys.stdout.flush()
if __name__=="__main__":
    t1=fr(3,11)
    for t,name in [(t1,'3/11'),(fr(-7,4),'-7/4')]:
        print('tau',name,'(*) n<=12, 20 x-orders:',star_check(t,12,20))
    print('control flip quad',star_check(t1,8,12,sign_quad=-1))
    print('control e^{-beta D}',star_check(t1,8,12,sign_exp=-1))
    hmax=int(sys.argv[1])
    for t,name in [(t1,'3/11'),(fr(-7,4),'-7/4'),(fr(13,5),'13/5')]:
        print('tau',name);extract(t,hmax)
