# Referee's own code (evidence only). Phi(alpha,Y;z) at rational numbers; no ring C, no import from FB_*.
# Phi = exp((z/2) G dd) exp(z^-1 sum_{m>=3} T_m(xi))|_{xi=0}, G=-H^-1, H=diag(q)+cJ, c=-Y/3, q_j=alpha_j+Y.
# Test 1: alpha=alpha*: Phi = Theta(-z/E) S(v,z), v=1/(Y+E/3). Test 2: y=0 point Y=-alpha_i: Phi = prod Theta(-z/a).
import sys,math,time
from fractions import Fraction as Fr
def bern(n):
    B=[Fr(0)]*(n+1);B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(Fr(math.comb(m+1,k))*B[k] for k in range(m))/(m+1)
    return B
B=bern(40)
def exps(l,N):
    r=[Fr(0)]*(N+1);r[0]=Fr(1)
    for n in range(1,N+1):
        r[n]=sum(k*l[k]*r[n-k] for k in range(1,n+1))/n
    return r
def Theta(c,N):  # Theta(c z)
    l=[Fr(0)]*(N+1)
    for m in range(1,(N+1)//2+1):
        if 2*m-1<=N: l[2*m-1]=B[2*m]/(2*m*(2*m-1))*Fr(c)**(2*m-1)
    return exps(l,N)
def mul(a,b,N):
    r=[Fr(0)]*(N+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j<=N: r[i+j]+=x*y
    return r
def rec(N):
    r=[Fr(1)]; rm=Fr(0)
    for m in range(0,N+2):
        nxt=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))/(27*(m+1))
        r.append(nxt)
    return r
def S(v,E,N):
    r=rec(N); out=[Fr(0)]*(N+1)
    for l in range(N+1):
        # r_l (vz)^l / prod_{j<=l+1}(1+jz/E)
        s=[Fr(0)]*(N+1); s[0]=Fr(1)
        for j in range(1,l+2):
            # divide by 1+(j/E) z
            c=Fr(j)/E
            t=[Fr(0)]*(N+1); t[0]=s[0]
            for n in range(1,N+1): t[n]=s[n]-c*t[n-1]
            s=t
        for n in range(N+1-l):
            out[n+l]+=r[l]*v**l*s[n]
    return out
def inv3(M):
    a=M
    det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    adj=[[0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            i1,i2=[k for k in range(3) if k!=j]
            j1,j2=[k for k in range(3) if k!=i]
            adj[i][j]=(-1)**(i+j)*(a[i1][j1]*a[i2][j2]-a[i1][j2]*a[i2][j1])
    return [[adj[i][j]/det for j in range(3)] for i in range(3)],det
def Phi(al,Y,N):
    q=[Fr(x)+Y for x in al]; c=-Fr(Y)/3
    H=[[ (q[i] if i==j else 0)+c for j in range(3)] for i in range(3)]
    Hi,det=inv3(H)
    e1=sum(al);e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2];e3=al[0]*al[1]*al[2]
    W=e1*Y*Y+2*e2*Y+3*e3
    assert det==W/3
    G=[[-Hi[i][j] for j in range(3)] for i in range(3)]
    # polys dict (a,b,c)->coef; grading k: T_{k+2} has weight k
    def Tm(m):
        d={}
        f=Fr(1,math.factorial(m))
        for j in range(3):
            e=[0,0,0];e[j]=m
            d[tuple(e)]=q[j]*f
        co=-3*Y*f/3**m
        for a in range(m+1):
            for b in range(m+1-a):
                cc=m-a-b
                k=(a,b,cc)
                d[k]=d.get(k,0)+co*math.factorial(m)*Fr(1,math.factorial(a)*math.factorial(b)*math.factorial(cc))
        return d
    K=2*N
    T={k:Tm(k+2) for k in range(1,K+1)}
    Es=[{(0,0,0):Fr(1)}]
    for k in range(1,K+1):
        cur={}
        for j in range(1,k+1):
            pref=Fr(j,k)
            for m1,c1 in T[j].items():
                for m2,c2 in Es[k-j].items():
                    mm=(m1[0]+m2[0],m1[1]+m2[1],m1[2]+m2[2])
                    cur[mm]=cur.get(mm,0)+pref*c1*c2
        Es.append(cur)
    memo={(0,0,0):Fr(1)}
    def mom(b):
        if b in memo: return memo[b]
        if sum(b)%2: return Fr(0)
        i=next(j for j in range(3) if b[j])
        bb=list(b);bb[i]-=1
        tot=Fr(0)
        for d in range(3):
            if bb[d]==0: continue
            n=bb[d]; b2=list(bb);b2[d]-=1
            tot+=n*G[i][d]*mom(tuple(b2))
        memo[b]=tot
        return tot
    out=[]
    for n in range(N+1):
        tot=Fr(0)
        for m,cf in Es[2*n].items():
            tot+=cf*mom(m)
        out.append(tot)
    return out
if __name__=="__main__":
    N=int(sys.argv[1]) if len(sys.argv)>1 else 4
    # Test 1
    for (e1,e2,Y) in [(Fr(3),Fr(7),Fr(-2,3)),(Fr(2),Fr(-5,3),Fr(5,7))]:
        E=e1+e2; t0=time.time()
        P=Phi([E/3]*3,Y,N)
        v=1/(Y+E/3)
        rhs=mul(Theta(-1/E,N),S(v,E,N),N)
        ctl=mul(Theta(1/E,N),S(v,E,N),N)
        print("T_F eps",e1,e2,"Y",Y,"Phi==Theta(-z/E)S:",[P[n]==rhs[n] for n in range(N+1)],"control Theta(+):",[P[n]==ctl[n] for n in range(N+1)],"%.1fs"%(time.time()-t0),flush=True)
    # Test 2
    al=[Fr(1),Fr(2),Fr(5)]
    for i in range(3):
        Pi=Phi(al,-al[i],min(N,4))
        a=[3*al[i]]+[al[j]-al[i] for j in range(3) if j!=i]
        M=min(N,4)
        r=[Fr(1)]+[Fr(0)]*M
        for x in a: r=mul(r,Theta(-1/x,M),M)
        print("y=0 colour",i,"Phi==prod Theta(-z/a):",[Pi[n]==r[n] for n in range(M+1)],flush=True)
