# Referee's own formal Wick expansion of the 3-dim oscillating integral. Tests of Stage B's algebraic row:
# (a) Phi(alpha,Y=-alpha_i) = h_{(3a_i,a_j-a_i,a_k-a_i)} (BS Normalisation, up to eps factors),
# (b) Phi(E/3,Y) = Theta(-z/E) S(v,z),
# Evidence only, exact fractions.
import sys, math, itertools
from fractions import Fraction as Fr
from tf import theta_coeffs, mul_ser, rec_r
from hod import B

def inv3(M):
    a,b,c=M[0]; d,e,f=M[1]; g,h,i=M[2]
    det=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    return [[(e*i-f*h)/det,(c*h-b*i)/det,(b*f-c*e)/det],
            [(f*g-d*i)/det,(a*i-c*g)/det,(c*d-a*f)/det],
            [(d*h-e*g)/det,(b*g-a*h)/det,(a*e-b*d)/det]],det

def Tm(alpha,Y,m):
    """T_m as dict monomial->coef"""
    out={}
    fm=math.factorial(m)
    for j in range(3):
        e=[0,0,0]; e[j]=m
        out[tuple(e)]=out.get(tuple(e),Fr(0))+(alpha[j]+Y)/fm
    # -3Y (sum xi/3)^m/m!
    c=-3*Y/fm/Fr(3)**m
    for i in range(m+1):
        for j in range(m+1-i):
            k=m-i-j
            mult=math.factorial(m)//(math.factorial(i)*math.factorial(j)*math.factorial(k))
            out[(i,j,k)]=out.get((i,j,k),Fr(0))+c*mult
    return out

def pmul(a,b):
    o={}
    for x,u in a.items():
        for y,v in b.items():
            k=(x[0]+y[0],x[1]+y[1],x[2]+y[2])
            o[k]=o.get(k,0)+u*v
    return o

def Phi(alpha,Y,N):
    alpha=[Fr(x) for x in alpha]; Y=Fr(Y)
    H=[[ (alpha[i]+Y if i==j else Fr(0))-Y/3 for j in range(3)] for i in range(3)]
    Hi,det=inv3(H)
    G=[[-Hi[i][j] for j in range(3)] for i in range(3)]
    Ts={m:Tm(alpha,Y,m) for m in range(3,2*N+3)}
    P=[{(0,0,0):Fr(1)}]
    for k in range(1,2*N+1):
        cur={}
        for m in range(3,k+3):
            q=pmul(Ts[m],P[k-(m-2)])
            for x,u in q.items():
                cur[x]=cur.get(x,Fr(0))+u*(m-2)
        P.append({x:u/k for x,u in cur.items()})
    memo={}
    def mom(b):
        if b in memo: return memo[b]
        if sum(b)%2: return Fr(0)
        if sum(b)==0: return Fr(1)
        i=next(k for k in range(3) if b[k]>0)
        bb=list(b); bb[i]-=1
        tot=Fr(0)
        for j in range(3):
            if bb[j]>0 and G[i][j]!=0:
                c=list(bb); c[j]-=1
                tot+=G[i][j]*bb[j]*mom(tuple(c))
        memo[b]=tot
        return tot
    out=[]
    for n in range(N+1):
        out.append(sum(u*mom(x) for x,u in P[2*n].items()))
    return out,det

def h_of(wts,N):
    r=[Fr(1)]+[Fr(0)]*N
    for a in wts:
        r=mul_ser(r,theta_coeffs(-1/Fr(a),N),N)
    return r

def Sser(v,E,N):
    r=rec_r(N+2)
    out=[Fr(0)]*(N+1)
    for l in range(N+1):
        ser=[Fr(0)]*(N+1); ser[0]=Fr(1)
        for j in range(1,l+2):
            for n in range(1,N+1): ser[n]-=Fr(j)/E*ser[n-1]
        for k in range(l,N+1):
            out[k]+=r[l]*v**l*ser[k-l]
    return out

if __name__=="__main__":
    N=int(sys.argv[1]) if len(sys.argv)>1 else 3
    # (a)
    al=[Fr(2,3),Fr(5,7),Fr(-3,11)]
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        wt=[3*al[i],al[j]-al[i],al[k]-al[i]]
        ph,det=Phi(al,-al[i],N)
        hh=h_of(wt,N); hbad=h_of([-x for x in wt],N)
        print("(a) colour",i,"Phi-h_wt:",[str(a-b) for a,b in zip(ph,hh)]," control h_{-wt} differs:",any(a!=b for a,b in zip(ph,hbad)), " z^1 coeff:",ph[1])
    # (b)
    for (e1,e2,Y) in ((3,7,Fr(5,4)),(1,2,Fr(-1,5)),(Fr(5,2),Fr(-3,7),Fr(2,9))):
        E=Fr(e1)+Fr(e2)
        ph,det=Phi([E/3]*3,Y,N)
        v=3/(E+3*Y)
        th=theta_coeffs(-1/E,N)
        rhs=mul_ser(th,Sser(v,E,N),N)
        rhsbad=mul_ser(theta_coeffs(1/E,N),Sser(v,E,N),N)
        print("(b)",e1,e2,Y,"Phi - Theta(-z/E)S:",[str(a-b) for a,b in zip(ph,rhs)]," control Theta(+z/E) differs:",any(a!=b for a,b in zip(ph,rhsbad)))
