# Part 1 over Q at arbitrary rational samples: the sqrt(W_j) enter M only through conjugation by the constant matrix diag(s_j).
# With R=rR D^-1, Psi=rPsi D^-1 (D=diag(s_j), theta D=0): M = D M_rat D^-1, M_rat computed from (rR,rPsi) over Q.  Also F_p for p=1 mod 4, Tonelli-Shanks.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from rows_p2 import *
from fractions import Fraction as Fr
import part1, cl_lib

def tonelli(x):
    p=type(x).p
    F=type(x)
    if x.v==0: return x
    if pow(x.v,(p-1)//2,p)!=1: return None
    q=p-1;s=0
    while q%2==0: q//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1: z+=1
    m=s;c=pow(z,q,p);t=pow(x.v,q,p);r=pow(x.v,(q+1)//2,p)
    while t!=1:
        i=0;tt=t
        while tt!=1:
            tt=tt*tt%p;i+=1
        b=pow(c,1<<(m-i-1),p)
        m=i;c=b*b%p;t=t*c%p;r=r*b%p
    return F(r)

def run_rat(alpha,eps,D,K):
    F=Fr
    P=P2Data(F,alpha,eps,D,K)
    t=P.tt
    rr=P.rows_rational()
    Psi=[[rr[a][j].zcoef(0) for j in range(3)] for a in range(3)]
    etai=mconst(t,P.etai_f);eta=P.eta;Hs=P.Hs
    Psii=minv(Psi)
    Rend=mmul(mmul(etai,rr),mmul(Psii,eta))
    B=minv(Rend)
    M,z0=defect_endo(B,Psi,eta,Hs,P.thetaQ)
    out={'M_rat=0':miszero(M),'comm no psi^0':z0}
    for nm,BB,HH in (('R^*',mmap(lambda a:a.zneg(),Rend),Hs),('R(-z)',mmap(lambda a:a.zneg(),B),Hs),('H*->-H*',B,mscal(Hs,F(-1)))):
        Mx,z=defect_endo(BB,Psi,eta,HH,P.thetaQ)
        nz=mzero_report(Mx)
        out[nm]='M_rat=0' if not nz else 'M_rat!=0, psi-orders %s, min y %d'%(sorted(set(e[3] for e in nz)),min(e[2] for e in nz))
    print('Q, sample',tuple(map(str,alpha)),tuple(map(str,eps)),'D',D,'K',K,out,flush=True)

for al,ep in [((1,2,5),(3,7)),((2,-3,8),(1,4)),((1,3,8),(Fr(5,2),Fr(-3,7))),((3,7,Fr(-11,2)),(Fr(5,3),Fr(-2,7))),((2,-3,7),(1,4)),((Fr(1,3),Fr(5,7),Fr(-2,9)),(Fr(11,4),Fr(13,6)))]:
    run_rat(al,ep,3,5)
run_rat((1,2,5),(3,7),5,7)
# sample that needs p=1 mod 4
cl_lib.fp_sqrt=tonelli
part1.fsqrt=lambda x:(rat_sqrt(x) if isinstance(x,Fr) else tonelli(x))
n=2**89+1
cnt=0
while cnt<2:
    n+=4
    if not is_prime(n): continue
    F=make_fp(n)
    al,ep=(1,3,8),(Fr(5,2),Fr(-3,7))
    P=P2Data(F,al,ep,1,1)
    if any(tonelli(P.W0(i)) is None for i in range(3)): continue
    cnt+=1
    part1.run(F,al,ep,3,5,'F_p p=1 mod 4, p=%d'%n)
