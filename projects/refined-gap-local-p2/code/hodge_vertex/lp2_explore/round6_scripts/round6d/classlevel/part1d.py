# M_0 identity of A2.8 tested on a non-flat R: M_0=-W_Psi+[Lambda,R^f_1] with R(-z), R^*; skewness of M_0; sample values of class-level numbers.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cl_lib import *
from rows_p2 import *
from fractions import Fraction as Fr
p=618970019642690137449563171
F=make_fp(p)
P=P2Data(F,(1,2,5),(3,7),3,5)
t=P.tt
roots=[fsqrt(P.W0(i)) for i in range(3)]
rr=P.rows_rational()
R=[[rr[a][j]*(F(1)/roots[j]) for j in range(3)] for a in range(3)]
Psi=[[R[a][j].zcoef(0) for j in range(3)] for a in range(3)]
etai=mconst(t,P.etai_f);eta=P.eta;Hs=P.Hs
Psii=minv(Psi)
Rend=mmul(mmul(etai,R),mmul(Psii,eta))
B=minv(Rend)
WPsi=mmul(Psii,mmap(P.thetaQ,Psi))
Lam=[[ (P.Ainv*P.yb(P.Y[i],t) if i==j else const(t,0)) for j in range(3)] for i in range(3)]
def f_frame(X):  # X End in flat coordinates -> Psi^-1 eta X eta^-1 Psi
    return mmul(mmul(mmul(Psii,eta),X),mmul(etai,Psi))
R1f=f_frame(mmap(lambda a:a.zcoef(1),Rend))
print('R^f_1 symmetric:',miszero(msub(R1f,mT(R1f))))
print('W_Psi skew:',miszero(madd(WPsi,mT(WPsi))))
comm=msub(mmul(Lam,R1f),mmul(R1f,Lam))
print('W_Psi=[Lambda,R^f_1] (flat R):',miszero(msub(WPsi,comm)))
# non-flat: B' = R(-psi)^... : legs of R(-z). M_0 should be -W_Psi-[Lambda,R^f_1]=-2 W_Psi
Bneg=mmap(lambda a:a.zneg(),B)
Mx,_=defect_endo(Bneg,Psi,eta,Hs,P.thetaQ)
M0=mmap(lambda a:a.zcoef(0),Mx)
print('R(-z): M_0=-W_Psi+[Lambda,(R^f_1 of R(-z))]=-2W_Psi:',miszero(madd(M0,mscal(WPsi,F(2)))))
print('R(-z): M_0 skew:',miszero(madd(M0,mT(M0))))
