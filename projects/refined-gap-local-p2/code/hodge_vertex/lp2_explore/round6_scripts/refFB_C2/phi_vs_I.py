# Stage B pressing point (1): the algebraic row Theta(z/e1)Theta(z/e2) Phi(alpha,Y;z) at the root Y_j(y) against
# the I-function row (J=I at P_j) times h_{wt_j}, generic alpha, tiny rational y. Evidence only.
from fractions import Fraction as Fr
from phi import Phi
from tf import theta_coeffs, mul_ser
from gen import Gen
al=(Fr(1,2),Fr(2),Fr(-5,3)); e1,e2=Fr(3),Fr(7)
ny=6; Z=4
G=Gen(al,e1,e2,ny,Z)
import sys
y=Fr(1,10**int(sys.argv[1]))
a=[Fr(x) for x in al]
e1a=sum(a)
for j in range(3):
    # Y_j series from the same recursion as Gen: recompute by Newton on the cubic at the rational y (exact to y^ny via iteration)
    Y=-a[j]
    for it in range(4):
        L=(Y+a[0])*(Y+a[1])*(Y+a[2])+27*y*Y**3
        dL=(Y+a[1])*(Y+a[2])+(Y+a[0])*(Y+a[2])+(Y+a[0])*(Y+a[1])+81*y*Y**2
        Y=Y-L/dL
        Y=Fr(round(Y*10**150),10**150)
    ph,det=Phi(al,Y,Z)
    pre=mul_ser(theta_coeffs(1/e1,Z),theta_coeffs(1/e2,Z),Z)
    alg=mul_ser(pre,ph,Z)
    # I-function row at y
    row=[sum(c*y**k for k,c in enumerate(G.rows[j][0][n])) for n in range(Z+1)]
    print("colour",j,"relative differences z^0..z^4:",[float((x-r)/(y**(ny-1))) for x,r in zip(alg,row)]," row z^1:",float(row[1]))
