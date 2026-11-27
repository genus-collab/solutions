from fractions import Fraction as Fr
from alg import *
from core import *
prec=40
B=Base(Fr(1),prec); r=rec_r(prec+2)
Sinf=LS({},prec); vp=B.v0; fact=Fr(1)
for l in range(prec):
    fact*=(l+1); Sinf=Sinf+(vp*(r[l]/fact)); vp=vp*B.v0
print('S_inf=-t_c through t_c^%d:'%(prec-1),all(Sinf.get(k)==(-1 if k==1 else 0) for k in range(prec)))
rr=rec_r(200)
print('r_k nonzero k<=200:',all(x!=0 for x in rr),'signs first 20:',''.join('+' if x>0 else '-' for x in rr[:20]))
