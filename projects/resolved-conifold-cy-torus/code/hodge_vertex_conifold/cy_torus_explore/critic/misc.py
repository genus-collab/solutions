import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cylib import *
from gsum import *
import gsum
import gen_q
# (1) non-CY control at generic Q: totals must not match the closed form
a0,a1,a2,e1,e2=Fr(7,5),Fr(-3,8),Fr(13,6),Fr(11,9),Fr(2)
w=weights(a0,a1,a2,e1,e2)
import random
rnd=random.Random(3)
while True:
    Q0=rnd.randrange(2,P)
    disc=(w.a0*w.a0+(2*w.a0*(w.a1+w.a2)+4*w.a1*w.a2)*Q0+(w.a1-w.a2)**2*Q0*Q0)%P
    if pow(disc,(P-1)//2,P)==1: break
gsum.build_rows=gen_q.rows_at_Q0(Q0)
for g in (2,3):
    S=Sum(w,g,'BA',4); tot,rows=S.run()
    cf=gen_q.closed_form(w,g,Q0,S.L); c0=constant_term(w,g)
    print('nonCY g',g,'closed form holds at t^0..t^3:',[lcoef(tot,k)==(cf[k]+(c0 if k==0 else 0))%P for k in range(4)])
# (2) alpha dependence of the constant map term at fixed (e1,e2)
e1,e2=Fr(3,5),Fr(17,4)
E=e1+e2
for g in (2,3,4):
    vals=[]
    for a1,a2 in ((Fr(1,3),Fr(2,7)),(Fr(5,9),Fr(-4,3)),(Fr(11,2),Fr(1,5)),(e1,e2)):
        a0=E-a1-a2
        if a0==0:
            vals.append(None); continue
        w=weights(a0,a1,a2,e1,e2)
        vals.append(constant_term(w,g))
    print('g',g,'F_g(0) at 3 alpha samples equal:',vals[0]==vals[1]==vals[2])
# (3) frame-B individual graph values at generic Q depend on alpha? g=2, compare two alpha at the same e, Q0
