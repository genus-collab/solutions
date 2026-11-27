from fractions import Fraction as Fr
import core, checkj, check0
from core import *
orig=core.Base
offs=lambda E,prec: orig(E,prec,fsec=[Fr(1),Fr(1,10)])
print('### colour j, off-shell section f=1+Delta/10')
checkj.Base=offs
checkj.run(Fr(-2,5),Fr(1),K=2,BW=3)
checkj.Base=orig
print('### colour j, r_2 perturbed')
rr=core.rec_r
checkj.rec_r=lambda N: rr(N,perturb={2:Fr(1,7)})
checkj.run(Fr(-2,5),Fr(1),K=3,BW=3)
checkj.rec_r=rr
print('### colour 0, off-shell section')
check0.Base=offs
check0.run(Fr(-2,5),Fr(1),K=2,BW=5,verbose=False)
check0.Base=orig
print('### colour 0, r_2 perturbed')
check0.rec_r=lambda N: rr(N,perturb={2:Fr(1,7)})
check0.run(Fr(-2,5),Fr(1),K=3,BW=5,verbose=False)
