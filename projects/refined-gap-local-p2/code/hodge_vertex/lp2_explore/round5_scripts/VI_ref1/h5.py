import os
import sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import faber
from fractions import Fraction as Fr
t=Fr(2,7); print('sigma Z_2 at 2/7 from beta', faber.beta(2,-t,1+t), 'Faber', t*(1+t)*faber.Z(2,t))
t=Fr(3,11); s=t*(1+t); z=faber.Z(5,t); print('h=5 tau=3/11', s**4*z==faber.beta(8,-t,1+t), s**4*z)
