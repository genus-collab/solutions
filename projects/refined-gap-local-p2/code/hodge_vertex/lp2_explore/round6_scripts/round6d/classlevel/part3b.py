import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from part3 import *
run(Fr,(Fr(-12),Fr(-6),Fr(3),Fr(1),Fr(216)),6,7,'Q (literal square sample), high order')
p=618970019642690137449562219
run(make_fp(p),(1,2,3,5,7),6,9,'F_p, high order, p=%d'%p)
