import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from barnes_check import *
th = Dd(1)/Dd(3)
def rep(lab, v, r):
    print(lab, "LHS re %s" % str(v[0])[:50])
    print("   target  %s" % str(r)[:50])
    print("   Re diff %.3E  Im %.3E" % (v[0]-r, v[1]), flush=True)
h = Dd(1)/Dd(120)
e = 5*th
v = lhs_sigma(2*th, 2*th, 2*th, 4*th, e, -Dd(1)/Dd(2), h, Dd(16))
rep("CStep6 (2/3,2/3,2/3,4/3,5/3), Re sigma=-1/2, h=1/120", v, (2*PI/Dd(3).sqrt())**3)
print("   formula RHS minus (2pi/sqrt3)^3: %.3E" % (rhs(2*th,2*th,2*th,4*th,e)-(2*PI/Dd(3).sqrt())**3), flush=True)
a,b,c,d = Dd("0.3"),Dd("0.45"),Dd("0.7"),Dd("0.2")
v = lhs_sigma(a,b,c,d,a+b+c-d+1, -Dd("0.15"), Dd(1)/Dd(80), Dd(16))
rep("generic (0.3,0.45,0.7,0.2), Re sigma=-0.15, h=1/80", v, rhs(a,b,c,d,a+b+c-d+1))
print("control G(1/2)^2-pi %.3E" % (rgamma(Dd("0.5"))**2-PI))
