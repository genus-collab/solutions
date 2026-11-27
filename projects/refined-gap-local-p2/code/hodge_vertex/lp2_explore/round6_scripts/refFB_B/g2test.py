import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from g3_ref import *
gl=graphs(2);print(len(gl))
gl3=graphs(3);print(len(gl3), sum(1 for x in gl3 if len(x[1])==6))
print(sorted(x[2] for x in gl))
def bs(e1,e2,y,F):
    X=1/(1+27*y);E=e1+e2;k=e1*e2
    f02=5*F**3/(8*X)+(5/(8*X)-Fr(1,2))*F**2+Fr(1,96)*(13*X+20/X-32)*F-(27*X**3-98*X**2+121*X-50)/(2160*X)
    f11=F**2/16+Fr(1,96)*(8-7*X)*F+(84*X**2-196*X+121)/5760
    f20=X/384*F-(21*X**2-29*X+8)/5760
    return E**4/k*f20+E**2*f11+k*f02
for e1,e2,t,F in [(Fr(3),Fr(7),Fr(1,2),Fr(1)),(Fr(1),Fr(1),Fr(-2,3),Fr(5,3)),(Fr(5,2),Fr(-3,7),Fr(3,11),Fr(-2,3))]:
    T=Total(e1,e2,t,F,5)
    v=T.total(2,gl)
    print(e1,e2,t,F,v.b==0,v.a==bs(e1,e2,t**3,F))
print("hodge g=2:",hodge_const(2,Fr(3),Fr(7)),Fr(100,640))
print("refined g=2,d=1:",refined(2,Fr(3),Fr(7),1)[1],"expected 48241/1680")
