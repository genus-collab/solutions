import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from part2 import *
S=Setup((1,2,5),(3,7),2,7)
for g,ins in [(1,((1,1),)),(2,((2,2),)),(2,((4,1),)),(2,((2,1),(2,2))),(2,((3,0),(1,2))),(1,((1,1),(1,2)))]:
    print(g,ins,'GP N_d (d=0,1,2):',[str(x) for x in S.gw(g,ins)],'graph sum y^0..2:',[str(x) for x in S.R.integral(g,list(ins))],'GP in y:',[str(x) for x in S.gw_y(g,ins)])
