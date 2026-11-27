from rows import *
c=setup((Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7)),2,6)
D=c.D
for i in range(3):
  for j in range(3):
    tot=[{} for _ in range(D+1)]
    for a in range(3):
      for b in range(3):
        Pjm=[{p:x*(-1)**p for p,x in c.Pt[j][b][d].items()} for d in range(D+1)]
        m=yzmul(c.Pt[i][a],Pjm,D,c.Zmax)
        for d in range(D+1):
          for p,x in m[d].items():
            tot[d][p]=tot[d].get(p,Fr(0))+c.etainv[a][b]*x
    print(i,j,[{p:str(x) for p,x in sorted(t.items()) if x!=0} for t in tot][:3])
print(c.etainv); print(c.eta)
print(c.r[0][1][0],c.r[0][2][0],c.r[0][1][1],c.r[0][2][1])
