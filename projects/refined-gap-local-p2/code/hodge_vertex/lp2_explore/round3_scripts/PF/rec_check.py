# check recurrence 27(n+1)r_{n+1} = 3(3n^2+3n+1) r_n - n^3 r_{n-1} against Wick numbers
from fractions import Fraction as F
import subprocess,sys
L=int(sys.argv[1])
out=subprocess.run(['py','wick.py',str(L)],capture_output=True,text=True).stdout.split('\n')
w=[F(x.split()[1]) for x in out if x.strip()]
r=[F(1)]
for n in range(L):
    rm=r[n-1] if n>=1 else F(0)
    r.append((3*(3*n*n+3*n+1)*r[n]-n**3*rm)/(27*(n+1)))
print("all equal up to l=%d:"%L, all(a==b for a,b in zip(w,r)), len(w))
print(r[30] if L>=30 else r[-1])
