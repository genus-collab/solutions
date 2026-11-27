# Test of thm:gau-R0 (i)-(ii) at first order in u:
#   Phi^can = Fs_0(tau zh) + u p_{1,4}(0) zh Fs_4(zh) + O(u^2),  tau = 1 - 4u/9 + ...
# so [u^1] ph_n = -(4/9) n [t^n]Fs_0 + p_{1,4}(0) [t^{n-1}]Fs_4  for every n, with ONE unknown p_{1,4}(0).
import importlib.util, sys, io, contextlib
from fractions import Fraction as Fr
spec=importlib.util.spec_from_file_location("chk","check_hpf.py")
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)
ph,Fs,N=m.ph,m.Fs,m.N
# n=1 fixes p14
p14=(ph[1].get(1,Fr(0))+Fr(4,9)*1*Fs[0][1])/Fs[4][0]
print("p_{1,4}(0) =",p14)
ok=[]
for n in range(1,N+1):
    lhs=ph[n].get(1,Fr(0))
    rhs=-Fr(4,9)*n*Fs[0][n]+p14*Fs[4][n-1]
    ok.append(lhs==rhs)
print("order-u identity for n=1..%d:"%N,ok)
print("Fs_4 =",Fs[4][:4])
# second check: mu(b,k) table
nu=[0,-1,3,2,1]
bad=[(b,k) for b in range(5) for k in range(0,40) if (4*k+nu[b]-nu[(b+k)%5])%5!=0 or (4*k+nu[b]-nu[(b+k)%5])<max(0,4*k-4)]
print("mu violations:",bad)
print("mu=0 cases:",[(b,k) for b in range(5) for k in range(0,10) if 4*k+nu[b]-nu[(b+k)%5]==0])
