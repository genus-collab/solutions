# Critic: re-run the row_audit per-graph test at fresh and at special CY samples near Q = 1 (frame B).
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
from cylib import *
from gsum import *

def run(name, a0, a1, a2, e1, gs=(2,3,4), frame='B'):
    a0,a1,a2,e1 = Fr(a0),Fr(a1),Fr(a2),Fr(e1)
    e2 = a0+a1+a2-e1
    w = weights(a0,a1,a2,e1,e2)
    for g in gs:
        t0=time.time()
        S = Sum(w,g,frame,6)
        tot, rows = S.run()
        nreg=nsing=0; dec=[]
        for (i,gv,A,percol,gtot) in rows:
            he_ = sum(sum(r) for r in A)>0
            for col,t in percol:
                pol=[k for k in range(-2*g,0) if lcoef(t,k)]
                if he_:
                    if pol: nsing+=1
                    else: nreg+=1
                elif pol: dec.append((col,pol))
        pred,Pg = predicted(w,g,S.L)
        gap = lcoef(pred,2-2*g)
        vt = [c for c in dec if c[0]==(1,)]
        vok = len(dec)==1 and dec[0][1]==[2-2*g] if gap else len(dec)==0
        tot_ok = all(lcoef(tot,k)==lcoef(pred,k) for k in range(1-2*g,4) if k!=0)
        print(name,(str(a0),str(a1),str(a2),str(e1),str(e2)),'g',g,frame,'edge colourings regular',nreg,'singular',nsing,'vertex polar only t^(2-2g):',vok,'total ok:',tot_ok,'%.1fs'%(time.time()-t0),flush=True)

if __name__=='__main__':
    run('fresh1', '7/5','-3/8','13/6','11/9')
    run('fresh2', '-5/4','19/7','2/3','17/5')
    run('u=0 (a1=E/2)', '1','3','2','5/2')        # E=6, a1=E/2
    run('e1=e2', '3/2','-1/3','5/7','(3/2-1/3+5/7)/2'.replace('(3/2-1/3+5/7)/2', str((Fr(3,2)+Fr(-1,3)+Fr(5,7))/2)))
    run('a1=a2', '2/3','5/4','5/4','1/7')
    run('a1=a2,e1=e2', '2','3/5','3/5',str((2+Fr(6,5))/2))
    run('P_2=0 (3 a2=E)', '1','3','2','7/5')       # E=6, a2=2
    run('A=a1 (B=a2... a0+a2=a1)', '1','3','2','1/9')
    run('a0=0, a1=e2 only', '0','5/3','2/7','5/3')
    run('a0=-a1-a2+small', '-3','1','5/2','7')     # E=1/2, large e's: e2 negative
    run('a1=-a2 (a0=E)', '5','7/3','-7/3','2')
    run('a0=a1', '4/3','4/3','1/5','3')
