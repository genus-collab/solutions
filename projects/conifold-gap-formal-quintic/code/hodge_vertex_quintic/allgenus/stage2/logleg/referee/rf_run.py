"""Referee run: independent TR omegas + falsifier (b) + (V1) step + controls.
usage: py -I -B rf_run.py PRIME_EXP KM g,n [g,n ...]      (PRIME = 2^PRIME_EXP - 1)"""
import sys, time
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
from fractions import Fraction
from rf_rec import Curve
from rf_build import build, closure
from rf_pair import leg_pow, leg_log, leg_log_ebasis, contract, Fderiv, modq, ratrec
from rf_tr import inv

def parts(a, smax, mmin=1):
    """non-decreasing tuples of a positive ints with sum <= smax"""
    if a == 0:
        yield (); return
    for m in range(mmin, smax // a + 1):
        for rest in parts(a - 1, smax - m, m):
            yield (m,) + rest

def main():
    e2, KM = int(sys.argv[1]), int(sys.argv[2])
    targets = [tuple(int(x) for x in s.split(',')) for s in sys.argv[3:]]
    P = 2 ** e2 - 1
    T0 = time.time()
    print(f"rf_run: P = 2^{e2}-1, KM = {KM}, targets {targets}")
    cv = Curve(P, KM)
    print(f"  curve set-up {time.time()-T0:.1f} s")
    LL, LE = leg_log(cv), leg_log_ebasis(cv)
    print(f"  [{'PASS' if LL == LE else 'FAIL'}] log functional: zeta-route == e-basis route Int_0^oo dp/(p-c)^(m+1)")
    order = closure(targets)
    W = build(cv, order)
    Lp = {}
    def lp(m):
        if m not in Lp: Lp[m] = leg_pow(cv, m)
        return Lp[m]
    vals = {}
    nsing = nbad = nmix = 0; nsel = 0
    for (h, n) in order:
        chi = 2 * h - 2 + n
        T = W[(h, n)]
        Tk = [T]
        for k in range(1, n + 1):
            Tk.append(contract(Tk[-1], LL, P))
        for k in range(n + 1):
            a = n - k
            if a == 0:
                val = Tk[n].get((), 0)
                vals[(h, 0, k, ())] = val
                print(f"  [{'PASS' if val == modq(Fderiv(h, n), P) else 'FAIL'}] ({h},{n}) <(-log p)^{n}> = {ratrec(val, P)}"
                      f"  Gap/F_h^({n})(1) = {Fderiv(h, n)}  frozen (-5)^(1-h) x = {Fraction(-5) ** (1 - h) * Fderiv(h, n)}")
                continue
            smax = 5 * chi + 5
            ns = nz = 0
            cache = {(): Tk[k]}
            for ms in parts(a, smax):
                # contract L-legs one at a time with prefix caching
                pref = ()
                for m in ms:
                    nxt = pref + (m,)
                    if nxt not in cache:
                        cache[nxt] = contract(cache[pref], lp(m), P)
                    pref = nxt
                val = cache[ms].get((), 0)
                vals[(h, a, k, ms)] = val
                if sum(ms) < 5 * chi:
                    ns += 1
                    if sum(ms) % 5: nsel += (val != 0)
                    if val: nz += 1; print(f"     NONZERO singular: ({h},{n}) L-legs {ms} log {k}: {ratrec(val, P)}")
            nsing += ns; nbad += nz; nmix += ns if k else 0
            print(f"  [{'PASS' if nz == 0 else 'FAIL'}] ({h},{n}) a={a} L-legs, k={k} log legs: {ns} singular multisets"
                  f" (Sum m < {5*chi}, all residues mod 5), non-zero {nz}")
    print(f"  TOTAL singular multisets {nsing} (mixed {nmix}), non-zero {nbad}, Z/5 selection violations {nsel}")
    return cv, W, vals, order

if __name__ == '__main__':
    main()
