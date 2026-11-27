"""Referee checks for RIGIDITY_PROOF.md (exact rationals, own code, no project code).
Sections: A s.1 identities mod P(U); B Lemma 5.3 root; C Lemma 5.5 limits (s-series);
D Lemma 5.6 Li; E power sums / frame constants; F (0.1), P_g, gap from Z_2;
G fresh GP localisation at d=1, g=1,2 (quintuple Hodge integrals) vs P_g on T_CY and BS display.
"""
from fractions import Fraction as Fr
import random, sys
random.seed(20261103)

def rnd():
    while True:
        v = Fr(random.randint(-9, 9), random.randint(1, 7))
        if v != 0:
            return v

def weights():
    while True:
        a0, a1, a2, e1 = rnd(), rnd(), rnd(), rnd()
        E = a0 + a1 + a2
        e2 = E - e1
        A, B = a0 + a2, a0 + a1
        if 0 not in (a1, a2, A, B, e1, e2, E, a1 - a2):
            return a0, a1, a2, e1, e2, A, B, E

# ---------- quadratic algebra Kw(Q)[U]/(P(U)) : elements (c0, c1) = c0 + c1 U ----------
class QA:
    def __init__(s, p2, p1, p0):  # P(U) = p2 U^2 + p1 U + p0
        s.p2, s.p1, s.p0 = p2, p1, p0
    def mul(s, x, y):
        c0 = x[0]*y[0]; c1 = x[0]*y[1] + x[1]*y[0]; c2 = x[1]*y[1]
        # U^2 = -(p1 U + p0)/p2
        return (c0 - c2*s.p0/s.p2, c1 - c2*s.p1/s.p2)
    def add(s, x, y): return (x[0]+y[0], x[1]+y[1])
    def sc(s, k, x): return (k*x[0], k*x[1])
    def inv(s, x):
        # conjugate: U -> Ssum - U with Ssum = -p1/p2 ; norm = x * xbar in Kw
        Ssum = -s.p1/s.p2
        xb = (x[0] + x[1]*Ssum, -x[1])
        n = s.mul(x, xb)
        assert n[1] == 0
        return s.sc(1/n[0], xb)
    def norm(s, x):
        Ssum = -s.p1/s.p2
        xb = (x[0] + x[1]*Ssum, -x[1])
        n = s.mul(x, xb); assert n[1] == 0
        return n[0]

def secA():
    ok = True
    for _ in range(25):
        a0, a1, a2, e1, e2, A, B, E = weights()
        Q = rnd()
        if Q == 1: continue
        tau = 1 - 1/Q
        qa = QA(B*Q, -(a0 + Q*(a1 - a2)), -a2)
        U = (Fr(0), Fr(1)); one = (Fr(1), Fr(0))
        # v = U/(1-U), w = QU/(1-QU)
        v = qa.mul(U, qa.inv(qa.add(one, qa.sc(-1, U))))
        w = qa.mul(qa.sc(Q, U), qa.inv(qa.add(one, qa.sc(-Q, U))))
        # W_1 = -a2 - A v - a1 w ; critical: must be 0
        W1 = qa.add(qa.add(qa.sc(-a2, one), qa.sc(-A, v)), qa.sc(-a1, w))
        ok &= (W1 == (0, 0))
        # L = -a1 w ; check tau L^2 - (E - a0 tau) L + a1 a2 = 0
        L = qa.sc(-a1, w)
        L2 = qa.mul(L, L)
        chk = qa.add(qa.add(qa.sc(tau, L2), qa.sc(-(E - a0*tau), L)), qa.sc(a1*a2, one))
        ok &= (chk == (0, 0))
        # W_2 = -A v(1+v) - a1 w(1+w)   [theta v = v + v^2]
        W2 = qa.add(qa.sc(-A, qa.mul(v, qa.add(one, v))), qa.sc(-a1, qa.mul(w, qa.add(one, w))))
        Delta_a = qa.sc(a1*A*e1*e2, W2)
        Delta_b = qa.sc(e1*e2, qa.add(qa.add(qa.sc(a0*a1*a2, one), qa.sc(2*a1*a2, L)), qa.sc(-E, L2)))
        ok &= (Delta_a == Delta_b)
        # (a1a2 - E L) norm = a1 a2 A B
        ok &= (qa.norm(qa.add(qa.sc(a1*a2, one), qa.sc(-E, L))) == a1*a2*A*B)
        # Delta_i = 1/eta(e_i,e_i), eta = (1/(ABe1e2))[[E/(a1a2),1],[1,-a0]], e_i=(H-L_j)/(L_i-L_j), L_j = E/tau - a0 - L_i
        Lj = qa.add(qa.sc(E/tau - a0, one), qa.sc(-1, L))
        dL = qa.add(L, qa.sc(-1, Lj))
        etaHH = -a0/(A*B*e1*e2); eta1H = 1/(A*B*e1*e2); eta11 = E/(a1*a2*A*B*e1*e2)
        num = qa.add(qa.add(qa.sc(etaHH, one), qa.sc(-2*eta1H, Lj)), qa.sc(eta11, qa.mul(Lj, Lj)))
        etaee = qa.mul(num, qa.inv(qa.mul(dL, dL)))
        ok &= (qa.mul(etaee, Delta_b) == one)
        # Disc identities
        Disc = a0**2 + (2*a0*(a1+a2) + 4*a1*a2)*Q + (a1-a2)**2*Q**2
        ok &= (Disc == (a0 + Q*(a1-a2))**2 + 4*a2*B*Q)
        ok &= (Disc == Q**2*((E - a0*tau)**2 - 4*a1*a2*tau))
        ok &= ((2*a0*(a1+a2) + 4*a1*a2)**2 - 4*(a1-a2)**2*a0**2 == 16*a1*a2*A*B)
        # P(1), P(1/Q)
        P = lambda u: B*Q*u*u - (a0 + Q*(a1-a2))*u - a2
        ok &= (P(Fr(1)) == A*(Q-1)) and (P(1/Q) == a1*(1-Q)/Q)
        # eta^{-1}
        det = eta11*etaHH - eta1H**2
        inv = [[etaHH/det, -eta1H/det], [-eta1H/det, eta11/det]]
        ok &= (inv == [[e1*e2*a0*a1*a2, e1*e2*a1*a2], [e1*e2*a1*a2, -e1*e2*E]])
    print("A  s.1 identities (W_1=0, L quadratic, (1.1), (1.4), norm, Delta=1/eta(e,e), Disc, (1.3), eta^-1):", ok)
    return ok

# ---------- B: Lemma 5.3, root Q_1 in K[[a0]] ----------
def ser_mul(x, y, N):
    r = [Fr(0)]*N
    for i, xi in enumerate(x):
        if xi == 0: continue
        for j, yj in enumerate(y):
            if i + j < N: r[i+j] += xi*yj
    return r

def secB():
    ok = True
    N = 11
    for _ in range(5):
        a1, a2 = rnd(), rnd()
        while a1 == a2: a2 = rnd()
        # Disc(Q) = a0^2 + (2a0(a1+a2) + 4a1a2) Q + (a1-a2)^2 Q^2, as series in a0 with Q = series
        c1 = [4*a1*a2, 2*(a1+a2)] + [Fr(0)]*(N-2)
        c0 = [Fr(0), Fr(0), Fr(1)] + [Fr(0)]*(N-3)
        c2 = [(a1-a2)**2] + [Fr(0)]*(N-1)
        Q = [Fr(0)]*N  # Newton: Q <- Q - Disc(Q)/Disc'(Q) in K[[a0]]
        for it in range(8):
            D = [c0[i] + 0 for i in range(N)]
            D = [D[i] + ser_mul(c1, Q, N)[i] + ser_mul(c2, ser_mul(Q, Q, N), N)[i] for i in range(N)]
            Dp = [c1[i] + 2*ser_mul(c2, Q, N)[i] for i in range(N)]
            # invert Dp (unit, Dp[0] = 4 a1 a2)
            inv = [Fr(0)]*N; inv[0] = 1/Dp[0]
            for n in range(1, N):
                inv[n] = -sum(Dp[k]*inv[n-k] for k in range(1, n+1))/Dp[0]
            corr = ser_mul(D, inv, N)
            Q = [Q[i] - corr[i] for i in range(N)]
        D = [c0[i] + ser_mul(c1, Q, N)[i] + ser_mul(c2, ser_mul(Q, Q, N), N)[i] for i in range(N)]
        ok &= all(d == 0 for d in D)
        ok &= (Q[0] == 0 and Q[1] == 0 and Q[2] == -1/(4*a1*a2))
        ok &= (Q[3] == (a1+a2)/(8*a1**2*a2**2))
        ok &= (Q[4] == -(5*a1**2 + 6*a1*a2 + 5*a2**2)/(64*a1**3*a2**3))
    print("B  Lemma 5.3: Disc(Q_1)=0 to O(a0^11), val(Q_1)=2, c2,c3,c4 as stated:", ok)
    return ok

# ---------- C: Lemma 5.5, s = 1/Q series ----------
def secC():
    ok = True
    N = 8
    for _ in range(6):
        a0, a1, a2, e1, e2, A, B, E = weights()
        # U_1(s): B U^2 - (a0 s + a1 - a2) U - a2 s = 0 ; Newton in K[[s]]
        lin = [a1 - a2, a0] + [Fr(0)]*(N-2)
        U = [(a1-a2)/B] + [Fr(0)]*(N-1)
        for it in range(6):
            F = [B*ser_mul(U, U, N)[i] - ser_mul(lin, U, N)[i] - (a2 if i == 1 else 0) for i in range(N)]
            Fp = [2*B*U[i] - lin[i] for i in range(N)]
            inv = [Fr(0)]*N; inv[0] = 1/Fp[0]
            for n in range(1, N): inv[n] = -sum(Fp[k]*inv[n-k] for k in range(1, n+1))/Fp[0]
            corr = ser_mul(F, inv, N); U = [U[i] - corr[i] for i in range(N)]
        F = [B*ser_mul(U, U, N)[i] - ser_mul(lin, U, N)[i] - (a2 if i == 1 else 0) for i in range(N)]
        ok &= all(f == 0 for f in F) and U[0] == (a1-a2)/B
        # V_2(s): B s V^2 - (a0 s + a1 - a2) V - a2 = 0
        V = [-a2/(a1-a2)] + [Fr(0)]*(N-1)
        for it in range(6):
            VV = ser_mul(V, V, N)
            F = [B*(VV[i-1] if i >= 1 else 0) - ser_mul(lin, V, N)[i] - (a2 if i == 0 else 0) for i in range(N)]
            Fp = [2*B*(V[i-1] if i >= 1 else 0) - lin[i] for i in range(N)]
            inv = [Fr(0)]*N; inv[0] = 1/Fp[0]
            for n in range(1, N): inv[n] = -sum(Fp[k]*inv[n-k] for k in range(1, n+1))/Fp[0]
            corr = ser_mul(F, inv, N); V = [V[i] - corr[i] for i in range(N)]
        VV = ser_mul(V, V, N)
        F = [B*(VV[i-1] if i >= 1 else 0) - ser_mul(lin, V, N)[i] - (a2 if i == 0 else 0) for i in range(N)]
        ok &= all(f == 0 for f in F) and V[0] == -a2/(a1-a2)
        # limits at s = 0 (values): colour 1: u = U/(1-U), w = QU/(1-QU) = -1/(1 - s/U) -> -1
        u1 = U[0]/(1-U[0]); w1 = Fr(-1)
        W2_1 = -A*u1*(1+u1) - a1*w1*(1+w1)
        L1 = -a1*w1
        ok &= (L1 == a1) and (W2_1 == -(a1-a2)*B/A)
        # colour 2: u = sV/(1-sV) -> 0 ; w = V/(1-V)
        u2 = Fr(0); w2 = V[0]/(1-V[0])
        W2_2 = -A*u2*(1+u2) - a1*w2*(1+w2)
        L2 = -a1*w2
        ok &= (L2 == a2) and (W2_2 == a2*(a1-a2)/a1)
        D1 = a1*A*e1*e2*W2_1; D2 = a1*A*e1*e2*W2_2
        ok &= (D1 == -a1*B*(a1-a2)*e1*e2) and (D2 == a2*A*(a1-a2)*e1*e2)
        ok &= (D1 == e1*e2*(L1-L2)*(a1*a2 - E*L1)) and (D2 == e1*e2*(L2-L1)*(a1*a2 - E*L2))
        # flop weights: e(wt'_-) = D1, e(wt'_+) = D2
        ok &= (D1 == (a2-a1)*a1*B*(-e1)*(-e2)) and (D2 == (a1-a2)*a2*A*(-e1)*(-e2))
        # 1 - U_1(0) = A/B, 1 - V_2(0) = a1/(a1-a2)
        ok &= (1 - U[0] == A/B) and (1 - V[0] == a1/(a1-a2))
    print("C  Lemma 5.5 / (1.5): U_1(s), V_2(s) in K[[s]], limits L_i, W_2, Delta_i, flop e(wt'):", ok)
    return ok

# ---------- D: Lemma 5.6 ----------
def polmul(p, q):
    r = [Fr(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i+j] += a*b
    return r
def secD():
    ok = True
    # A_n via recursion A_{m+1} = Q[(1-Q)A_m' + (m+1)A_m]
    A = [[Fr(0), Fr(1)]]  # A_0 = Q
    for m in range(0, 8):
        Am = A[-1]
        dA = [i*Am[i] for i in range(1, len(Am))] + [Fr(0)]
        t1 = polmul([Fr(1), Fr(-1)], dA)
        t2 = [(m+1)*c for c in Am] + [Fr(0)]*(len(t1)-len(Am))
        s = [t1[i] + t2[i] for i in range(len(t1))]
        A.append([Fr(0)] + s)
    for n in range(1, 8):
        An = A[n]
        while An[-1] == 0: An = An[:-1]
        ok &= (len(An) - 1 == n) and (An[0] == 0) and (An[-1] == 1)
    ok &= (A[1][:2] == [0, 1]) and (A[2][:3] == [0, 1, 1]) and (A[3][:4] == [0, 1, 4, 1]) and (A[4][:5] == [0, 1, 11, 11, 1])
    def Li(n, Q):
        num = sum(c*Q**i for i, c in enumerate(A[n])); return num/(1-Q)**(n+1)
    for n in range(1, 8):
        for _ in range(5):
            Q = rnd()
            if Q == 1: continue
            ok &= (Li(n, 1/Q) == (-1)**(n+1)*Li(n, Q))
            # series check
            S = sum(Fr(d)**n * Fr(1, 3)**d for d in range(1, 120))
            ok &= abs(Li(n, Fr(1, 3)) - S) < Fr(1, 10**20)
    ok &= (Li(0, Fr(5)) == -1 - Li(0, Fr(1, 5)))
    # expansion of Li_{-n}(e^t) at t=0: polar part (-1)^{n+1} n! t^{-n-1}
    # e^t/(1-e^t) = -1/t - 1/2 - sum B_{2k} t^{2k-1}/(2k)!  (Bernoulli)
    Bern = [Fr(1), Fr(-1, 2), Fr(1, 6), Fr(0), Fr(-1, 30), Fr(0), Fr(1, 42), Fr(0), Fr(-1, 30), Fr(0), Fr(5, 66), Fr(0), Fr(-691, 2730), Fr(0), Fr(7, 6)]
    from math import factorial
    # Laurent coefficients of Li_0(e^t) = -1/(1-e^{-t}) ... use e^t/(1-e^t) = -1 - 1/(e^t - 1); 1/(e^t-1) = sum B_k t^{k-1}/k!
    lo = {}
    for k in range(0, 15):
        lo[k-1] = lo.get(k-1, Fr(0)) - Bern[k]/factorial(k)
    lo[0] -= 1
    def deriv(f):
        return {k-1: k*v for k, v in f.items() if k != 0}
    f = dict(lo)
    for n in range(1, 7):
        f = deriv(f)
        neg = {k: v for k, v in f.items() if k < 0 and v != 0}
        ok &= (neg == {-n-1: Fr((-1)**(n+1)*factorial(n))})
        if n == 1: ok &= (f[0] == Fr(-1, 12))
        if n == 3: ok &= (f[0] == Fr(1, 120))
    print("D  Lemma 5.6: A_n, degrees, Li(1/Q) symmetry, series, polar part (-1)^(n+1) n! t^(-n-1), constants:", ok)
    return ok

# ---------- E: frame constants (power sums) ----------
def secE():
    ok = True
    for _ in range(10):
        a0, a1, a2, e1, e2, A, B, E = weights()
        for k in (1, 2, 3):
            ps = lambda v: sum(Fr(x)**(1-2*k) for x in v)
            nw = (a1, A, -e1, -e2)
            ok &= (ps(nw) + ps((E, -a1, -A)) == ps((E, -e1, -e2)))
            ok &= (ps((a2-a1, -A, B)) + ps(nw) == ps((a2-a1, a1, B, -e1, -e2)))
            ok &= (ps((a2, -a1, a1-a2)) + ps(nw) == ps((a1-a2, a2, A, -e1, -e2)))
        # a0 = 0 case
        a0z = Fr(0); Ez = a1 + a2; Az = a2; e2z = Ez - e1
        for k in (1, 2, 3):
            ps = lambda v: sum(Fr(x)**(1-2*k) for x in v)
            ok &= (ps((a1, Az, -e1, -e2z)) - ps((Ez, -e1, -e2z)) == ps((a1, a2, -Ez)))
    print("E  power sums: h_nw h_(E,-a1,-A) = h_w ; flop constants ; a0=0: h_nw h_w^-1 = h_(a1,a2,-E):", ok)
    return ok

# ---------- F: (0.1), P_g, gap_2 from Z_2 ----------
def secF():
    ok = True
    # S(y) = y/(2 sinh(y/2)) = 1/sum_{k} y^{2k}/(4^k (2k+1)!)
    from math import factorial
    N = 5
    den = [Fr(1, 4**k*factorial(2*k+1)) for k in range(N)]  # in y^2
    s = [Fr(0)]*N; s[0] = Fr(1)
    for n in range(1, N): s[n] = -sum(den[k]*s[n-k] for k in range(1, n+1))
    ok &= (s[:4] == [1, Fr(-1, 24), Fr(7, 5760), Fr(-31, 967680)])
    def P(g, e1, e2): return sum(s[i]*s[g-i]*e1**(2*i)*e2**(2*(g-i)) for i in range(g+1))
    for _ in range(5):
        e1, e2 = rnd(), rnd()
        ok &= (P(2, e1, e2) == (7*e1**4 + 10*e1**2*e2**2 + 7*e2**4)/5760)
        ok &= (6*P(3, e1, e2)/(e1*e2) == -(31*e1**6 + 49*e1**4*e2**2 + 49*e1**2*e2**4 + 31*e2**6)/(161280*e1*e2))
        # beta_n = (n-1)! [t^n] 1/(4 sinh sinh) = (n-1)! P_{(n+2)/2}/(e1 e2)
        # Theorem VI at g=2 via Z_2 = -(24 + 28 sg + 7 sg^2)/5760, sg = 1/sigma, sigma = -e1 e2/E^2, Z[2,0] = E^-3 Z_2
        E = e1 + e2
        if E == 0: continue
        sigma = -e1*e2/E**2; sg = 1/sigma
        Z2 = -(24 + 28*sg + 7*sg**2)/5760
        Z20 = Z2/E**3
        ok &= ((-e1*e2*E**3)*Z20 == P(2, e1, e2)/(e1*e2))
    # BS display degree 4 on T_CY equals P_2
    for _ in range(5):
        e1, e2 = rnd(), rnd()
        ep, em = (e1+e2)/2, (e1-e2)/2
        ok &= ((3*ep**4 + 8*ep**2*em**2 + 3*em**4)/720 == P(2, e1, e2))
        ok &= (-(ep**2 + em**2)/12 == P(1, e1, e2))
    print("F  S(y) coefficients, P_2, P_3, beta_2, beta_4, gap_2 from Z_2 (Theorem VI at g=2), BS display vs P_g:", ok)
    return ok

# ---------- G: fresh GP localisation, d = 1 ----------
# e1 e2 GW_{g,1} = sum_{g+ + g- = g} V_+(g+) V_-(g-), V_p(0)=1,
# V_p(h) = Int_{Mbar_{h,1}} prod_a Lambda^vee_h(w_a) / (w_1 (w_1 - psi)),  w = wt_p, w_1 = tangent weight of P^1.
# Hodge integrals on Mbar_{1,1}: psi = lambda_1 = 1/24.
# On Mbar_{2,1} (i = power of lambda_1, j = power of lambda_2, k = power of psi):
H2 = {(0, 0, 4): Fr(1, 1152), (1, 0, 3): Fr(1, 480), (0, 1, 2): Fr(7, 5760), (2, 0, 2): Fr(7, 2880),
      (1, 1, 1): Fr(1, 2880), (3, 0, 1): Fr(1, 1440), (0, 2, 0): Fr(0), (2, 1, 0): Fr(0), (4, 0, 0): Fr(0)}

def V1(w):
    # genus 1: prod (w_a - lambda_1) = prod w - lambda_1 * e4(w); 1/(w1 (w1 - psi)) = (1/w1^2)(1 + psi/w1 + ...)
    pw = Fr(1)
    for x in w: pw *= x
    e4 = sum(pw/x for x in w)
    w1 = w[0]
    return (pw*Fr(1, 24)/w1 - e4*Fr(1, 24))/w1**2

def V2(w):
    # prod_a (w_a^2 - lambda_1 w_a + lambda_2): collect coefficients of lambda_1^i lambda_2^j
    poly = {(0, 0): Fr(1)}
    for x in w:
        new = {}
        for (i, j), c in poly.items():
            for (di, dj, f) in ((0, 0, x*x), (1, 0, -x), (0, 1, Fr(1))):
                new[(i+di, j+dj)] = new.get((i+di, j+dj), Fr(0)) + c*f
        poly = new
    w1 = w[0]
    tot = Fr(0)
    for (i, j), c in poly.items():
        k = 4 - i - 2*j
        if k < 0: continue
        tot += c*H2[(i, j, k)]/w1**(k+2)   # 1/(w1 (w1 - psi)) = sum psi^k / w1^(k+2)
    return tot

def secG():
    ok_cy = True; ok_bs = True; fail_off = True
    for _ in range(8):
        a0, a1, a2, e1, e2, A, B, E = weights()
        wp = (a0, a1, a2, -e1, -e2); wm = (-a0, B, A, -e1, -e2)
        g1 = V1(wp) + V1(wm)
        g2 = V2(wp) + V1(wp)*V1(wm) + V2(wm)
        s = [Fr(1), Fr(-1, 24), Fr(7, 5760)]
        P1 = -(e1**2 + e2**2)/24; P2 = (7*e1**4 + 10*e1**2*e2**2 + 7*e2**4)/5760
        ok_cy &= (g1 == P1) and (g2 == P2)
    # off-CY: general 5 weights, compare with BS display (line ~1694) and check g2 != P2 generically
    for _ in range(8):
        a0, a1, a2, e1, e2 = rnd(), rnd(), rnd(), rnd(), rnd()
        A, B = a0 + a2, a0 + a1
        if 0 in (a1, a2, A, B, e1, e2, a0): continue
        wp = (a0, a1, a2, -e1, -e2); wm = (-a0, B, A, -e1, -e2)
        g1 = V1(wp) + V1(wm); g2 = V2(wp) + V1(wp)*V1(wm) + V2(wm)
        ep, em = (e1+e2)/2, (e1-e2)/2; c1K = e1 + e2 - a0 - a1 - a2
        bs1 = -(ep**2 + em**2)/12 + c1K*ep/12
        bs2 = (3*ep**4 + 8*ep**2*em**2 + 3*em**4)/720 + c1K*ep*(2*a1*a2 + (a1+a2)*a0 + a0**2 - 7*(a0+a1+a2)*ep - 14*em**2)/1440
        ok_bs &= (g1 == bs1) and (g2 == bs2)
        P2 = (7*e1**4 + 10*e1**2*e2**2 + 7*e2**4)/5760
        if c1K != 0: fail_off &= (g2 != P2)
        # Lemma 5.7 structure at d=1: g1 = e1 e2/12 - (e1+e2)(a0+a1+a2)/24
        ok_bs &= (g1 == e1*e2/12 - (e1+e2)*(a0+a1+a2)/24)
    print("G  fresh d=1 localisation: g=1,2 equal P_1, P_2 on T_CY:", ok_cy,
          "| equal BS line-1694 display and Lemma 5.7 shape on T_max:", ok_bs, "| g=2 differs from P_2 off T_CY:", fail_off)
    return ok_cy and ok_bs and fail_off

if __name__ == "__main__":
    res = [secA(), secB(), secC(), secD(), secE(), secF(), secG()]
    print("ALL:", all(res))
