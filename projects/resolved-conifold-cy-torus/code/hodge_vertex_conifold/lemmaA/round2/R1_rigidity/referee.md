# Lemma A round 2: referee R1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- T1 (closed form of the polar-order relations): proof re-derived step by step (s.2-s.7); independent KP propagation (refkp.pm) plus direct leaf-curve sum agree with the closed form, 0 mismatches. Exact: K = 5 (2/7), K = 8 (-5/8, 7/4), K = 12 (-5/8, 9/4). Mod 1e9+7: K = 10, 18, 22 at 12 framings including special ones. Mod 998244353: K = 20 (3/11, -13/7).
- The referee's propagation of all Y(T) at tau = 2/7, level <= 6 is identical to the route's kp_point.pl (E_0 formulation of Route A's Theorem D). The full system is consistent for arbitrary one-leaf data at generic tau (0 failures in more than 300000 overdetermination checks), so T1's linear form is canonical.
- Flow theorem / free energy on the slice (Theorem 6.1), checked on random paths: direct Sum_T Y(T)P^T/Aut T equals Res Psi log(W u/u(phi^-1)) to orders 11-15 at -5/8, 7/4, -2/5, -3/10, 5/2. phi_t(0) = 0, phi_t(1) = 1 and the inverse are verified. A perturbed formula fails.
- Branch-point identity d[u^-j]_+(1+tau) = j Pi_j and N(1+tau) = 0: re-derived by hand.
- Curve identity on the leaf curve: re-derived, including the residue identity via Lambda = e^(t mu) and the Pfaff transformation (DLMF 15.8.1). Numerically, Kflow - K_t is mu-independent with constant log H_1 - tau^2 t (orders <= 15). The residue identity holds for n <= 8, t^12, 4 framings (416/416). W = alpha holds with zeta built independently from Theorem 11.1's s*.
- Special-framing lemma (Lemma 8.1): re-derived; row spaces equal at all 176 framings with q <= 2G-2 for G = 3..8, using the referee's own two-leaf values.
- Forcing theorem (Theorem 9.1): proof re-derived. Independent rank test that uses neither T1 nor Lemma 8.1 (direct relation_k plus two-leaf regularity): all 83 tuples of the hypothesis range for G <= 11 are forced, with all one-leaf values y_D..y_{2G-3} vanishing at tau0. An exact BigRat case (G,D,q,k0) = (6,3,5,2) is also forced.
- Zero count (Proposition 10.1): Z1 brute force equals the formula for D <= 14 (Route A's Z(D)). Z2 by Moebius counting equals direct enumeration. The count for G <= 1000 fails only at (G,D) = (4..7, 3). Smallest margin for G >= 8 is 1, at (9,3). Bounds (B1) and (B2) re-derived by hand (constant 0.037602 > 0.0376) and never violated numerically.
- (OL) holds for every G >= 8, and Lemma A holds for every h given the cited inputs: the induction (Theorem A, even D by T1, odd D by forcing plus count, Corollary E, base cases h <= 4 data and h = 5, 6, 7 Theorem F) and the degree bound deg Q <= G-2 from brief 2(a),(b) were checked.
- Leaf-curve normalisation b_j = A_j alpha^-j/c^2 = (-1)^j (Pi_j^2/j) t^(j-1) H_j/H_1^j: c-scaling, string elimination and zeta cancellation confirmed. Agrees with Route D's T_k weights, its Lemma 5.2 cancellation and its residual -(1-x)^2 x^2/(12c^4) Xi^c(2).
- (I1) needs only genera < G: no stable graph of genus G with an edge has a genus-G vertex, confirmed by the dilaton-degree / tree argument.
- Global nullities reproduced by an independent implementation. Window G, TLR: 4, 7, 6, 1, 0 (G = 3..7). Window G-1, TLR: 2, 3, 0, 0 (G = 3..6). At G = 3, 4: TLR+CF 1, 2; full 1, 0; full+CF 0, 0.
- R1's remark that Route A's TLR nullities 11, 11, 7 included c = 2 is confirmed: the referee obtains 11, 11, 7 with c = 2 and 7, 6, 1 without.
- TR audit clean: no topological recursion, loop equations, BKMP, P2 or TR scripts in R1's work or in the referee's checks.

## issues

### [1]

#### description

The text says 'Kflow_t - K_t = log C'. With C := alpha e^(tau^2 t)/H_1, the constant is log(W/C) = log(alpha/C) = log H_1(T) - tau^2 t. Verified numerically (ref_flow.pl, orders <= 15). Only mu-independence is used, so T1 is unaffected.

#### location

report s.7, Corollary 7.2

#### severity

cosmetic

#### suggested_fix

Write Kflow_t - K_t = log(alpha/C) = log H_1(1-e^-t) - tau^2 t, as R1_checks.wls (C3) already does.

### [2]

#### description

W(t) = alpha(t) on the leaf curve is used without proof. It is needed for Pcal_t = Q_t(w) + alpha w and for the constant C.

#### location

report s.2 and s.7(a)

#### severity

minor

#### suggested_fix

Add: alpha = 1 + Sum_{j>=2} Pi_j P_j alpha^j is the sum rule Sum_{j>=1} (-1)^j (Pi_j^2/j) zeta^j H_j = -t (leaf_structure_proof.md Cor. 11.2(a)); the root W = 1 + O(t) is unique. Verified in ref_flow.pl with zeta built independently from s*.

### [3]

#### description

Typo: 'Xi_0(S+{1,n}) = (s+n)!/(s+n+1-m)' is missing a factorial.

#### location

report s.1, input (I4)

#### severity

cosmetic

#### suggested_fix

Xi_0(S+{1,n}) = (s+n)!/(s+n+1-m)!.

### [4]

#### description

The argument is complete, but it also rests on computer-verified base cases: Theorem F for h = 5, 6, 7 (rigorous modular ranks; Route E's engine also covers h = 5, 6, 7) and exact Hodge data for h = 3, 4. It is also conditional on (H) = Brini-Schuler through Route D. The plain label 'proved-modulo-literature' hides the computer-verified part.

#### location

claims 'Lemma A holds for every h' and Corollary R (status field)

#### severity

minor

#### suggested_fix

Label it 'proved modulo the cited inputs (Route A Thm D/Cor E/Thm A, Route D (I1) with (H), Theorem 11.1, brief 2(a),(b)) and the computer-verified base cases h <= 7'.

### [5]

#### description

(I1) needs only genera < G because a stable graph of genus G with an edge has no genus-G vertex. I checked this: dilaton leaves have psi-degree >= 2, so a genus-0 vertex with e edge half-edges and k leaves needs k <= e-3, hence e >= 3; a tree with one genus-G vertex and all other vertices of degree >= 3 contradicts the degree sum. The step is correct but is not spelled out anywhere in the route's text.

#### location

input (I1) as used in s.1 and s.11 (Route D Lemma 5.3 / brief 2(e))

#### severity

minor

#### suggested_fix

Add this two-line argument where (I1) is introduced.

### [6]

#### description

The raw output file begins with stray diff lines ('33,42d32 < G=31 ...').

#### location

out/count_zeros2_G3to1000_raw.txt

#### severity

cosmetic

#### suggested_fix

Regenerate the file from a clean run.

### [7]

#### description

Route A's TLR nullities 11, 11, 7 at G = 4, 5, 6 include c = 2 as an unknown, although Route A's TLR(G) statement sets v_2 = 0. The referee's ref_global_c.pl reproduces 11, 11, 7 with c = 2. Without c = 2 the nullities are 7, 6, 1 (window G) and 3, 0, 0 (window G-1, using the Hurwitz vanishing). R1's remark is correct.

#### location

Route A record (A_mv/report.md s.3.7), as remarked in R1 s.12

#### severity

minor

#### suggested_fix

Correct Route A's record. Optionally restate TLR in window G-1, where it already holds at G = 5, 6.


## notes

The referee's scripts are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\referee\: refkp.pm, ref_t1.pl, ref_dump.pl, ref_flow.pl, ref_res.pl, ref_special.pl, ref_forcing.pl, ref_zeros.pl, ref_global.pl and ref_global_c.pl. Outputs are in referee\out\.

R1_checks.wls was not run (no wolframscript). The main session may run it: its checks (C1)-(C6) duplicate checks I did independently in Perl.

Modular runs used p = 1000000007 and 998244353, with exact BigRat spot checks.

The pointwise-scan claim for G <= 12 was reproduced only in part: all forcing tuples for G <= 11, and the outside-range tuples for G <= 8.

## refuted


## route

R1: rigidity of the KP system (TLR and variants), round 2 of the attack on Lemma A. Claims checked: T1, flow theorem, branch-point identity, curve identity, special-framing lemma, forcing theorem, zero count, (OL) for G >= 8, Lemma A in all genera, global nullities

## verdict

sound

## report_text

# Referee report: Route R1, round 2 (rigidity of the KP system; T1, forcing theorem, Lemma A in all genera)

**Referee folder.** G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\referee\ holds the scripts; raw outputs are in referee\out\.

**Working rules followed.**
- No topological recursion of any kind.
- No wolframscript.
- No report file was written.
- All referee code was written from scratch. It does not copy the route's kpmod.pm or kp_point.pl, nor its E_0 construction.
- Arithmetic is exact (Math::BigRat) or modular, with p = 1000000007 and p = 998244353. The route used 2^31-1.

## 0. Verdict: SOUND

The main line holds, given the project's established inputs:
- T1 (closed form of the polar-order relations);
- the flow theorem and the free-energy formula on the slice;
- the branch-point identity;
- the curve identity on the leaf curve;
- the special-framing lemma (Lemma 8.1);
- the forcing theorem (Theorem 9.1);
- the zero count (Proposition 10.1);
- (OL) for every G >= 8, and hence Lemma A (M_irr(h)) for every h and Corollary R.

The inputs it rests on:
- Route A: Theorem D, Corollary E, Theorem A and Theorem F;
- Route D: (I1), that is (H) = Brini-Schuler + Lemma 5.3 + (V);
- leaf_structure_proof Theorem 11.1;
- brief 2(a),(b);
- the exact data for h <= 4.

How it was checked:
- I re-derived every step of s.1-s.11 by hand and found no gap.
- Every computational claim I tested was confirmed by an independent re-implementation.
- Each re-implementation also passed a sensitivity test: a deliberately wrong variant fails.

Nothing was refuted. The issues found (s.4) are minor or cosmetic:
- a misnamed constant in Cor. 7.2;
- W = alpha is used without its one-line proof;
- a factorial typo;
- status labelling of the final claim;
- one housekeeping point;
- one correction that belongs in Route A's record, which R1 itself flagged and which I confirm.

## 1. Re-derivations by hand

**Lemma 2.1 (level bookkeeping).**
- Write g_k(P) := d_1 d_k F_0(P) for the genus-0 coefficients of Y_P(z) = z - Sum_k g_k z^-k. Each genus-0 factor g_k z^-a w^-b (a+b = k+1) at t-order r'' raises the level by r''+k+1 >= 2.
- At equal level, the dependencies run in one direction only:
  - G(r) depends on V(r);
  - V(r) depends on U(r) and on G(r') with r' <= r-1, because P_j = O(t^(j-1)) and j >= 2;
  - U(r) depends on G(r') with r' < r.
- There are no cycles. So the formal solution is unique, and on every coefficient of level <= 2G-3 it equals the true one. In particular relation_k = [t^k]Phi for k <= 2G-3.
- The (1,n) entries of (H) are tautological. The identity (Y(w)-Y(z)) G = V(z)-V(w) is equivalent to G_1n = G_n1 = V_n together with the anti-diagonal recursion G_{m+1,n} = G_{m,n+1} + Sum_k g_k (G_{m-k,n} - G_{m,n-k}).

**Section 3.**
- Lemma 3.1: expand in |l| > |l_w|. Since [A_- hatV]_+ = 0, one gets -[A hatV]_+ + hatV A_+ = [A_+ hatV]_-.
- (S'), (E') and (P') re-derived. The identities used:
  - theta = -w d_w = h d_l;
  - z'/z = -1/h;
  - d_t Pcal|_l = d_t T|_z - (W'/W) h dPcal;
  - b = N - (W'/W) a.

**Lemma 4.1 (branch point).**
- [g]_+(l0) = Res g/(l-l0), with g = (u^-j)' = j u^-j (l-1-tau)/(l(l-1)). The factor (l-1-tau) cancels, and the residue is j binom(j tau+j-1, j-1) = j Pi_j.
- N(1+tau) = -d_t Sum_j Pi_j P_j W^j + W' = 0, by the genus-0 equation.

**Theorem 5.1 (flow).**
- phi_t(1) = 1 follows from binom(j tau+j, j) - binom(j tau+j-1, j) = Pi_j and W = 1 + Sum Pi_j P_j W^j.
- Using h = -(l+tau) - sigma/(l-1-tau), [poly/(l-c)]_- = poly(c)/(l-c) and N(1+tau) = 0, one gets [d_t Pcal]_+ = hN + W'(l+tau) - nu. Hence d_t phi = -hN, a polynomial.
- In part (iii), both (S') and (E') check out for hatU = [Psi o phi]_- and hatV = [Vfull]_- + delta:
  - the pole terms at l = 1+tau cancel;
  - N delta is a polynomial.
- The t-evolution is unique, because [a.]_- is invertible with a = W + O(t). This gives hatU_t = [Psi o phi_t]_-.

**Theorem 6.1.**
- d_t psi = -v_t(psi).
- The substitution mu = phi_t(l) is legitimate for formal residues at infinity. The perturbation is a t-series of polynomials, and the log term is handled by d log.
- (W'/W + v_t/h) a = -b and Res(bR) = 0 give Phi'.
- Kflow is a t-series of polynomials, because psi_t fixes 0 and 1.

**Theorem 7.1 (curve identity).**
- (a) Needs W = alpha (issue 2). This follows from the sum rule Sum_j A_j = -c^2.
- (b) Lagrange: [w^n]g = (1/n) Res f^-n dmu = (C^n/n) Res u_t^n.
- (c) The residue identity, re-derived:
  - substitute Lambda = e^(t mu), then Lambda = 1 + (e^t-1)s;
  - expand (1+(e^t-1)s)^(-n tau-1) binomially and take the residue at infinity of (1-1/s)^(n tau) s^(m-n);
  - this gives (-1)^(n-1) t^(n-1) e^(-n sigma t) Pi_n 2F1(n+n tau, -n tau; n; 1-e^t);
  - Pfaff (DLMF 15.8.1, 2F1(a,b;c;z) = (1-z)^(-b) 2F1(c-a,b;c;z/(z-1))) turns this into e^(n tau t) H_n(T), exactly as stated.
  - The passage from real tau in (0,1) and small t to the formal identity is valid: both sides are t-series with coefficients polynomial in tau.
  - Hence [w^n](g+Q_t) = (Pi_n H_n (-1)^(n-1)/(n t)) (C^n t^n e^(-n tau^2 t) - zeta^n), which vanishes for C = alpha e^(tau^2 t)/H_1.
- (d) [t^m] Z = eps_m (l-1-tau)/(l(l-1)) is a polynomial that is O(1/l), so it vanishes and eps = 0.
- Consequence: Kflow - K_t = log(W/C) = log H_1(T) - tau^2 t, which is independent of mu. The text writes "log C" (issue 1); only mu-independence is used.

**Corollary 7.2 and the s-form.**
- [l^-1] u^N l^k = (-1)^d binom(N tau, d), and [l^-1] u^N (l-1)^k = (-1)^d binom(N tau + k, d) = (-1)^d binom(N(1+tau)+d-1, d).
- Res Psi s/(1-sl)^2 = -(1/s) Psi'(1/s).
- h(1/s) = -(1-s)/(s(1-(1+tau)s)) and h((1+s)/s) = -(1+s)/(s(1-tau s)).
- The s-form follows. It is pure algebra, valid for every k.

**Lemma 8.1.**
- At tau0 = -k0/q, Lagrange for n = qm gives the factor (1-y)^(m k0).
- Sum_r (-1)^r binom(m k0, r)(qm-r) y^r = qm (1-y)^((m-1)k0) (1-y)^(k0-1)(1-(1+tau0)y).
- The u1-triangularity is correct, since c_{m,m-1} = 1.
- The conditions with m >= 2 follow from those with m = 1, because (m-1)k0 + q - 1 < qm.
- The coefficients used, f_N with N <= 2G-2, involve only y_L with L <= 2G-3.
- The true two-leaf values vanish at tau0 for q|n: Xi_G(m,n) is a polynomial in 1/sigma, and sigma0 != 0.

**Theorem 9.1 (forcing).**
- q R(1/l) = l X(l) with X = A r(l) + B r(l+1). This identity holds; I checked both terms.
- The even coefficients l^-2, ..., l^-(2G-2) of X vanish. Reason: R_0 = Res Psi = 0, and R_k = 0 for even 2 <= k <= 2G-4 by (I1)+T1, since kappa_k = B_k/(k k!) != 0.
- E = N/(l(l^2-1))^A with deg N <= 3A-2G < 0, so E = 0 and X is odd.
- Oddness matches the pole orders at -1 and +1, so ord_0 pi >= A-B.
- From (l+1)^B s1(l) + (l-1)^B s1(-l) = 0 we get (l-1)^B | s1, while deg s1 <= 2B-D-2 < B since B <= D+1. Hence s1 = 0.
- pi_1(1-l) = -(A/B) pi_1(l) applied twice gives pi_1 = (A/B)^2 pi_1, and A != B. So pi_1 = 0 and p = 0.
- Pi_{D+1}(tau0) != 0, since q >= D+2.
- The argument actually proves more: F = 0 mod y^(2G-1) at tau0, i.e. every y_L(tau0) = 0 for D <= L <= 2G-3. My computation confirms this.

**Proposition 10.1 (zero count).**
- Z1 multiplicities from Corollary E:
  - a = 2: tau-order D-1, so sigma-order (D-1)/2;
  - a >= 3: order floor(D/(a-1)) - [a | D+1]; here the tau-order equals the sigma-order since 2 tau0 + 1 != 0, and there are phi(a)/2 sigma-values.
- Denominators: Z1 points have a <= D+1, Z2 points have q >= D+2, so the two sets are disjoint. Within Z2, the condition 2k0 < q avoids double counting k0 <-> q-k0.
- (B1): checked by hand, including the losses from 2k0 < q for k0 = 3, 4 at D = 3, 5.
- (B2): 2 Phi(n) - 1 >= (2 - pi^2/6) n^2, and the constant is (2-pi^2/6)/4 (1/4+1/9+1/16) = 0.037602 > 0.0376.
- The split at D+1 = (7G-6)/24 is correct. For G > 300 the second branch has D >= 87 > 80.

**Induction (s.11).**
- deg Q <= G-2 follows from brief 2(a),(b) applied to R = (D+1) in the range. Note (J-1)...(J-N) contains the factor 0 iff D <= 2G-3.
- Base cases: G <= 4 from exact data; G = 5, 6, 7 from Theorem F (and independently from Route E's engine).
- For G >= 8:
  - D = 1 by Theorem A;
  - even D by T1 at k = D (equivalent to parity README (f));
  - odd D by (OL);
  - then Corollary E gives all T at level D.
- The structure is correct and uniform in G.

## 2. Independent computations (referee code)

| script | what it does | range | result |
|---|---|---|---|
| refkp.pm | Own KP propagation (see note below the table). | all runs below | Overdetermined system consistent every time, e.g. 145145 checks at K = 18 (7 framings), 141108 at K = 22; 0 failures. |
| ref_dump.pl vs kp_point.pl (DUMPY) | Every Y(T) as a linear form in y_1..y_6, compared with the route's E_0-based propagation. | tau = 2/7, exact, 29 multisets | Identical. |
| ref_t1.pl | T1. relation_k computed directly as [t^k] Sum_T prod P_t Y(T)/Aut T, with P_t = (-1)^t (Pi_t/t) t^(t-1) H_t/H_1^t built from the 2F1 series. Compared with the closed form; kappa_k taken from the series of log((e^x-1)/x). | Exact: K = 5 (2/7); K = 8 (-5/8, 7/4); K = 12 (-5/8, 9/4). Mod 1e9+7: K = 10 (3 framings); K = 18 (-5/8, 7/4, 1/3, -11/13, 4, -1/3, -2/5); K = 22 (-7/9, 5/2). Mod 998244353: K = 20 (3/11, -13/7). | 0 mismatches: 441 exact and 4336 modular coefficients. |
| ref_flow.pl | Theorem 6.1 on random paths: direct Sum_T Y(T)P(t)^T/Aut T against Res Psi Kflow, with phi, psi = phi^-1 and log u(psi) built as bivariate (t, mu) series. Also checks phi(0) = 0, phi(1) = 1 and psi/mu, (psi-1)/(mu-1) exact. On the leaf curve: W equals zeta H_1/t, with zeta computed independently from Theorem 11.1's s*; Kflow - K_t is mu-independent and equals log H_1 - tau^2 t. | K = 6, 11 (-5/8); 11 (7/4); 10 (-2/5); 15 (-3/10); 14 (5/2); 2 random paths each, plus the leaf curve | 0 mismatches. Perturbation 1+tau -> 1+2tau gives 18 mismatching orders. |
| ref_res.pl | Residue identity of s.7(c). | n <= 8, t^0..t^12, framings -5/8, 7/4, -2/5, 1/3 | 416/416. Perturbed exponent: 131 mismatches. |
| ref_special.pl | Lemma 8.1: row space of my own two-leaf regularity conditions against the p-conditions [y^N] p, q <= N <= 2G-2. | G = 3..8, all 176 framings with q <= 2G-2 | 0 differences. |
| ref_forcing.pl | Forcing theorem without T1, Lemma 8.1 or the odd-function argument. Unknowns y_D..y_{2G-3}. Rows: the directly computed relation_k = 0 (k <= 2G-3), plus two-leaf regularity Y(m,n)(tau0) = 0 for q \| n, m+n <= 2G-1. | Every tuple in the hypothesis range for G <= 11: 83 tuples, the same set as the route's t_forcing.pl. Also all tuples outside the range for G <= 8, and one exact (BigRat) case (G,D,q,k0) = (6,3,5,2). | y_D forced in 83/83, and all y_D..y_{2G-3} forced in 83/83. Outside the range, 248 of 359 forced. The exact case is forced. |
| ref_zeros.pl | Z1 by brute force over partitions (D <= 14) against the formula. Z2 by Moebius counting and by direct enumeration (G <= 60). Full count and the bounds (B1), (B2). | G <= 1000 | Z1 brute force = formula (Route A's 0,1,2,5,6,11,13,19,22,29). Moebius = direct. Failures only at (G,D) = (4..7, 3), with Z1+Z2 = 2,3,3,5. Smallest margin for G >= 8 is 1, at (9,3). No violation of (B1) or (B2). G = 800, D = 3 margin 618, as the route says. |
| ref_global.pl, ref_global_c.pl | Global nullities: polynomial ansatz, sampling, left-kernel elimination. | G = 3..7 | Window G, TLR: 4, 7, 6, 1, 0 (G = 3..7). Window G-1, TLR: 2, 3, 0, 0 (G = 3..6). At G = 3, 4: TLR+CF 1, 2; full 1, 0; full+CF 0, 0. All as claimed. With c = 2 added as an unknown: 11, 11, 7 at G = 4, 5, 6, which are exactly Route A's numbers (see issue 7). |

Note on refkp.pm. It writes the linearised dispersionless Hirota identity at every point of the slice in anti-diagonal form:
- G_{1n} = G_{n1} = V_n;
- G_{m+1,n} = G_{m,n+1} + Sum_k g_k (G_{m-k,n} - G_{m,n-k});
- <V_n>_R = (J_R + n) Y(R+n) (string);
- the genus-0 coefficients <g_k>_S = sigma Pi_k Pi_S (s+k)!/(s+k+1-m)! are taken directly from Xi_0 = Gamma(J)/Gamma(J-r+3), with no use of W.

Re-runs of the route's own scripts:
- kp_point.pl gives KP consistency rank 0 at (D, Lmax, tau) = (1, 7, 2/7) and (1, 9, -3/11). So the system is consistent for arbitrary one-leaf data at generic tau, and T1's linear form is canonical rather than an artefact of the propagation scheme.
- t_forcing.pl 8 and 11: 18 and 83 tuples, 0 failures.
- t_special.pl 6: 31 framings, 0 differences.
- t_res_ut.pl: all agree.
- The existing ode_path / ode_path2 outputs were inspected.

## 3. Inputs and normalisations

**Leaf curve.**
- Xi^c_h(R) = c^-(3h-3+2r) Xi_h(R) (Route D, consistent with my own scaling count c^(3-2r) times c^(-3h)). This gives a factor c^-2 per leaf.
- Eliminating the string term gives Sum_n e^n/n! (J)_n c^-2n = (1 - e/c^2)^-J with e = A_1 + c^2, i.e. alpha^-J with alpha = -A_1/c^2 = zeta H_1/t.
- C_j = c^2 Pi_j^2/j and H_j = 2F1(-j tau, -j tau; j; T), since j x/c = -j tau.
- So b_j = (-1)^j (Pi_j^2/j) t^(j-1) H_j/H_1^j, and zeta cancels. This is exact.
- Cross-checks with Route D's independent data:
  - the leading order gives the T_k weights;
  - [t^1] rho_j = 0 reproduces Lemma 5.2 (the D = k-1 terms cancel), and T1 correspondingly has beta_{k,k-1} = 0;
  - b_2^(3) = -(1-x)^2 x^2/(12c^6) reproduces Route D's residual -(1-x)^2 x^2/(12c^4) Xi^c(2) (check_polar_relations.out, 4 values of x);
  - kappa_1 Pi_2 = (1+4 sigma)/2 is D5.1.

**(I1) needs only genera < G.**
- A stable graph of genus G with an edge has no genus-G vertex.
- Reason: dilaton leaves carry psi-degree >= 2, so a genus-0 vertex with e edge half-edges and k leaves needs k <= e-3, hence e >= 3. A tree with one genus-G vertex and all other vertices of degree >= 3 violates the degree sum.
- So all edge graphs involve only Xi_h with h < G and are regular by Theorem R's counting. (H) then leaves only the pole t^(2-2G). This confirms Route D's Lemma 5.3 premise (issue 5: worth writing out).

**Theorem D at every point of the slice.**
- My anti-diagonal form and the route's E_0 form give identical propagations.
- The level and excess restrictions are respected (Lemma 2.1).

**Genus 0.**
- d_1 d_n F_0 = sigma Pi_n W^(n+1)/(n+1) (plane-forest count, checked against the Catalan case).
- Y_0(1/u(l)) = l + tau. Checked by Lagrange: [u^n](1/y) = -sigma Pi_n/(n+1), with [u^0] = -tau.

**Literature.**
- Pfaff transformation: DLMF 15.8.1, used in the correct form.
- Lagrange inversion: standard.
- No other literature enters beyond the project inputs, which the route lists correctly.

## 4. Issues (none affects correctness)

1. **(cosmetic) s.7, Corollary 7.2.**
   - The text says "Kflow_t - K_t = log C". The constant is log(W/C) = log(alpha/C) = log H_1(T) - tau^2 t; I verified this numerically.
   - R1_checks.wls (C3) has it right. Only mu-independence matters.
2. **(minor) s.2 and s.7(a).**
   - "W(t) = alpha(t) on the leaf curve" is used without proof.
   - It is one line: alpha = 1 + Sum_{j>=2} Pi_j P_j alpha^j is the sum rule Sum_j (-1)^j (Pi_j^2/j) zeta^j H_j = -t (leaf_structure_proof Corollary 11.2(a)). Add it.
   - Verified by ref_flow.pl with zeta built from s*.
3. **(cosmetic) s.1, (I4).** Should read Xi_0(S+{1,n}) = (s+n)!/(s+n+1-m)! (the factorial is missing).
4. **(minor) Status of "Lemma A for all h" / Corollary R.**
   - The proof is complete, but it also rests on computer-verified base cases: Theorem F for h = 5, 6, 7 (rigorous modular ranks; Route E's engine also covers h = 5, 6, 7) and exact computations for h = 3, 4.
   - It is conditional on (H) = Brini-Schuler through Route D, as the route says.
   - Label it "proved modulo the cited inputs and the computer-verified base cases h <= 7".
5. **(minor, exposition of an input)** Record the "no genus-G vertex in a graph with an edge" argument of s.3 above. It is what makes (I1) need only genera < G.
6. **(cosmetic)** out/count_zeros2_G3to1000_raw.txt begins with stray diff lines ("33,42d32 ...").
7. **(minor; Route A's record, flagged by R1 and confirmed here)**
   - Route A's TLR nullities 11, 11, 7 at G = 4, 5, 6 include c = 2 as an unknown, although Route A's statement of TLR(G) imposes v_2 = 0.
   - With c = 2 excluded the nullities are 7, 6, 1 (window G) and 3, 0, 0 (window G-1).
   - So the natural form of TLR, with the Hurwitz window G-1, already holds at G = 5, 6.

## 5. TR audit

Clean.
- The R1 scripts contain no TR, loop equation, spectral-curve recursion, BKMP, P2 or TR script. The only grep hit is mvlib.pm's comment "No topological recursion anywhere".
- The uniformising coordinate u(l) = (l-1)^tau l^(-1-tau) is a genus-0 change of variables of the dispersionless Lax function, which the brief allows as algebra. It is applied to the Hodge CohFT, where TR would even be allowed.
- The only inputs about the conifold theory are (H) (Brini-Schuler, degeneration) and the leaf closed forms (no TR).
- My own scripts use none of these either.

## 6. Obstructions and simpler path

**The stated obstructions are real.**
- The uniform count fails exactly at G = 4..7, D = 3, so base cases are needed below genus 8; Theorem F supplies them.
- TLR(G) is unproved but no longer needed.

**Simpler path.** I found nothing shorter than T1 + Lemma 8.1 + odd-function argument + counting. Two remarks:
- The forcing theorem actually shows F = 0 mod y^(2G-1) at each admissible tau0.
- Many tuples outside 3(q-k0) < 2G are also forced (248/359 for G <= 8). A sharper analysis of E = N/(l(l^2-1))^A when 3A >= 2G might remove the need for Theorem F at G = 5..7. Optional.

## 7. Limitations

- R1_checks.wls was not run (no wolframscript). Its checks (C1)-(C6) are covered by my Perl checks: residue identity, branch point and phi(1) = 1, curve identity, T1, forcing, zero count. The main session may still run it.
- Modular results are evidence, not proof. Two primes and exact spot checks agree, and the main statements have complete proofs.
- The pointwise-scan claim (G <= 12) was reproduced only in part: the forcing tuples for G <= 11 and the outside tuples for G <= 8.

## 8. Files (referee\)

- refkp.pm
- ref_t1.pl
- ref_dump.pl
- ref_flow.pl
- ref_res.pl
- ref_special.pl
- ref_forcing.pl
- ref_zeros.pl
- ref_global.pl
- ref_global_c.pl
- out\ (all outputs listed in s.2)
