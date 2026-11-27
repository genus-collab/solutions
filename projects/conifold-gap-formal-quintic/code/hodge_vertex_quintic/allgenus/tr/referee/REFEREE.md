# Referee report on allgenus/tr (TR for the model curve x = p^4 + v/p, y = p) -- 2026-11-13

Labels as in the project: PROVED / VERIFIED / CONJECTURAL / REFUTED. AI referee; no human has read this. Referee code: rf_*.py in this
folder, written from scratch (it does not import any tr_*.py; project modules fq_fp, fq3_mvertex imported read-only; data read:
natural_form/lg/lg_symbol_rows_Q.json). Every script writes rf_<name>.log ending in "ALL CHECKS PASSED". Runs: one or two processes,
longest single run 6.5 min (rf_t8, F_5 at the special point).

## 0. Verdict: SOUND, with minor corrections

Every central claim of TR_MODEL.md was reproduced by an independent implementation along a different computational route, and
extended: the dictionary TR <-> frozen theory (including the normalisation Lam^(g-1) kappa^(-D), Lam = -5 kappa) for kappa = 1, 4, -1
against MVertex6 (the worker compared kappa = 1 only); R_TR = M_0(z/kappa) exactly to z^10 (worker z^8); Theorem G exactly over Q at
random points of H_{0;1,r}, r = 1..5, mod p at genus 4 and 5, and on a curve with TWO simple poles (a test of the locality mechanism
the worker did not run); F_4, F_5 at the special point; the u -> 0 regularity and limits; all quoted exact correlators, incl.
<t5(H)>_3 = 53/1209600. The literature claims were checked against the PDFs. Nothing was refuted. Corrections (s.6) are minor: author
names in files (hard rule), two loose or incomplete remarks ("x = -W flips psi -> -psi"; the metric-sign remark on (2-7)), wording
("any point"), and a small gap in the n = 0 identification (Phi-pairing for all d) that is closed here (s.3.3).

## 1. What was read

* Worker: TR_MODEL.md (249 lines), all 12 tr_*.py, all 8 logs.
* arXiv:1509.06954v2 (PDF pp. 1-3, 6-13, 16-19, 24-27): (2-7), (2-26)-(2-35), Thm 5.1 with conditions (1)-(4), Prop. 5.2, Thm 8.1,
  proof of Thm 8.1 (p_Dub = p - a_1/(n+1), d = (n-2)/n). The worker's transcription in TR_MODEL.md s.1 is accurate.
* arXiv:2502.18115v2 (PDF pp. 1-21): Thm 1.1, Cor. 1.2 and (1.3), s.2.1 hypotheses (compact Sigma, meromorphic dx, dy, simple
  ramification), kernel K_i (p. 6) and F_g (1.1) = (3.1) (identical to arXiv:1509.06954 (2-28) and to the referee's engine), the
  point sets of p. 13 (b-points: singular points of y with Res y dx != 0 or of x with Res x dy != 0; here b = 0 (pole of x) and
  b = infinity), Thm 3.2 and proof, s.3.2 (x = z + 1/z), Cor. 3.3, Prop. 3.4.
* Context: LG_MODEL.md, UNIFORM_REDUCTION.md s.6.3-6.5, the plan (one-point law, Y^C).

## 2. Claim by claim

| # | worker claim (status) | referee | how |
|---|---|---|---|
| 1 | arXiv:1509.06954 Thm 8.1 applies to p^4 + v/p, v != 0 (PROVED, quoted) | confirmed | PDF read; conditions (1)-(4) of Thm 5.1 hold (d_lambda p has no zeros over finite lambda; distinct critical values 5 p_c^4) |
| 2 | R_TR(z) = M_0(z/kappa), row = branch point (VERIFIED z^8; all orders argued) | confirmed (z^10 exact) | rf_t5_rmat.log, own local coordinate, Q(zeta5); transposed reading fails from z^2 |
| 3 | (2-32) gives the frozen unit and the dilaton leaf (VERIFIED) | confirmed | rf_t5_rmat.log: same four colour-0 coefficients -3/20, -9/800, 11571/80000, 2415609/6400000 |
| 4 | omega = Lam^(g-1) Sum kappa^(-D) <..>^frozen prod d xi, Lam = -5 kappa (VERIFIED) | confirmed, extended | re-derived (s.3.1); rf_t4_frozen.log: kappa = 1, 4, -1 x two primes, 6495 colour correlators each (1660 non-zero), 0 mismatches |
| 5 | exact correlators (<H>_1 = -1/12, ..., F_3 = 1/25200, genus-2/3 one-point lists) (VERIFIED) | confirmed | rf_t4_frozen.log: 18 flat values reconstructed identically at 2^61-1 and 2^89-1 |
| 6 | Theorem G (PROVED modulo arXiv:2502.18115 Cor. 1.2) | confirmed | statement and hypotheses checked in the PDF; extraction, degree count, pole evaluation, inversion re-derived (s.3.2) and recomputed (rf_t3_swapformula.log) |
| 7 | F_2..F_5(x = W) = B_2g/(2g(2g-2)); c_2..c_5 (VERIFIED) | confirmed | rf_t4 (F_2, F_3 vs MVertex6), rf_t8_fg_special.log (F_4 both primes, F_5 at 2^61-1) |
| 8 | F_g depends only on the residue on H_{0;1,r}; primary-sector lemma g >= 2 (PROVED mod quoted) | confirmed | rf_t2_thmG.log exact over Q, r = 1..5; rf_t7_extra.log: primary sector at (2,3), (3,2) |
| 9 | printed closed forms of arXiv:2502.18115 s.3.2, Cor. 3.3, Prop. 3.4 have the wrong sign (VERIFIED) | confirmed | p. 12 prints +B_2g/(2g(2g-2)) for x = z + 1/z, p. 19 derives B_2g/(2g(2-2g)) from "(-1)^m (m-1)! z^m"; correct is (-1)^(m-1) (dz/dx = -z^2); referee TR on x = p + 1/p gives -1/240 exactly |
| 10 | homogeneity expansion (PROVED) | confirmed | re-derived (s.3.4) |
| 11 | u -> 0: c_m = 0 for |m| < 5(2g-2+n), non-zero limit (VERIFIED 13 cases) | confirmed for 8 cases | rf_t6_v0.log, both primes, all quoted limits; control: one branch point alone does not cancel |
| 12 | limits = 4-spin correlators (CONJECTURAL) | not checked | label appropriate |
| 13 | <H>_g + <t1(1)>_g = 0 at g = 2, 3; 1/8 at g = 1 (VERIFIED) | confirmed | values in rf_t4_frozen.log |
| 14 | the plan's linear one-point law fails for <t5(H)>_3 (REFUTED) | confirmed | 53/1209600 reproduced; the law gives 1/604800 in absolute value; robust to signs and to z -> mu z (s.5) |

## 3. Independent derivations

### 3.1 Normalisation (claims 2-4)
With x - x(c_j) = -w_j^2/2 ((2-31) convention) and x''(c_j) = 5 kappa c_j^2: w_j = sqrt(Lam) c_j (p - c_j) + .., Lam = -5 kappa.
(2-30): Delta_j^(1/2) = sqrt(Lam) c_j; (2-34): xi^j = -1/(sqrt(Lam) c_j (p - c_j)) = xit^j/sqrt(Lam). The frozen TFT has
Deltahat_j^(1/2) = zeta^(4j) = c_j (MVertex6). So e_j = sqrt(Lam) et_j, and as tensors on span(et_j) the TR TFT is Lam^(g-1) x frozen
with metric /Lam; a CohFT rescaled by Lam^(g-1) with metric /Lam is again a CohFT (vertices Lam^(g_v-1), edges Lam, total Lam^(g-1)),
and the R-matrix (graph-sum) action commutes with this. R_TR(z) = M_0(z/kappa) multiplies the degree-D part by kappa^(-D) (R-coefficients of
order k come with psi^k, edge terms of order a+b+1 with a class of degree a+b+1, leaves T_e with kappa_(e-1)). Units: Sum_k
Delta_k^(-1/2) e_k = Sum_k et_k/c_k = the frozen unit. Hence frozen = Lam^(1-g) kappa^D x (coefficient of prod d xit), as stated.
For kappa = -1 this is 5^(1-g) (-1)^D: the (-1)^D is psi -> -psi only up to the overall sign (-1)^(3g-3+n) (wording, s.6).

### 3.2 Theorem G (claim 6)
(1.3) of arXiv:2502.18115 with y = z and no logarithmic singularity: S(t) - 1 = Sum_{k>=1} t^(2k)/(4^k (2k+1)!), the exponent is
Sum_k u^(2k+1) hbar^(2k) x^(2k)/(4^k (2k+1)!), and [hbar^(2g-1)][u^m] of exp(..)/(u hbar) keeps j factors with Sum k_i = g and
m = 2g-1+j, giving the worker's I_g (re-derived by hand). At infinity every term has degree <= 1 - (2g-2) r - 2g <= -5 (checked).
At the pole: 1/x' = -(z^2/v)(1 + O(z^2)), D^n z = (-1)^n n! z^(n+1) v^(-n) (1 + O(z^2)), x^(2k) = (2k)! v z^(-2k-1) (1 + O(z^(2k+1))):
every term is z^(-1)(1 + O(z)), so only leading coefficients enter, giving v^(2-2g) A_g; the binomial resummation and the inversion
[s^1] artanh(s)^(1-2g) = (1-2g) [w^(2g)] (w coth w) = (1-2g) 4^g B_2g/(2g)! were re-derived: A_g = -B_2g/(2g). The b-points of
p. 13 for x = P(z) + v/(z-a), y = z are exactly a and infinity (Res_a x dy = v, Res_inf x dy = -v, Res_inf y dx = v).
Status: PROVED modulo arXiv:2502.18115 Cor. 1.2 (itself resting on the x-y duality theorems quoted there), as the worker says.

### 3.3 The n = 0 identification in all genera (gap closed)
F_g^TR = (1/(2-2g)) Sum Res Phi omega_{g,1} equals Lam^(g-1) kappa^(3-3g) F_g^frozen in ALL genera iff
P_{i,d} := Sum_j Res_{c_j} Phi d xit^i_d = -delta_{d,1}/c_i for every d (the worker checked d <= 3). Proof: integrating by parts at
each branch point, P_{i,d} = -Sum_branch Res xit^i_d y dx = Res_0 + Res_inf of the rational form xit^i_d y dx, y dx = kappa (p^4 - 1/p) dp.
At p = 0: Res = -kappa xit^i_d(0), which is 1/c_i^2 (-kappa) for d = 0 and 0 for d >= 1 (1/x' = p^2/(kappa(p^5-1)) vanishes at 0). At
infinity xit^i_d = O(p^(-1-4d)): for d = 0 the residue is kappa c_i^3 = kappa c_i^(-2) (cancels), for d = 1 it is -1/c_i, for d >= 2 zero.
PROVED; VERIFIED d <= 10, kappa = 1, 4, -1, both primes (rf_t7_extra.log).

### 3.4 Homogeneity (claim 10)
x_u(k p) = k^4 x_1(p) and y(k p) = k y(p) for u = k^5, so phi_k^* (y dx_u) = k^5 y dx_1, B invariant, and TR scaling gives
omega^(u)(p) = Sum_m c_m k^(5(2-2g-n) + |m|) prod dp/p^(m+1) = Sum_m c_m u^(2-2g-n+|m|/5) ..: the worker's statement.

## 4. Independent computations (referee engine rf_tr.py: local coordinate x(a+t) - x(a) = l2 zeta^2 with t(zeta) by series
inversion, deck map zeta -> -zeta, all branch points, no symmetry; the worker used t and the series sigma(t) with Z/5 symmetry)

* rf_t1_regress.log: Airy curve = psi-class intersection numbers for (0,3)..(3,1) and dilaton equations, Q and 2^61-1; Gaussian
  curve x = p + 1/p exact over Q: F_2, F_3, F_4 = -1/240, 1/1008, -1/1440 = B_2g/(2g(2g-2)) (sign anchor for claim 9).
* rf_t2_thmG.log (exact over Q, rational branch points): F_2 at r = 1..5 (two points at r = 4), F_3 at r = 1..4, pole moved to
  b = 3/2, -2/7: all equal B_2g/(2g(2g-2)) v^(2-2g). Sum Res x omega_{g,1} = 0 (g = 2, 3): F_g is invariant under y -> y + const.
  Mod 2^61-1 and 2^89-1: F_4 v^6 = -1/1440 at r = 2, 3, 4, F_5 v^8 = 1/1056 at r = 2.
* rf_t3_swapformula.log: A_g by the finite sum = by series arithmetic = -B_2g/(2g), g = 2..14; the full integrand of (1.3) (no leading-term
  shortcut) for five polynomial parts, g = 2..5: Res_inf = 0 and F_g = B_2g/(2g(2g-2)) v^(2-2g); TWO simple poles,
  x = z + v1/z + v2/(z - b): TR F_2, F_3 == formula (1.3) == B_2g/(2g(2g-2)) (v1^(2-2g) + v2^(2-2g)) at both primes; control: the
  one-pole value alone does not match. This confirms the locality mechanism behind Theorem G beyond the cases the worker tested.
* rf_t4_frozen.log: TR on x = kappa W, kappa = 1, 4, -1, both primes, converted with the referee's xi basis (partial fractions,
  triangular elimination, zero remainders) against MVertex6(M_0): (0,3), (1,1), (0,4), (1,2), (2,1), (1,3), (2,2), (3,1), every
  colour correlator with D <= 6: 6495 per (kappa, prime), 1660 non-zero, 0 mismatches (the worker's 4665 non-zero minus its 3005
  of (0,5) is exactly 1660); F_2 = 1/1200, F_3 = 1/25200 from all three kappa. Flat values: the 18 of TR_MODEL.md s.3.
* rf_t5_rmat.log: claims 2, 3 exact over Q(zeta5) through z^10 for kappa = 1, 4, -1.
* rf_t6_v0.log: claim 11 for (0,3), (1,1), (0,4), (1,2), (2,1), (1,3), (2,2), (3,1), both primes, with the quoted limits
  -5/8, 2079/128, 3132675/1024 and the (0,3), (1,2) monomials.
* rf_t7_extra.log: s.3.3; primary sector: at (2,3) only <H,H,H>_2 = -1/50 = d_u^3(c_2 u^-2), at (3,2) only <H,H>_3 = 1/1260 =
  d_u^2(c_3 u^-4) survive among all flat primary correlators (both primes); genus 1, n = 3: only <H,H,H>_1 = -1/6 (information).
* rf_t8_fg_special.log: special point, F_4 = -1/1440 (both primes), F_5 = 1/1056 (2^61-1); (-5)^(1-g) F_g = c_4, c_5 on record.

## 5. Attempts to refute; controls

* Tests can fail: (C1) x = -W read with 5^(1-g) and no kappa^D: 5-145 mismatches per (g,n); (C2) Lam = +5 instead of -5: every
  non-zero entry fails at g = 0 and g = 2 (factor (-1)^(g-1)), g = 1 unaffected as it must be; (C3) MVertex6 built from M_0(-z):
  fails at (1,1), (2,1) (cannot fail at (0,3), D = 0 only). The transposed reading of (2-31) fails from z^2. In rf_t6 the
  cancellation at low |m| needs all five branch points. In rf_t3 the one-pole value is rejected on the two-pole curve.
* Hypotheses of quoted theorems: no logarithmic singularity, compact P^1, simple ramification, y unramified: all hold. The
  common pole of x and y at infinity is a b-point in the sense of arXiv:2502.18115 p. 13; its residue vanishes by the degree count
  (also recomputed). arXiv:1509.06954 Thm 8.1: semisimple points, here distinct critical values.
* One-point law (claim 14): the law ties <t_d(H^a)>_3 to [z^-d] Y^C_a; the five leading cases match in absolute value with
  ratio 1 per power of z, so a rescaling z -> mu z needs |mu| = 1, while <t5(H)>_3 / <t4(H^2)>_3 = -53/8 against the law's -mu/4.
  The refutation of the LINEAR law is therefore robust to sign and scale conventions; it is not a statement about the plan's
  non-linear V(q) or jet-corrected versions (the worker says so).
* "Regular" vs "zero": the worker keeps them apart (non-zero limits listed). Low-order evidence is labelled VERIFIED/CONJECTURAL.
* Not reproduced: (0,5) against MVertex6, the (1,4), (2,3), (3,2), (4,1) expansions at infinity, the exact Q(zeta) TR run (the
  referee's special-point values are mod two primes with identical rational reconstruction).

## 6. Corrections requested (minor)

1. Hard rule "no person's name in any file; cite by arXiv number": TR_MODEL.md cites papers by their authors' surnames (lines 26,
   45, 62, 65, 157, 172, 192, 235, 242), and so do the docstring and log header of the worker's formula script (its lines 2, 111;
   log line 1; the file name itself carries a surname). Replace by arXiv:2502.18115, arXiv:2211.12259, arXiv:1509.06954.
   (Standard names of mathematical objects, e.g. the Airy curve, are not at issue; author citations are. The referee's files
   contain none.)
2. TR_MODEL.md s.0 item 2, s.2.4 and the summary: "x = -W ... flips psi -> -psi" is exact only up to (-1)^(3g-3+n); state the
   kappa^(-D) formula (already correct in s.2.4) as the claim.
3. Add s.3.3 above (or a pointer) to s.4.1/s.0 item 3: the frozen F_g = (-5)^(1-g) F_g^TR in all genera needs the Phi-pairing
   for all d, not only d <= 3.
4. s.2.3 "PROVED modulo that standard uniqueness": name the inputs (homogeneity of R_TR: arXiv:1509.06954 Prop. 5.2 under Thm 8.1;
   M_0 satisfies (5/4)[C_(n+1), Z] = (V - n) C_n: LG_MODEL.md s.3-4). A direct route also exists: xit^i x' = -(1/c_i) Sum_k c_i^k p^(2-k),
   so the Laplace transform of d xi^i is an exact combination of the stationary-phase series J_k, k = -2..2, of LG_MODEL.md s.1.
5. s.0 item 3 "any point" of H_{0;1,r}: read "any point where x is simply ramified" (as in the statement of Theorem G); F_g is then
   constant on that dense connected set and extends.
6. s.1, remark on (2-7): "the two CohFTs differ by Omega_g -> (-1)^(g-1) Omega_g" is incomplete. Combining the verified dictionary
   with LG_MODEL.md s.5.4 (frozen = Hurwitz CohFT with eta x 20, z -> -4z) gives TR(lambda) = (-1)^(g-1+D) x (Hurwitz CohFT in the
   LG_MODEL convention): metric sign AND z -> -z (s.2.3). For F_g (D = 3g-3) the two signs cancel, consistent with
   F_2^TR(lambda, v = 4) = -1/3840 = F_2^Hurwitz(a_5 = 4) of LG_MODEL.md s.5.4. Harmless (the dictionary is fixed by computation).

## 7. Hard rules (worker)

Own folder only: no project file outside allgenus/ was modified after 17:30 today (directory scan). Python stdlib only, py -I -B
style, no wolframscript seen. Files < 250 lines (TR_MODEL.md 249, tr_field.py 248). ASCII only (byte scan: 0 non-ASCII). Labels
present and honest. Author names: see s.6.1 (violation, minor). Run times <= 140 s.

## 8. Files (this folder)

rf_ser.py (fields Q, Q(zeta5), F_p; Laurent series), rf_tr.py (TR engine), rf_curves.py (curves), rf_xi.py (xi basis, conversion),
rf_t1_regress.py, rf_t2_thmG.py, rf_t3_swapformula.py, rf_t4_frozen.py, rf_t5_rmat.py, rf_t6_v0.py, rf_t7_extra.py, rf_t8_fg_special.py
and their logs rf_t*.log.
