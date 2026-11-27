# Referee report on the gauge identity (G) (task key `gauge`)

Date 2026-11-13. Folder `allgenus/gauge/referee/`. AI referee; no human has read this. Labels: PROVED / VERIFIED / CONJECTURAL /
REFUTED for statements made here; for the worker's claims the verdicts are confirmed / weakened / refuted / unverifiable.
Every log quoted was produced today by the scripts in this folder (py -I -B, stdlib only, exact Q or mod 2^89-1, 2^61-1).

## 0. Verdict: SOUND WITH MINOR CORRECTIONS

* Every mathematical claim of `GAUGE_IDENTITY.md` and `GAUGE_IDENTITY_HAT.md` that is labelled PROVED was re-derived on paper (s.2);
  no gap was found beyond the stated dependence on Theorem U (UNIFORM_REDUCTION s.3, AI-audited only) and, for the actual R-matrix,
  on (H). Every VERIFIED claim was re-derived by OWN code that imports no `gi_*` module, mostly by a different route (s.3); the one
  exception is the linear-dressing regression of HAT s.5 (statuses re-read from the record, Fourier dictionary confirmed exactly).
* New evidence (not in the worker's files): (i) the regular gauge computed by a route independent of Theorem U's construction
  (u-recursion with unknown initial value C(z) = P(0, z) plus linear algebra over Q): the regularity equations have EXACTLY the
  solution space C -> C c(z^5) below z^24, and the canonical C equals the project gauge, both frames; (ii) the actual G equals the
  canonical gauge EXACTLY OVER Q on the referee's rational table (397 Fourier components, zhat <= 16, u <= 11); (iii) mod p with M_0 taken
  from the LG-model recursion (no quintic input) instead of the quintic rows; (iv) negative controls the worker did not run.
* Corrections (all minor, s.4): the claimed frame dependence of the leaf bookkeeping is spurious (the value is -(g-1) c_g/15 in EVERY
  frame); the production-row statement of GAUGE_IDENTITY s.0 item 5 is an inference, not a test; (G*) had no negative control; the
  float diagnosis is a reconstruction (the first script is not preserved), strongly supported; an author-initials citation tag in a heading; a sign clash in
  the use of `c_5`; the one-point dimension caveat is stated for vertex genus h = g only.

## 1. What was read and re-run

Read in full: both notes, all nine `gi_*.py`, all 19 logs, the worker's scratch `diag_gstar_orig.py`; background UNIFORM_REDUCTION.md,
LG_MODEL.md, NATURAL_FORM.md s.0-s.6, the plan; project modules `ur_gauge.py`, `fq_conifold.py`, `fq_series.py`, `fq_sser.EdgeLeaf`,
`natural_form/referee/ref_hat.py`, `ref_rows.py`, `ref_cross.py`, `natural_form/extra/nat_kregular_linear.{py,log}`.
Re-runs of the worker's own scripts from COPIES in the referee scratch folder (nothing in `gauge/` was touched): `gi_diag_old.py`
gives a byte-identical log (119 identities, 73 float failures, first four (1,0,2) (1,2,0) (1,2,10) (1,2,15)); `gi_hat_edge.py ref 61
1/3 -2/7 24 18` gives an identical log apart from timings. All 19 worker logs end in ALL CHECKS PASSED / DIAGNOSIS CONFIRMED.

## 2. Paper check of the proofs (all confirmed)

* Lemma V: the residue of x^m exp(-(x^4/4 + v/x)/z) dx gives I^C_m (index j = m + 4i + 1); (i) and (ii) re-derived; Y^C solves the model
  system (the b = 1 step uses (ii) at m = 0: I_3 = v I_-2); Y^C(0) = -e_1 (only I_-1 has a v^0 term). VERIFIED term by term (rf_misc).
* Lemma 1: det((n+1) z - B_-1) = ((n+1) z)^5 since B_-1 = r E_12 is nilpotent; a unit of Q((z)). Correct; no monodromy needed.
* Lemma Q: the induction ord y^(n+1)_a = ord b_a + ord y^(n)_(a+1) + 1 uses c_(a+1) = c_a - 1 + 5 delta_(a,1) (checked for all a) and
  exactness of leading terms (integrands have valuation >= 0; b_1 meets y_2 of order (4n+4)/5 >= 4). Correct.
* Theorem G: P Y~ lies in L[[u]] (at each u-order a power series in z times a polynomial in 1/z), solves the quintic system, has
  u^0-term P(0,z)(-e_1) = -P_11(0,z) e_1 by (R1); Lemma 1 concludes. Correct. Remark: the whole content is (R1) + Lemma 1; note that
  EVERY regular solution with a z-independent initial value is a series in 1/z (e.g. value e_0), so "pure 1/z" does not single out
  Y^C; its initial value -e_1 (image of the log-monodromy, re-derived: W_2(0) = z) does. The notes use the value, so no error.
* Proposition S: z S' = At^T S - S At needs A^T eta = eta A, which is b_0 = b_2, b_3 = b_4 (re-derived entrywise); T = eta^-1 S is a
  regular self-gauge of At with T_0 invertible; Theorem U uniqueness (for (At, At), identity included) gives T = t(z) 1; entry (2,1) at
  u = 0 gives t = c(z) c(-z). Correct. eta^-1 = eta/25, so the forms P^T eta P and P eta^-1 P^T agree (checked).
* (G*), the u-expansion (G)_n, (G*)_n, (R1), (R2), (R3a-c), and the explicit u^3, u^4 row-2 relations (including (1 - P_22(0,z))/(24 z^4)
  from (eta_4)_2 = (theta_4)_2 = -1/(24 z^4)): all re-derived. Correct.
* Diagonal: prod_a p_(0,a) constant from the z^1 diagonal; p_(0,0)^5 = prod_(j<4) (b_j/bt_j)^(4-j) = v/(u varpi_0^5) using
  v'^5 = v/(u (1-u)^5); hence p_(0,0) = tau^(-1/4)/varpi_0. Re-derived; also g_0 g_3 = 1 (needed for HAT s.2) re-derived.
* Hatted: [zhat^k] Ghat_(b+k,b) = u^mu(b,k) tau^((1-c_b)/4) p_(k,b); Lambda eta Lambda = s^6 eta (c_a + c_(3-a) = 2); homogeneity of Y^C;
  the regrouping bijection (needs m = a - 1 mod 5, checked); closed form of phi (signs (-1)^(m+i+j) = 1, phase zeta^(i(m-2)));
  Lemma L (difference quotient; the sign (-1)^m from Q'(-psi)); (E) from h = 0; (T) from [x^-1] (Q' f)_0 = -t_c, unit row zeta^i/xhat,
  [x^-1] f_i = -u tau^(-3/2) zeta^(-i)/5; leading term T'_2 = u p_(1,4)(u). All re-derived. Correct.
* Conditional leaf bookkeeping (HAT s.5): the algebra is right; kappa = -1 from the dilaton equation is right; status CONJECTURAL
  is the correct label. The "(E) is not sufficient" example (string equation, support rule c = 1 - m) is correct.
* Actual R-matrix (s.6): given (H), P^act = P^can c^act(z^5), c^act symplectic, hence exp(a_1 z^5 + a_3 z^15 + ...): z^10, z^20
  forced. Correct. Not circular: the referee rows fix the free WKB coefficients (u^(4n/5) of ph_n, 5 | n) by the QRR values at u = 1
  (ref_rows.wkb), i.e. at large radius, while the canonical gauge is normalised at u = 0.

## 3. Independent computations (own code; no gi_* import)

| script / log | route | result |
|---|---|---|
| `rf_core.py` | own Frobenius solution of the PF operator (5-term recursion on monomials), own log obstruction for g, own b, bt, tau, v, Y^C, Y^C_Q | gamma = 7/10; t_c, varpi_0, v agree with the worker's numbers in both frames; prod b = prod bt |
| `rf_frob_run.py` -> `rf_frob_run_24_ab_{0_0,1o3_m2o7}.log` | NEW ROUTE: gauge equation as a u-recursion with unknown C(z) = P(0,z) (85 graded unknowns, (R1) imposed), regularity as 1574 linear equations over Q, structural completeness ("reach") | homogeneous solution space below z^24 = exactly the c(z^5) directions (free: C_11 at z^0, z^5, z^10, z^15, z^20; only the top slot (1,2,24) undetermined); canonical C determined for k < 24 and EQUAL to ur_gauge's P(0,z)/P_11(0,z) (84 slots); 1171 slots of P^(n), n <= 22, agree. Both frames. VERIFIED |
| `rf_gi.py` -> `rf_gi_ab_*.log` (K = 20, mod u^30) | ur_gauge P (read-only), certified by DIRECT SUBSTITUTION into the gauge equation built from own b, bt (2940 coefficients, stored and canonical) | Theorem G for the canonical AND the stored gauge (P^s Y^C = c(z) Y^C_Q, 1214 slots each; 127 positive-z slots non-zero for the stored one); (G*) canonical (417) and stored-twisted (557) pass; CONTROLS stored-untwisted and canonical-with-c fail (425 of 557); Prop S for both normalisations, z^0..z^20 mod u^29; (R1)-(R3); Lemma Q valuations n <= 35; p_(0,0) formula. Both frames. VERIFIED |
| same | float diagnosis with the z-POWER in (-1)**m | old bound: 107 of 202 (frame 0,0) / 115 of 197 float failures, max 8.9e-16; sharp bound: 357 of 560 / 361 of 555; the exact versions: 0 failures. VERIFIED |
| `rf_actual.py` -> `rf_actual_ab_*.log` | (A) EXACT over Q: referee table `G_fourier_NZ16_K12.json` vs canonical gauge hatted and Fourier-transformed here (dictionary Xhat(m,m') = (1/5) eps_a eps_b [Ghat(u,0)^-1 Ghat]_ab, a = m-1, b = -m'-1, derived in s.2); (B) HatRows mod p with M_0 from the LG-MODEL recursion | (A) all 397 components equal, same support; control: stored normalisation differs in 169. (B) 2250 graded slots (k <= 24, N < 18, including those below u^mu(b,k)), every off-grade entry zero, p89 and p61, frames (0,0) and (1/3,-2/7): all equal; model rows = referee rows at u^0 through zhat^24; control: stored differs in 603/607 slots from (5,0,0,4). VERIFIED |
| `rf_hat.py` -> `rf_hat_ab_*.log` (zhat <= 16, u < 12) | FLAT frame, exact Q, no zeta, no EdgeLeaf: E' = Psi0^T Ehat Psi0, Ehat = [eta^-1 - Ghat(-psi)^T eta^-1 Ghat(-psi')]/(psi+psi'); a-priori degree bound (5N+4)/4 for completeness | hatted Ghat symplectic; hatted (G) on 1275 slots; (E) 880 coefficients zero; controls: no (-1)^m 60 non-zero, stored normalisation 50 non-zero; Lemma L with random covectors 120 identities; (T) through u^11 with leading coefficient -1/30 IN BOTH FRAMES; (T) control: the stored normalisation fails from u^5 (= c_5 u^4 times O(u), as predicted). VERIFIED |
| `rf_misc.py` -> `rf_misc.log` | term-level checks | Lemma V (i), (ii), model system, Y^C(0) = -e_1; closed form of phi m <= 80; 2 c_10 - c_5^2 value; pi_8 canonical = stored + c_5 pi_3. VERIFIED |

## 4. Attempts to refute; issues found

1. Leaf bookkeeping "frame dependence" (HAT s.4, s.5; worker claim on (T)): with gamma_t + w_1 = gamma_g = 7/10 in every conifold
   frame (UNIFORM_REDUCTION s.5), p_(1,4)(0)/2 = 2/3 - gamma_t - w_1 = -1/30 and -(g-1)(2 gamma_t + 2 w_1 - 4/3) c_g = -(g-1) c_g/15 in EVERY
   frame (rf_hat: leading coefficient -1/30 in frame (1/3,-2/7)). The formulas are right; "frame dependence" and "-u^2/30 in frame
   (0,0)" understate them. UNIFORM_REDUCTION s.6.2 already says the bookkeeping is frame independent. Pi(T') itself is frame dependent
   from u^3 on (through 1/varpi_0). WEAKENS the wording only.
2. GAUGE_IDENTITY s.0 item 5 / HAT s.1: "the production R' of genus 2-5 IS the canonical gauge at the orders used" is tested here only on
   the REFEREE rows (zhat <= 24, u < 18); on the production rows the worker's own test covers zhat <= 6, u < 6. The statement follows
   from the recorded agreement of production and referee rows (NATURAL_FORM s.1 item 1, zhat <= 16, u <= 11), not from this task.
3. (G*) was checked without a negative control. Supplied here: the stored gauge fails untwisted (425 of 557) and passes twisted by
   its c(z); the canonical gauge fails with the wrong twist. The positive result stands.
4. The first gi_relations.py is not preserved (no copy in the worker's scratch either). The diagnosis is therefore a reconstruction.
   It is strongly supported: the reconstruction reproduces the four slots of the brief in the brief's order; the exact identity holds
   in all slots (proof + exact check); the float mechanism is reproduced by independent code (z-power sign). Confirmed with this caveat.
5. Hard rule "no person's name": an author-initials citation tag (two capitals and a year) appears in the GAUGE_IDENTITY.md s.6 heading and in the worker's
   claims; should be an arXiv number. ASCII, < 250 lines, stdlib, py -I -B, writes only in `gauge/` and its scratch: complied
   (no project file outside `allgenus/` modified after 15:00 except `genus5_cluster/RESUME_FRAME2_P61.md`, not attributable to this task).
6. Notation: GAUGE_IDENTITY s.4 writes `2 c_10 - c_5^2` with c = P^s_11(0,z) (so c_5 = -43397/807975000), while s.5, s.6 use
   c_5 = +43397/807975000. Harmless (the value quoted is right, rf_misc), but should be unified.
7. HAT s.5 dimension caveat: the bound d <= (5g-5)/2 < 3g-2 concerns a leaf on the genus-g vertex (vertex graph). For leaves on a vertex
   of genus h < g inside an edge graph the comparison is with 3h-2; the conclusion "no one-point law for all d" is unaffected.
8. Looked for and NOT found: tests that cannot fail (the central checks (G), (G*), S, actual = canonical, (E), (T) now all have a
   failing control; Lemma L and Lemma V are identities tested on random or exact data); "regular" read as "zero"; low-order
   evidence sold as all-order (all-order statements are proved on paper, the computations are labelled VERIFIED with windows); quoted
   theorems with unchecked hypotheses (Theorem U hypotheses hold for both systems in both frames: units except a simple pole of b_1,
   equal products, checked); tuned tests (the two test changes of the second version, integer sign and the exact (S') prediction,
   are corrections of demonstrable errors and are re-derived here).
9. Scope reminder (agrees with the worker): (G), (E), (T) are statements about the gauge and the edge/leaf factors only. They do not
   prove Conjecture N; the model-side singular structure (non-linear descendant extension, multi-point laws, genus-1 anomaly) is open.
   Everything is conditional on Theorem U; the new route of s.3 gives independent numerical support for its existence and uniqueness
   through z^23 in two frames. (H) and c^act = 1 beyond z^20 remain CONJECTURAL.

## 5. Claim-by-claim

| worker claim | verdict |
|---|---|
| (G*) failure = float bug, not truncation, not a wrong statement | confirmed (reconstruction caveat, s.4 item 4) |
| Theorem G for every regular gauge; canonical = unique one with P Y^C = Y^C_Q; c in Q[[z^5]] | confirmed (PROVED given Theorem U; VERIFIED incl. a non-canonical gauge) |
| Lemma Q (no log, exact valuations, [z^-1] = -t_c e_0) | confirmed |
| Proposition S and the stored defect (c(z)c(-z) - 1) eta^-1 from z^10 | confirmed |
| (G*), (R1)-(R3), u^3, u^4 relations | confirmed |
| diagonal p_(0,0) = tau^(-1/4)/varpi_0 | confirmed |
| actual gauge = canonical in zhat^k u^N, k <= 24, N < 18; c^act = 1 + O(z^25) | confirmed (two further routes, one exact over Q) |
| z^10, z^20 forced by symplecticity given (H) | confirmed |
| (H) and c^act = 1 to all orders | CONJECTURAL, correctly labelled |
| closed form of phi; phi(tau_m(Hhat^c)) formula | confirmed |
| Lemma L | confirmed |
| (E) | confirmed (edge-level statement only) |
| (T) | confirmed; "frame dependence" wording weakened (s.4 item 1) |
| conditional one-point law with kappa = -1 | correctly CONJECTURAL; frame-dependence wording weakened |
| no one-point law for all d | confirmed (vertex-genus caveat, s.4 item 7) |
| (E) not sufficient | confirmed |
| universality of Theorem G | confirmed (gauge identity part; needs Theorem U for the pair) |
| regression of the linear dressings (4 of 4; one converse exception) | confirmed (statuses re-read from nat_kregular_linear.log; dictionary confirmed by rf_actual (A)) |
