# Referee report on allgenus/model (MODEL_SINGULAR.md, MODEL_SINGULAR_2.md) -- 2026-11-13

Adversarial referee, AI only; no human has read this. Labels: PROVED / VERIFIED / CONJECTURAL / REFUTED. Everything in this folder
(`rf_*.py`, `*.log`) is new; no worker file was modified; no `ms_*` module is imported by any `rf_*` script. Project modules imported
read-only: `fq3_graphs` (stable graphs) and `natural_form/theory/ur_gauge` (Theorem-U recursion, flat Lemma-K check only). The worker's
logs were read as data (rf_compare.py parses `ms_singlist_full.log`). One worker script was re-run unchanged (`ms_edgeK.py`, item 7).

## 1. Verdict in brief

**Sound, with minor corrections.** Every numerical claim I re-derived by an independent route agrees exactly. Both refutations are real.
The Model Lemma (ML) is supported by all data, including two ranges beyond the worker's ((0,7), (1,6)); it stays CONJECTURAL. Lemma K is
correct and elementary. Theorem R is a correct CONDITIONAL reduction: Conjecture N follows from (ML) for all (h, n), (G) in canonical
normalisation, symplecticity, Theorem U and Corollary F (hence (H)). It is not a proof of Conjecture N. Corrections needed:
(a) Lemma K's stated hypothesis "phi_j -> 0" is false in the hatted frame where it is applied; the conclusion survives (item 6);
(b) "rules out ... any jet expression" holds only for genus-independent V (item 3); (c) "no resonance" should read "log-free" (item 5);
(d) hard rules: author names in MODEL_SINGULAR_2.md s.8, a 987-line log (s.4); (e) naming clash: "K" and "Lemma K" already denote the
frame K = (1, uH, H^2, H^3, H^4) and an older failed statement in UNIFORM_REDUCTION.md and the plan; rename (e.g. K_vc).

## 2. The referee's independent engine (rf_core.py, rf_gsum.py, rf_flat.py, rf_tables.py)

Different route from the worker's flat "charge" engine: the normalised CANONICAL frame with an explicit primitive 5th root of unity in
F_p (p = 2^61-1 and 2^89-1, both = 1 mod 5); colour sums by variable elimination; own psi-intersection (WK) numbers (string, dilaton, Virasoro recursion;
<tau_4>_2 = 1/1152, <tau_7>_3 = 1/82944, <tau_10>_4 = 1/7962624, <tau_2^3>_2 = 7/240); own R-matrix from a recursion DERIVED here, not the
stored rows. With g^a := N_a v^(-c_a/5) H^a the product is the group algebra of Z/5; the normalised idempotents give
(Psi^-1 H^a)_i = N_a zeta^(i(a+1)), sqrt(Delta_i) = zeta^(4i), U = -(5/4) diag(zeta^i); S = Psi R e^(U/z) solving z dS/dv = (H *) S with
R_k(v) = r_k v^(-4k/5) gives (r_(k+1))_ij = -(4/5)((V + k) r_k)_ij/(zeta^j - zeta^i) (i != j), (k+1)(r_(k+1))_ii = -Sum_j V_ij (r_(k+1))_ji,
V_ij = (1/5) Sum_a mu_a zeta^((i-j)(a+1)), mu_a = (c_a - 1)/4. R is symplectic through z^12 at both primes. The sign of z was NOT
assumed: rf_regress.log runs R(z) and R(-z) against record numbers only (F_2 = 1/1200, F_3 = 1/25200, <tau_1(1)>_1 = 5/24, <H^n>_h for
(0,3),(0,4),(1,1..3),(2,1..2),(3,1)); R(z) passes 11/11 at both primes, R(-z) fails 6. Exact values by CRT of the two primes and
rational reconstruction. rf_flat.py recomputes S, Pi_i, kappa(X,k) and K-membership (S symplectic on the log-free columns through
z^-15; Pi_i(kappa) = 0 for all i; PiC and the list of basis insertions in K agree with the worker's s.1). Tables: (0,3..7), (1,1..6),
(2,0..4), (3,0..3), (4,0..1) at 2^61-1; all but (1,6) also at 2^89-1 ((3,3) at 2^89-1 via rf_chunk.py, graph list saved once).

## 3. Item-by-item

1. **Engine (worker: VERIFIED). Confirmed.** rf_compare.log: the one-point table h = 1..4 agrees value for value with ms_onept.log; ALL
   986 entries of the worker's singular-correlator list (ms_singlist_full.log, h <= 4, n <= 4) agree exactly over Q (two primes),
   0 mismatches. This includes the 316 (3,3) entries, 16 of them of psi-deficit 7-9, which the worker could check
   only through string/dilaton (the table produced by the first attempt's process). F_4 = 1/180000, <H^n>_h reproduced (rf_misc.log).
2. **One-point law (worker: REFUTED as stated). Confirmed, exact over Q.** <tau_5(H)>_3 = 53/1209600 against PiC F_3' = -1/604800; at
   genus 4 the remainders of tau_5(H), tau_6(1), tau_7(H^4) are 19/100800000 x (1, -1, -1/2) = 19/100800000 x Pi_5(X). On K the law holds
   for every singular one-point function at h = 2, 3, 4.
3. **No genus-independent V (worker: REFUTED). Confirmed, with a scope correction.** Hessians reproduced exactly: t1(H) t2(H^3): -23/12
   (h=2) vs -13/48 (h=3); t2(H) t1(H^3): -11/6 vs -31/60; t3(H) t0(H^3): -3/2 vs -11/60 (all four correlators singular). The argument
   holds even for non-homogeneous V (the coefficient of the single monomial v^e must match at both genera). It excludes
   genus-independent V, including V built from jets of the genus-0 solution. It does NOT exclude genus-mixing structures (an
   hbar-dependent V, or second-difference/Toda corrections that couple F_h to lower genera); the claim "rules out ... any jet
   expression" overstates this. The D-directions remain unclassified, as the worker says.
4. **Model Lemma (worker: CONJECTURAL, VERIFIED on data). Confirmed and extended.** rf_ml.py enumerates EVERY K-generator tuple with
   v-exponent e < 0 (pruned recursion, all k), with referee generators and tables. Counts per (h, n) are identical to the worker's
   ((0,4) 6, (0,5) 37, (0,6) 152, (1,2) 3, (1,3) 19, (1,4) 83, (1,5) 313, (2,1) 2, (2,2) 17, (2,3) 90, (2,4) 376, (3,1) 4, (3,2) 41,
   (3,3) 246, (4,1) 6; total 1395). NONE is non-zero, at both primes (rf_ml_A.log; (3,3): rf_ml_C.log, earlier one-prime run rf_ml_B.log). (ML1)
   holds in every case. CONTROLS, absent from the worker's ms_ktest: the wrong-sign calibration kappa^- = Sum_j (-1)^j S_j X z^(k-j) gives
   476 non-zero singular correlators out of 1395, bare tau_k(X) gives 456: the vanishing is not automatic. NEW RANGE, never computed
   before (rf_ml2.py; tuples with at most 4 kappas, the tau_0(H) legs contracted first): (0,7) 191 singular tuples, 0 non-zero (both
   primes; control 51 non-zero); (1,6) 376, 0 non-zero (2^61-1 only; control 131). Status: CONJECTURAL for general (h, n).
5. **H^2 tower (worker: PROVED via arXiv:2211.12259 Thm 2). Accepted, with a precision.** I read Theorem 2 and its proof (pp. 7-10 of
   the arXiv PDF): at the special point, with the calibration (1.24), log D restricted to t^alpha_d = delta^(alpha,N) delta_(d,0)
   (alpha = 2..N) is Sum eps^(2g-2)/n! RHM_(g;k_1+1..k_n+1) prod t^1_(k_i)/(k_i+1)!; all descendants of t^N = a_N vanish except t^N_0,
   i.e. the slice K restricted to the t^1 tower. Euler: V - 5B + (n + B) = 2 - 2h gives e = V >= 1 (checked). The worker says the
   column S(z)H^2 has "no resonance"; in fact the homogeneous equation allows one constant: S_1 H^2 may contain c H (v^0), but then
   S_2 H^2 contains -c log v H^2. The column beta = 1 of (1.24) is a residue without logarithm, so c = 0 and the dictionary holds. The
   constants (metric x 20, z -> -4z, sign of the metric) do not affect exponents. PROVED given the quoted theorem. Data agree:
   181 non-zero pure kappa(H^2) correlators (h <= 3, n <= 3), all with e >= 1; <kappa(H^2,4)>_1 = v^3/32 (rf_misc.log).
6. **Lemma K (worker: PROVED). Confirmed; elementary.** The proof is correct (shorter: Pi_m is the m = 0 functional of w^-m phi(w),
   which satisfies the same hypothesis). Gap: the hypothesis "phi_j -> 0 u-adically" FAILS in the hatted frame (phi_(i,m) = O(u)
   exactly, no decay: GAUGE_IDENTITY_HAT s.2), which is where ms_edgeK applies it. The conclusion survives because every hatted edge is
   a polynomial in (psi, psi') at each u-order (Theorem D), so all sums are finite. In the flat frame phi_j = O(u^(j - floor(j/5))).
   **Independent check (rf_lemmaK_*.log; flat frame, exact over Q, P from ur_gauge to z^22, u^9; a different route from the worker's
   hatted mod-p rows):** with P_11(0,z) normalised to 1: symplectic through z^22; all 780 reliable positive-power coefficients of P Y^C
   vanish; Pi_m(first slot of the edge) = 0 on all 4515 reliable coefficients (m <= 6; 742 sums of >= 2 terms; a random covector is
   non-zero on 868); leaf: Pi_0(-z R^-1 1) = -t_c through u^8 (normalisation Pi_0(-z) = -v), Pi_1..Pi_6 = 0. Controls: the stored
   ur_gauge normalisation is not symplectic from z^10, violates (G) (102 coefficients) and the edge statement (627 non-zero), and moves
   Pi_0 of the leaf at u^5; P exp(z^5/7) is symplectic but violates (G) and the edge statement (627 non-zero) and gives the leaf a Pi_5
   component exactly (1/7) t_c, as the worker states. (A first run omitted the z^0 term Y_1 = -1 in the (G) sum and showed 129
   spurious non-zeros; corrected, documented in rf_lemmaK.py.)
7. **Edge-K on the actual R' tables (worker: VERIFIED). Confirmed.** Re-running ms_edgeK.py ref 61 1/3 -2/7 24 18 reproduces its log
   line by line (timings aside). With item 6 there are now two routes (flat/Q and hatted/F_p).
8. **Theorem R (worker: PROVED, conditional). Accepted as a conditional reduction.** Checked: (a) every vertex of a graph with an edge
   carries an edge slot, in K_0 by Lemma K; (b) K = C tau_0(H) + K_0 and Pi_0 is linear, so a vertex with >= 1 K_0 slot expands into
   correlators of K-generators with >= 1 kappa, regular by (ML2); the vertex graph gives F_g(v + Pi_0(T)) = F_g(t_c) + regular, with
   t_c - v = (13/90) u^2 + ...; (c) the kappa-decomposition of edge slots and leaves has power-series coefficients (S_j Y has the v-power
   x(j, Y) = (c_Y - c_(Y+j) + 4j)/5 >= 1 for j >= 1, Y != H); (d) u-adic convergence: in the flat frame each leaf of psi-degree >= 2 raises
   e by >= 3/5 and the slot degrees at a vertex with m leaves are at most 3h - 3 + n - m, so each u-order receives finitely many terms;
   unstable vertices vanish because T = O(z^2). I found no logical gap beyond item 6. The hypotheses carry the weight: (ML) for ALL (h, n)
   is open (proved only on the H^2 tower, and for the pure H-sector via the gap F_g, which the TR worker proves modulo a quoted
   theorem); (G), the canonical normalisation, (H), and "actual R' = canonical gauge" are VERIFIED to finite order only.
   Conjecture N therefore remains CONJECTURAL; Theorem R is a reduction, not a proof.
9. **Primary-sector lemma (worker: VERIFIED, consequence of ML). Correct**: w(a,0) = (4 - c_a)/5 > 0, so every primary correlator at
   h >= 1 has e < 0, and tau_0(X) = kappa(X,0) for X != H. Confirmed in the computed range.
10. **Delta_g(m) (worker: Delta_4(5) = 0, VERIFIED). Confirmed exactly** (rf_delta.log): leaf -1/8400000, loop 0, bridges +1/8400000;
   Delta_3(5) = -1/640, Delta_2(5) = 0; m = 1, 3, 7 vanish by the charge rule. I re-derived the formula (first order of exp(c z^m), m odd:
   leaf c z^(m+1) 1, edge c Sum_j (-1)^j w^(m-1-j) z^j eta^-1). So genus 4 cannot see the z^5 normalisation at u^-2.

## 4. Hard rules

* MODEL_SINGULAR_2.md s.8 names the four authors of arXiv:2211.12259: violates "no person's name in any file; cite by arXiv number".
  (Names of mathematical objects such as Barnes G or Toda are terms, not citations.)
* ms_singlist_full.log has 987 lines (rule: pieces of fewer than 250 lines); all other worker files comply and are ASCII.
* No write outside the worker's folders was found; caches are in the worker's scratch folder.
* Labels: acceptable. The summary sentence on Theorem R should keep "if (ML) holds for all (h, n)" attached.

## 5. Files

rf_core.py (WK, model data, R recursion), rf_gsum.py (graph sum), rf_flat.py (S, Pi, kappa, CRT), rf_tables.py (flat tables),
rf_regress.py/.log (sign fixed by record numbers), rf_run_amp.py, rf_chunk.py (amplitude runs), rf_compare.py/.log (one-point table,
the worker's singular list, Hessians; rf_compare_pre33.log is the run before (3,3) existed), rf_ml.py, rf_ml2.py, rf_ml_A.log,
rf_ml_B.log, rf_ml_C.log ((ML) with controls; new range), rf_lemmaK.py, rf_lemmaK_{canon,stored,c5}.log (flat Lemma K), rf_delta.py/.log,
rf_misc.py/.log. Amplitude caches: scratch folder allgenus_model_ref/ (JSON, mod p).
