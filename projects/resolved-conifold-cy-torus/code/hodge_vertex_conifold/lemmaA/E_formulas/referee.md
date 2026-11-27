# Route E referee report (round 1)

## confirmed

- NO TR: no part of the route, its engine, its Mathematica cross-checks or its data uses topological recursion, spectral curves, EO, loop equations, BKMP, the P2 dictionary or 'mod P2' results. The regularity consequence goes through Corollary R (leaf_structure_proof.md s.13.2), which states it is TR-free (not refereed here).
- Engine (computer-verified): the referee engine (RefWK/RefHodge, written from scratch) is kappa-free, keeps the dilaton leaves as explicit points, runs DVV on the smallest exponent and builds Omega directly at numeric tau. Its component polynomials agree coefficient by coefficient with every author engine file used: (h,r) = (1,1-3), (2,2-3), (3,2-3), (4,2-3), (5,1-7), (6,1-3), 2 566 213 coefficients, 0 mismatches. The referee engine also matches all Mathematica tables in lemmaA/data, which the author did not use: psi, mv and q tables, including the full s-polynomial, 2532 comparisons. It matches hvc_multiprobe_g4b (52 values), Faber, Getzler-Pandharipande, the lambda_g formula, FP Thm 2 and Dijkgraaf (g <= 6), and the one-part Marino-Vafa formula (derived from LLZ, normalisation fixed on the tables) through genus 7. The author's engine passes Harer-Zagier chi(M_{h,r}) (h <= 6), (G) and one-part MV (h <= 6); hodge_known.pl, consistency.pl, check61_engine.pl and verifyH.pl were re-run with 0 failures.
- Prop. W (forgetful recursion with K, proved): re-derived, correct (kappa_0 = 2h-3+n on the target); the implementation matches.
- (6.1) at h = 5, 6 (computer-verified): re-run, and checked independently for h = 2..7 with redundant grids. It holds identically at genus 7 too.
- M_irr(5) and M_irr(6) (computer-verified), independently: 44/44 and 96/96 irreducible points with every component zero; 22/22 and 42/42 controls with ALL components non-zero; Conjecture-M lattice (parts >= 0) zero at 716 points (r <= 6) and 668 points (r <= 5). Also M_irr(2,3,4) exhaustively.
- Genus 7 (computer-verified, partial in the report): confirmed and completed. All 194/194 points are zero in all components, so M_irr(7) holds. The author's ch_7 run uses HODGE_NCH=4 (checked in its environment) and agrees (190/194 at 19:21).
- Binomial-support equivalence, phi_a formula and the excess-Chern reformulation (proved): re-derived. NB: the excess-Chern form concerns the Hurwitz component only.
- Kappa reduction for r = 2 (proved): re-derived.
- J-only impossibility and top j-degree = WK numbers (proved): re-derived.
- SN => Lemma A (proved): correct as an implication (see issue 2 on its content).
- NR for r = 2 (computer-verified): an independent implementation is consistent with a UNIQUE solution at genus 1..7 (author: h <= 5), with 93 conditions per sigma-power at genus 7. A_{k,n'} = O(sigma^{2k-n'}) holds; A_11 = e^{-beta}-1 through L^7; the Gaussian leading A_kk hold for k <= 7.
- SN (conjectural) evidence, with the author's ORD=2 normalisation: reproduced and extended to genus 7 on 12 slices (r = 2: m = 2,3,4,5,6,9; r = 3: (2,2),(2,3),(3,3),(1,2),(0,3); r = 4: (2,2,2)), and to genus 5 on (2,2,2,2). The fitted C_s equal the author's stored slices (m = 2, 4, 5, 9). The slices are polynomial in m with degree <= 10 and vanish at m <= s (115 sequences, 6 redundant points each).
- Closed formula 0(e) for Xi^(h-1)(j1,j2) (conjectural): holds for h = 2..7; genus 7 is a parameter-free prediction.
- Mumford decomposition c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum n! ch_n (proved): re-derived. X_a = -X_b != 0 in the Lemma-A range was checked for h = 3, 4, 5, 6 (author: h = 3).
- Hurwitz extra layer (computer-verified): exhaustive for h = 2..5 (3, 7, 15, 30 zeros at D = 2h-1; non-zero at every point with D = 2h-2 or D = 2h). Extended to genus 6 (55/56 points, r <= 10) and genus 7 (89 points, r <= 8).
- Only the lambda_g component is divisible by a line J = l, r = 2 (computer-verified): confirmed for ALL complex l (gcd over Q) at h = 3, 6, 7.

## issues

### [1]

#### description

SN is stated with only [l^0]C_s = delta_{s0}. Every test imposes an extra normalisation that the summary never states: [l^{2h}]C_s = 0 for s > N_h = 2h-3+r (the 'ORD=2 convention'; the stored slices with m >= 4 were made this way). Without it the C_s are not unique, because the products for s >= N_h all equal 1. So 'the' C_s of the report are the normalised ones: their Q[sigma] coefficients, their polynomiality in m with zeros at m <= s, and the Borel kernel K_m. The number of test conditions per genus (h per power of sigma) is the same either way. The stronger, natural normalisation C_s = O(l^{2s}) is FALSE. It is the default ORD=1 of the author's own slice_fitJ.pl: 'ORD=1 perl slice_fitJ.pl 4 2' and '... 5 3' print 'genus 2: ANSATZ FAILS'. The referee's t_sn.pl with ORD1=1 fails the same way (sn_m4_ORD1.txt).

#### location

Summary 0(c) and the claim 'Conjecture SN'; slice_fitJ.pl (default ORD=1); slice_fitR.pl header; stored slice_m4..m16.txt

#### severity

minor

#### suggested_fix

State SN with '[l^{2h}] C_s = 0 whenever s > 2h-3+r', i.e. C_s = O(l^{2 ceil((s+3-r)/2)}). Say that every statement about C_s refers to these normalised series. Change the default of slice_fitJ.pl to ORD=2.

### [2]

#### description

The implication is correct, but at fixed genus the SN shape is EQUIVALENT to Lemma A on the slice. Take f(j) = sigma^h Xi_h(m,j) with deg_j f <= N+h and 1 <= M <= N. (a) Every SN term contains P(j) = prod_{i=0}^{N-M}(j-i), so P divides f. (b) Conversely, if P divides f, then f/P is always representable. The products E_s = prod_{l=s+1}^{M-1}(M+j-l) times b_i = [l^{2i}]B^j (degree exactly i) span every degree <= h+M-2. The top coefficient of f is forced by string: sigma^h <tau_0^{r-1} tau_{3h-3+r}>_h = sigma^h/(24^h h!). Hence, up to integrality in sigma, the only content of SN beyond Lemma A is the h-independence of the C_s, and the implication never uses it. SN, and likewise NR for r = 2, is a strengthening of Lemma A, not a reduction. Only an all-genus argument can profit from it; a genus-by-genus attack gains nothing.

#### location

Summary 0(c): Proposition 'SN => Lemma A', the headline 'THE STRUCTURE THAT MAKES LEMMA A EVIDENT', and s.4 'SN (or NR) is the precise identity'

#### severity

minor

#### suggested_fix

Record the fixed-genus equivalence. Reword the headline: SN/NR are all-genus strengthenings whose only extra content is the universality of the series.

### [3]

#### description

No derivation appears anywhere in the report or the scripts; PsiNewton.pl and fitNR.pl only fit. The referee confirmed the values numerically for k <= 7, and checked Dijkgraaf's formula itself against WK numbers for g <= 6. The status 'proved' is still unsupported.

#### location

Summary 0(d): 'leading A_kk = (-sigma l^2/8)^k/(k!(2k+1)) ... forced by Dijkgraaf's two-point function: proved given NR'

#### severity

minor

#### suggested_fix

Downgrade to 'computer-verified, k <= 7', or write out the derivation from the top j-degree part of NR and Dijkgraaf's formula.

### [4]

#### description

The genus-7 statement is partial: 161 of 194 points at submission, with the 33 points of r >= 6 open. The referee's independent engine completed it. All 194 irreducible points (r <= 11) have all eight components zero. All 76 computed controls (every control with r <= 11; only {2^12} not computed) are non-zero in every component. 1015 Conjecture-M lattice points (r <= 5) vanish. So M_irr(7) holds, and with Corollary R the regularity theorem holds through genus 7. The author's own run had 190/194 at 19:21, with no violation.

#### location

Summary 0(b) and the genus-7 claim (mirr_g7.txt)

#### severity

minor

#### suggested_fix

Update the genus-7 claim: M_irr(7) is computer-verified in full (referee data in E_formulas\referee, tally_final.txt).

### [5]

#### description

The sum is empty for M = 0 (r = 1, or m = (0)). Yet (6.1) and the dilaton slice are its s = 0 term. This is why slice_m0.txt reads 'ANSATZ FAILS', while the report says m = (0) reproduces the dilaton equation. HFun.pm inserts C_0(0) = S1S2 by hand.

#### location

SN formula Sum_{s=0}^{M-1}; slice_m0.txt

#### severity

cosmetic

#### suggested_fix

Use the range 0 <= s <= max(M-1, 0).

### [6]

#### description

The header announces check '(c) lambda_g component = b_h (J-1)...(J-N)', but the script never implements it. (G) is proved in the brief, and the referee verified it on the author's engine (h <= 6, at points in and out of range) and on the referee engine.

#### location

consistency.pl header

#### severity

cosmetic

#### suggested_fix

Implement (c), or remove it from the header.

### [7]

#### description

The script tests only integer l in [-3, 3h+2]. The referee's t_lines.pl computes a gcd over Q, which covers all complex l. For h = 3, 6, 7 no component other than lambda_g is divisible by any line J = l. The lambda_g component is divisible exactly for l = 1..2h-1. So the claim holds in the strong form.

#### location

explore_lines.pl / claim 0(f)(1)

#### severity

cosmetic

#### suggested_fix

State the claim for all complex l, citing the gcd computation.

### [8]

#### description

The X_b piece of the Mumford decomposition needs ch_n for odd n <= 2h-3 (ch_9 at genus 6), i.e. HODGE_NCH >= h-1. This is fine at h = 3 (the author's range). Extending to h >= 5 with the default NCH gives false failures; the referee hit 41 at h = 6 before adding ch_9.

#### location

comp_hm1.pl / Hodge.pm HODGE_NCH

#### severity

cosmetic

#### suggested_fix

Set NCH from h in comp_hm1.pl, or document the requirement.

### [9]

#### description

At 17:25 the author's mirr_points.pl 7 (PID 8048) and fitNR.pl 3 4 0 (PID 10187) were still running, and the runaway check_lemmaA.pl jobs (3721, 3931) had ended. The referee stopped only processes it had started itself. These were three t_moments.pl runs writing more than 1 GB into the Drive folder (their partial outputs deleted), duplicate or superseded lattice jobs, and its last background jobs at submission. Its large full-format files were moved to the session scratchpad: <user home>\AppData\Local\Temp\claude\G--My-Drive-Projects-Old-refGW\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\routeE_referee_bigfiles.

#### location

Housekeeping

#### severity

cosmetic

#### suggested_fix

None needed.


## notes

REFEREE REPORT ON ROUTE E (closed formulas and structure for r >= 2), the attack on Lemma A = M_irr(h). Adversarial referee, 2026-10-31.

I did not write referee.md, because subagents may not write report files; this text is the report. Scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\E_formulas\referee\ (index in section 6). No wolframscript was run, and no topological recursion of any kind was used.

0. VERDICT: SOUND, minor issues only.

No claim with status proved, proved-modulo-literature or computer-verified was refuted. I redid every computation with an engine I wrote from scratch, which shares no code and no algorithmic path with the author's (section 2), and pushed every check one genus further. New facts, all without TR:
 * (6.1) holds identically in (j, tau) at genus 7 (author: h <= 6).
 * The closed formula 0(e) for Xi^(h-1)_h(j1,j2) holds at genus 7. This is a parameter-free prediction: the author found the formula from h <= 5 and checked h = 6.
 * NR (r = 2) is consistent, with a unique solution, at genus 6 and 7 (author: h <= 5). A_11 = e^{-beta}-1 and the Gaussian leading coefficients of A_kk hold through L^7.
 * SN, with the author's normalisation (issue 1), holds through genus 7 on 12 slices: r = 2 with m = 2, 3, 4, 5, 6, 9; r = 3 with (2,2), (2,3), (3,3), (1,2), (0,3); r = 4 with (2,2,2). It holds through genus 5 for r = 5, (2,2,2,2); the genus-6/7 steps of that slice were not completed.
 * M_irr(h) holds, checked independently, for h = 2, 3, 4, 5, 6 AND 7. Every irreducible point vanishes in every component; at genus 7 that is 194/194 points, r <= 11. Every control point D = 2h-2 is non-zero in every component; at genus 7, 76 of the 77 controls were computed. Hence, with Corollary R (leaf_structure_proof.md s.13.2), the regularity theorem holds through genus 7, conditional only on what that corollary assumes.
 * The Hurwitz extra layer, the Mumford cancellation X_a = -X_b != 0, and "only lambda_g is divisible by a line J = l" all persist up to h = 7.

Issues:
 (1) The SN statement omits the normalisation used in all tests, and the natural normalisation C_s = O(l^{2s}) is FALSE. It fails at genus 2 for M >= 4, also with the author's own slice_fitJ.pl at its default.
 (2) SN => Lemma A is correct, but genus by genus the SN shape is EQUIVALENT to Lemma A on the slice. SN and NR are strengthenings of Lemma A, not reductions.
 (3) "Proved given NR" for the Gaussian moments has no argument.
 (7) The genus-7 claim is partial; it is completed here.
 (4)-(6), (8), (9) are cosmetic or housekeeping.
Nothing depends on topological recursion.

1. NO TOPOLOGICAL RECURSION — checked.
 * The author's Perl engine (WK.pm, Hodge.pm, HodgePt.pm) uses only DVV, Mumford's GRR formula, the Arbarello-Cornalba kappa calculus and the forgetful recursion.
 * The Mathematica cross-checks are Hodge-side: E_crosscheck_g5*.wls, and E_crosscheck_g4r2.wls via hvc_mp_header.m -> hvc_lib.m hodgeNormU4, which is the Mumford recursion. So are the data used (hvc_multiprobe*_out.txt, lemmaA/data).
 * A grep over all route files finds no spectral curve, EO, loop equations, BKMP, dictionary P2, "mod P2", status_tr*, referee1_tr*, hvc_tr_tau or hvc_moment_test.
 * proof_status_mixed_scale.md s.0 mentions TR confirmations of (6.1) at genus 5/6; the route does not use them.
 * The regularity consequence uses Corollary R, which states that it uses neither TR nor P2 (not refereed here).
 * My own engine is TR-free by construction.

2. MY INDEPENDENT ENGINE (referee/RefWK.pm, RefHodge.pm, RefParse.pm; exact Math::BigRat). It differs from the author's at every step:
 * WK numbers: string, dilaton, then DVV applied to the SMALLEST exponent (the author applies it to the largest).
 * Hodge integrals I(g, mu, e) = Int ch^mu psi^e use NO kappa classes at all. The kappa_l term of Mumford's formula becomes an extra marked point: Int kappa_l X = Int_{n+1} psi_{n+1}^{l+1} X. String and dilaton remove exponents 0 and 1; they are valid because ch(E) is pulled back.
 * lambda_k comes from Newton's identities, and the Bernoulli numbers are computed. Omega = prod_w Lambda^vee(w)/w^h is built at numeric tau straight from the brief's definition, then interpolated in v = -1/sigma using h+1 values, with a check at an (h+2)-th value.
 * The dilaton leaves stay explicit marked points with 1/m!, as in the brief's definition. So K, Prop. W and HodgePt are never used. Components come from tau-independent moments over leaf-exponent multisets.

Validation of my engine (all exact, 0 failures):
 * WK numbers: <tau_{3g-2}>_g, <tau_2 tau_3>_2, <tau_4^2>_3.
 * Lambda-class formulas: the lambda_g formula for g <= 6; Faber's formula, and lambda_{g-1}^3 = 2 x Faber, for g = 2..6; Getzler-Pandharipande for several psi-vectors, g = 2..6; Faber-Pandharipande Thm 2 (the formula the report cites), g <= 6, k <= 5; Dijkgraaf's two-point function, g <= 6.
 * The lemmaA/data tables, which the author did not use: psi_g1..g4, mv_g1..g4 (l <= 3), and q_g2..g4 (r <= 2, full s-polynomial), 2532 comparisons in total.
 * hvc_multiprobe_g4b: 52 comparisons.
 * The one-part Marino-Vafa formula, genus 1..7 (statement below).
 * Harer-Zagier chi(M_{h,r}) for h = 3..6, and (G) for h = 3..6.

Checks on the AUTHOR's engine with closed formulas the author did not use (0 failures): Harer-Zagier for h = 2..6, r = 0..3; (G) at points in and out of the range; one-part MV for h = 2..6.

The one-part MV formula. From LLZ Thm 1 with mu = (d): chi_hook = (-1)^k, kappa_hook = d(d-1-2k), V_hook = 1/([d][d-k-1]![k]!), and the q-binomial theorem give
   Int_{Mbar_{g,1}} Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau)/(1 - d psi) = d^{-2} [l^{2g}] prod_{a=1}^{d-1} f((d tau + a) l) / prod_{a=1}^{d} f(a l),
with f(y) = sin(y/2)/(y/2). The normalisation was fixed on mv_g1..4_l1 (96 comparisons).

Harer-Zagier: Xi^(0)_h(0^r) = chi(M_{h,r}) = B_{2h}/(2h(2h-2)) prod_{i<r}(2-2h-i).

Engine agreement: my component polynomials coincide with the author's files coefficient by coefficient for (h,r) = (1,1..3), (2,2..3), (3,2..3), (4,2..3), (5,1..7), (6,1..3). That is 2 566 213 coefficients with 0 mismatches, including engine_h5_r7.txt (2 016 585 coefficients).

3. CLAIM BY CLAIM

3.1 Engine (computer-verified): CONFIRMED.
 * The author's own checks pass on re-run: hodge_known.pl, consistency.pl (h <= 6), check61_engine.pl and verifyH.pl.
 * Normalisations: Xi^c = c^{-(3h-3+2r)} Xi^{c=1}, with tau = -x/(1+x) and sigma = -x/c^2.
 * Dilaton leaves carry 1/m!.
 * Omega is symmetric in the three weights, and every formula used is invariant under tau <-> -1-tau.

3.2 Prop. W (proved): CONFIRMED by re-derivation. kappa_b = pi^* kappa_b + psi_n^b and K = pi^*K (1 - psi_n). psi_i^a = pi^* psi_i^a + D_in pi^* psi_i^{a-1}, with psi_n D_in = 0, D_in D_jn = 0 and pi_* D_in = 1. pi_* psi_n^c = kappa_{c-1}, and kappa_0 = 2h-2+(n-1) on the target. The D-terms survive only when A is empty and a_n = 0.

3.3 (6.1) at h = 5, 6: CONFIRMED and EXTENDED to h = 7 (t_61.pl, with h+3 tau values and 3h+2 j values). The author's h = 6 check uses exactly 7 tau values, which is enough but leaves no redundancy.

3.4 M_irr(5) and M_irr(6): CONFIRMED independently.
 * Genus 5: 44/44 irreducible points zero; 22/22 controls with all components non-zero; 716 Conjecture-M lattice points (r <= 6) zero.
 * Genus 6: 96/96 zero; 42/42 controls; 668 lattice points (r <= 5) zero.
 * Genus 2, 3, 4 are also exhaustive: 1, 6 and 18 points.

3.5 Genus 7: COMPLETED. All 194 points vanish in all 8 components; the counts by r = 1..11 are 11, 30, 41, 38, 29, 19, 12, 7, 4, 2, 1. 76 controls are non-zero in all components. 1015 Conjecture-M lattice points (r <= 5) vanish. The author's run uses HODGE_NCH=4, checked in its environment; at 19:21 it had 190 of 194 points zero, with no violation.

3.6 Binomial support (proved), 3.7 kappa reduction for r = 2 (proved) and 3.8 J-only impossibility plus top degree (proved): CONFIRMED by re-derivation. NB: the excess-Chern reformulation 1.2 concerns only the Hurwitz component.

3.9 SN (conjectural): evidence REPRODUCED and EXTENDED, with my own fitting code, as listed in section 0. The fitted C_s equal the author's stored slices for m = 2, 4, 5 and 9 (64 + 76 + 111 coefficients for m = 4, 5, 9). The slices are polynomial in m of degree <= 10 and vanish at m <= s (115 sequences, 6 redundant points each), once C_0(0) = S1S2 is inserted.

3.10 SN => Lemma A (proved): CONFIRMED, subject to issues 2 and 4.

3.11 NR (computer-verified): CONFIRMED and EXTENDED with an independent implementation. The fit is consistent with a unique solution at every genus 1..7; at genus 7 the full grid gives 93 conditions per power of sigma. A_{k,n'} = O(sigma^{2k-n'}) holds. A_11 = e^{-beta} - 1 through L^7. The leading A_kk = (-sigma/8)^k/(k!(2k+1)) for k <= 7.

3.12 Closed formula 0(e) (conjectural): holds for h = 2..7.

3.13 Mumford decomposition (proved): CONFIRMED by re-derivation, and numerically for h = 3..6, with X_a = -X_b != 0 in the range. Pitfall: X_b needs ch_n for odd n up to 2h-3, i.e. ch_9 at genus 6.

3.14 Hurwitz extra layer (computer-verified): CONFIRMED exhaustively for h = 2..5. Extended to genus 6 (55 of the 56 points, r <= 10) and genus 7 (89 points, r <= 8).

3.15 Lines J = l (computer-verified): CONFIRMED for all complex l. The gcd over Q is constant for k < h; for k = h its roots are exactly l = 1..2h-1.

4. ISSUES (details are in the issues field)
 (1) The SN normalisation is missing from the statement, and the natural one is false.
 (2) SN is equivalent to Lemma A at fixed genus, so it is a strengthening, not a reduction.
 (3) "Proved given NR" is unsupported.
 (4) The SN range Sum_{s=0}^{M-1} is empty for M = 0.
 (5) consistency.pl advertises check (c) but does not implement it.
 (6) explore_lines.pl tests only integer l.
 (7) The genus-7 claim was partial; completed here.
 (8) Pitfall: comp_hm1.pl needs NCH >= h-1 beyond h = 4.
 (9) Housekeeping: I stopped only my own jobs. These were three runs writing more than 1 GB into Drive (partial outputs deleted), duplicate jobs, and my last background jobs at submission. My large files were moved to <user home>\AppData\Local\Temp\claude\G--My-Drive-Projects-Old-refGW\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\routeE_referee_bigfiles.

5. ARE THE OBSTRUCTIONS REAL? IS THERE A SIMPLER PATH?
 * "No mechanism for genus-independent coefficients": real, and sharper than stated. By issue 2, universality is the ONLY content of SN/NR beyond Lemma A, so only an all-genus argument can profit from them. They are credible targets now, verified at genus <= 7 on independent data. (6.1) and 0(e) are parameter-free statements verified one genus beyond their discovery range.
 * Excess intersection: this covers Xi^(0) only.
 * "(6.1) unproved for h >= 7": it is now verified at h = 7, but still unproved.
 * The Marino-Vafa p_1-shift remark: not checked by me. It is plausible, but should be quantified before it is used to discard the MV route.
 * Simpler path: none found. The equivalent forms of section 1 are genuinely equivalent. The natural rigorous route for r = 1 is (6.1) through MV with partitions (j, 1^k); my genus-7 check of the one-part MV formula shows those normalisations are under control. For r >= 2 the NR Gaussian structure points to a Wick/R-matrix mechanism.
 * Practical: my moment-based algorithm makes genus-by-genus verification cheap: genus 7, r = 8 took 45 minutes in Perl. M_irr(8) is feasible.

6. FILES (all in referee\)
 * Engine: RefWK.pm, RefHodge.pm, RefParse.pm.
 * Tests: t_basic.pl, t_fp.pl, t_dijkgraaf.pl, t_tables.pl, t_g4b.pl, t_mv1.pl, t_author_closed.pl, t_61.pl, t_e.pl, t_nr.pl, t_sn.pl, t_slicepoly.pl, t_mumford.pl, t_lines.pl, t_compact_check.pl.
 * Data: t_moments.pl produces ref_engine_*.txt and ref_mset_*.txt.
 * Lattice: t_lattice*.pl produces lat*.txt; t_tally.pl summarises them in tally_final.txt.
 * Logs: cmp_h5_r6.txt and cmp_h5_r7.txt (engine comparisons); sn_m*.txt (SN slices); t_nr_out.txt (NR).

## refuted


## route

E: closed formulas and structure for r >= 2 (attack on Lemma A = M_irr(h))

## verdict

sound

