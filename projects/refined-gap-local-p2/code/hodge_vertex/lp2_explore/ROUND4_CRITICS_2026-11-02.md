# Round 4 critics and referees on LP2_REGULARITY_PROOF_draft.md (workflow wf_bab37015-382, 2026-11-02)

Main session re-derived by hand the z^0 pairing identity in S4 of critic 1 and the W_i erratum; everything else is as reported.

## thmR: proved_modulo_cited
Open:
- Lemma A for 5 <= h <= g-1 is proved only modulo literature (2D Toda/mKP, MV formula; citation tightening listed in lemmaA/round2/round2_status.md gaps). Its statement must be restricted to non-empty R. Not re-examined here.
- Base chain write-up: R3 Corollary B step (i) (formal content of R2 Lemma 3 step 2, F3 applied three times plus F2) is still marked 'P mod writing out'. R2 Prop 5, Theorem 6 and Prop 8, and hence every factor input, rest on it and on BS C1-C3.
- Gap corollary: needs M_irr(g) (Lemma A at h = g); Prop 13.2_P2 (graph sum = F_g^cf via R2 Prop 11, with Cf and C5 = BS line ~4056); and the match of t_c (normalised t_c = Delta_c + O(Delta_c^2)) with BS's t_cf in cfasym. Only then is the single-vertex pole (-3)^{g-1} gap_g^cf the polar part of F_g^cf, which is BS rmk:cfgap.
- Documentation: state the string and lam_0 identities with their explicit E^{-2} factors, and the Lemma D0 proof, in the consolidated document. LSP s.13.1 and PSM s.3 still state M_irr(h) up to h = g. Graphs with an edge need only h <= g-1.

Referee: sound_with_minor_fixes. Minor:
- Notation clash. The report writes P^_i for the normalised row at colour i and also writes 'P^_0 = S(v_0,z)'. In R3 s.0 / Part II, S is component 0 of the colour-0 row, not the whole row. The report should say explicitly that every dilaton leaf Dil_i is built from component 0 (the P0 entry) of P^_i. That is also why the sign of c1 does not matter for the dilatons and Delta_i at the regular colours: Lemma G only covers the pairings, so it cannot be cited for the dilatons. Edges are covered by Lemma G; dilatons are covered by 'Dil uses P0 only', which should be stated.
- (Agen) The scaling Xi^{(E)}(R) = E^{-(3h-3+2r)} Xi^{(1)}(R;tau) holds only for the normalised class Omega_w = prod_a Lambda^vee_h(w_a)/w_a^h, the LSP Prop 13.2 convention. I re-derived it by weight counting: leaves have weight -1, T has weight +1, and the normalised Hodge factor has weight 0. With the unnormalised product the exponent would be 3-2r. The convention should be stated. It does not affect vanishing, only the constant in the corollary.
- (Agen) and the corollary use the conifold's CY relation c = 1 + x (LSP Appendix A). This relation is what makes (E,-e1,-e2)/E = (1,tau,-1-tau) agree with the conifold (c,-x,-1)/c, and what makes c = E at e2 = 1, so that 3e1e2E^3/(-xc^3) = -3. The report lists 'Scale: c -> E' as a free substitution. It should cite c = 1+x explicitly.
- (D0) F3 cites R2 Cor 9 (Form). That item is a statement about [z^k](v h_w S), polynomial in v of degree <= k+1. The fact actually used, that [z^k]S is a polynomial in v of degree <= k, follows directly from the definition of S (sum over l <= k). The citation should be corrected.
- The identification of the dilaton leaf as q_t = -psi S(v_0,-psi), with d = 1 as the shift, is read off the code (graphValue / graphs.m line 54). It is not derived from Givental's T(z) = z(1 - R^{-1}(z)1) together with symplecticity. This is harmless for Theorem R_P2, which is a statement about the decoration sum as defined. The theorem should say so, and the GW meaning of Dil should be deferred to the open Prop 13.2_P2.
- The checks.py string and lam_0 tests run on the genus-0 trivial CohFT, while the identities are used at h >= 2 with Omega_w. This is evidence only and not load-bearing, since the all-genus proofs are complete. It would be better to label the test as a control of the algebra, not of the Hodge statement.
- The statement of M_irr(h) should exclude the empty multiset explicitly: Xi_h(empty) = Z[h,0] is the gap and is non-zero. The report notes this restriction in the open steps, but the theorem statement and (iii) should carry it, with the empty R' handled by C_P2.

Referee: sound_with_minor_fixes. Minor:
- (Agen) The scaling exponent Xi^{(E)}(R) = E^{-(3h-3+2r)} Xi^{(1)}(R) depends on how Omega_w is normalised. With the unnormalised product Lambda^vee(w_1)Lambda^vee(w_2)Lambda^vee(w_3), rescaling psi -> psi/E gives E^{3-2r}: each lam leaf contributes E^{-1}, each T contributes E^{+1}, and the degree-d part of Omega contributes E^{3h-d}. The quoted exponent needs an extra E^{-3h}, i.e. division by (w_1 w_2 w_3)^h. Since it is a nonzero factor, the vanishing (iii) is unaffected. The consolidated text should still state which normalisation of Omega_w is used. The single-vertex corollary avoids the issue by setting e2 = 1, but that also needs the conifold scale c to equal E = 1+x (the CY condition). This should be said explicitly, not left inside 'at e2 = 1'.
- checks_run: the string/lam_0 test in genus 0 only exercises the (2-2h) term at h = 0. I re-ran both identities at genus 1 (D vs D+2 control). The consolidated document should cite a test at h >= 1, or simply rely on the proof, which is complete.
- (D0): the vanishing of the constant term for k >= 1 rests on Z-infinity, i.e. S_inf(v_0) = -t_c with t_c the conifold period normalised by C1cf = theta t_c. The report cites this correctly but should say that it is exactly the place where the normalisation of t_c (and hence rho(0) = 1 and the pole coefficient) enters.
- Normalisation clash: R2 Cor 9 uses lam_j = (1 - j psi/E)^{-1}, while R3 and the report use lam_j = 1/(E - j psi). Spans are unaffected, but the beta_j coefficients, and the statements sigma_0^(0) = E lam_1 and Sum beta_j/j = E delta_k0, are only correct in the R3 normalisation. State this once.
- (iv): the counting display '>= Sum (2-2h_v+D_v) + 2#{genus-0 colour-0} + 0' should say explicitly that regular-colour vertices, genus-1 colour-0 vertices and all non-allocated edge/dilaton orders contribute >= 0. It should also say that the allocated D_v at a rescaled vertex is bounded by the sum of the t_c-orders of the factors feeding it (split for rescaled-rescaled edges, wholly to v for loops). The argument is correct; the write-up is terse.
- Open items are correctly flagged and none is hidden: Lemma A for h >= 5 holds only modulo literature; the base chain R3 Cor B step (i) is 'P mod writing out'; the gap identification needs M_irr(g), Prop 13.2_P2 and matching t_c with t_cf. The outcome 'proved_modulo_cited' is appropriate.

## lemmaA: proved_modulo_cited
Open:
- Interface for the Theorem R_P2 port (item h, not Lemma A itself): lp2_lib.m uses unnormalised rows with vertex factor sqD^(2g-2+n+m+nd). Write the vertex as Delta_0^{h-1} d^n Fpot_h(q_t) with rows normalised so that the dilaton shift is exactly -psi (r_0[0] = 1). The same rho must normalise the leaves and give rho^2 Delta_0 = 3 e1 e2 E^3/t_c^2. Restate L1 for Dil_0 as rho q_t - E in span(lam_1..lam_{k+1}) for k >= 1 (from Z-infinity), with [t^0] rho q_t = q_0 exactly (round 3 item (d): 'to be restated').
- The literature for R3 has not been pinned (conifold critic gap, plan P2): theorem numbers for the 2D Toda / mKP equations of content-product tau functions, or the referee's three-line fermionic derivation written in, plus LLZ / OP theorem numbers. This carries over unchanged to local P2.
- R1 is not admissible under the local-P2 all-orders ruling (Q5), because its base h <= 7 and its zero count are computer-assisted. Cite R3 as the proof of record, and R2 (Hodge-side TR, allowed) as the backup.
- Gap corollary: statement (d) gives the vertex value (-3)^{g-1} gap_g. Its identification with the local-P2 GW polar part (BS rmk:cfgap) still needs Prop 13.2_P2 (local-P2 frame-B formula on T_F), which is open and outside this task.
- Mode convention: the transfer holds for mode B only, where the leaf scale E equals the weight w_1. In mode A (w -> -w) the leaf scale is not a weight, and Lemma A does not apply verbatim.

Referee: sound_with_minor_fixes. Minor:
- Step 2(iii) marks as proved that the Taylor coefficients in Theorem R_P2 step (ii) are exactly Xi^P2_h(R). That is the port of Theorem R, interface item (h) in ROUND3_STATUS, and it is open. It also has a rank-3 feature with no conifold counterpart: H2 gives sigma^(-1) = (0, E, -E^2), a t_c^{-1} term in components 1 and 2. The step should be labelled conditional on (h) and on the restated L1/F3 for Dil_0, or moved out of the Lemma A_P2 proof. Statements (a)-(c) do not depend on it.
- Open step 1 and Step 6 do not agree. Step 6 calls rho^2 Delta_0 = 3 e1 e2 E^3/t_c^2 proved, while the open step says the same rho still has to normalise the leaves. ROUND3_PROOFS_H2_edges_PF.md line 14 already defines rho = -E v_0/t_c, states rho^2 Delta_0 = 3e1e2E^3/t_c^2, and uses this same rho in H2 (sigma = rho P0hat(-psi)). That part of the open step is therefore discharged. What really remains open is the restatement of Dil_0 (item (d)) and rho(0) = 1, which depends on the cited normalisation t_c = Delta_c + O(Delta_c^2).
- Step 3 cites brief s.2(b) for Xi_h(R) in Q[sigma, 1/sigma]. The direct source is proof_status_mixed_scale.md s.1.1(e): Xi_h(R; tau) = Sum_{k=0}^{h} (-1/sigma)^k Xi^(k)_h(R) with tau-independent Xi^(k). That makes the specialisation at tau = -e1/E, the diagonal included, immediate. The remark about tau being transcendental is redundant.
- (b), extension to parts >= 1. This is correct only because the string identity of proof_status s.1.2 gives Xi(1^p) = prod_{i<p} (J_empty + i) Xi(empty) = 0, since J_empty = 0, and because parts equal to 0 are excluded: Xi(0) = (2-2h) Fpot_h(q_0) != 0. The statement should say explicitly that R contains no part 0. The scale-c form Xi(1,R) = J_R c^{-2} Xi(R) is consistent with the scale law (checked by degree count).
- (d) is stated only at e2 = 1. Give the homogeneous form in (e1, e2), or note that it follows from the scale law. Route B-c rests on conifold Prop 13.2, whose Teleman input was read only through summaries (conifold critic, gap list). Record (d) as conditional on that, or on R2 Cor. 7 (Hodge-side TR, allowed).
- Step 4's claim that lp2_lib.m's sq and etaF carry (e1 e2)^{-1} and E^{-3} was not checked by this referee. It is not load-bearing.

Referee: sound_with_minor_fixes. Minor:
- Step 2(iii) says the Taylor coefficients in Theorem R_P2 step (ii) are 'exactly' Xi^P2_h(R). That depends on the row-normalisation interface (r_0[0]=1, one rho normalising both the leaves and Delta_0, with [t^0] rho q_t = q_0 exactly), and the report itself lists that interface as open step 1. Statement (a) should say the identification holds once that interface is settled. Statement (b), a pure Hodge-integral identity, does not depend on it.
- R3 divides by Pi_2 = 2tau+1, which vanishes on the diagonal tau = -1/2 (e1 = e2), and by Pi_k, which vanishes at tau = -a/k. Step 4 is right that this does no harm: R3's induction hypotheses and conclusions are identities in Q(tau), and Xi_h lies a priori in Q[sigma, 1/sigma]. It should say explicitly that every inductive hypothesis is used as an identity in tau and never at a specialised tau. Otherwise the diagonal case looks like it needs a division by zero.
- The step that matches the MV formula to Xi_h (hbar = sigma x^2, q_k = x sigma Pi_k p_k, Route A normalisation (1.1)) is recorded in the conifold notes as 'normalisation exact', backed by 2170 + 996 data comparisons, and is listed under round2_status s.5.1 'not re-verified line by line'. The all-orders ruling Q5 requires a derivation. The report mentions this only as a 'shared input'; it belongs in open_steps next to the R3 literature-pinning item.
- In checks_run, 'stored local-P2 single-vertex pole' is compared with README / summary values that are total polar parts. They agree with the single-vertex graph at leading order: edge graphs are subleading, and the genus-4 per-graph line for graph {4} shows the gap followed by zeros. The wording should be made precise. This is evidence only and does not affect the proof.
- Step 5 uses m = 0 in vertexOptions. That is correct for F_g, because there are no external markings, but the report should say why. The Taylor 'legs' of Theorem R step (ii) appear only after the stable-graph sum at q_t has been formed, and they do not enter the stability count.

## Author notes
FILE: G:/My Drive/Projects/Old/refGW/hodge_vertex/lp2_explore/LP2_REGULARITY_PROOF_draft.md (about 15,700 words; not yet refereed as a whole). Chunks and assembly script: <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round4/author/

SECTIONS
- 0 Summary: Theorem LC_P2, Theorem R_P2, Corollary SV_P2, conditional Corollary G (not claimed), what is cited, what is not claimed.
- 1 Conventions: leaves, t_c and ρ (Lemma 1.2), rows defined by the closed form (Def 1.3), factors and code bookkeeping (Lemma 1.4), Hodge side, modes.
- 2 Rows on T_F: Theorem A (PF), F1-F5, Lemma 2.4, normalisation, continuity and closed form (Thm 2.7), code status, Cor 2.8 (Lemma CL), Cor 2.9, symplecticity (Prop 2.10), Givental-Teleman on T_F (Prop 2.11).
- 3 Singular colour: F1, F2, L1, Z-infinity, K = 1, Lemma D0, closure (Prop 3.1, not used).
- 4 Lemma H2 and Lemma G.
- 5 Residue lemma and Theorem E00.
- 6 Theorem M and Cor 6.1 (LC_P2 source table).
- 7 Hodge side: dictionary, string/λ_0 identities, why λ_0 does not occur, range Lemma 7.1, Lemma A_P2, proof of R_P2, single vertex and gap, audit.
- 8 Status table, remaining inputs, Prop 13.2_P2 and gap as final step.
- 9 BS v2 errata (14 items).
- Appendix: computer checks with paths.

ROUTINE-BUT-UNWRITTEN ITEMS NOW WRITTEN
- R2 Lemma 3 step 2 is Lemma 2.4, using explicit amplitudes (zθ+α_j)I[h] = I[q_j h + z∂_j h + zθh], with no "mod total derivatives". This discharges R3 Cor B step (i).
- Restated Dil_0 L1/F3 is Lemma D0.

UNSURE POINTS AND INCONSISTENCIES BETWEEN ROUNDS
1. lem:modrel (BS proof is a sketch). It is a real input of K = 1, not only of the gap. Prop 2.10 step 3 compares with eq:RHb at c1 = A, and that needs d2 = e_1(α*)E(y), i.e. E = (D-X)/2. Round 1 listed it (N4) only for the frame-B skeleton; rounds 2-3 folded it into C3 without comment. Flagged in s.8.3 as the weakest allowed input.
2. Lemma 1.2(c), t_c = t_cf, is new here and unrefereed. Round 4 listed it as open. I derive it from BS line 5145, eq:pfinf and the exponents (0,0) of L2b. I did not read eq:pfinf itself, so its operator being L2b∘θ rests on a round-2 referee's remark. A referee should read it.
3. R3 Theorem A proof says Delta_i^{-1/2} commutes with P. It does not, since it depends on y. Fixed: it is absorbed via vS.
4. The H2 and E00 inputs are as stated in R3, re-derived by hand (Steps 1-3, E00(i),(iv), M-b); no new computation was run.
5. Normalisation and citation clashes resolved:
   - λ_j in R2 vs R3;
   - D0's citation of (Form);
   - the code's v_i vs ours;
   - e_k(α) vs weights;
   - round-4 lemmaA open step 1 vs its Step 6 (ρ part discharged);
   - s.7.5(ii) now supplies the Taylor-coefficient claim marked conditional in lemmaA Step 2(iii);
   - LSP/PSM still assume M_irr up to h = g for edge graphs, while only h <= g-1 is needed.
6. (VI) source. Round 4 offered routes B-a and B-c. I record B-a (Hodge-side TR, allowed) as the proof, because B-c rests on Teleman read through summaries. BS prop:cfgraphsum (eq:betagiv, 5122) is related but not used, and I did not reconcile its graph sum with our single-vertex Z[g,0].
7. MV-to-Xi_h matching in Lemma A (Route A normalisation) is backed by data comparisons only. I list it as open as a written derivation under Q5; this is stricter than the conifold status "proved modulo literature". R1 is marked inadmissible.
8. The round-1 frame-B skeleton has no verifier, and it quotes the wrong W_i sign from BS. Whether the code's dilaton leaf -ψS is the Givental translation of R' is underived; the scoping note's "no analogue of Lemma 13.3" is unresolved. So Prop 13.2_P2 and Corollary G are not claimed.
9. Line numbers differ across rounds. I checked the labels I cite against main.tex by grep, but lines marked "~" (3382, 4056, 4370, 4245, 5320, 5399) are carried over or approximate.
10. Symbol reuse remains: P for the conjugated operator and the superpotential term vs rows; σ for the shift operator in Theorem A vs the rescaled row. All are local, with no logical effect, but they should be cleaned up before circulation.
11. A failed write during assembly emptied the draft (my own file from this session). I rebuilt it from the chunks. s.0 was retyped with three small wording changes: labels on Lemma A and its h = 2 status, a softened coefficient field in LC_P2, and "Lemma A(h), h <= g" in SV_P2.

## Critic 1
Read end to end. Theorem LC_P2 and the conditional Theorem R_P2 (with M_irr(h) as a hypothesis) are structurally complete modulo BS C1-C3; I found no step that fails.

What I re-derived by hand, all agreeing:
- Lemma 1.2; Lemmas 4.1-4.3; H2 Steps 1-4, including sigma_0^(1).
- Lemma R and E00(iv); (M-b), (M-d).
- Drop-out, and the eta values from eq:etaflat.
- Lemma 2.4, Thm 2.7(e) step 2, Prop 2.6 steps 2-3, Z-infinity.
- (S_P2), (H0), (7.2), Lemma 7.1, and steps (i)-(iv) of the counting.
- The g = 2 display.

I did not re-derive the F1-F4 calculus in full.

No TR leak found: Hodge-side TR enters only through route R2, used for (VI) and as the Lemma A backup.

The draft is not yet complete as an unconditional all-orders result:
- Lemma A still rests on finite-order or unwritten items (S1, S2).
- (VI) has two finite-order steps (S3). One is fixable in three lines, given in S3.
- lem:modrel is an open sub-step inside a cited input (S4), and I propose a way to remove it.
- Prop 13.2_P2 and the gap corollary are open (S5).
- s.8.2's 'Nothing else' and the [C] label on Lemma A_P2(b) overstate the position and should be corrected.

M1-M12 are local fixes. M1 and M2 touch load-bearing steps, but both are true with one added line.

The Appendix's claim that no proof step uses computer evidence is false for S1, S3 and for proof (ii) in S2, until those items are written.

Serious:
- S1. FINITE-ORDER INPUT IN LEMMA A: matching the MV series with Xi_h (s.7.4 input (iii), s.8.2 item 2). Under ruling Q5 the matching has no written derivation: only 'normalisation exact' plus 2170 + 996 data comparisons. The draft admits this ('open as a written derivation'), but s.8.1 still labels Lemma A_P2(b) [C], and s.0.3 says 'proved modulo literature'. Theorem R_P2 is formally safe because M_irr is a hypothesis. The unconditional all-genus claim is not safe. Fix: derive lemmaA/A_mv/report.md (1.1), line 179, from LLZ Thm 1 for all g. Then relabel, or until then mark (b) as [C + O].
- S2. PROOF OF RECORD OF LEMMA A NOT PINNED, AND THE INPUT LIST IS INCOMPLETE (s.7.4(b), s.8.2 'Nothing else'). R3 Theorem L has two proofs and the draft does not choose one. Proof (ii) needs every u_k-linear member P_k, whose form 'was fixed by direct checks for k <= 4' (R3 report.md line 206; round2_status s.5.2 item 2 says 'Add the derivation of P_k'). As written that is finite-order input for an all-k step. Proof (i) needs only P_2 (fermionic derivation in R3 referee 1 s.1.1) but also needs Route A Theorem B (I5, cut-and-join) or Cor. E. The draft's list (i)-(iii) omits four things: (I4) Route A Prop 2.1 / Cor 2.2 / Lemma 3.1; (I5); the empty-profile table of R3 Lemma 2.3, which R3 itself asks a referee to re-derive; and corrections (a)-(d) of round2_status s.1.1. Critic plan item P1 (one written all-genus proof, RESUME.md step 1) also has no file in lemmaA/ or lemmaA/round2/. Fix: declare proof (i) as the record, or write the bilinear-identity derivation of P_k for all k; extend the list in s.8.2; and write P1.
- S3. VERTEX IDENTITY (VI), FOR COROLLARY SV_P2 AND HENCE G ONLY (s.7.6). There are two finite-order steps. (a) The conversion of R2 Cor 7 into (VI) is 'checked exactly for h <= 12 [E]', and no proof is written. It has a three-line all-genus proof, which I did by hand. With S(y) = (y/2)/sinh(y/2) (R2 report line 121), e_i = eps_i/E and sigma = -e1 e2, one gets S(e1 l)S(e2 l) = e1 e2 l^2/(4 sinh sinh). So sigma^h Int K Omega = -e1 e2 beta_{2h-2}(e) = sigma beta_{2h-2}(e). Using (7.2) at r = 0 and beta_{2h-2}(eps) = E^{2h-2} beta_{2h-2}(e), this gives (-eps1 eps2 E^3)^{h-1} Z[h,0] = beta_{2h-2}(eps). This also needs the PSM s.1.1(a) identification of Int K Omega with Fpot(q_0) at c = 1. Insert this and demote the Perl check to evidence. (b) round2_status s.5.1 item 6 says the ABDKS/Eynard normalisation behind R2 was 'fixed against TR-free data for g <= 9'. s.8.2 item 3 hides this and the ABDKS / Eynard 2014 / BKMP-for-C^3 inputs under 'its own literature'. Either derive that normalisation for all g or list it as open. Route (B-c) carries the Teleman caveat instead.
- S4. WEAKEST CITED STEP IN THE MAIN CHAIN: lem:modrel E = (D-X)/2 (s.8.3 remark, s.9 item 11). The chain Prop 2.10 -> K = 1 -> H2 Step 3, E00, Theorem M depends on d2 = e_1 E(y). BS gives only a sketch ('Proof (sketch)' at BS 2697), so the label 'P mod C1-C3' has an [O] sub-step inside it. It can be removed on T_F using only (Symplecticity) and the Zagier-Zinger relation A^2 D = X. BS's rows of eq:RHb and Definition 1.3 at c1 = A differ only by (1/D)(d2 - e_1 E) R_1 in row 2. At z^0 the pairing of the closed form equals 1 algebraically: I checked eps E^3/9 - 6 eps Delta_c Y^3 - 3 eps E Y^2 = 3 eps E (Y+E/3)^2 using Delta_c Y^3 = Y^3 - (Y+E/3)^3. So (Symplecticity) at i = j, z^0 forces 2 eta^12 (1/D)(d2 - e_1 E) (R_1^{z^0})^2 = 0, hence d2 = e_1 E. This hand argument should be refereed. It needs eq:RHb and (Symplecticity) to be derived in BS independently of the second relation; eq:RHb uses the definitions in eq:ADEFXdef, so this looks fine, but it needs confirming.
- S5. GAP COROLLARY NOT PROVED (s.8.3, in scope as the final step). Proposition 13.2_P2 is [O] with five listed sub-items and no verifier, so Corollary G is conditional and G_g is not identified with any GW potential. The draft states this correctly. It is the one in-scope deliverable that is entirely missing. Theorems LC_P2 and R_P2 do not depend on it, but they are therefore statements about the decoration sum only.

Minor:
- M1. HYPOTHESIS NOT MET: Fpot outside its domain (s.1.5 vs s.7.1 (H0), s.7.5(i)-(ii), SV_P2). Fpot_h is defined only for q = -psi + O(psi^2), but it is applied at rho q_t, whose psi-coefficient is -rho = -1 - t_c/18 + ..., and at q_0 + delta with delta having a psi-linear part. The Taylor expansion and (7.1) with lambda = rho need Fpot defined t_c-adically. Fix: define it for q = -a psi + O(psi^2), a in 1 + t_c C[[t_c]], and prove (H0) by the binomial series Sum_p (1-a)^p (2h-2+m)_p/p! = a^{-(2h-2+m)} (checked by hand).
- M2. Theorem A's rigour paragraph (s.2.1) gives only the z-adic statement and says it 'replaces the V-adic framing', but Cor 2.9 (Z-infinity), which is load-bearing for L1 and H2, uses 'coefficientwise in v, as identities of rational functions in z'. That is the V-adic version. It is true because P e_n is a finite combination of e_{n-3..n} with e_n in V^n Q(z). State it explicitly.
- M3. False statement, not load-bearing (s.6 opening): v_{j,k} is called 'a vector of rational functions of w'. At j = 1, 2, v_j(0) != 0, so [t_c^k] S(v_j, -w) is an infinite series in w; Cor 2.8's rationality bullet uses v_0 = O(t_c) and covers colour 0 only. (M-a) to (M-c) work with formal series in w, so replace 'rational' by 'formal power series'.
- M4. Cor 2.8, first bullet, is false as written but its conclusion holds. Delta_0^{±1/2} is listed as analytic at the conifold, though it is proportional to v_0^{∓1}. The rows are called polynomial in v_i, though Lemma 4.1 has 1/v and 1/v^2, so they are Laurent polynomials. The conclusion 'meromorphic germ' is fine.
- M5. Lemma 1.4 and s.7.3 cite the wrong code. graphs.m:46 vertexOptions uses pure-psi weights wk with total degree = dim. The Hodge weights come from vertexOptionsH / hodgeNormU / usubOf in hv_lib3.m (hv_lib4.m hodgeNormU4 at genus 4), and s.1.5 cites only conifold hvc_lib.m. lp2_percolour.wls also uses its own loop, not graphValue. So '[P by reading]' does not currently establish that dec[[v,4]] = Int Omega_w psi^a / nd!. Re-read and re-cite. I confirmed only that hCoeffs[wB,-1] = h_w^{-1}.
- M6. Circular definition: s.0.1 says the decorations 'are those of s.1.4' and that (0.1) is 'as defined in s.1.4', but s.1.4 lists factors, not decorations. Make (0.1) the definition. State the range of a, the ordered d-tuples, and the finiteness of the sum by dimension (dim >= 2 nd; psi-degree <= dim because Omega_w has mixed degree). This finiteness is what makes 'analytic' meaningful in Theorem R_P2.
- M7. Non-factorised edges: (1.2) and step (i) are written for half-edge series f_r, but edges are tensors, and loops and rescaled-rescaled edges are not products. Say 'by multilinearity, leaf by leaf (E00(i))'.
- M8. Specialisation of eps: Prop 2.6 and Thm 2.7 fix generic alpha over BS's generic weights, while Theorem LC_P2 and s.7.7(f) claim every complex eps1 eps2 E != 0 with no argument. Add one line: each z-coefficient of Pi_ij - delta_ij is rational in eps with denominators only powers of eps1 eps2 E, so the identity specialises. This is the only place BS enters.
- M9. Status table s.8.1: H2, E00 and Theorem M are labelled 'P', but H2 Step 3, E00 (K = 1) and (M-d) rest on Prop 2.10, so they are P mod C1-C3. Lemma A_P2(b): see S1.
- M10. Dangling references. Theta(w, beta) in F5 and Thm 2.7 is never defined in the draft (s.1.3 defines Theta(x) only); cite BS about line 4278, whose conventions match. Round-1 N3 and N8 (s.8.3, s.7.7(c)) and the 'round-1 skeleton' exist only in journal wf_e928f771-9e3, not in ROUND1_STATUS, and the body gives no pointer. The round-4 reports (thmR, lemmaA) and their referees are named only in the header comment, with no journal path. R3's own Theorem L omits 'non-empty'; the critic's correction is cited, which is fine.
- M11. Notation traps. s.1.1 says psi = -z, but in (1.1) and Ehat(z, w) the first argument z IS psi (rows at -z). m_v means the dilaton count in s.0.1 and markings in Lemma 1.4. 'Lemma H2' and 'Proposition H2' are both used. Cor 2.9's clause 'within the class of Lemma 2.5(iii), this fixes c(z) = 1' is vestigial; c = 1 comes from Thm 2.7(e).
- M12. Lemma 1.2(c) is flagged in the draft as new and unrefereed. I checked it: eq:pfinf = L2b o theta, and BS 5145 gives t_cf = Delta_c + O(Delta_c^2). It looks right and is used only for Corollary G.

## Critic 2
Consistency verdict: no numerical or normalisation inconsistency. Every convention I tested in the draft matches the round 1-3 files and the code. All findings are statement, definition, citation or labelling issues, fixable in place, and none invalidates the proof chain. All round-2 and round-3 referee fixes in the status notes are carried over.

By-hand comparison against code and files, all agree:
- eta from etainvF: eta^00 = e1e2E^3/9, eta^11 = -3e1e2E, eta^12 = -3e1e2, eta^22 = 0; also BS line about 5397.
- Delta_i^{1/2} = sqDsym, including the z^0 term of cPP_z9.
- Y_i, v_i, W, and rows P1, P2 versus buildRmatrix.
- c2, d2 versus eq:RHb.
- h_w versus hCoeffs (hv_lib3.m line 75) and LSP (C8), with Upsilon at BS about 4281.
- t_c, C1cf versus tcSetup.
- eq:pfinf = L2b o theta.
- Telescoping of Ed, and exponents of sqD.
- Omega_w versus hodgeClass.
- Round-citation labels R2:x, R3:x exist and match.
- Sampled BS labels and lines match. eq:DEXrel is at 2702, inside the quoted range.
- Errata 1-2 re-derived (eq:Wadm gives +; determinant 243 W).

Spot checks, exact (py, fractions), scripts in <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round4/consistency/ (check1_row.py, check2.py, check3_hodge.py, check4.py, lser.py):

1. Definition 1.3 closed form versus stored rows (check1_row.py). E^2 sqrt(3) sq P_0 = h_w(z) v S(v,z) equals cluster_genus4/cPP_z9.m at z^0..z^9, at (7,3), (5,2), (4,9), at two rational L. Control h_w -> 1/h_w fails.

2. E00 and H2 versus stored edges (check2.py, at (7,3), (4,9)).
- Ehat_00 from (1.1) on the closed form equals the stored leaf table of edges_sing_out_7_3.txt, s^0..s^3, with zero residual at a sample (psi, psi').
- E00(iv): -(e1e2E^4/27) lam_1 lam_1' = -70000/9 lam_1 lam_1' at (7,3).
- E00(i) direct minus theorem = 0 for n = 2..5.
- rho = 1 + t_c/18 + ..., sigma^(-1) = (0,E,-E^2), sigma^(0) = (E lam_1, (E^2/9) lam_1, 0).
- K=1 residual zero to prec 10; H2(b) leaf supports hold for k=0..4.
- Z-infinity: S_inf = -t_c to t_c^12, and Sum beta_j/j = E delta_k0.
- Lemma 4.1 and Prop 3.1 hold on and off shell.
- By hand against stored: row 1 at s^-1 and s^0 (10/3, -20/9).

3. Theorem M versus stored rows (check4.py).
- Orthogonality z^0..z^6, and (M-b) for k=0..2, b=0..2, at both regular colours.
- chi_1, chi_2 at (7,3) equal the stored w^0..w^3 (for example chi_1 w^2: 350/27 + (490/27) omega).
- Stored Phat0_j(z;0) z^0..z^3 matches at both colours.
- Off-shell psi-linear coefficient E Delta_c C1cf^2 theta(f)/t_c gives 0, -1, 7/6, -403/486 at (7,3), Q=1/10; the first values match the stored ones.

4. Hodge side at genus 2 from known Hodge integrals (check3_hodge.py).
- (-e1e2E^3) Z[2,0] = beta_2 exactly at four weights; (-3) beta_2 reproduces the logged -83/6480 and -9773/327680.
- Xi({2}) = 0 (M_irr(2)), Xi({1}) = Xi({1,1}) = 0, and Xi(3) != 0.
- Xi(1,3) = 3/E^2 Xi(3); Xi(0) = -(2/E^2) Z; Xi(0,3) = 0; scale law (7.2) with diff 0.
- Controls: mode A (w -> -w) and wrong leaf scale both give Xi({2}) != 0, confirming s.1.6.

Scope limits: these are finite-order evidence only and were not used as proof steps. I did not re-referee the R3/Lemma A literature chain or Prop 13.2_P2. BS line references outside the sampled labels were not all checked.

Most important edits, in order: S1 (Fpot domain), S2-S3 (rationality overclaims), S6 (z versus psi), S4 (cite hv_lib3.m), S5 (table labels).

Serious:
- S1. Definition does not cover its use: LP2_REGULARITY_PROOF_draft.md s.1.5 defines Fpot_h(q) only for q = -psi + O(psi^2), where T = O(psi^2) makes the m-sum finite. s.7.1 (H0)/(7.1) and s.7.5(i)-(ii) and the proof of Cor SV_P2 evaluate it at rho q_t, whose psi-linear coefficient is -rho != -1 (rho = 1 + t_c/18 + ...), and at lambda q. There T has a psi-linear term and the sum over m is not finite, so 'Fpot_h(rho q_t)' is undefined as written. The same gap sits in LSP s.13.1 (leaf_structure_proof.md), which says only 'dilaton-shifted coordinates'. Fix: define Fpot_h on q = -u psi + O(psi^2), with u a unit in 1 + t_c C[[t_c]], by the Taylor series at q_0, Fpot_h(q_0 + delta) := Sum_r (1/r!) d^r Fpot_h(q_0)[delta^r], t_c-adically convergent. Prove (H0) for such lambda from the Euler identity dFpot(q)[q] = (2-2h) Fpot(q). At q_0 this identity is exactly (L0_P2) + (C_P2): E^2 Xi(0) - E^2 Xi(1) = (2-2h) Z. I checked those two exactly at genus 2.
- S2. False statement at s.6, opening paragraph: 'v_{j,k} := [t_c^k] P^^{(j)}, a vector of rational functions of w'. At colours 1, 2, v_j(0) != 0, so S(v_j, -w) is a non-terminating l-sum with poles at every z = -E/j. Its coefficients are only formal power series in w; my fixed-z evaluation there does not terminate and fails, while the z-series version passes. R3 Part IV has no such claim; the draft introduced it. Fix: 'a vector of formal power series in w'. (M-a)-(M-c) survive because they are taken [w^b]-wise.
- S3. Overclaim at s.2.6 Cor 2.8 (inherited from R3:Lemma CL), labelled [P]. (a) 'polynomial in v_i times functions analytic at the conifold (... Delta_i^{+-1/2})': Delta_0^{1/2} has a simple pole, and by Lemma 4.1 P^_1, P^_2 carry v^{-1}, so the correct word is Laurent polynomial or meromorphic factors. (b) 'each [t_c^k] coefficient is a rational function of z' holds only at colour 0 and only in frame B; h_w is not rational, and colours 1, 2 fail as in S2. s.3.1 and Step 3 of H2 use only the colour-0 frame-B case, so the proof is unaffected. Restate the corollary accordingly.
- S4. Lemma 1.4 cites the wrong code for the vertex, and it is labelled [P, by reading]. It (and s.7.3, and table s.8.1) cites graphs.m lines 46-54 vertexOptions. That function has weight wk[...] (line 51), the pure psi Witten-Kontsevich integral, with no Hodge class. The actual vertex is in hv_lib3.m: hodgeClass line 57, hodgeNormU line 72, vertexOptionsH line 119, prepareGraphH line 131, evalPointRows line 143 (plus hCoeffs and rowsTimes). lp2_percolour.wls line 41 loads hv_lib3.m and line 105 uses prepareGraphH. Also, lp2_percolour.wls does not call graphValue; it has its own product at about lines 116-118 of the same form. The content of Lemma 1.4 is correct (I re-derived the telescoping and the sqD count), and the omega, sign and normalisation agree with hv_lib3.m lines 2-3 and 54-60. But the reading proof must cite these lines, and s.1.5 should cite hv_lib3.m rather than only hvc_lib.m, which is the conifold file.
- S5. Status table s.8.1 misstates dependencies. 'Lemmas 4.1-4.3, Prop H2, Lemma G', 'Lemma R, Theorem E00' and 'Theorem M' are labelled plain P. But H2 Step 3 (constants of sigma_2), E00 and (M-d) use K=1 or orthogonality, i.e. Prop 2.10, which is P mod (C1)-(C3). Per the header labels these must read 'P mod C' (Lemmas 4.1-4.3, Lemma R and H2 Steps 1-2 stay P). The text of s.0.2 and Cor 6.1 is right; only the table is inconsistent, and it is the document's accounting.
- S6. Sign trap in the variable z. s.1.1 says 'Edge variable: psi = -z', i.e. z is the ROW argument. But (1.1), Lemma 1.4 and (M-d) use z as the EDGE variable, P_ij(z,w) = -<P^(i)(-z), P^(j)(-w)>, so that z = +psi: E00 has tau(psi) = P^(-psi), and edges_sing.pl has lam_j = 1/(E - j z). Read uniformly, (1.1) gives the wrong sign of the leaf argument. R3 Part IV has the same clash. Fix: write (1.1) as Ehat_ij(psi, w) with P_ij(psi, w) := -<P^(i)(-psi), P^(j)(-w)>, and state 'row argument z = -psi'.

Minor:
- Lemma 1.4 proof says 'm_v = 0 (no markings)' (the code's gr["m"] = legs), while s.0.1 and (0.1) use m_v for the number of dilaton insertions (as in LSP (13.1)). Rename the code's m to 'legs'.
- s.0.1 says the decorations are 'those of s.1.4', but s.1.4 only refers back to (0.1); nowhere are they actually defined. Define them once: a >= 0, ordered dilaton tuple d >= 2, weight 1/m_v!, degree constraint. LSP (13.1) can serve as the model.
- (M-e) asserts Ehat_21(z,w) = Ehat_12(w,z) but the proof given does not derive it; it needs orthogonality as in (M-d). Likewise Ehat_jj symmetry is never stated, yet T(Gamma,c) needs orientation independence of every edge. One line from Prop 2.10 suffices: P_ij(-w,w) = -delta_ij Delta_i.
- s.6 Remarks and the table 'chi_j not determined in closed form' are stale. (i) P^_2^(j)|_{t_c=0} = 0 identically, because Definition 1.3 has the factor Delta_c and the bracket is analytic at c1 = C1cf; the order-t_c^-1 part of orthogonality with sigma^(-1) forces the same thing. Verified z^0..z^3. So the 'rank-2 relation' is a one-liner, not conjectural, and the inconsistent third entry of the stored c-vector is explained because it multiplies a zero row. (ii) (M-b) at k=0 already gives chi_j(w) = -<(E, E^2/9, 0), v_{j,0}(-w)>/(E+w), which reproduces the stored chi_1, chi_2 at (7,3) exactly (see overall).
- s.3.3 clashes with s.3: 'v := v_0' throughout s.3-5, but Prop 3.1 uses v for the row vector. The QDE 'z d/dt_c v = (A_0 - L_i C1i) v' uses an undefined L_i (code's Y_i), while Lemma 4.2 defines L_1 := C1i Y, which would double-count C1i. Write (A_0 - C1i Y_i). I verified the row-2 identity itself exactly, on and off shell.
- Label collisions: F1-F5 (stationary phase, s.2.2) versus F1-F3 (Thm LC), with 'F2' used for both in s.2.3-2.7 and s.3.1, s.8.2. L1 (leaf statement) versus L_1 (Lemma 4.2). C1-C5 versus C1cf/C1i. Rename the stationary-phase items, e.g. SP1-SP5.
- Symbol overloads. T: T_F, T = q_0 + psi, T = z theta + q in Theorem A, and the T of Prop 2.10. sigma: shift operator in Theorem A versus the rescaled row versus sigma_tau. P: conjugated operator, y^{-1/3} e^{Sum x/3}, P_ij, rows. N: N(z), the number operator, N(z,w) in Lemma 1.4, N in Theorem M. Lambda(gamma) versus Lambda_m. E versus BS's E(y). Lemma numbering is also out of order (Lemma 2.1 comes after Cor 2.2). None causes an error at its local use.
- s.1.1: lambda_0 is admitted but (LAM1), (LAM2) are false or undefined at j = 0; restrict to j != 0.
- s.1.3: S^[m] needs a convention at l + 1 - m = -1 (l=0, m=2), where the upper limit is -1. The product must be read as 1; Lemma 4.1 passes exactly with this reading.
- s.1.4 Remark (signs) reasons with code factors. In the normalised form (0.1) the invariance is trivial: P^ = Delta^{1/2} P^B is invariant, and T(Gamma,c) contains no square root at all (S, C1i, d2 are single valued). Say so; s.7.7(b) would simplify.
- s.2.2-2.7 and the table say 'P mod F1-F5', although F1-F5 are written out and proved (R3 labelled them [P]; s.2.2 [S]). 'mod' suggests an open input and conflicts with s.8.2 'Nothing else'.
- s.0.3 first sentence restricts the R3 proof to h >= 5, then says R3 covers every h >= 2 (as s.7.4 and Theorem L say). Drop the restriction. s.7.4 and s.8.2 also omit the proved chain inputs L3 (dilaton calculus, A_mv s.2, s.3.1) and do not say which proof of Theorem L is of record. Proof (ii) needs all P_k and no Theorem B; proof (i) needs P_2 + Theorem B or Cor E (round2_status.md s.1.1, s.3 table L1-L5).
- s.0.3, s.7.6: 'BS Prop 3.1' is the numbering of arXiv:2410.00118, not a BS v2 label. The header's [C] rule requires label + line: prop:ResConifold, main.tex line 1708. See also prop:cfgraphsum, line 6366, which is not cited and is relevant to (VI).
- s.1.1 'the code uses the symbol v_i for 1+3t_i': this is not in lp2_lib.m, edges_sing.pl or hrow.pl (grep count 0). The wording comes from a round-3 referee (journal wf_74ed0422-50f). Say 'a referee script', or drop it.
- s.3.3 pole 3 e1 e2 E^2 (XX(0)+2/9)/t_c versus ROUND1_STATUS 'e1e2E^2(XX(0)+2/9)': not a contradiction. Round 1 is per 1/s, with t_c ~ 3s. My exact off-shell run confirms the draft's form (-630 at (7,3), Q=1/10; -9126/5 at (4,9)). Worth a parenthetical.

