# Lemma A round 2: critic (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## dependency_chain

Corollary R (genus g >= 2, frame-B Hodge-vertex graph sum of the resolved conifold x C^2 at t = 0): graphs with an edge are regular; the single-vertex graph has the single pole gap_g t^{2-2g}; F4 holds.

A. Graphs with an edge:
- A1 Theorem LC (F1-F3, L1-L3): proved modulo literature (Stirling remainder, DLMF); refereed.
- A2 Theorem R (i)-(iv): homogeneity, leaf expansion, string reduction, counting. Proved.
- A3 Theorem C and the string/dilaton reduction to M_irr: proved.
- A4 No contributing genus-g vertex in a genus-g graph with an edge: proved (dilaton degree >= 2 gives genus-0 vertices >= 3 half-edges; tree count).
- A5 Lemma A for 2 <= h <= g-1: proved modulo literature.

B. Single-vertex graph, any one of:
- (B-a) Theorem R (v) + Lemma A(g) + R2 Cor. 7 (gap formula = F4 at n = 0; Hodge-side TR, allowed). Needs no Prop. 13.2 and no (H).
- (B-b) Theorem R (v) + Lemma A(g), compared with F_g via Prop. 13.2 and (H).
- (B-c) Prop. 13.2 + (H) + A1-A5, without Lemma A(g).

C. Interpretation (graph sum = F_g):
- C1 Prop. 13.2 / Lemma 13.3: proved modulo Teleman, Givental, PPZ, Graber-Pandharipande, Coates-Givental, Mumford, Faber-Pandharipande. Refereed twice and revised; Teleman read only via summaries.
- C2 (H) = Brini-Schuler Prop. 3.1 (degeneration; no TR): proved in the literature.

Lemma A via R3, the proof of record:
- L1 MV formula (LLZ 2003 Thm 1; OP 2004): literature; normalisation exact.
- L2 2D Toda first-time-set (mKP) equations for tau_n = Z_MV(e^{n tau x} o p) (Jimbo-Miwa, Ueno-Takasaki, Takasaki, Orlov-Scherbin, Harnad-Orlov): literature; form derived from the fermionic bilinear identity and verified exactly.
- L3 Route A Prop. 2.1, Cor. 2.2, Lemma 3.1 (dilaton calculus): proved.
- L4 Route A Thm B (proof (i) only): proved modulo MV cut-and-join.
- L5 Theorem M, Cor. I_2, Thm M_P, Cor. I_P, Theorem L: proved modulo L1-L3.

Alternatives for Lemma A:
- R2: ABDKS Prop. 3.1 + Eynard 2014 (+ BKMP for C^3); non-inductive; Hodge-side TR.
- R1: Route A Thm D / Cor. E / Thm A / Thm F; Route D (I1) with (H), Prop. 13.2, Theorem LC, Thm 11.1; lambda_g; computer-assisted base h <= 7 and zero count G <= 1000.

No input uses TR for the conifold x C^2 theory.

## gaps

- Statement of Lemma A: brief s.1 and leaf_structure_proof.md s.13 literally include R = empty, but Xi_h(empty) = Int K Omega != 0. Restrict to non-empty R.
- R3 literature: the 2D Toda / mKP citations are loose (Okounkov math/0004128 s.2.7 is unconfirmed). Give theorem numbers (Jimbo-Miwa 1983, Ueno-Takasaki 1984, Orlov-Scherbin 2001, Harnad-Orlov 2015), or write in the referee's three-line fermionic derivation.
- MV formula and MV cut-and-join: give theorem numbers (LLZ Thm 1; OP). Their normalisation is already fixed exactly by the data.
- R2 literature: the proofs of ABDKS Prop. 3.1 / Thm 2.8 and of Eynard 2014 were not re-checked, only their statements and hypotheses. Eynard was taken from ABDKS (49)-(50), not compared with the original; (49) has an index slip; the Hodge form rests on BKMP for C^3 (Chen; Zhou). TR1 is an alternative route to P2.
- Prop. 13.2 (frame-B formula): Teleman was read only through tool summaries, and the CohFT axioms for the non-compact Z go through the twisted P^1 theory. It is needed only to identify the graph sum with F_g (routes B-b, B-c), not for the regularity of the graph terms.
- Theorem LC: its analytic input, Stirling's series with remainder (DLMF), is cited, not proved.
- R3 write-up: G^(g)_0 is defined only up to an additive constant; state that v >= 1 for empty-U blocks; state the single exponential sector; the odd-parity sampling claim needs n+3 points; fix the I_2 step (d)4 wording; add the derivation of P_k; correct the remark that genus 7 did not finish; separate the dependencies of proof (i).
- R2 write-up: non-empty R; the h = 0 term is a convention; Theorem 3 in multiplied-out form; the engine count is 1051, and the script counted no-data cases as OK; fix V5 and V20; the genus-7 table is not new data; literature keys.
- R1 write-up: the Cor. 7.2 constant; the justification of W = alpha; the (I4) factorial; the (I1) argument; the Lagrange step of s.7(b); label the zero count computer-assisted; scan labels; the scope of Thm 6.1. Its dependency list must include Prop. 13.2, Theorem LC and Theorems A-C of the conifold note.
- R4 write-up: the G = 3 claim is refuted (with Xi(2) = 0, no Hodge data is needed for h <= 20); K1' evidence; the cosection observation is heuristic; N1-N4 are ansatz failures in finite ranges; regenerate the lam1_rec outputs; write out H_2.
- Route A record: the TLR nullities 11, 11, 7 included c = 2. Correct values: 7, 6, 1 (window G) and 3, 0, 0 (window G-1).
- The three proofs are independent in mechanism but share the MV formula. R3 and R1 use it directly; R2 probably through BKMP for C^3 behind Eynard/ABDKS (to be confirmed).
- F4 at n = 0 in all genera without Prop. 13.2 and (H) rests on R2's Cor. 7, which uses TR on the Hodge side (allowed). A TR-free proof of (6.1) at j = 0 is still missing; the TR-free alternative is Prop. 13.2 + (H).
- Documentation is stale: leaf_structure_proof s.0.4, s.13.2, G1; proof_status s.0, s.2, s.3; conifold_regularity_note s.3a, s.4, s.7b; proof_M_residue s.5, s.6, s.9; brief s.2(g); also the memory and README. The R2 referee2.md has not been saved.
- No logical gap, no uniformity gap, and no TR on the conifold x C^2 theory was found in R1-R4.

## lemmaA_status

proved-modulo-literature-all-genera

## plan

- P1 (essential, low): write and referee a single all-genus proof of Lemma A, with R3 proof (ii) as the record and proof (i) / the Cor. E shortcut as variants, including all referee corrections and the non-empty statement. Restate Corollary R as proved modulo literature in every genus, and add s.13.4 (Lemma A only below the genus; F4 at n = 0 in all genera).
- P2 (essential, low-medium): pin the literature. MV (LLZ, OP); 2D Toda / content-product (Jimbo-Miwa, Ueno-Takasaki, Orlov-Scherbin, Harnad-Orlov), or include the fermionic derivation; ABDKS Prop. 3.1 / Thm 2.8 in the published version; Eynard 2014 and the ABDKS (49) slip; BKMP for C^3; Teleman read directly; Brini-Schuler in the published version.
- P3 (high, medium): integration critic across both round-2 workflows (R1-R4 and TR1-TR2) once TR1's second referee is in. Fix the presentation order: R3, then R2 or TR2 + TR1, then R1.
- P4 (medium, low): record the new all-genus results in the notes: (6.1), the gap formula and F4 at n = 0; Lemmas B and D; the vanishing direction of Conjecture M for all h and r; K1; T1 and the slice free-energy formula.
- P5 (optional, medium-hard): a TR-free proof of (6.1) and of the gap formula, from R3's M_P members used as a recursion off the range, or from E-operators / Orlov-Scherbin closed formulas.
- P6 (low): bookkeeping. Document fixes, the Route A TLR record, regeneration of R2/R4 outputs, saving R2's referee2.md, then the memory, README and tarball.
- P7 (optional, low): the Mathematica confirmations; exact Theorem F certificates for G = 8..20; DVV-free MV checks of M_irr(5..7).
- P8 (research, optional): multi-leaf closed forms (SN, NR) from R2's Theorem 1; merge formulas from R1's slice formula; R4's components k >= 2.

## proved

- Lemma A (precise form): for every h >= 2 and every NON-EMPTY multiset R of integers >= 2 with D(R) <= 2h-3, Xi_h(R) = 0 identically in tau. Proved modulo literature, three ways. (1) R3 Theorem L: MV formula + 2D Toda first-time-set equations + Route A dilaton calculus, plus Thm B or Cor. E for proof (i). (2) R2 Cor. 4, non-inductive: ABDKS Prop. 3.1 + Eynard 2014; TR on the Hodge CohFT only. (3) R1: KP + conifold polar relations + forcing, with base h <= 7 computer-assisted. Each has two 'sound' referees.
- R3 Theorem M (2-insertion relation), for all G, eps, R. It is the w^{n+r} dilaton-point coefficient of the first mKP equation 2 d_{p2}(F-F1) + d_{p1}^2(F+F1) + (d_{p1}(F-F1))^2 = 0, with F = log Z_MV and F1(p) = F(e^{tau x} o p). Proved modulo literature; 426 exact Q[tau] identities at genus 0-7.
- R3 Corollary I_2: under (H1) Lemma A in genera < G and (H2) Xi_G(T) = 0 for D(T) < D, Xi_G(R+2) = 0 for every R with D(R) = D-1. Proved.
- R3 Theorem M_P / Corollary I_P: each member P_k = D_{t_k} + p_k(-D~) gives Xi_G(R+k) = 0 under (H1),(H2), with target coefficient -(k-1) sigma^G tau Pi_k delta. Proved modulo the 2D Toda bilinear identity.
- R3: the lambda_g limit of Theorem M. Proved.
- R2 Theorem 1: closed n-point formula for the TR correlators of C_tau. Proved modulo ABDKS arXiv:2312.16950v3 Prop. 3.1; in formal form it holds for every tau with sigma != 0.
- R2 Theorem 2 (Lemma D): [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n. Proved modulo Theorem 1.
- R2 Theorem 3 (P2 = Lemma B): omega_{g,n} = sigma^{1-g} Sum Int K Omega prod psi^k d xi_k, for all g, n and sigma != 0. Proved modulo Eynard 2014 as stated in ABDKS (49)-(50).
- R2 Corollary 6: formula (6.1) for all h >= 1, all j and all tau. Proved modulo literature.
- R2 Corollary 7 (gap formula): sigma^h Int K Omega = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l). Critic's consequence: (-x c^3)^{h-1} Z[h,0] = (2h-3)! P_h(x,1)/x = gap_h, i.e. F4 at n = 0 in every genus; the conversion was checked exactly for h <= 12.
- R2: Theorem 5 (divided differences), the moment lemma, the two-point difference form and F1. Proved.
- R1 T1: closed form of the polar-order relations, relation_k = kappa_k Sum_L y_L (-1)^{k-L}[(1+tau)binom((L+1)tau,k-L) - tau binom((L+1)(1+tau)+k-L-1,k-L)]. Proved modulo Route A Thm D and Theorem 11.1.
- R1: flow theorem and free energy on the slice P_1 = 0. Proved modulo Thm D, as formal identities graded by level.
- R1: branch-point identity, and the curve identity FC with C = W e^{tau^2 t}/H_1. Proved.
- R1: special-framing Lemma 8.1 and forcing Theorem 9.1. Proved modulo Thm D, T1 and Route D (I1).
- R1 zero count: Z1(D) + Z2(G,D) >= G-1 for G >= 8 and odd D. Proved, computer-assisted; for G >= 8 the Hurwitz input is not needed.
- R4 Theorem K1: Xi^(1)_h(j) = Int c(T(-log)) lambda_1/(1 - j psi) = 0 for h >= 2 and 2 <= j <= 2h-2. K1' also holds. Both proved modulo Route B Thms B, C and Mumford's ch_1 formula. Props. M1-M3 and Claim E proved.
- Brief s.2(f) parity: T_D(G) = Pi_{D+1} kappa_D Xi_G(D+1) with kappa_D = [y^D](tau y + log((e^y-1)/y)), i.e. kappa_1 = tau + 1/2 and kappa_D = B_D/(D D!), zero for odd D >= 3. Re-derived; exact polynomial identity in tau for D <= 14.
- Brief s.3: Lemma A for all h <=> (OL) for every G >= 3. The induction uses exactly the hypotheses of Cor. E and Theorem D5. Confirmed.
- Corollary R at genus g needs Lemma A only for h <= g-1, given Prop. 13.2 and (H). Confirmed: dilaton leaves have psi-degree >= 2, so genus-0 vertices need >= 3 half-edges, so a genus-g graph with an edge has no genus-g vertex.
- Corollary R in every genus g >= 2 (graphs with an edge regular at t = 0; single-vertex pole gap_g t^{2-2g}; F4). Proved modulo literature via the chain of report s.3.

## requested_checks

- Mathematica (new, seconds): wolframscript -file "G:/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/round2/critic/critic_gap_check.wls". Expected: h = 2..10 all True, then 'ALL TRUE: True'. This is F4 at n = 0 in all genera from R2's Cor. 7, checked against rc_lib.m's Pg / cyPolar.
- Mathematica (minutes), from the hodge_vertex_conifold folder: wolframscript -file lemmaA/round2/R3_toda/R3_thmM_check.wls. Expected: every line 0, 'ALL ZERO: True', and a non-zero control. This is Theorem M symbolic in tau on the q-tables, G = 2-4.
- Mathematica (optional, heavy): wolframscript -file "G:/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/round2/R1_rigidity/R1_checks.wls". Expected: (C1), (C2) True; (C3) True three times; (C4) True at tau = 2/7, -3/5, 3; (C5) {tuples, 0}; (C6) 0 failures for 8 <= G <= 60 and 4 for G = 4..7. Already duplicated by two Perl referee implementations.
- Mathematica (optional, TR line): extend HVC\hvc_tr_tau.wls at tau = 2, 3 to genus 4-5, n <= 3, and compare with 'perl xy_npoint.pl TAU G N JMAX' in round2\R2_onepoint. Also: TR2_recursion\tr2_closed_check.wls <gmax> (expect all differences 0), and TR1's request (hvc_tr_tau.wls at a non-integer tau, genus 3-4).
- No longer needed: R4's requested Xi^(1)_10 check against C1, now implied by (6.1) (R2) together with the exact C1 = -[sigma^{h-1}](6.1) identity for h <= 14.
- Main session (bookkeeping): save R2's second referee report as round2\R2_onepoint\referee2.md, and collect TR1's second referee before the integration critic.

## summary

Lemma A (M_irr(h)) now holds in all genera, proved modulo literature by three independent arguments. Each has two referee reports with the verdict "sound" and no major issue.
- R3 is the proof of record. It uses the modified-KP equations of the 2D Toda lattice of the Marino-Vafa (MV) content-product tau function at the dilaton point. It uses no TR, no conifold input and no computer input.
- R2 uses the x-y swap for log TR on the Hodge curve C_tau (ABDKS) plus Eynard's kappa-class formula. It needs no induction, and it applies TR only to the Hodge CohFT, which is allowed.
- R1 uses KP along the leaf curve, forcing at special framings and the conifold polar relations. Its cases h <= 7 are computer-assisted.

Hence Corollary R holds in every genus, modulo literature. The inputs are Theorem LC, Theorem R (i)-(v), Theorem C and Lemma A, plus Prop. 13.2 and Brini-Schuler (H) for the identification with F_g. The dependency chain is in s.3 of the report.

I checked the three main-session claims myself, and all hold.
- **Parity, brief s.2(f).** Re-derived by Lagrange inversion. It is exact as a polynomial identity in tau for D <= 14 (147 evaluations).
- **Reduction of Lemma A to (OL), brief s.3.** The induction runs with exactly the hypotheses of Cor. E and Theorem D5.
- **"Corollary R at genus g needs Lemma A only up to g-1".** Correct, but the justification needs fixing. Dilaton leaves have psi-degree >= 2, so genus-0 vertices need at least 3 half-edges, and a tree containing a genus-g vertex is then impossible. "Unstable" is the wrong reason.

One new observation: R2's gap formula (Corollary 7) is the vertex identity F4 at n = 0 in every genus. I checked the conversion exactly up to h = 12. So the single-vertex statement of Corollary R follows from Lemma A(g) without Prop. 13.2 and (H).

No logical or uniformity gap was found. What remains:
- pinning down the literature, mainly the 2D Toda citations, ABDKS/Eynard and Teleman;
- the referees' write-up corrections, including the precise statement of Lemma A, which must say non-empty R;
- updating the documents;
- optional Mathematica confirmations.

My scripts are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\critic\.

## report_text

# Lemma A, round 2: completeness critic for routes R1-R4, and status (2026-11-01)

The main session saves this as HVC\lemmaA\round2\round2_status.md, where HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold.

**Scope**
- Routes R1_rigidity, R2_onepoint, R3_toda and R4_geometry (workflow wf_e42146b3-114).
- Their four reports and seven referee reports: two each for R1, R2 and R3, one for R4. R2's second report was read from the structured results.
- TR1_dictionary and TR2_recursion belong to the other round-2 workflow. I did not assess them. They appear only where they supply alternative inputs.

**What I read**
- round2\brief.md, the four report.md files and the saved referee reports.
- lemmaA\parity\README.txt and kappa_parity.pl.
- leaf_proof\leaf_structure_proof.md: s.0, s.13.1-13.3 and s.14.
- A_mv\report.md s.2-3: Prop. 2.1, Cor. 2.2, Theorems A, B, D, Cor. E and Theorem F.
- D_localization\report.md s.5: Lemmas 5.2, 5.3 and Theorem D5.
- conifold_regularity_note.md: F4 and s.3a (Theorem C).
- proof_status_mixed_scale.md s.1.2 and proof_M_residue.md s.5.
- hvc_lib.m (vertexOptionsH4) and rc_lib.m (Pg, cyPolar).
- lemmaA\RESUME.md, status_2026-10-31.md and round1_status.md.

**My own computations**
All are exact Perl (Math::BigRat); no wolframscript was run. Files are in HVC\lemmaA\round2\critic\:
- kappa_poly_check.pl, with output kappa_poly_check_D14.out: the parity identity as an identity of polynomials in tau, D <= 14. 147 evaluations, 0 mismatches.
- gap_identity_check.pl, with output gap_identity_check_h12.out: the conversion of R2's Cor. 7 into F4 at n = 0. 44 comparisons, 0 mismatches.
- critic_gap_check.wls: the same check against rc_lib.m's own definitions, for the main session to run.
- I also re-ran kappa_parity.pl; every value was OK.

**Rules**
- From brief s.0: no TR for the resolved conifold x C^2 theory.
- Allowed: TR for the Hodge CohFT (curve C_tau, BKMP for C^3, P2), and KP and 2D Toda for the MV series.
- The user's reply relayed with this run ("yes, It's an allowed input") agrees with these rules.
- So R2's Hodge-side TR and R3's 2D Toda input are both admissible.
- No route applies TR to the conifold x C^2 theory.

## 0. Bottom line

1. **Lemma A = M_irr(h) holds for every h >= 2, proved modulo literature.** There are three independent arguments. Each has two "sound" referee reports and no major issue.
   - **R3, the proof of record.** It uses the modified-KP (first-time-set) equations of the 2D Toda lattice tau_n = Z_MV(e^{n tau x} o p) at the dilaton point. It needs no TR, no conifold input, no computer input and no data.
   - **R2.** It uses the x-y swap for log TR on the Hodge curve C_tau (ABDKS) and Eynard's kappa-class formula, which gives P2. It is non-inductive, and its TR is applied only to the Hodge CohFT.
   - **R1.** It integrates KP (Route A Thm D) along the leaf curve, forces zeros at special framings and uses the conifold polar relations. The cases h <= 7 are computer-assisted.
2. **Corollary R holds in every genus, modulo the literature in s.3.** This is the all-genus regularity of the frame-B graph sum of the resolved conifold x C^2 at the conifold point.
3. **New all-genus results survived the referees.**
   - R3: Theorem M.
   - R2: the closed n-point formula on C_tau, Lemma D (integer tau), P2, (6.1) for all h, and the gap formula.
   - R1: T1, the flow theorem and the forcing theorem.
   - R4: Theorem K1.
   - R2's gap formula is the vertex identity F4 at n = 0 in every genus (s.4.4).
4. **The three main-session claims are correct**, with one corrected justification (s.4.3):
   - the parity reduction of brief s.2(f);
   - the reduction of Lemma A to (OL) in brief s.3;
   - "Corollary R at genus g needs Lemma A only up to genus g-1".
5. **No logical gap.** What remains is pinning literature statements, the referees' write-up corrections, documentation, and optional confirmations (s.5-s.7).

## 1. Rigorously proved claims (only those that survived the referees)

**Conventions** (brief s.1)
- c = 1, w = (1, tau, -1-tau), sigma = tau(1+tau), lam_j = 1/(1 - j psi), T = -psi^2/(1-psi).
- Pi_j = binom(j(1+tau)-1, j-1), D(R) = Sum(j_i - 1), J = Sum j_i, r = #R, delta_g(T) = D(T) - 2g + 2.
- Xi_h(R) = Int_{Mbar_{h,r}} K Omega prod 1/(1 - j_i psi_i), with K = exp(-Sum kappa_k/k).

**Precise statement of Lemma A.** This corrects brief s.1 and leaf_structure_proof.md s.13. Both literally include R = empty, but Xi_h(empty) = Int K Omega != 0.
- For every h >= 2 and every NON-EMPTY multiset R of integers >= 2 with D(R) <= 2h-3, Xi_h(R) = 0 identically in tau.
- By Xi_h(1 + R') = J_{R'} Xi_h(R') and Xi_h(1^p) = 0, this is equivalent to the same statement for non-empty R with parts >= 1.

### 1.1 R3 (2D Toda; both referees sound)

**Theorem M (2-insertion relation).** For all G >= 0, eps in {0,1}, n = 2G + eps and every multiset R of integers >= 2:

    -2 Pi_2 Sum_{2g+k=n+1, k>=1} sigma^g tau^k (delta_g(R+2))_k X_g(R+2)/k!
      + [eps=0] 2 sigma^{G+1} B_G(R)
      + Sum_{2g+k=n, k>=1} sigma^{g+1} tau^k C^g_k(R)/k!
      + Sum_{S subset R} Sum_{Sum(2g_i+k_i)=n+2} Phi^{g1}_{k1}(R_S) Phi^{g2}_{k2}(R-R_S) = 0.

- B, C, Phi and the empty-profile values are as in R3 s.4.
- The relation is the w^{n+r} coefficient, with w -> 0 the dilaton point, of the q_R coefficient at x^n of the first mKP equation
  2 d_{p2}(F - F1) + d_{p1}^2(F + F1) + (d_{p1}(F - F1))^2 = 0, where F = log Z_MV and F1(p) = F(e^{tau x} o p).
- **Status: proved modulo literature.** The inputs are:
  - the MV formula;
  - the first-time-set equation P_2 of the 2D Toda lattice at p' = q^rho. Referee 1 derived it from the fermionic bilinear identity on |0>(x)|1>, which fixes tau_1 = Z(e^{+tau x} o p);
  - Route A's dilaton-point calculus (Prop. 2.1, Cor. 2.2, Lemma 3.1; proved and refereed in round 1).
- **Checks.**
  - Referee 1: 426 exact identities in Q[tau] at genus 0-7, including genus 7 at R = (13), (14).
  - Referee 2: 592 transcription checks.
  - The author's 374 + 182 + 60 + 14.
- **Corrections to carry over:**
  - (a) G^(g)_0 lies in w^{2g-2}Q[w] only up to an additive constant.
  - (b) Empty-profile functions with v = 0 are more singular than their x-order. This is harmless, because an empty-U block consists of d_{p1}'s only, so v = |B| >= 1.
  - (c) State explicitly that every term lies in the single sector e^{u(J+2)}, because log Z_MV contains all unstable genus-0 terms.
  - (d) Odd n has tau-degree <= 2G+3, so the 2G+3-point sampling proves only the even relations. Only those are used, and the referees' exact Q[tau] checks cover both parities.

**Corollary I_2.** Assume:
- (H1) Lemma A holds in genera 2..G-1;
- (H2) Xi_G(T) = 0 for every non-empty T with parts >= 2 and D(T) < D, where 1 <= D <= 2G-3.

Then Xi_G(R+2) = 0 for every R (parts >= 2, possibly empty) with D(R) = D-1.
- **Status: proved.**
- Mechanically confirmed for G = 2..12 with all R (6259 relations), and for G = 13..21 with r <= 3.
- Wording fix: in step (d)4 the partner is Phi^0_1(R - T_i).

**Theorem M_P / Corollary I_P.** Let k >= 2 and take the member P_k = D_{t_k} + p_k(-D~) of the 2D Toda bilinear identity at equal second times (charges (0,1)).
- Under (H1) and (H2) it gives Xi_G(R+k) = 0 for every R with D(R) = D-k+1.
- The target coefficient is -(k-1) sigma^G tau Pi_k delta_G(R+k).
- **Status: proved modulo the 2D Toda bilinear identity.**
- Exact checks:
  - in Q[tau]: P_3 for G <= 5, P_4 for G <= 4, P_5 at G = 2;
  - in bilinear form on Z_MV: P_2 through x^120, and P_3..P_6 through x^60.

**Theorem L (Lemma A for every h >= 2).** Induction on h, and inside genus h on D = 1..2h-3, with exactly (H1) and (H2).
- Proof (ii): apply Corollary I_{P_k} at R = T - k, for any part k of T.
- Proof (i): Corollary I_2, then Route A Theorem B (cut-and-join), by induction on r_T and then min T.
- Shortcut variant: Corollary I_2 at R = (D), then Route A Cor. E at T = (D,2). Its factor (D+1) Pi_{D+1}/(Pi_D Pi_2) is non-zero.
- **Status: proved modulo literature.**
  - Proof (ii) needs only the MV formula, the 2D Toda identity (all P_k) and the dilaton calculus.
  - Proof (i) needs P_2, plus Theorem B (MV cut-and-join, Liu-Liu-Zhou (LLZ)) or Cor. E (KP).
- There is no circularity: Theorems M, M_P and B carry no hypothesis, and the base case h = 2 needs none.

**Other R3 claims.**
- The lambda_g limit of Theorem M is proved.
- The first mKP equation for Z_MV follows from the cited 2D Toda structure. It is also verified exactly, on Z_MV and on generic content-product tau functions of degree <= 9.
- Route A Theorem B is re-verified: exactly at genus 2-6, and mod p at genus 7-10.
- Line (2) (the bosonised Eq and the genus-0 CY Toda equation) is computer-verified. The "no termwise analysis" statement is a remark.
- Line (3) is conjectural and not needed.

### 1.2 R2 (x-y swap on C_tau; both referees sound)

**Theorem 1 (closed n-point formula on C_tau).** For integer tau >= 1 and 2g-2+n > 0:

    [x^mu]F_{g,n} = ((-1)^{n-1}/prod mu_i) [hb^{2g-2+n}] [v_1^-1]...[v_n^-1] ( prod_i Q_{mu_i}(v_i) Sum_{n-cycles s} prod_i 1/(v_i - v_{s(i)} + (mu_i+mu_{s(i)}) hb/2) ),
    Q_mu(v) = prod_{m in M(e2 mu)} (v + 1 + m hb/e2) / prod_{m in M(tau mu)} (v + m hb/tau).

- **Status: proved modulo ABDKS**, arXiv:2312.16950v3, Prop. 3.1 (eq. (17)) with Rem. 2.6, 3.2 and 3.4.
- The hypotheses were checked by both referees:
  - dX has one simple zero, at v = tau;
  - v = 0 and v = -1 are LogTR-vital, and v = oo is not;
  - the dual TR is trivial;
  - the standing assumptions of Thm 2.8 hold.
- Referee 1: in a formal form, the formula holds for every tau with sigma != 0.
- Checks:
  - 3058 values (mod p) against the TR-free engine, for n <= 7, g <= 7, integer and non-integer tau;
  - 403 + 316 values against exact Eynard-Orantin (EO) recursion, at tau = 1 and tau = -1/2.

**Theorem 2 (Lemma D).** [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n.
- This holds for integer tau >= 1, and in the formal form for every tau != 0, -1.
- **Status: proved modulo Theorem 1.**
- Wording: R_mu is a polynomial in c, Laurent in hb with hb-exponent >= -1, homogeneous of degree |mu|.

**Theorem 3 (P2 = Lemma B).** For all g, n and every tau with sigma != 0:
- omega_{g,n} = sigma^{1-g} Sum_k Int K Omega prod psi_i^{k_i} d xi_{k_i};
- equivalently, [x^mu]F_{g,n} = sigma^{1-g} prod Pi_{mu_i} Xi_g(mu).

**Status: proved modulo Eynard, CNTP 8 (2014)**, in the form of ABDKS (49)-(50).
- ABDKS (49) has an index slip, xi_k -> xi_{k+2}. R2 uses the consistent form.
- ABDKS attribute the plain Hodge form (44) to BKMP for C^3 (Bouchard-Marino; Zhou; Chen).
- The normalisation was confirmed against TR-free data at genus 2-9.
- TR1 offers an alternative proof via Zhou's BKMP theorem and EO symmetry.

**Moment lemma.** Proved.

**Corollary 4 (Lemma A for all h, non-inductive).** For every h >= 2 and every non-empty R of positive integers with D(R) <= 2h-3, Xi_h(R) = 0 identically in tau.
- Proof: Theorems 2 and 3 at tau = 1, 2, 3, ..., then polynomiality in tau.
- **Status: proved modulo literature.**

**Theorem 5 (divided differences).** Proved.

**Corollary 6 ((6.1) for all h >= 1, all j, all tau).**

    sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) [lambda^{2h}] S(e1 lambda)^{1+j e2} S(e2 lambda)^{1+j e1}.

**Status: proved modulo literature.** The h = 0 term of the generating-function form is a convention.

**Corollary 7 (gap formula).** For h >= 2:
- sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [l^{2h}] S(e1 l)S(e2 l);
- Xi_h(1+R') = |R'| Xi_h(R');
- Xi_h(0) = (2-2h) Int K Omega.

**Status: proved modulo literature.** It gives F4 at n = 0 in every genus (s.4.4).

**Also proved:** the two-point difference form and F1 (Murnaghan-Nakayama).

**Corrections**
- The engine total is 1051, not 1072. V11 had no data, and the summary lines counted no-data as OK.
- V5 was one sigma-value short of an identity. Referee 2 fixed this componentwise at genus 7-9.
- V20 only tests the consistency of the code.
- The genus-7 table is identical to the existing RefHodge table, so it is not new data.
- M_irr(7) is complete.
- Literature keys: Eynard's paper is [Eyn14] in ABDKS's list.

### 1.3 R1 (KP rigidity; both referees sound)

**T1 (closed form of the polar-order relations).** Assume Lemma A in genera < G. Then for 1 <= k <= 2G-3,

    relation_k := [t^k] Sum_T prod b_t Xi_G(T)/|Aut T|
               = kappa_k Sum_{L<=k} y_L (-1)^{k-L} [(1+tau) binom((L+1)tau, k-L) - tau binom((L+1)(1+tau)+k-L-1, k-L)],

where y_L = Pi_{L+1} Xi_G(L+1). On the leaf curve, b_j = (-1)^j (Pi_j^2/j) t^{j-1} H_j/H_1^j with H_j = 2F1(-j tau, -j tau; j; 1-e^{-t}).

**Status: proved modulo Route A Thm D** (applied at every point of the slice P_1 = 0) **and Theorem 11.1.**

**Flow theorem (Thm 5.1) and free energy on the slice (Thm 6.1).** Proved modulo Thm D, as identities of formal series graded by level. For the true F_G they hold at levels <= 2G-3.

**Branch-point identity (Lemma 4.1).** Proved.

**Curve identity FC (Thm 7.1).** Proved, with C = W e^{tau^2 t}/H_1.
- W = alpha by the sum rule (Theorem 11.1 / Cor. 11.2(a)).
- In Cor. 7.2 the constant is log H_1 - tau^2 t.
- The Lagrange-inversion step of s.7(b) should be written out.

**Special-framing Lemma 8.1.** Proved modulo Thm D (R empty).

**Forcing Theorem 9.1.** Under the hypotheses of (OL), Xi_G(D+1) vanishes at tau0 = -k0/q whenever:
- gcd(k0,q) = 1,
- k0 <= D+1,
- 2k0 < q,
- q >= D+2,
- 3(q-k0) < 2G.

**Status: proved modulo T1, Lemma 8.1 and Route D's (I1).**

**Zero count (Prop. 10.1).** Z1(D) + Z2(G,D) >= G-1 for all G >= 8 and odd 3 <= D <= 2G-3.
- **Status: proved, computer-assisted.** The count is exact for G <= 1000, the bounds are analytic beyond, and (B1) holds for all G and D.
- The margin is >= 1, so for G >= 8 the Hurwitz input brief 2(b) is not needed. The lambda_g input 2(a) is.

**(OL) for G >= 8, hence Lemma A for all h (the third proof).** **Status: proved modulo the cited inputs and the computer-verified base h <= 7.**
- Cited inputs:
  - Route A: Thm D, Cor. E, Thm A;
  - Route D (I1), with (H), Prop. 13.2, Theorem LC and Thm 11.1, Lemma 5.3 and (V);
  - lambda_g.
- Base cases:
  - Theorem F at G = 5, 6, 7, exact over Q;
  - h = 3, 4 from exact data or rank certificates.

**Computer-verified**
- The global nullities in window G-1.
- Route A's "TLR nullity 11, 11, 7" included c = 2. The correct values are 7, 6, 1 in window G and 3, 0, 0 in window G-1.
- Pointwise-scan claims are evidence only: mod-p forcing can give false positives.

**Refuted:** FC1 (constant alpha), the deformed-curve ansatz, and the guessed zeta.

### 1.4 R4 (geometry; one referee, sound)

**Theorem K1.** For every h >= 2 and every integer 2 <= j <= 2h-2,

    Xi^(1)_h(j) = Int_{Mbar_{h,1}} c(T(-log)) lambda_1/(1 - j psi) = 0.

- **Status: proved modulo Route B Thms B, C** (ELSV; Okounkov math/0004128 eq. (10); Hurwitz genus 0) **and Mumford's ch_1 formula.**
- It uses no KP, no (6.1) and no TR. Checked for h <= 22.

**K1' (boundary value).** Proved modulo the same inputs. The evidence is hbdry_rec.pl or ref_rec_A.pl, not delta_bdry.pl.

**The supporting results**
- M1: proved.
- M2: proved modulo Route B.
- M3: proved modulo M2 and Claim E.
- Claim E: proved.
- The K2 reduction and the reformulations are proved, with B - 1 in the Sheffer form.
- H_2 is open and now moot.

**Computer-verified**
- C1 for h <= 22. It equals minus the sigma^{h-1} part of (6.1) for h <= 14, so at those h it now follows from R2.
- M3' at j = 2g.
- KP plus Theorem H plus K1 plus Xi(2) = 0. By the referee's correction, genus 3 and 4 then close without Hodge data, so Lemma A for h <= 20 needs no tables.
- (6.1) at genus 7-9.

**Downgraded**
- The cosection observation is heuristic.
- N1-N4 are failures of specific ansaetze in finite ranges.

**Refuted:** "at G = 3 one unknown remains, Xi^(2)_3(4)".

### 1.5 Data beyond genus 7 (computer-verified; not needed)

- Genus 8 (r <= 3, 1737 zeros) and genus 9 (r = 1), TR-free engine (R2 referee 2).
- 16 top-level profiles at genus 7-10, MV Schur side mod p (R3 referee 2).

## 2. Status of Lemma A for all h

**Proved modulo literature in all genera.**

| | R3 (record) | R2 | R1 |
|---|---|---|---|
| mechanism | 2D Toda / mKP members at the dilaton point | log-TR x-y swap on C_tau, then P2 | KP along the leaf curve; forcing at tau = -k0/q; conifold polar relations |
| literature | MV formula; 2D Toda for content-product tau functions | ABDKS Prop. 3.1 / Thm 2.8; Eynard 2014 (and BKMP for C^3) | MV; KP for MV and Takasaki-Takebe; Brini-Schuler (H); Teleman etc. (Prop. 13.2) |
| project inputs | Route A dilaton calculus (proof (i) also Thm B or Cor. E) | none | Route A Thm D, Cor. E, Thm A, Thm F; Route D (I1); Thm LC, 11.1; lambda_g |
| induction | genus and level | none | genus and level |
| computer | none | none | Thm F (G = 5-7), zero count to G = 1000, base h <= 4 |
| TR | none | Hodge CohFT only (allowed) | none |
| referees | sound, sound | sound, sound | sound, sound |

**Recommendation**
- R3 proof (ii) as the proof of record: it has the fewest and most standard inputs.
- R2 second. It is independent of the dilaton calculus and of the integrable hierarchy, and it also gives (6.1) and the gap formula.
- R1 third.

**Common input.** The MV formula enters R3 and R1 directly, and probably R2 as well: P2 goes through Eynard's formula, whose Hodge form ABDKS attribute to BKMP for C^3, and that theorem is proved from MV (to be confirmed). So the three proofs are independent in their mechanism, not in their geometry.

## 3. Dependency chain of the all-genus theorem (Corollary R)

**Statement** (leaf_structure_proof.md s.13.2). For every g >= 2:
- every coloured stable graph of genus g with an edge, in the frame-B graph sum (13.1), is regular at t = 0 after summing over decorations;
- the single-vertex graph has the single pole gap_g t^{2-2g}, with gap_g = (2g-3)! P_g(x,1)/x;
- F4 holds in genus g.

**A. Graphs with an edge**

| | input | source | status |
|---|---|---|---|
| A1 | Theorem LC (F1-F3, L1-L3, all orders) | leaf_structure_proof.md s.3-10 | proved modulo literature. The one analytic input is Prop. 3.15: Lemma 3.8 (proved) and Stirling's remainder (DLMF 5.11.1, 5.11(ii), cited). Refereed sound. |
| A2 | Theorem R (i)-(iv): homogeneity, leaf expansion, string reduction, counting | proof_status s.3; leaf_structure_proof s.13.1 | proved |
| A3 | Theorem C (Z[h,n] = 0 for n >= 1); string/dilaton reduction to M_irr | conifold_regularity_note s.3a; proof_status s.1.2 | proved |
| A4 | no contributing genus-g vertex in a genus-g graph with an edge | s.4.3; Route D Lemma 5.3 | proved |
| A5 | Lemma A for 2 <= h <= g-1 | s.2 | proved modulo literature |

**B. The single-vertex graph**

Any one of three routes gives the single-pole statement:
- **(B-a)** Theorem R (v) with M_irr(g) gives the pole (-x c^3)^{g-1} Fpot_g(q_0) t^{2-2g}. R2 Cor. 7 (s.4.4) then evaluates it as gap_g. This route needs no Prop. 13.2 and no (H).
- **(B-b)** Theorem R (v) with M_irr(g), plus comparison with F_g via Prop. 13.2 and (H). This is the route of the current text.
- **(B-c)** Prop. 13.2 and (H) with A1-A5, without M_irr(g) (s.4.3).

**C. Interpretation (graph sum = GW potential)**

| | input | source | status |
|---|---|---|---|
| C1 | Prop. 13.2 and Lemma 13.3 (frame-B formula) | leaf_structure_proof s.13.3 | proved modulo Teleman 2012, Givental 2001, PPZ 2015, GP 1999, CG 2007, Mumford 1983 and FP 2000. Two referees (minor revision), revised. Unverified: Teleman was read via summaries; the CohFT axioms for non-compact Z go through the twisted P^1 theory. |
| C2 | (H): F_g - F_g(0) = P_g Li_{3-2g}(Q)/(eps1 eps2) | Brini-Schuler arXiv:2410.00118, Prop. 3.1 = Thm A1 (degeneration, Mumford, BP; no TR) | proved in the literature. The statement was checked, and localisation was confirmed to (4,1), (3,2), (2,3). |

**L. Lemma A via R3**

| | input | source | status |
|---|---|---|---|
| L1 | MV formula | LLZ JDG 65 (2003) Thm 1; OP GT 8 (2004) | literature; normalisation exact (2170 + 996 comparisons) |
| L2 | first-time-set (mKP) equations of tau_n = Z_MV(e^{n tau x} o p): P_2, and all P_k for proof (ii) | Jimbo-Miwa 1983; Ueno-Takasaki 1984; Takasaki 1984; Orlov-Scherbin 2001; Harnad-Orlov 2015 | literature. Derived by R3 referee 1 from the fermionic bilinear identity. Verified generically (degree <= 9) and on Z_MV. |
| L3 | dilaton calculus: Prop. 2.1, Cor. 2.2, Lemma 3.1 | A_mv s.2, s.3.1 | proved (refereed). Degree bound re-checked at genus 1-8. |
| L4 | Route A Thm B (proof (i) only) | A_mv s.3.3 | proved modulo the MV cut-and-join |
| L5 | Thm M, I_2, M_P, I_P, Theorem L | R3 s.2-5 | proved modulo L1-L3 (and L4) |

**Assembly.** Corollary R in every genus = A1-A5 + L1-L5 + one of (B-a), (B-b), (B-c). The identification with the conifold x C^2 potential uses C1-C2. Nothing applies TR to the conifold x C^2 theory. The only TR is R2's Cor. 7 in route (B-a), which is Hodge-side and allowed; (B-b) avoids it.

## 4. Checks of the main-session claims

### 4.1 The parity reduction (brief s.2(f)): confirmed

**Step 1: substitute Cor. E into T_D.** Take T_D(G) = Sum_{D(R)=D} (1/|Aut R|) prod((-1)^{j_i} Pi_{j_i}^2/j_i) Xi_G(R), with R irreducible (parts >= 2). After substituting Cor. E,

    T_D(G) = Pi_{D+1} kappa_D Xi_G(D+1),
    kappa_D = Sum_T (1/|Aut T|) (J-1)...(J-r+1) prod (-1)^t Pi_t/t.

**Step 2: a generating series.** Since (J-1)!/(J-r)! = (D+r-1)!/D!, we get kappa_D = (1/D)[y^D](1 - alpha)^{-D}, where alpha = Sum_{t>=2} (-1)^t (Pi_t/t) y^{t-1}.

**Step 3: Lagrange inversion.**
- It gives kappa_D = [y^D] log(w/y), where w is the inverse of y = w - Sum_t (-1)^t (Pi_t/t) w^t = -log B_{1+tau}(-w).
- Here log B_s(z) = Sum_k binom(sk-1,k-1) z^k/k, and B = 1 + z B^s.
- Hence w = e^{tau y}(e^y - 1).

**Result.**
- kappa_D = [y^D](tau y + log((e^y-1)/y)).
- kappa_1 = tau + 1/2, and kappa_D = B_D/(D D!) for D >= 2. This is tau-independent, non-zero for even D, and zero for odd D >= 3.
- D = 1 reproduces T_1 = ((1+4 sigma)/2) Xi(2), which is D5.1. It is also the leading coefficient of R1's T1.

**Exact checks**
- kappa_parity.pl re-run: all OK.
- New: critic\kappa_poly_check.pl. Since kappa_D has tau-degree <= D, checking D+3 values of tau is a proof for each D <= 14. Result: 147 evaluations, 0 mismatches.

The reduction is no longer needed by R3 or R2; R1 still uses it.

### 4.2 The reduction of Lemma A to (OL) (brief s.3): confirmed

Let IH(G,D) be: Lemma A(h) for 2 <= h < G, and Xi_G(T) = 0 for every non-empty T with parts >= 2 and D(T) < D.

**The two cited results, against IH(G,D)**
- **Cor. E** (A_mv s.3.5). Its hypothesis is exactly IH(G,D).
- **Theorem D5 at k = D.**
  - Its hypotheses are (H) at order D, and Lemma A at levels <= D-2 in genera <= G. Both are implied by IH(G,D).
  - Its framework is Prop. 13.2, Theorem LC, Lemma 5.2 (via Thm 11.1), Lemma 5.3 and (V).

**The induction**
- G = 2: Lemma A(2) is Xi_2(2) = 0, which is (c).
- D = 1: (c).
- Even D: (f) gives Pi_{D+1} kappa_D Xi_G(D+1) = 0 with kappa_D != 0. Then Cor. E covers all of level D.
- Odd D >= 3, including the top level 2G-3: (OL), then Cor. E.
- Conversely, Lemma A(G) contains the conclusion of (OL).

Hence Lemma A for all h <=> (OL) for every G >= 3. The reduction inherits KP and the conifold inputs. It is the basis of R1 and of the computer-assisted ranges.

### 4.3 Lemma A only below the genus (leaf_structure_proof s.13): confirmed, with a corrected justification

**(i) Genus-0 vertices need >= 3 half-edges.**
- Dilaton leaves in (13.1) have psi-degree >= 2: vertexOptionsH4 enforces Min >= 2, and the PPZ translation has no z^0 or z^1 term.
- At a genus-0 vertex with e half-edges and m leaves, dimension requires e + m - 3 >= 2m. So e >= m + 3 >= 3.

**(ii) No genus-g vertex.**
- If a genus-g graph with an edge had a genus-g vertex, it would be a tree, and all its other vertices would have genus 0 and hence degree >= 3.
- Then 2(V-1) >= 1 + 3(V-1), which is impossible.
- So Theorem R (i)-(iv) needs M_irr(h) only for h <= g-1.

**(iii) The single-vertex graph.**
- Prop. 13.2 and (H) give: the total is F_g, whose singular part is gap_g t^{2-2g}.
- So the single-vertex graph carries that pole, without M_irr(g).
- Its other polar coefficients then vanish. This is Route D's (I1), which R1 uses at genus G given Lemma A(< G); there is no circularity.

**Wording fix.** status_2026-10-31.md says "a genus-0 one is unstable". That is wrong: such vertices are stable, but their integrals vanish by (i).

### 4.4 New: R2's Corollary 7 is F4 at n = 0 in every genus

F4 is (-x c^3)^{h-1} Z[h,0] = gap_h, with Z[h,0] = Fpot_h(q_0). Until now it was exact only for h = 2, 3, 4. conifold_regularity_note s.3a calls it "the one quantity the argument leaves free".

**The conversion**
- Take c = 1 + x and tau = -x/c. Then e1 = x/c, e2 = 1/c and sigma = -x/c^2.
- P_h is homogeneous of degree 2h.
- Route D's normalisation gives Z[h,0] = c^{-(3h-3)} Int K Omega.

Cor. 7 then gives

    (-x c^3)^{h-1} Z[h,0] = (2h-3)! P_h(x,1)/x,

which is the t^{2-2h} coefficient of P_h Li_{3-2h}(e^t)/(e1 e2), i.e. gap_h in the sense of rc_lib.m (Pg, cyPolar).

**Checks**
- critic\gap_identity_check.pl: h <= 12 at 4 values of x, 44 comparisons, 0 mismatches.
- The normalisation agrees with the exact F4 at h = 2, 3, 4.
- critic_gap_check.wls repeats the check with rc_lib.m's own code.

**Consequences**
- Route (B-a) of s.3.
- An all-genus consistency check of Prop. 13.2's normalisation against (H).
- Route (B-a) uses Hodge-side TR (allowed). A TR-free proof of (6.1) at j = 0 would remove that (plan P5).

## 5. Gaps

### 5.1 Literature and normalisations not re-verified line by line

1. **R3, the 2D Toda input.**
   - The citations are loose; "Okounkov math/0004128 s.2.7" is unconfirmed.
   - Give theorem numbers in Jimbo-Miwa, Ueno-Takasaki, Orlov-Scherbin and Harnad-Orlov.
   - Or write the three-line fermionic derivation (R3 referee 1 s.1.1) into the proof.
2. **MV formula and cut-and-join.** Give theorem numbers (LLZ; OP).
3. **R2.**
   - For ABDKS Thm 2.8 / Prop. 3.1, the statements were checked but the proof was not re-checked.
   - Eynard 2014 was taken from ABDKS (49)-(50), which quote it without proof; it was not compared with the original. Note the (49) index slip.
   - BKMP for C^3 (Chen; Zhou) lies behind (44). Cite theorem numbers.
   - TR1 gives an alternative route to P2.
4. **R1.** KP for MV and the Takasaki-Takebe dispersionless Hirota equations. Full Hirota checks: genus 1-2 (round 1) and G = 7-10 (R2 referee 2).
5. **Corollary R.**
   - Stirling's remainder, cited (Theorem LC).
   - Teleman, read only via summaries, and the CohFT axioms for the non-compact Z (Prop. 13.2).
   - Brini-Schuler's own inputs are cited, not re-proved. App. A has misprints.
   - lambda_g: Faber-Pandharipande, Ann. of Math. 157 (2003).
6. **Normalisations.**
   - MV: exact.
   - P_2: derived from the fermionic identity; the wrong sign fails.
   - ABDKS/Eynard: fixed against TR-free data for g <= 9.
   - s.4.4: consistent with F4 at h <= 4. Confirm with critic_gap_check.wls.

### 5.2 Logical and write-up gaps

No logical gap was found in any main line. Write-up items:

1. **Lemma A statement.** Must say "non-empty R" (brief s.1; leaf_structure_proof s.13).
2. **R3.**
   - Items (a)-(d) of s.1.1.
   - The wording of I_2, and of the summary sentence on genus-G terms.
   - Add the derivation of P_k.
   - The line (2) no-go is only a remark.
   - Genus 7 did finish: g7_points_A.out exists and exact genus-7 checks exist.
   - Separate the dependency sets of proof (i) via Thm B and via Cor. E.
3. **R2.**
   - Non-empty R.
   - The h = 0 convention.
   - Theorem 3 in the multiplied-out form.
   - The counts and the no-data script flaw.
   - V5 and V20.
   - The genus-7 table.
   - M_irr(7) is complete.
   - Literature keys.
4. **R1.**
   - The Cor. 7.2 constant.
   - W = alpha.
   - The (I4) factorial.
   - The two-line (I1) argument (s.4.3).
   - The s.7(b) Lagrange step.
   - Label the zero count "computer-assisted", with (B1) stated in general.
   - Scan labels, Thm 6.1 scope and notation.
   - A clean count_zeros2 raw file.
   - "2(b) not needed for G >= 8".
   - List as dependencies the frame-B premise, Theorem LC, Theorems A-C of the conifold note and the regularity of the colour-2 single vertex.
5. **R4.**
   - The G = 3 claim.
   - The K1' evidence.
   - Cosection status.
   - N1-N4 wording.
   - Regenerate lam1_rec_h8.txt and lam1_rec_h16.txt.
   - Define phi_a.
   - Write out H_2.
   - B - 1 in the Sheffer form.
   - Section references.
   - Give the constrained rigidity its own name, not TLR.
6. **Route A record.** The TLR nullities 11, 11, 7 included c = 2.

### 5.3 Uniformity

- R3 is uniform; its mechanical checks only confirm.
- R2 is non-inductive.
- R1's finite parts are computer-assisted, with analytic bounds beyond.
- No gap.

### 5.4 TR audit

- R3, R1 and R4 are TR-free (referee audits).
- R2 uses TR only for the Hodge CohFT (allowed; confirmed again with this run). It does not use input (H).
- DVV appears only in the Hodge engines (Witten-Kontsevich numbers), and only for checks.
- R1's genus-0 Y(z) is genus-0 algebra.
- Theorem LC's Beta integral is not a mirror curve.

### 5.5 Bookkeeping

1. Save R2's second referee report as R2_onepoint\referee2.md.
2. Update every place that still says "conditional on M_irr(h), h >= 5" or "(6.1) conjectural":
   - leaf_structure_proof s.0.4, s.13.2 and G1;
   - proof_status s.0, s.2 (Lemmas B and D) and s.3;
   - conifold_regularity_note s.3a, s.4 and s.7b;
   - proof_M_residue s.5, s.6 and s.9;
   - brief s.2(g).
3. Then update memory, the README and the tarball.

## 6. Ranked plan for round 3

1. **P1 (essential, low).** Write and referee one proof of Lemma A for all h.
   - Main proof: R3 proof (ii), with proof (i) and the Cor. E shortcut as variants.
   - Include all referee corrections, the non-empty statement and the fermionic derivation of P_k.
   - Restate Corollary R as proved modulo literature in every genus, with the chain of s.3.
   - Add s.13.4: Lemma A only below the genus (s.4.3), and F4 at n = 0 in all genera (s.4.4).
2. **P2 (essential, low-medium).** Literature pinning:
   - MV (LLZ, OP);
   - 2D Toda / content-product (Jimbo-Miwa, Ueno-Takasaki, Orlov-Scherbin, Harnad-Orlov);
   - ABDKS Prop. 3.1 / Thm 2.8 (published version);
   - Eynard 2014, and the index slip in ABDKS (49);
   - BKMP for C^3 (Chen; Zhou);
   - Teleman, read directly;
   - the published Brini-Schuler.
3. **P3 (high, medium).** Integration critic across both workflows (R1-R4, TR1-TR2), after TR1's second referee. Recommended presentation order: R3, then R2 (or TR2 + TR1), then R1.
4. **P4 (medium, low).** Record the new all-genus results:
   - (6.1) and the gap formula;
   - Lemmas B and D (R2; TR2 for all tau);
   - the vanishing direction of Conjecture M for all h and r;
   - K1;
   - T1 and the free energy on the slice.
5. **P5 (optional, medium-hard).** A TR-free proof of (6.1) and of the gap formula, via R3's M_P as a recursion off the range, or via E-operators and Orlov-Scherbin. This would make F4 at n = 0 TR-free without Prop. 13.2 and (H).
6. **P6 (low).** Bookkeeping (s.5.2, s.5.5).
7. **P7 (optional, low).** The Mathematica confirmations of s.7, and exact certificates: Theorem F for G = 8..20, and DVV-free checks of M_irr(5..7).
8. **P8 (research, optional).**
   - Multi-leaf closed forms (SN, NR) from R2's Theorem 1.
   - Merge formulas from R1's slice formula.
   - R4's components k >= 2.

## 7. Consolidated Mathematica checks for the main session

None of these is needed for a proof. Run at most 2 kernels, and read each script first.

1. **critic_gap_check.wls** (new; runs in seconds).
   - Path: HVC\lemmaA\round2\critic\critic_gap_check.wls.
   - Run: wolframscript -file "G:/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/round2/critic/critic_gap_check.wls".
   - Expected: h = 2..10 all True, then "ALL TRUE: True".
   - Purpose: F4 at n = 0 in all genera, checked with rc_lib.m's own Pg and cyPolar.
2. **R3_thmM_check.wls** (minutes).
   - Path: HVC\lemmaA\round2\R3_toda\R3_thmM_check.wls.
   - Run from the hodge_vertex_conifold folder: wolframscript -file lemmaA/round2/R3_toda/R3_thmM_check.wls.
   - Expected: every line 0, then "ALL ZERO: True", then a non-zero control.
   - Purpose: Theorem M symbolic in tau, on the Mathematica q-tables (G = 2-4).
3. **R1_checks.wls** (optional; heavy).
   - Path: HVC\lemmaA\round2\R1_rigidity\R1_checks.wls.
   - Expected:
     - (C1), (C2): True;
     - (C3): True three times;
     - (C4): True at tau = 2/7, -3/5 and 3;
     - (C5): {tuples, 0};
     - (C6): 0 failures for 8 <= G <= 60, and 4 failures for G = 4..7.
   - These are already duplicated by two Perl referee implementations.
4. **TR line** (optional).
   - Extend HVC\hvc_tr_tau.wls at tau = 2, 3 to genus 4-5, n <= 3, and compare with perl xy_npoint.pl TAU G N JMAX (R2).
   - TR2_recursion\tr2_closed_check.wls gmax: all differences 0.
   - TR1's request: hvc_tr_tau.wls at a non-integer tau, genus 3-4.
5. **No longer needed:** R4's Xi^(1)_10 check. C1 equals minus the sigma^{h-1} part of (6.1) for h <= 14, and (6.1) is now proved.

## 8. Files written by the critic

All in HVC\lemmaA\round2\critic\:
- kappa_poly_check.pl and kappa_poly_check_D14.out;
- gap_identity_check.pl and gap_identity_check_h12.out;
- critic_gap_check.wls.

No other file was modified.

