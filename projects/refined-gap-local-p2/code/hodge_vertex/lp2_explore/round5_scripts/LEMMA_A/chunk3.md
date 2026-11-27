---

## 8. Complements (proved, not used in s.3-s.6)

These state the whole matching of the MV series with the Hodge side. They live over Q(tau): the series below are literally undefined at zeros of Pi_k, and must not be evaluated there. Specialised tau is reached only through s.7.1.

### 8.1 All-profile, all-genus matching [P modulo (MV)]

In Q(tau)((x))[[s, P_2, P_3, ...]] put u = s/(1+s), x' = x/(1+s), p0_1 = s/(sigma x), p0_k = P_k/(sigma x Pi_k) (k >= 2), and

    FF(x,s;P) := log Sum_nu e^{-u|nu|} e^{(tau+1/2) kappa_nu x'/2} prod_box (2 sinh(h x'/2))^{-1} s_nu(p0).

**Proposition 8.1.** FF = Sum_{g>=0} Sum_R hbar^{g-1} QQ_g(s;R) P_R/|Aut R|, over multisets R of integers >= 2, where

    QQ_g(s;R) := (1+s)^{-(2g-2+r)} e^{-uJ} Sum'_{a>=0} M_g(R + 1^a) (u e^{-u})^a/a!,

and the prime omits (R,a) = (empty, 0). In particular only the powers x^{2g-2} occur, and sigma^{1-g} [x^{2g-2}] d_{P_{j_1}}...d_{P_{j_r}} FF |_{P=0} = QQ_g(s;R), with no Aut factor.

*Proof.* The map x' -> x/(1+s), p_1 -> e^{-u} s/(sigma x), p_k -> e^{-uk} P_k/(sigma x Pi_k) is a well-defined continuous ring map: p_1 lands in (s), p_k in (P), and each s-coefficient of x'^n = x^n (1+s)^{-n} is a monomial in x. Homogeneity gives s_nu(e^{-uk} p0_k) = e^{-u|nu|} s_nu(p0), so this map sends F to FF. Apply it to (2.1) at mu = R + 1^a, using 1/(1+s) = 1-u and s(1-u) = u. All x and sigma powers cancel except hbar^{g-1}. QED

### 8.2 Finite formula [P modulo (MV)]

For h >= 1, 2h-2+r > 0, R a non-empty multiset of integers >= 2, N = 3h-3+r, in Q(tau):

    Xi_h(R) = ( |Aut R| / (sigma^{h+r-1} Pi_R) ) Sum_{a=0}^{N} (-1)^a sigma^{-a} L_{N-a}^{(2h-2+r+a)}(-J-a) [x^{2h-2+r+a}] [p_R p_1^a] F.

This is (2.1) inserted into the Laguerre form of Lemma 3.2. For R empty and h >= 2 the same holds for a >= 1, plus the a = 0 term L_N^{(2h-2)}(0) <Omega_h>, with <Omega_h> := Int_{Mbar_h} Omega_h (s.8.5).

### 8.3 Reformulation of Lemma A [P]

Theorem L is equivalent to: for each non-empty multiset R of integers >= 2, the series Sum_h hbar^{h-1} sigma^{1-h} [x^{2h-2}] d_{P_R} FF |_{P=0, s=-1}, taken coefficientwise as polynomials in s, is a Laurent polynomial in hbar of degree <= floor(D(R)/2). Indeed Xi_h(R) = 0 for all h with D(R) <= 2h-3 means that only h - 1 <= floor(D(R)/2) survive.

### 8.4 Table of QQ_g [P modulo (MV)]

Using Lagrange inversion for u = y e^u:

| (g, R) | QQ_g(s;R) | value at s = -1 |
|---|---|---|
| r >= 1, 2g-2+r > 0 | Q_g(s;R), polynomial of degree <= N | Xi_g(R) |
| g = 1, R = (j) | (j - 1 - 1/sigma + s)/24 | (j - 2 - 1/sigma)/24 |
| g = 0, r >= 3 | polynomial depending only on (J, r) | (J-1)...(J-r+3) |
| (0, (j,k)) | 1/(j+k) | 1/(j+k) |
| (0, (j)) | (1+s)/j^2 - s/(j(j+1)) | 1/(j(j+1)) |
| g >= 2, R empty | Q_g(s;empty) - (1+s)^{-(2g-2)} <Omega_g> | none: pole of order 2g-2 |
| (1, empty) | (1/24) [ log(1+s) + (1+1/sigma)(1/(1+s) - 1) ] | none: log and simple pole |
| (0, empty) | s + s^2/4 + s^3/(6(1+s)) | none: simple pole |

So QQ_g has a value at s = -1 for every non-empty R, and for R empty it has none in any genus. For R empty and g >= 2, the value Xi_g(empty) is recovered only after restoring the a = 0 term, as in s.8.2.

Derivations. (0,(j,k)): M_0 = (J+a)^{a-1}, so G = e^{Ju}/J. (0,(j)): Lemma 3.3. (0,empty): M_0 = a^{a-3}, (y d/dy)^2 G = u, G(0) = 0, so G = u - 3u^2/4 + u^3/6. (1,(j)): Omega_1 = 1 - (1+1/sigma) lambda_1 on Mbar_{1,n}, as -1/tau + 1/(1+tau) = -1/sigma and lambda_1^2 = 0, with <tau_1>_1 = <lambda_1>_1 = <tau_0 tau_2>_1 = 1/24; k >= 2 is excluded by degree. (1,empty): y d/dy G = y G_(1) = s Q_1(s;(1)) = s (s - 1/sigma)/24 by Lemma 3.2, and y d/dy = s(1+s)^2 d/ds, G(s=0) = 0. (0, r >= 3): A_m(X) = Sum_a (X+a)^{a+m} y^a/a! satisfies (X + theta) A_{m-1} = A_m, which gives Q^[m] = X Q^[m-1] + m s Q^[m-1] + s(1+s) d_s Q^[m-1] with Q^[-1] = 1/X, so Q^[m](-1) = (X-m) Q^[m-1](-1).

### 8.5 The empty profile (not part of Lemma A)

- **<Omega_g> from the MV series alone.** [P modulo (MV)] Lemma 3.2 at r = 0, J = 0 and the degree bound at m = N+1 give Sum_{a<=m} (M_a/a!) (-1)^{m-a} L_{m-a}^{(2h-3+a)}(-a) = 0. The a = 0 coefficient is L_{N+1}^{(2h-3)}(0) = binom(N+2h-2, N+1) != 0. So M_0 = <Omega_h> is a finite combination of MV coefficients M_a, a >= 1.
- **Closed form.** [P modulo C: Mumford, Faber-Pandharipande] <Omega_g> = -A_g sigma^{1-g}, A_g = Int lambda_{g-2} lambda_{g-1} lambda_g. On w_1+w_2+w_3 = 0 the degree 3g-3 part of prod Lambda^vee(w_a) integrates to (-1)^{g-1} [ e_3(w) Int lambda_{g-1}^3 + (Sum_{a != b} w_a^2 w_b) A_g ]: the (3,0,0) terms die by lambda_g^2 = 0. Mumford gives lambda_{g-1}^2 = 2 lambda_{g-2} lambda_g, so Int lambda_{g-1}^3 = 2 A_g, and Sum_{a != b} w_a^2 w_b = -3 e_3(w). Hence Int prod Lambda^vee(w_a) = (-1)^g A_g e_3(w); at w = (1, tau, -1-tau), e_3 = -sigma, this gives the stated value. The closed form of A_g (as quoted in OP s.0.6) is not needed here.
- **Xi_h(empty).** Xi_h(empty) = Q_h(-1;empty) is the gap. Non-vanishing as an element of Q[tau, 1/sigma]: [C: R2 Cor. 7, sigma^h Xi_h(empty) = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l), e1 = -tau, e2 = 1+tau, S(y) = (y/2)/sinh(y/2); Hodge-side TR, allowed]. This non-vanishing, **unlike Theorem L, does not specialise**: at h = 2 one has Xi_2(empty) = -(24 + 28 s' + 7 s'^2)/5760 with s' = 1/sigma, which has zeros at sigma = -0.804..., -0.363... (complex tau). Local P2 statements allowing complex weights must take this into account.
- **Part 0.** Xi_h(U + 0) = (D(U) - 2h + 2) Xi_h(U) is s.7.2. This identity needs no citation.

### 8.6 Variants, not used

- Proof (i) of R3: P_2 alone (2-insertion, Theorem M), then Route A Theorem B (MV cut-and-join) or Route A Cor. E (KP). More inputs and no gain.
- R2 (lemmaA/round2/R2_onepoint): ABDKS x-y swap plus Eynard 2014 on the Hodge curve C_tau, non-inductive. Allowed as Hodge-side TR. It stays a backup. Its literature and its normalisation "fixed against TR-free data for g <= 9" are unpinned, so it does not meet the all-orders ruling as written (critic 1 S3(b) of round 4).
- R1: not admissible under the all-orders ruling, since its base h <= 7 and zero count are computer-assisted.

---

## 9. Literature: what was and was not confirmed

| item | used for | status | confirmed | not confirmed |
|---|---|---|---|---|
| Liu-Liu-Zhou, J. Differential Geom. 65 (2003) 289-340; arXiv:math/0306434v3 | **the input** | [C] | Introduction read directly (author and referees of round 5): C_{g,mu} with -sqrt(-1)^{|mu|+l}/|Aut mu|; the l < 3 extension; eqs. (1), (2), (3); kappa_nu. Theorems 1, 2 are the cut-and-join for R and C, so "LLZ Thm 1" for the formula is wrong by arXiv numbering. Prop. 2.2 (tau = 0) agrees. Journal data confirmed on Project Euclid. | Journal text not opened, so "eq. (1)" is confirmed for arXiv v3 only (dated 31 Dec 2004, after publication). LLZ's proof not re-checked. |
| Okounkov-Pandharipande, Geom. Topol. 8 (2004) 675-699; arXiv:math/0307209v3 | second derivation of (2.1); independent published proof of the formula | [C], not load-bearing | Theorem 1 (0.8), (0.1), (0.5), s.0.6, (1.24), (1.25), (1.28) read directly | The printed 1/|Aut mu| in (0.5), (0.8) is inconsistent at repeated parts (s.2.2) |
| String and dilaton with pulled-back classes | Lemma 3.2 | textbook | Mechanism stated in Lemma 3.2 | No theorem number pinned (e.g. Faber-Pandharipande, Invent. Math. 139 (2000), not opened for one) |
| Genus-0 multinomial formula | M_0, s.1 | textbook | Stated in s.1 | No number pinned |
| Row-hook lemma | Prop 2.1(a) | textbook, sketch in place | LLZ (8)-(9) point to Macdonald pp. 10-11 | Macdonald I.1 Ex. 1 from memory |
| Principal specialisation of Schur functions | Prop 2.1(b) | textbook | Exponent computation done in place | Macdonald I.3 Ex. 2 and Stanley EC2 Cor. 7.21.3, both from memory |
| Frobenius; Murnaghan-Nakayama | s.4 (F1) | textbook | | No number pinned |
| Lagrange inversion | Lemma 3.3, s.8.4 | textbook | | |
| Jimbo-Miwa, Publ. RIMS 19 (1983) 943-1001 | corroboration of (4.2), not an input | not an input | Read in scan by the ThmL author: Thm 1.1-1.2; Appendix 1 p. 996, table "1st Modified KP hierarchy". Exact linear algebra: P_2..P_5 lie in the span of the transcribed entries | Equation number of the bilinear identity not looked up; no referee re-read the table |
| Harnad-Orlov, CMP 338 (2015); arXiv:1407.7800v9 | corroboration only | not an input | (2.2), (2.8)-(2.10), (5.1)-(5.4) read. (2.2) as printed has exponent e^{-xi(dt,z)}, which appears to be a sign slip: it fails exactly at level >= 2 on a generic content-product tau, while e^{+xi} holds and coincides with (4.1) | Journal version not compared |
| Okounkov, Math. Res. Lett. 7 (2000) 447-453 | no longer needed | not an input | Abstract read | "s.2.7" of R3 unconfirmed |
| Orlov-Scherbin TMP 128 (2001) 906-926; Ueno-Takasaki ASPM 4 (1984) 1-95; Takasaki ASPM 4 (1984) 139-163; Takebe LMP 21 (1991) 77-84 | no longer needed | not an input | Bibliographic data from Harnad-Orlov's list only | Texts and theorem numbers not read |
| Mumford relation; Faber-Pandharipande A_g | s.8.5 only | [C], outside Lemma A | A_g closed form as quoted in OP s.0.6 | FP not opened |
| Arbarello-Cornalba pushforward | only to rewrite the m-sum as K = exp(-Sum kappa_k/k); the m-sum is the definition of record here | textbook, not used | Rederived: 1 - s w/(1-w) = (1 - (1+s) w)/(1-w) gives K_s = exp(Sum ((1+s)^k - 1) kappa_k/k) | |

Web sources read in round 5: arXiv:math/0306434, projecteuclid.org/euclid.jdg/1090511689, arXiv:math/0307209, arXiv:1407.7800.

---

## 10. Status table and open points

### 10.1 Status

| result | status | rests on |
|---|---|---|
| Prop 2.1, (2.1) real form | P mod C | LLZ (1)-(3), textbook specialisation; all signs derived |
| (2.2)-(2.4) | P | Prop 2.1, Pi_k != 0 in Q(tau) |
| Lemma 3.1, Lemma 3.4, Cor 3.5, Remark 3.6 | P | Lemmas 3.2, 3.3 |
| Lemma 3.2 | P mod textbook | string, dilaton |
| Lemma 3.3 | P mod textbook | Lagrange inversion |
| (F1)-(F6), P_k for all k, (4.2), (4.3) | P mod textbook | Murnaghan-Nakayama, Frobenius |
| Theorem M_k (5.1), all k_0, n, R | P mod (MV) | s.2-s.4 |
| Theorem L | P mod (MV) | Theorem M_k, (S); no hypothesis beyond (H1), (H2) |
| Cor L1, L2, L3 | P mod (MV) | Theorem L; L3 by weight count |
| s.7.4, s.7.5 transfers of the vanishing | P mod (MV) | Cor L1-L3 |
| Prop 8.1, s.8.2, s.8.3, s.8.4 | P mod (MV) | not used in Theorem L |
| s.8.5, closed form of <Omega_g> | P mod C | Mumford; outside Lemma A |
| s.8.5, non-vanishing of Xi_h(empty) | C | R2 Cor. 7, Hodge-side TR; outside Lemma A |

**Precise claim.** Theorem L, for non-empty R, is proved to all orders, modulo the single cited theorem (MV) and the textbook facts of s.0. No step of s.2-s.6 is open.

### 10.2 Open points

- **O1 (refereeing).** This assembled document has not had its own referee. Its two source reports had two referees each, and s.2-s.6 were re-derived by the writer. Notation was renamed during assembly (dictionary in s.1), and doc_check.py (s.12) was written against the final formulas to guard against transcription errors, but that is evidence only.
- **O2 (citation details; none is mathematical).** Journal numbering of LLZ. Macdonald and Stanley exercise numbers, from memory. No pinned theorem number for string/dilaton with pulled-back classes, for Frobenius or for Murnaghan-Nakayama. Harnad-Orlov's sign in the journal version.
- **O3 (input not re-proved).** The proof of (MV) in LLZ was not re-checked. (MV) has a second published proof (OP), with the caveat on the printed 1/|Aut mu|.
- **O4 (editing not done, outside this task).** The other files still carry the items of s.11.3.
- **O5 (outside Lemma A, unchanged by this document).**
  - Empty profile: the gap formula and F4 / (VI) still rest on R2 Cor. 7 (Hodge-side TR) or on the frame-B route. No TR-free proof of the gap formula exists yet (plan P5). Critic 1 S3 of round 4 is not addressed here.
  - The identification of Theorem R's Taylor coefficients with Xi_h(R), conifold (Theorem R) and local P2 (interface item (h), restated Dil_0).
  - Prop. 13.2 and Prop. 13.2_P2 (frame-B), lem:modrel (critic 1 S4), the gap corollary (S5).

---

## 11. Disposition of critic and referee items

### 11.1 Items now discharged

| item | how |
|---|---|
| Round 4 critic 1 S1 (MV-Xi matching only "normalisation exact" with data comparisons) | Prop 2.1 and Lemma 3.2, all genera, all profiles; no constant from data |
| Round 4 critic 1 S2 (no proof of record; P_k fixed by checks for k <= 4; input list incomplete; P1 missing) | Record = R3 proof (ii) rewritten here; s.4 derives P_k for all k; input list in s.0; this file is P1 |
| Round 4 critic 2 minor on s.7.4, s.8.2 (input list, proof of record) | s.0 |
| Round 4 lemmaA referee: MV matching belongs among open steps | Closed by s.2 |
| Round 4 lemmaA referee: hypotheses never used at specialised tau; direct source of the Laurent form; no part 0 | s.6 Uniformity; s.7.1; Cor L2 |
| round2_status s.5.2 item 1 (non-empty R) | In the statement |
| round2_status s.1.1 (a) (additive constant of G_empty) | s.2.3: the a = 0 term is explicit; R3 Lemma 2.1 needs the correction, which is not a new finding |
| round2_status s.1.1 (b) (v = 0 functions) | Moot: every block has U non-empty |
| round2_status s.1.1 (c) (single sector) | s.5, by a purely formal argument |
| round2_status s.1.1 (d) (odd n sampling) | (5.1) is proved for all n; only n = 2G is used; no sampling claim |
| round2_status s.5.2 item 2 (P_k, I_2 wording, dependencies) | s.4; the I_2 wording does not arise; proof (i) dependencies are not used |
| R3 Lemma 2.3 (empty-profile table, "referee should re-derive") | Remark 3.6, corollary of Lemma 3.4, not needed separately |

Item 2 of LP2 draft s.8.2: its first bullet (2D Toda input) is replaced by s.4, its second bullet (LLZ/OP numbers) is settled against arXiv v3 only (O2), and its third bullet (matching) is closed by s.2.

### 11.2 Round-5 referee fixes applied here

From the MVmatch referees:
- False sentence on poles of QQ_g corrected: s.8.4 says that every empty-profile QQ_g lacks a value at s = -1, with pole order 2g-2 at g >= 2.
- Quantifier of the reformulation fixed: s.8.3.
- Overclaim of "discharged" removed: s.11.1 and O5 state exactly what is and is not closed.
- LLZ Prop. 2.1 not cited for the hook form: Prop 2.1(a) uses (LLZ-3) and its own derivation.
- <Omega_h> shown to be determined by MV data: s.8.5.
- E3 given as a cross-reference to round2_status s.1.1(a).
- Citation limits recorded: s.9.
- Scope note: Prop 8.1 and s.8.2 hold over Q(tau) only, stated at the head of s.8.

From the ThmL referees:
- Non-vanishing of Xi_h(empty) moved out of the theorem to s.8.5, with its source and with the warning that it does not specialise (exact h = 2 formula, complex zeros).
- The (2-2h) identity proved: s.7.2.
- The vacuous sub-case "T_B != R" of the induction deleted: s.6.
- Definition of G^(g)_U and the empty-profile part of F_g made explicit in all genera: s.2.3.
- Principal specialisation added to the textbook list, and the hook identity written in, not only cited.
- Notation cleaned: T no longer means both a multiset and a leaf class, D no longer means level, y d/dy and Hirota derivative, u no longer means both tree variable and Hirota shift, tau no longer means both framing and tau function.
- X_0 flagged as not a Xi: s.1.
- (S) stated for g >= 2, the only use.
- Harnad-Orlov worded as "appears to be" a slip: s.9.
- The label "P mod (MV)" restricted to the vanishing, and not extended to item (h): s.7.5, O5.

### 11.3 Errata and edits to make in other files (not done here)

- **E1.** lemmaA/A_mv/report.md line 177 misquotes LLZ: prefactor -(sqrt(-1))^l instead of -(sqrt(-1))^{|mu|+l}, and V_nu in [m] = q^{m/2} - q^{-m/2} instead of the printed 2 sin. The two slips cancel only together, because prod 1/[h] = (-i)^{|nu|} V_nu is exactly the rescaling p_k -> (-i)^k p_k. Its (1.1) is correct.
- **E1'.** "LLZ Thm 1" should read "LLZ eq. (1)": A_mv/report.md line 8, R3 report line 250, round2_status lines 26 and 460, LP2 draft s.7.4 and s.8.
- **E2.** Route A Cor. 2.3: "+ (unstable constant)" is not constant in s at fixed x. It is -Sum_{g>=2} A_g x'^{2g-2} plus the genus-0, 1 empty terms of s.8.4. "Xi is the genus-by-genus value at s = -1" fails for R empty in every genus: at g = 0, 1 there is no value, and at g >= 2 it needs the a = 0 term restored. "P_1 := s" means that parts equal to 1 are not separate variables.
- **E3.** R3 Lemma 2.1: G^(g)_0 lies in w^{2g-2}Q[w] only after adding the constant <Omega_g> (round2_status s.1.1(a)). Harmless there.
- **E4.** OP (0.5), (0.8): printed 1/|Aut mu|, s.2.2.
- **E5.** LLZ Prop. 2.1 as printed has 2^l (the text extraction shows it as 2l) where (LLZ-3) gives one factor 2 per box.
- **E6.** R3 report lines 11, 34, 246: statement lacks "non-empty".
- **LP2 draft.** s.0.3, s.7.4(b) with its input list, s.8.1 label of Lemma A_P2(b) (to "P mod (MV), LEMMA_A_PROOF.md", after O1), s.8.2 item 2 (replace by a pointer here), s.8.2 "Nothing else" and Appendix.
- **Conifold notes.** leaf_structure_proof s.0.4, s.13, G1; proof_status; conifold_regularity_note; proof_M_residue (round2_status s.5.5); RESUME.md and memory.

---

## 12. Evidence [E] (no proof step uses any of this)

All exact unless stated. Paths are under <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round5/ (session scratchpad; copy them into the project if they are to be kept).

| check | range | result | path |
|---|---|---|---|
| Fresh check of this document's formulas: Prop 2.1 Schur side, Laguerre form, (5.1), target coefficient | tau = 2/7, -5/3; G = 2, 3; k_0 <= 4; n = 2G, 2G+1; 48 relations each | 0 failures; target response exact 24/24; flipping c(odd,k) or beta -> -beta breaks 13/24 (the rest are terms that vanish anyway). X_2(3) = 1367/186624, X_3(5) = 205027/47029248 at tau = 2/7, equal to the reports' values | LEMMA_A/doc_check.py, doc_check_t27.out, doc_check_tm53.out |
| (2.1) against an independent Hodge engine (DVV, Mumford), all new code | author: g <= 4, 150 pairs (g, mu), 1387 comparisons; referee 2: g <= 5 at l <= 2, g <= 4 at l = 3, 595 comparisons, several tau including zeros of Pi | 0 mismatches; 9 kinds of control all fail | MVmatch/, MVmatch_ref2/ |
| LLZ (1)-(3) literally in lambda over Q(i) | 90 cases | match; the mixed misquote fails 66/90 | MVmatch/check3.py |
| Prop 8.1 directly from the Schur sum | g <= 3 (author), g <= 6 empty profile (referee) | 0 mismatches; x' = x, no e^{-u|nu|}, s -> -s, no Pi all fail | MVmatch/, MVmatch_ref2/testC*.py |
| Lemma 3.2 and operators | genus 0-4, non-integer leaves | 0 failures | MVmatch/, MVmatch_ref2/testB*.py |
| s.8.2 against the definition, h = 2, 3, 4; Schur-side zeros at genus 4-7 | Xi_4(7) = 193375/214990848 at tau = 2/7 non-zero, as expected out of range | 0 differences | MVmatch/check4.py, MVmatch_ref2/testD*.py |
| P_2..P_7 on generic content-product tau; P_k <= 10 at low weight | degree 8, 9 | 0 non-zero; swapped order, tau_{-1}, wrong c_(k) non-zero | ThmL/check1_toda.py, ThmL_ref2/t6_pk.py |
| P_k on Z_MV as Laurent series in x | tau = 2/7, -5/3; k <= 6, weight <= 8, order 30 | 114 coefficients, 0 failures | ThmL/check6_mvbilinear.py |
| (5.1), MV-side values | G <= 4 (author, 1638 relations), G <= 6 mod 2^61-1 (referee 1, 218 + 196 relations) | 0 failures | ThmL/check3_relation.py, ThmL_ref1/ref_m41.py |
| (5.1), Hodge-side values, MV not assumed | G <= 5, including tau = -1/3, a zero of Pi_3 | 186 relations, 0 failures; 5 kinds of control break it | ThmL_ref2/t2_*.py, t3_run.py |
| Target coefficient for every part k_0, including repeated | | 202/202, 124/124, 32/32 | ThmL/, ThmL_ref2/t7_resp.py |
| Mechanical run of the s.6 logic over all terms | all T for G <= 5 | 53.4 million terms, 0 survivors, one target per case | ThmL/check4_logic.py |
| Theorem L directly on the Hodge side | h = 3 (r <= 4), h = 4 (r <= 5, all in-range T), h = 5 (r <= 3), h = 6 (r <= 2), at tau = 2/7, -5/3, -1/2, 13/9 | 383 tests, all zero; just outside the range all 228 controls non-zero | ThmL_ref2/t1_thmL.py |
| Earlier rounds | genus-7 engines (194 zeros), R3 referees' mechanical logic G <= 12, Mathematica q-tables | see round2_status s.1 | lemmaA/round2/ |

Limits: the g = 6 direct Hodge comparison of (2.1) did not finish. The G = 6-8 mechanical logic run with r <= 2 did not finish. Neither is claimed.
