# Round 6 status, 2026-11-02

File: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md. Started by the round-6d fix step so that the mark in the conifold files has a target; extended in round 6e into the full status of round 6. Proof of record: FRAMEB_PROOF.md ("FB"), with LP2_REGULARITY_PROOF_draft.md ("the draft") and hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md. Conifold: hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md s.13. BS is BS_v2_2026-11-02/main.tex, not edited.

"P mod X" means proved in these files given X, and nothing beyond X. **All referee passes of round 6 are by AI referees, not human ones. The proof has not been read by a human referee.** Finite-order computations are evidence only and are never proof steps.

## 1. What is proved, and modulo what

Scope: local P2 x C^2 on the fibrewise torus T_F (eps1 eps2 E != 0), per colouring. Nothing is claimed off T_F.

**Corollary G (refined conifold gap, BS rmk:cfgap).** For g >= 2, GW^cf_g := Cf(F_g^{T_F}) = (-3)^{g-1} beta_{2g-2}(eps1,eps2)/t_cf^{2g-2} + O(1), where the O(1) term is holomorphic near t_cf = 0 and its value is not determined. Corollary G_{k,h} is the same statement for the coefficients F_{k,h}. Cf is the reading of BS's substitution adopted as a definition (FB s.0.2, Definition (Cf)); it is a substitution on the polynomial in F, not an analytic continuation of the series.

The chain, with the status of each link:

| Statement | Status |
|---|---|
| Theorem FB_A (CohFT at generic alpha, Teleman's R, its constants, graph sum) | P mod eq:J=I (CCIT15 Thm 31) and the genus-0 QDE (derived in FB s.1.2) for (A2), (A3) and the row form of the graph sum; Teleman, PPZ, GP, virtual-class axioms, Mumford, FP Thm 4 (FB s.0.2) |
| Theorem FB_B (ring map to T_F; replaces draft Prop 2.11) | P mod Theorem FB_A and draft statements, part (b) without eq:J=I |
| Theorem FB_C (C1), (C2) (frame change Tbar_g[𝖥] = G_g[𝖥]) | P mod draft Prop 2.10, Mumford as in FP, PPZ; does not use eq:J=I |
| Proposition 13.2_P2 (Cf(F_g^{T_F}) = G_g[F_cf]) | P mod Stages A, B, draft Props 2.10, 2.12, PPZ; logarithmic form of Lemma (t_cf) |
| Lemma (t_cf) (t_cf = t_c) | P mod DLMF 15.6.7, Barnes's second lemma (primary source not opened), Gauss's values of psi |
| Theorem LC_P2 | P (draft), statement about the decoration sum |
| Theorem R_P2, Corollary SV_P2, Theorem VI | P mod (MV) = LLZ eq. (1) and textbook facts, statements about the decoration sum; through LEMMA_A_PROOF.md Theorem L |
| Corollaries G, G_{k,h} | The sum of the above, with Cor 2.8 and Theorem FB_B(e) |

End-to-end inputs of Corollary G (FB s.0.3): eq:J=I (CCIT15 Thm 31) and the genus-0 QDE (derived in FB s.1.2 from Pandharipande math/9806133v1 Lemma 1 and A1.4); Teleman's Theorem (2), Prop 8.7(ii), 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss); PPZ Defs 2.2, 2.5, 2.13, Prop 2.14; Graber-Pandharipande; the fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist; Mumford as in FP, and FP Theorem 4; DLMF 15.6.7, Barnes's second lemma, Gauss's digamma values; (MV) with the textbook facts of FB s.1.3. The proofs of these cited theorems were not re-derived in any pass.

**Proviso.** Every statement of the chain that rests on eq:J=I also carries the Q1 proviso of s.3.1 below: A2.1, A2.4-A2.6, A2.9, A3, parts (A2), (A3) and the row form of the graph sum of Theorem FB_A; Theorem FB_B except (b); Proposition 13.2_P2; Corollaries G, G_{k,h}. Not affected: A1, A2.3, A2.8, A3.2, Theorem FB_C (C1), (C2), Lemma (t_cf), and Theorems LC_P2, R_P2, Corollary SV_P2, Theorem VI as statements about the decoration sum (dependency reading of the current text; no referee has checked this partition).

**Conifold.** Corollary R (leaf_structure_proof.md s.13.2) is proved in every genus as a statement about the decoration sum, mod (MV) and textbook facts. Its identification with F_g is Proposition 13.2, P mod the literature inputs listed there; see the last section of this file and lemmaA/RESUME.md. Item O1 does not concern it.

## 2. Rounds 6a-6e

Verdict "sound" abbreviates "sound with minor fixes, no major item". Reports: round6_scripts/round6a_result.json, round6b_result.json, round6c_result.json, round6d_result.json; the detailed record of 6c-6e is the marker ROUND6C-REFEREES in s.0.4 of the draft.

| Round | Work | Pass | Verdict |
|---|---|---|---|
| 6a | LEMMA_A_PROOF.md | L_logic, L_comp (independent computation) | sound, both |
| 6a | Draft s.0-s.2; s.3-s.6; s.7.0-s.7.7; s.7.8 | D_rows; D_core; D_hodge; D_VI | sound, all four |
| 6b | Prover reports FB_A, FB_B, FB_C | refFB_A, refFB_A2 (Stage A); refFB_B (conditional on Stage A); refFB_C; refFB_C2 (conditional on Stages A, B) | sound, all five |
| 6c | FRAMEB_PROOF.md assembled; draft integrated; A2.8 rewritten after refA | refA (Stage A, A2.8 in its earlier form) | sound |
| 6c | | refBC (Stages B, C) | sound |
| 6c | | refS2 (s.2.1-2.5 as used by FB:A3.5, B2, B4) | sound |
| 6c | | refGap (Corollary G, restricted as its report says) | sound |
| 6c | | critic (assembled text) | sound |
| 6c | Conifold Prop 13.2 and s.13 | refCon | **major_gap**, one item: sign in the Prop 13.2(a) field-of-definition paragraph. Repaired in round 6d |
| 6d | A2.8 as rewritten | ref28hand (by hand); ref28num (finite order, evidence only) | sound, both |
| 6d | Text written after round 6c | refNew | sound |
| 6d | Repaired rank-2 paragraph | refCon2 | sound |
| 6d | Citations | lit | **major_gap**, one item at citation level: BS's source for eq:J=I covers T_F only (item O1) |
| 6d | Forward cold read of FRAMEB_PROOF.md | coldF | sound |
| 6d | Backward read of the chain of Corollary G | coldB | sound |
| 6e | Source of eq:J=I at T_max, its application to Z, QDE derivation | verify:ccit (by hand; exact check through d = 4 as evidence only) | sound; CCIT15 Thm 31 adopted as source of record |

No pass covers the proofs of the cited theorems. The journal print of CCIT15 was not read; numbering is from arXiv v5.

## 3. Items open

### 3.1 Open item: ruling on the T_max use (item O1)

Item O1 (FRAMEB_PROOF.md, head; report lit, round 6d; citation repaired in round 6e by report verify:ccit). BS 2658 cites eq:J=I to MR2486673 = Coates, arXiv:0810.2200v1, which covers T_F only. Stage A cannot be moved to T_F, because the Q = 0 product is not semisimple there. The source of record at T_max is eq:J=I (CCIT15 Thm 31): Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912, arXiv:1310.4163v5, Theorem 31 (FB s.1.2). The genus-0 QDE at T_max is derived in FB s.1.2 from Pandharipande's Lemma 1 and A1.4. The T_max use is outside the literal wording of ruling Q1, which names statements on T_F, and awaits the user's confirmation. It does not concern the conifold Proposition 13.2.

### 3.2 Not proved (FB s.0.3)

- Anything off T_F.
- The value of the regular part of GW^cf_g.
- BS's F_{k,h} in Q[X^{+-1},F]. Only F_g^{T_F} in K[y^{+-1},Delta_c^{-1}][F] is proved. Open for all g, independent of Corollary G.
- (Orbifold regularity) of BS def:dirint.
- The holomorphic anomaly equations.
- Uniqueness in BS def:dirint.
- The proofs of (MV), of Teleman's Theorem (2), of Harer stability, Looijenga and Madsen-Weiss, and of eq:J=I (CCIT15 Thm 31), which stay cited.

### 3.3 Journal requests of the critic (round 6c) that remain

- One consolidated paper replacing the three linked files, with one notation, no round bookkeeping, no code-reading Lemma 1.4, and Corollary G with all hypotheses in one place; a displayed genus-2 example; a deposit of scripts separate from the proof (critic 25).
- The relation of Cf to a change of polarisation of the Givental quantisation (critic 26). Not addressed.
- Evidence at g >= 5. No full total has been computed at g >= 5. What exists is partial (FB s.A.4): g = 5 per colouring for graphs with a vertex of genus >= 2 at three eps pairs; g = 6 for about 1000 of 5226 such graphs at one eps pair, with an unfinished log; the single-vertex polar part through g = 11 at (3,7). All of it is mod p or at finite order, and it proves nothing.
- A printed source for the form of Barnes's second lemma used in C Step 6 (critic 23; FB s.1.2 says this should be done before submission).

Requests closed in rounds 6c-6e: a pinned theorem for eq:J=I at T_max (critic 21, now CCIT15 Thm 31); Teleman's Prop 7.7 as an input (critic 22, removed by the direct derivation in A2.8); comparison with Lho-Pandharipande (critic 22, FB s.1.2, comparison only); a self-contained row 2 of eq:RHb (critic 24, A2.1).

### 3.4 Housekeeping

lemmaA/RESUME.md lists the stale conifold notes (proof_M_residue, README.txt, outline.md, memory notes, hvc_lib.m). Teleman's journal numbering has not been compared with arXiv v3.

## 4. BS errata

The full text and the reasons are FB s.5.1; none is applied to main.tex. Headlines:

1. Lines 1277, 2446 versus eq:gtH (3799), eq:F2ref (4217): GW_g printed with d > 0, while the graph sum includes d = 0.
2. Line 3537: f_j^* must be the constant covector at Q = 0.
3. eq:normNLR (4327) with the printed Asym: holds with a constant c, c^2 = -1.
4. Lines 4054-4056: regularity on the diagonal asserted in one clause; the localisation needed is not written (proved another way, B5).
5. Lines 4196-4198, 4345-4347, 3983, 4056: ansatz, "by continuity", rational tori and closure are not needed (replaced by B2, B4, B6).
6. Lines 5163-5167: Cf(f)(tau) = f(-1/(3 tau)) is false as written.
7. Line 5145: t_cf needs the branch log y = log(-y) + i pi.
8. prop:algindABC (2854): attribution to BFGW Lem 4.2 (the independence is in their s.4.1).
9. Lines 4071-4076: regularity at y = 0 is argued for coefficients, while what is given is the total, so Q[X^{+-1},F] does not follow as written; no counterexample is known.
10. About lines 3382-3384: printed sign of W_i on T_F (draft s.9 item 1). The other items of draft s.9 stand.
11. Lines 927-929: the expansions of Pi_gamma are normalisation statements, up to a constant prefactor BS does not state.
12. Lines 2633, 2658 (and the commented pointer at 2659): MR2486673 gives eq:J=I on T_F only. The correct citation at T_max is MR3414388, Theorem 31.

