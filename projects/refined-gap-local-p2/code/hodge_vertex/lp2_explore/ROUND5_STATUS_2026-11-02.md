# Round 5 status (local P2 all-genus programme), 2026-11-02

Workflow run wf_9d7fd94e-8d7, 13 agents, none failed. Full output: round5_scripts/round5_result.json. Scripts of the agents: round5_scripts/ (downloaded literature not copied). Main-session check: round5_scripts/main_session/check_main.py.

## Deliverables

- hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md: new; single proof of record of Lemma A (Theorem L), about 9,700 words.
- LP2_REGULARITY_PROOF_draft.md: edited in place, about 23,500 words. Backup: LP2_REGULARITY_PROOF_draft_round4.md.

## Results and verdicts

| Item | Outcome | Referees | Major items |
|---|---|---|---|
| MVmatch: MV series = Xi_h, all genera and profiles, no fitted constant | proved mod LLZ eq. (1) | 2, sound with minor fixes | 0 |
| ThmL: Theorem L, proof of record = R3 proof (ii), P_k derived for all k | proved mod (MV) | 2, sound with minor fixes | 0 |
| modrel: symplecticity from Theorem A only (Prop 2.10); both relations of lem:modrel (Prop 2.12) | proved, no citation | 1, sound with minor fixes | 0 |
| VI: (VI) in all genera (Theorem VI) | proved mod LLZ eqs. (1), (5), boson-fermion | 1, sound with minor fixes | 0 |
| Final critic | unconditional_mod_literature = true, in the sense below | - | 3 serious items |

## Status now

- Lemma A (Theorem L): P mod (MV) and textbook facts, for non-empty R, all h >= 2, both theories. No TR, no cited KP theorem, no computer step.
- Theorem LC_P2: P, no BS input.
- Theorem R_P2 (graphs with an edge, decoration sum G_g): P mod (MV).
- Corollary SV_P2 (polar part (-3)^{g-1} beta_{2g-2} t_c^{2-2g}): P mod (MV), LLZ (5), boson-fermion.
- BS genus-0 inputs (C1)-(C3) are not used by the three statements above. They are needed only to connect G_g to GW.
- Proposition 13.2_P2 (G_g = GW^cf_g) and Corollary G (the gap): OPEN. Not claimed.
- Scope: fibrewise torus T_F only, per colouring.

## Main-session checks

By hand:
- LEMMA_A s.4: (F2), (F3), (F5), (F6), P_2 at lowest degree directly from Schur functions.
- s.3: Lemma 3.1, Lemma 3.4 step.
- s.5: c(|B|,k), w-order count, single-sector argument.
- s.6: the whole induction (count, case split, target coefficient).
- Theorem VI: uniqueness step (w_n = C tau^{-1} (1+tau)^{n+1}, killed by tau -> -1-tau).
- Z_2 = -(24 + 28 s + 7 s^2)/5760, s = 1/sigma, from the known genus-2 intersection numbers, no MV. It agrees with sigma Z_2 = beta_2.

Exact computation, evidence only (check_main.py, fresh code, tau = 2/7, from Schur functions through LEMMA_A s.8.2):
- Xi_2((2)), Xi_3((2)), Xi_3((3)), Xi_3((4)), Xi_4((2)) = 0.
- Controls outside the range are non-zero: Xi_2((3)) = 1367/186624, Xi_3((5)) = 205027/47029248.
- sigma^{h-1} Xi_h(empty) = beta_{2h-2}(-tau, 1+tau) exactly at h = 2, 3, 4 (with A_g).
- A^2 D = X, 2E(y) = D - X through y^23.

Not re-derived by me: right side of (E1) in Theorem VI beyond the structure and the check above; Theorem A; H2 Steps 1-4 in this round; s.2.2-2.4 of the draft.

## Critic's serious items (all open)

S1. Prop 13.2_P2 is the one missing in-scope deliverable. s.8.3 needs four changes:
 (a) Thm 2.7(e) gives z-only factor exactly 1, so the frame change should be Mumford + PPZ Prop 2.14, as in conifold Prop 13.2(c),(e); less of an obstacle than s.8.3 says.
 (b) New: Prop 2.6 step 5 rests on BS (Normalisation) (line 3497), asserted without proof in BS; it fixes higher-genus constants by constant maps, so it is not on the allowed genus-0 list. Port conifold Prop 13.2(c) at generic alpha.
 (c) New: "if delta divides Q, some phi_{k,m} blows up" in Prop 2.11 needs algebraic independence of A, F over functions rational in y (prop:algindABC or BFGW Lem 4.2). s.0.3 wrongly says (C4) only makes Cf well defined.
 (d) New: Cf on coloured terms (these involve t = y^{1/3}): A, F stay independent over C(t).
 Also open: BS line 4056; Dil = Givental translation check. The critic's items (a)-(d) are its own reading, unverified.

S2. Neither assembled document has had a referee of its own. Labels should not change until they have.

S3. The conifold upgrade (Corollary R unconditional mod (MV), and genus 3-4 no longer computer-assisted) needs written fixes in leaf_structure_proof.md:
 (a) s.13.1(i) F_h(lambda q) needs the domain fix (draft s.1.5 and (H0), c -> E).
 (b) s.13, proof_status s.1.2, the note: add "non-empty" to M_irr.
 (c) Teleman citations wrong against arXiv:0712.0160v3: classification is Prop 8.7(ii), Hodge twist (8.6) + Prop 8.7(i), R-ambiguity Prop 8.4; add a field-embedding line.
 (d) Restate input (5) and step (v) of Theorem R: the single pole can come from Theorem VI.

## Minor items carried forward

- LEMMA_A_PROOF.md stale against the draft: O5 and s.8.5 (gap formula, TR-free), s.7.5 and O5 (Taylor identification, modrel). Sync.
- Draft labels R_P2, Lemma A_P2(b) to "P mod (MV)" after referee.
- Theorem VI Step 2(b): replace boson-fermion citation with LEMMA_A s.4 (4.2) at k = 2; Steps 1, 3 can cite LEMMA_A Prop 2.1, Lemma 3.2.
- Draft s.7.5(ii) last bullet vestigial.
- Xi_h(empty) != 0 holds in Q(tau) but does not specialise: Xi_2(empty) vanishes at sigma = -0.804, -0.363 (complex tau).
- Unconfirmed citations: JDG numbering of LLZ; Macdonald I.3 Ex. 2; Stanley EC2 7.21.3; Graber-Pandharipande, Mumford, DLMF 5.11.1 not opened. LLZ's proof of (MV) not re-checked.
- OP (0.5), (0.8): printed 1/|Aut mu| inconsistent at repeated parts (mu = (1,1)); OP (1.24), (1.28) have no such factor.
- LLZ Prop 2.1: printed 2^l should read 2^{|nu|}.
