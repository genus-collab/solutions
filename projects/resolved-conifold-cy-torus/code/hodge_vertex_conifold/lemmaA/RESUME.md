# RESUME NOTE: Lemma A / all-genus regularity, state at wrap-up (2026-10-31, about 21:40); resume point updated 2026-11-02 (round 6)

## RESUME POINT, 2026-11-02 (round 6; supersedes everything below where it differs)
- Lemma A is now lemmaA/LEMMA_A_PROOF.md (Theorem L; s.7.4 for the conifold convention), P modulo (MV) and textbook facts. Frozen copy refereed in round 6: LEMMA_A_PROOF_round5.md (two referees, sound with minor fixes, no major item; their minor items are about labels, stale pointers and citations in that file).
- Corollary R (leaf_proof/leaf_structure_proof.md s.13.2) is proved in every genus, as a statement about the decoration sum, modulo (MV) and textbook facts; its identification with F_g is Proposition 13.2 [repaired in round 6d; passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md], P modulo the literature inputs listed there; for the value of the pole also the further inputs of Theorem VI (LP2 draft s.7.8: lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2) at k = 2, with the Frobenius formula, and (MV) only at mu = 1^a; LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; two referee passes (round 5; round 6, D_VI)). Genus 3 and 4 are no longer computer-assisted. Graphs with an edge need Lemma A only for h <= g - 1.
- Round-5 critic item S3 applied to the conifold notes (backups *_before_round6.md): s.13.0 of leaf_structure_proof.md ports Fpot_h on a domain and (H0) (S3(a)); "non-empty" in M_irr (S3(b)); Teleman citations checked against the arXiv:0712.0160v3 text, field embedding added (S3(c)); input (5), step (v) of Theorem R restated via Theorem VI, [BS] Prop. 3.1 as a check (S3(d)).
- One referee pass (round 6c): s.13.0, s.13.3 (5). The field-of-definition paragraph of Proposition 13.2 had one pass in its round-6c form, which found the sign error (the family is at u = -s H). It was then rewritten with the direct derivation of FRAMEB_PROOF.md A2.8, with Teleman's Prop 7.7 only in a remark: [repaired in round 6d; passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md].
- The repaired paragraph of Proposition 13.2 has had its round-6d referee pass (refCon2; hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md). AI pass only; no human referee has read it.
- Next: proof_M_residue (round2_status s.5.5), README.txt, outline.md and the memory notes still carry the old status; hvc_lib.m line ~181 still calls the closed form a prediction; Teleman's journal numbering not compared with arXiv v3.
- Local-P2 status and open items (Prop 13.2_P2, gap corollary): hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md (full round-6 status; ROUND5_STATUS_2026-11-02.md for the rounds before). Round 6e: eq:J=I at generic alpha on T_max is now cited to CCIT15 Thm 31 (FRAMEB_PROOF.md s.1.2, item O1), so Prop 13.2_P2 and the gap corollary are P mod eq:J=I (CCIT15 Thm 31) and the other inputs listed there. The T_max use is outside the literal wording of ruling Q1, which names statements on T_F. The user confirmed it on 2026-11-02, accepting CCIT15 Theorem 31 at generic weights under that ruling.

Read first: lemmaA/status_2026-10-31.md (audited status), lemmaA/round1_status.md (round-1 synthesis, plus s.6 with today's Mathematica
checks), lemmaA/round2/brief.md (the rules and the known results given to the round-2 agents).

## Rules (user)
- No topological recursion for the resolved conifold x C^2 theory itself. TR for the Hodge CohFT (curve C_tau, BKMP for C^3, P2, Lemma D) is
  allowed, and the KP / 2D Toda hierarchy of the Marino-Vafa series is an allowed input.
- The frame-B formula is Givental-Teleman reconstruction (user, 2026-10-31): call it standard, with citations pending, not an open premise.
- Agents run no wolframscript; the main session runs the Mathematica checks, at most about 2 kernels at a time. No python (use Perl/BigRat).
  Never kill all Wolfram processes by name; stop only our own jobs, by PID.
- Generic numerical checks must use non-special framings: in the E_crosscheck*.wls scripts tau = -x/(1+x), so take x < 0 (e.g. x = -2/9,
  tau = 2/7).

## In flight at wrap-up (Claude Code session 37f24460-3601-4c85-b123-5156e19eb475)
Journals: <user home>/.claude/projects/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/subagents/workflows/<run>/journal.jsonl
Scripts:  <user home>/.claude/projects/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/workflows/scripts/<name>-<run>.js
- wf_e42146b3-114 (lemmaa-all-genus-step): routes R1-R4 done (report.md saved in round2/R1_rigidity, R2_onepoint, R3_toda, R4_geometry).
  FINISHED. All referees SOUND with no major issue (R1, R2, R3 two each; R4 one). Critic saved as round2/round2_status.md: Lemma A
  proved modulo literature in all genera (R3 = proof of record), hence Corollary R in every genus; its gaps and plan are in s.5-s.6.
- wf_b8c3aa8f-c61 (lemmaa-tr-hodge-route): routes TR1, TR2 done (report.md saved in round2/TR1_dictionary, TR2_recursion). Referees done
  and saved, all SOUND: TR2 (both), TR1 (first). PAUSED by the user on 2026-11-01 (about 00:45), then RESUMED and FINISHED: referee2 TR1 = SOUND, saved as TR1_dictionary/referee2.md: the workflow was stopped while referee2 TR1
  was running, and its three Perl jobs (r2_big.pl, TR1_dictionary/referee2) were stopped by PID. Partial outputs are in
  TR1_dictionary/referee2/runs/save_tau*_S9.txt. To resume in the same session: Workflow with scriptPath = the lemmaa-tr-hodge-route script
  and resumeFromRunId = wf_b8c3aa8f-c61; the finished agents replay from cache and only referee2 TR1 re-runs. In a new session: start a new
  workflow that runs only the referee2 TR1 prompt from that script. TR1's second referee is not needed for the round-2 conclusion: R3 is the
  proof of record. No critic in this workflow.
- wf_de5bc7e9-d4c (frame-b-givental-teleman-paragraph): FINISHED; inserted as Proposition 13.2 / Lemma 13.3 (see step 3 below).
- Background Perl job still running: PID 13332 = fitNR.pl 3 4 0 (a Route E fit from round 1); it can be left or stopped by PID.
  mirr_points.pl 7 FINISHED: the author's engine now also has all 194 genus-7 points zero, and all 77 controls non-zero (E_formulas/mirr_g7.txt).

## How to collect results
- Same session still open: the completion notifications arrive; save each new referee result with
    perl lemmaA/tools/save_r2.pl <journal> lemmaA/round2 "referee R1=R1_rigidity/referee.md" "referee2 R1=R1_rigidity/referee2.md" ...
  (it refuses to overwrite existing files). lemmaA/tools/dumpres.pl dumps one journal line as markdown.
- New session: the journals above are plain files, and every finished agent's return value is in a {"type":"result"} line (labels in the
  matching {"type":"started"} line). Save them with save_r2.pl. Workflow resume (resumeFromRunId) works only inside the original session, so
  re-run any agent that never finished with a new workflow; its prompt is in the script file.

## Next steps (in order)
1. All referees are done. Follow the critic's plan (round2/round2_status.md s.6): P1 one written
   all-genus proof with R3 as the record (non-empty R); P2 pin the literature; P3 integration critic across both workflows; then run an integration critic across
   both round-2 workflows, focusing on whichever claims survive (independent lines: R3 = 2D Toda; R1 = KP + conifold polar relations;
   R2 / TR2 + TR1 = ABDKS x-y swap + P2).
2. Mathematica checks requested by the round-2 writers (read each script first; at most 2 kernels):
   round2/R1_rigidity/R1_checks.wls (C1-C6), round2/R3_toda/R3_thmM_check.wls, round2/TR2_recursion/tr2_closed_check.wls, TR1's request
   (hvc_tr_tau.wls at a non-integer tau, genus 3-4), R2's request (hvc_tr_tau.wls at tau = 2, 3, genus 4-5), R4's Xi^(1)_10 check against
   its C1 (no script yet). Exact requests: requested_checks in each report.md.
3. DONE (2026-10-31, about 23:40): the frame-B paragraph is in leaf_structure_proof.md as Proposition 13.2 and Lemma 13.3 (s.13.3),
   with gap (G3) replaced and companion lines in s.12.1 and (G2). Workflow wf_de5bc7e9-d4c: writer, two referees (both "minor revision";
   the major point, that Givental 2001 covers compact targets only, is fixed by resting the proof on Teleman plus localisation at Q = 0), reviser.
   Pre-edit copy: leaf_proof/referee_scripts/leaf_structure_proof_before_frameB.md. Still to do: Appendix A notation (wt_+-, nw, Nb, Kw, sig,
   PF, km, Hodge_v, R_gen); the F0closed comment in hvc_lib.m (line ~177) calls F_g(0) a prediction, but it is now proved for all g >= 2;
   misprints in Brini-Schuler App. A to report to the authors (s^2 d_s^2 W = kappa Delta should read Delta = kappa W_2; the normalisation display).
   Unverified points: Teleman was read through tool summaries [2026-11-02: arXiv v3 text read, citations corrected]; the CohFT axioms for non-compact Z go through the twisted P^1 theory.
4. Doc fixes, after the referees:
   - DONE 2026-11-02: M_irr only below the genus: leaf_structure_proof.md s.0.4, s.13.1, s.13.2 and G1; proof_status_mixed_scale.md s.3 and Addendum 2 (now Addendum 3);
     conifold_regularity_note.md s.4 and s.7b.
   - DONE 2026-11-02, differently: s.13.3(5) now rests on Theorem VI; Brini-Schuler Prop. 3.1 is cited as an independent check.
   - Record round 2's outcome (K1; whichever all-genus route survives) in the notes, README and memory.
5. Refresh the tarball (hodge_vertex_conifold.tar.gz in the refGW folder; files over 50 MB are left out, see below).

## Today's main-session checks (all passed)
- Genus-5 Mathematica cross-checks at tau = 2/7 (11 values) and tau = 3/2 (5 values), controls included.
- Genus-4 two-leaf and genus-3 three-leaf symbolic comparison (869 coefficients, 0 mismatches).
(round1_status.md s.6; E_formulas/E_crosscheck_g5_M1_out.txt, E_crosscheck_g5_M2_out.txt, compare_g4r2_out.txt.)
