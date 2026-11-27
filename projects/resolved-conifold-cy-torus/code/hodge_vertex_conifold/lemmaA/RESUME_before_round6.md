# RESUME NOTE: Lemma A / all-genus regularity, state at wrap-up (2026-10-31, about 21:40)

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
   Unverified points: Teleman was read through tool summaries; the CohFT axioms for non-compact Z go through the twisted P^1 theory.
4. Doc fixes, after the referees:
   - M_irr only below the genus: leaf_structure_proof.md s.0.4, s.13.1, s.13.2 and G1; proof_status_mixed_scale.md s.3 and Addendum 2;
     conifold_regularity_note.md s.4 and s.7b.
   - Cite Brini-Schuler Prop. 3.1 for (H) at s.13.3(5).
   - Record round 2's outcome (K1; whichever all-genus route survives) in the notes, README and memory.
5. Refresh the tarball (hodge_vertex_conifold.tar.gz in the refGW folder; files over 50 MB are left out, see below).

## Today's main-session checks (all passed)
- Genus-5 Mathematica cross-checks at tau = 2/7 (11 values) and tau = 3/2 (5 values), controls included.
- Genus-4 two-leaf and genus-3 three-leaf symbolic comparison (869 coefficients, 0 mismatches).
(round1_status.md s.6; E_formulas/E_crosscheck_g5_M1_out.txt, E_crosscheck_g5_M2_out.txt, compare_g4r2_out.txt.)
