# Lemma A (M_irr(h) for all h): shared brief for the first attack (2026-10-31)

Project folder: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold (call it HVC below). Work files go in HVC\lemmaA\<route>\.

## 0. The rule that overrides everything
The user said: **"no TR at all"**. Do NOT use topological recursion, spectral curves, Eynard-Orantin invariants, loop equations,
the remodelling/BKMP correspondence, or the "dictionary P2" in any argument, and do not use any result that was proved only
"mod P2" (proof_status_mixed_scale.md s.1.6, s.1.7(b),(c), s.1.8 and the parts of s.2 about Lemmas B and D), nor the TR scripts
(review_checks/status_tr*.pl, referee1_tr*.pl, hvc_tr_tau.wls, hvc_moment_test.wls). Not even as a cross-check.
Allowed (cite precisely, with statement and normalisation): the Marino-Vafa formula for one-partition triple Hodge integrals
(Liu-Liu-Zhou 2003, Okounkov-Pandharipande 2004) and its cut-and-join equation; ELSV; Hurwitz theory (Burnside/character
formulas, cut-and-join, shifted symmetric functions, completed cycles); GRR, Mumford's relations, Chiodo's formula, r-th roots
and polynomiality in r (Janda-Pandharipande-Pixton-Zvonkine, Clader-Janda), Porteous and Fulton's excess intersection; virtual
localization relations (Graber-Pandharipande, Faber-Pandharipande); Lagrange inversion; generating functions and exact algebra.

## 1. The statement
Hodge class (mode B): Omega = prod_{a=1}^{3} Lambda^vee_h(w_a)/w_a^h, w = (c, -x, -1), c = 1 + x (Lambda^vee_h(u) = Sum_i (-1)^i lambda_i u^(h-i)).
Leaves lam_j(psi) = 1/(c - j psi); dilaton leaves T_c(psi) = -psi^2/(c - psi). For a multiset R = {j_1..j_r} of integers:
    Xi_h(R) := Sum_{m>=0} (1/m!) Int_{Mbar_{h,r+m}} Omega prod_{i=1}^r lam_{j_i}(psi_i) prod_{k=1}^m T_c(psi_{r+k})
(a finite sum: T_c has psi-degree >= 2). D(R) := Sum (j_i - 1), J := Sum j_i.
**Lemma A = M_irr(h): for every h >= 2 and every multiset R of integers >= 2 with D(R) <= 2h - 3: Xi_h(R) = 0 identically in x.**
(Then r <= 2h - 3.) Known: true for h = 2 (proved by hand), h = 3 and h = 4 (exact Hodge-integral computations). Needed: all h
(it is the only remaining input of the all-genus regularity theorem of leaf_proof/leaf_structure_proof.md s.13).
It is equivalent to the vanishing direction of "Conjecture M": Xi_h(R) = 0 for all integers j_i >= 0 with 1 <= J and D <= 2h - 3
(string and dilaton relations reduce to parts >= 2). Non-vanishing outside the range is PROVED (s.1.3 below), so only vanishing is open.

Normalisation c = 1 (used in the notes): weights (1, tau, -1-tau) with tau = -x/(1+x), sigma := tau(1+tau) = -x/c^2, leaves
1/(1 - j psi), T = -psi^2/(1 - psi); Xi^(c-norm) = c^-(3h-3+2r) Xi^(c=1). The Hodge class is then the Marino-Vafa class
Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau) up to the factor (-sigma)^-h, and the leaves are exactly the Marino-Vafa leaves
1/(1 - mu_i psi_i) with mu_i = j_i: Xi is a one-partition CY-vertex Hodge integral EXCEPT for the dilaton leaves.

## 2. What is proved without TR (read these sections; do not re-derive them, cite them)
proof_status_mixed_scale.md s.1.1-1.5 and s.1.7(a); proof_M_residue.md s.1-9 (Lemma 1, Lemma 2, Lemma 3, Propositions 3, 4, 6,
Theorem 7, the conjectural closed formula (6.1) with its equivalent forms s.6a, the exploratory s.8); conifold_regularity_note.md
s.3a (Theorem C) and s.7. Summary:
 (K)  kappa form: Sum_m 1/m! pi_*(prod T_c) = K := exp(-Sum_k kappa_k/k) (c = 1), so Xi_h(R) = Int_{Mbar_{h,r}} Omega K prod lam_{j_i}.
      Interpolation: K_s := exp(Sum_k ((1+s)^k - 1) kappa_k/k) comes from the dilaton leaves s psi^2/(1 - psi) = s(lam_1 - 1 - psi);
      s = -1 gives K, s = 0 gives no dilaton leaves (pure Marino-Vafa integrals).
 (L)  Lambda^vee(1) K = c(T(-log)) := c(-R pi_* omega_log^{-1}), the log tangent bundle of Mbar_{h,r} (GRR).
 (S)  c(T(-log)) prod 1/(1 - j_i psi_i) = c(-R pi_* S_D) prod p_{j_i}(psi_i), S_D = omega_log^{-1}(Sum j_i x_i), p_j(psi) = prod_{k=1}^{j-1}(1 - k psi);
      rank(-R pi_* S_D) = 3h - 3 - D, and the range D <= 2h - 3 is exactly deg S_D < 0 on smooth curves. If S_D had no sections on any
      stable curve the integral would vanish by degree (Lemma 3); jump loci (rational tails etc.) are the obstruction.
 (C)  components: with alpha_i the Chern roots of E and c = 1: Xi_h(R; tau) = Sum_{k=0}^{h} (-1/sigma)^k Xi^(k)_h(R),
      Xi^(k)_h(R) := Int c(T(-log)) e_k(alpha + alpha^2) prod lam_{j_i}, tau-INDEPENDENT. k = 0: the single class Lambda^vee(1)
      (sigma -> infinity, "Hurwitz component"); e_1 = lambda_1; k = h: c(T(-log)) e_h = lambda_h K (Mumford), the lambda_g case.
      M_irr(h) <=> Xi^(k)_h(R) = 0 for all 0 <= k <= h in the range.
 (SD) string/dilaton: Xi(1, R) = J_R Xi(R), Xi(0, R) = (D_R - 2h + 2) Xi(R) (c = 1); reduction to parts >= 2.
 (G)  lambda_g component, all h, r, PROVED: Xi^(h)_h(R) = b_h (J-1)(J-2)...(J-N), N = 2h-3+r, b_h = [z^(2h)] (z/2)/sin(z/2); it vanishes in the range.
      Hence sigma^h Xi_h(R) is a polynomial in sigma of degree <= h with zero constant term in the range; h independent linear
      conditions on Xi^(0..h-1) suffice (s.1.7(a)).
 (H)  Hurwitz component proved for j = 2 (all h >= 2) and j = 3 (h >= 3) (localized top Chern class + dimension of meromorphic strata);
      lambda_1-component proved for j = 2. Open: every other (j, k) for r = 1, and everything with r >= 2.
 (6.1) conjecture (r = 1, all tau), exact for h <= 4 against Hodge data:
      sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^(2h)] S(e1 lambda)^(1 + j e2) S(e2 lambda)^(1 + j e1),
      e1 = -tau, e2 = 1 + tau, S(y) = (y/2)/sinh(y/2); Hurwitz limit via Lagrange inversion in proof_M_residue.md s.6a(iii).
 No closed formula is known for r >= 2 (s.8 of proof_M_residue.md lists two failed guesses).

## 3. Exact data (Hodge side, computed by Mathematica with Mumford's formula; no TR) -- use them to fix every normalisation
 READ HVC\lemmaA\data\WARNING_READ_FIRST.txt FIRST: tables written before 12:01 were wrong; all were regenerated at 12:01.
 HVC\lemmaA\data\psi_g<g>.txt   g | n | {a_1..a_n} | <Omega psi^a>_g      (c-normalisation, x symbolic; g = 1, 2: n <= 5; g = 3: n <= 4; g = 4: n <= 3)
 HVC\lemmaA\data\mv_g<g>_l<l>.txt  <Omega prod lam_{m_i}(psi_i)>_g as a polynomial in m1..ml (no dilaton leaves), g <= 4, l <= 3
 HVC\lemmaA\data\q_g<g>_r<r>.txt   Q = Sum_k s^k/k! <Omega prod lam_{j_i} Tt^k>, Tt = psi^2/(c - psi): polynomial in s, j1..jr, g = 2, 3, 4,
                                r = 1, 2 (s = -1 gives Xi_g(j1..jr); s = 0 gives the plain Marino-Vafa-type integral)
 HVC\hvc_multiprobe_out.txt (genus 2-3 lattice tests and symbolic r = 2 polynomials in jj[1], jj[2], x), hvc_multiprobe_g4_out.txt
 (genus 4, r = 1 symbolic, r = 2 lattice), hvc_multiprobe_g4b_out.txt (genus 4, r = 3..5 points), hvc_which_class_out.txt,
 check_closed_formula_r1.pl ((6.1) against the genus <= 4 data), hvc_mp_header.m (the Mathematica definitions).
 Values are Mathematica InputForm (rational functions of x; ^ for powers). Perl with Math::BigRat is available (no python on this
 machine); a small parser for these expressions is easy (or evaluate at rational x and compare numbers exactly).

## 4. Working rules
 - Do not run wolframscript (licence: about two kernels, used by the main session). If a Mathematica computation is needed, write the
   script under lemmaA\<route>\ and list it in requested_checks; the main session runs it.
 - Perl + Math::BigRat for exact checks; keep scripts and outputs in lemmaA\<route>\. Write files with the Write tool (long shell
   heredocs have failed on this machine). ASCII only. Never kill processes. Do not edit files outside lemmaA\<route>\.
 - Every claim gets a status: proved (complete argument written out), proved-modulo-literature (cite the exact statement used and
   check its normalisation against the data), computer-verified (state the range), conjectural, or refuted.
 - Aim: a complete rigorous proof of Lemma A for all h. If that is out of reach, the strongest rigorous partial results (for example
   one component for all R and all h, or all components for r = 1) and a precise reduction of what remains, with the exact obstruction.
