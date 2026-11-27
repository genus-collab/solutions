# DIGEST_L_2: s.1, what is known (part 1): the conifold gap, the quintic, the formal quintic

Continuation of `DIGEST_L.md` (reading codes [A], [H], [P], [abs], [bib] defined there). Each claim: content; pinned source; status in
that source; how I read it; LIT check (`stage5/lit/LIT*.md`); use in the note. Bibliography keys are those of s.3 (`DIGEST_L_5.md`).

## 1.1 The conifold gap: origin (physics)

* **L1. Leading order.** The order `t_D^(2-2g)` of the leading conifold singularity of `F_g` goes back to BCOV (hep-th/9309140), as
  attributed by HKQ s.1, the sentence after (1.2) ("the order of the leading term was established in [22]"). Physics. Read: HKQ [P];
  BCOV [abs] only. Use: introduction, one clause, cited through HKQ's attribution (key Bershadsky:1993cx).
* **L2. Universal coefficient.** Ghoshal-Vafa (hep-th/9506122; Nucl. Phys. B453 (1995) 121-128): the non-critical `c = 1` string at the
  self-dual radius is equivalent to the topological string near a conifold point; its Penner free energy gives the universal behaviour
  of the topological partition function near the conifold (abstract). HKQ s.1 credits [34] = Ghoshal-Vafa with the coefficient of the
  leading term. Physics. Read [abs]. Key Ghoshal:1995wm.
* **L3. The gap.** Huang-Klemm (hep-th/0605195v2; JHEP 0709:054 (2007)), s.2.3 "Fixing the ambiguity" (p. 10-11): in Seiberg-Witten
  theory, the absence of the terms `1/a_D^k`, `k = 0..2g-3`, in `F_D^(g)` ("the gap"); the `2g-2` gap conditions fix the ambiguity;
  "the gap occurs also at the conifold in the local `O(-3) -> P^2` geometry"; extension to global models deferred to HKQ. Physics,
  observed. Read [P]. New key HK07.
* **L4. HKQ, the compact case.** Huang-Klemm-Quackenbush (hep-th/0612125v1; Lect. Notes Phys. 757 (2009) 45-102):
  * eq. (1.2): `F^(g)_conifold = (-1)^(g-1) B_2g / (2g(2g-2) t_D^(2g-2)) + O(t_D^0)` "in special local coordinates `t_D`";
  * s.2.2, eq. (2.29): a string-theoretic explanation (one light hypermultiplet): "This is precisely the gap condition";
  * s.3.5 "Expansions around the conifold point" and eq. (3.79): for the quintic "we expect this gap structure to be present at higher
    genus", in coordinates `t^_D`; it would impose `2g-2` conditions on the ambiguity; s.3.5: the gap is independent of the choice of
    polarisation;
  * s.3.6 ("Fixing the holomorphic ambiguity: a summary of results": the count of conditions), s.5.1 and the abstract: with the
    constant maps, orbifold regularity and the Castelnuovo bound, the quintic is determined up to genus 51 (s.5.1: "We claim that this
    bound is g <= 51").
  Status: physics conjecture, used as a boundary condition. Read [P] (signs of (1.2) as printed, minus signs partly lost in my
  extraction; the form above is the standard one and is confirmed by the same constant in CGYZ (L8)). New key HKQ09.

## 1.2 The quintic threefold (mathematics)

* **L5. GJR 2018: the gap is the open one of four predictions.** Guo-Janda-Ruan, arXiv:1812.11908v1 (only version; no journal
  reference on the abstract page). Abstract: the B-model predictions (i) Yamaguchi-Yau finite generation, (ii) holomorphic anomaly
  equations, (iii) orbifold regularity, (iv) the conifold gap; "This article is devoted to proving first three conjectures." Remark 1.6:
  the gap would determine the upper `2g-2` coefficients of `f_g`, reducing the initial conditions to `floor((2g-2)/5)`; "It is natural to
  speculate that our approach also implies the conifold gap condition. We leave this for future research." Remark 1.7: the parallel
  approach of Chang-Guo-Li-Li-Liu. Read [A], [P]. New key GJR18.
* **L6. Chang-Guo-Li.** arXiv:1809.11058v1 ("Polynomial structure of Gromov-Witten potential of quintic 3-folds via NMSP"; published as
  "Polynomial structure of Gromov-Witten potential of quintic 3-folds", Ann. of Math. (2) 194 (2021), no. 3, 585-645 [bib: CGYZ
  [NMSP2], Brini [21]]): finite generation `F_g(Q(q))` in the period ring (as summarised in CGYZ s.1.1). arXiv:1810.00394v2 ("BCOV's
  Feynman rule of quintic 3-folds"; journal version per CGYZ [NMSP3]: "BCOV Feynman rule of the quintic threefold", Duke Math. J. 175
  (2026), no. 4, 607-675, to be confirmed): the BCOV Feynman rule via an R-matrix graph sum from NMSP; genus 1 and 2 potentials; the two
  Yamaguchi-Yau equations (abstract). N-Mixed-Spin-P fields: arXiv:1809.08806v2, Geom. Topol. 25 (2021) 775-811. Genus one via MSP:
  arXiv:1711.10118v1, IMRN 2020, no. 19, 6421-6462. None of these claims the gap. Read [abs] only. New keys CGL21, CGL26, CGLL21, CGLZ20.
* **L7. Liu-Ruan.** arXiv:2210.13411v2: the Castelnuovo bound for the quintic, and the Gopakumar-Vafa invariants at the bound; physicists
  then compute all GW invariants up to genus 53 "provided that the conifold gap condition holds" (abstract). Gap assumed. Read [abs].
  New key LR22.
* **L8. Chang-Guo-You-Zhang 2026: the gap for the quintic is proved.** arXiv:2609.13921v2 (v1 2026-09-12, v2 2026-10-31; preprint, no
  journal reference). NEW PRIOR ART, not in any project source (the project's LIT task searched only (QRR), formal stationary phase and
  the chart mechanism). Read [H] s.1 in full, parts of s.3-s.4.
  * Definitions 1.1-1.2: period ring `R` generated by `I_0^(+-1)`, `(D tau)^(+-1)` over `Q(q)`; an admissible pair `(omega_0, omega_1)`
    at the conifold and the differential homomorphism `iota` ("conifold expansion"); Example 1.3 for the quintic: `Delta = 1 - 5^5 q`,
    `omega_0 = 1 + O(Delta^3)`, `omega_1 = Delta + (7/10) Delta^2 + O(Delta^3)`, `tau_c = omega_1/omega_0`; Remark 1.4: every
    `omega_0' = omega_0 + c_1 omega_1 + c_2 omega_2` is admissible.
  * Conjecture 1.5 (quintic conifold gap): for any admissible pair and `g >= 2`,
    `F_g^c(tau_c) - (-1)^(g-1) B_2g / (5^(g-1) 2g(2g-2) tau_c^(2g-2))` lies in `C[[tau_c]]`.
  * Definition 1.8 (conifold-regular CohFT; the "conifold gap property" with a constant `kappa`); Theorem 1 (abstract gap theorem:
    a conifold-regular master CohFT `Omega^M = R T_R(Omega^(X,tau) (+) Omega^pt)`, `R` from a regular local-global QDE system
    (Definition 1.7), and a translated point theory with the gap property imply the gap for `X`).
  * Theorem 2: for the quintic (MSP `[0,1]`-theory) the conditions hold, so "the conifold gap property holds for the quintic
    threefold", every `g >= 2`, any admissible `omega_0`; also `F_g(Q(q)) in R`. Extensions (s.5): Fermat hypersurfaces of degree 6, 8,
    10, local `P^2`, FJRW theories.
  * Mechanism: Proposition 3.1, a factorisation `R = R^reg R^prin` with `R^reg` conifold regular; "The factorization below forces every
    non-leading `R^prin`-graph contribution to vanish" (s.3, opening); the gap of the cone-vertex theory (Theorem 5, via Theorem 4 and a
    master space for the resolved conifold) is transferred to the quintic.
  Status there: proved (unrefereed preprint). New key CGYZ26.
* **L9. The formal quintic in CGYZ.** CGYZ s.1: "the approach of Guo-Janda-Ruan gives a computational verification for `g <= 5`, via the
  formal quintic and the log GLSM localization formula under development [GJR18] (see also Lei's lecture notes [Lei25])". Remark 1.9:
  that formula "reduces the conjecture to two key statements: the conifold gap for the formal quintic, and the conifold regularity of its
  correlators with special insertions. While we expect our method to prove the conifold gap for the formal quintic as well, a
  corresponding formal theory may not exist for general Calabi-Yau threefolds." So, per CGYZ v2, the all-genus gap for the formal
  quintic is not proved there, nor (to their knowledge) elsewhere. Read [H]. The `g <= 5` verification is known to me only through
  CGYZ (Lei's notes not read).

## 1.3 Other targets

* **L10. Resolved conifold.** The gap holds trivially: it follows from the multiple-cover formula (Brini arXiv:2509.19298v3, s.1.1
  footnote 1; CGYZ s.1: "the resolved-conifold gap, which follows from the multiple-cover formula of Faber and Pandharipande [FP00]").
  Read [H] both.
* **L11. Local `P^2`.** Brini, arXiv:2509.19298v3 (v3 2026-10-13; comments: "Accepted version, to appear on DMJ"): Conjecture (1.5)
  `GW_g^CF = 3^(g-1) B_2g / (2g(2g-2)) t_CF^(2-2g) + O(1)` with `t/(e^t - 1) = Sum B_n t^n/n!`; Theorem (= Corollary 4.4): the Conifold
  Gap Conjecture holds for `K_P^2` in all genera; method: a conifold mean field model (ensemble of repulsive particles), with the
  Eynard-Orantin recursion for its correlators (v3 comments), then BCOV direct integration. s.1.1: "the Conifold Gap does not emerge as a
  general consequence of the R-matrix quantisation formalism for Gromov-Witten theory". Earlier check to genus 7: Coates-Iritani,
  arXiv:1804.03292v2, Kyoto J. Math. 61 (2021), no. 3, 543-706 (as stated in Brini s.1.2). Read [H], [abs]. New keys Br25, CI21.
* **L12. Topological recursion.** Fang-Chen-Chen-Zong, arXiv:2608.11960v1: a conifold gap theorem for the topological recursion of
  toric mirror curves under a generic one-node degeneration; one universal polar term of the conifold-polarised free energy, the rest
  jointly holomorphic (abstract). The note's model curve `x^4 + v/x` is not a toric mirror curve and degenerates by five branch points
  merging into the simple pole (NATURAL_FORM s.0), so it is outside the stated scope; the note does not use it. Read [abs]. Key FCCZ26.
* **L13. Not to be cited.** arXiv:2610.38473 (a physics computation to genus 83) was withdrawn (v2, 2026-10-30). Read [abs].

## 1.4 The formal quintic

* **L14. Definition.** Lho-Pandharipande, arXiv:1702.06096v4 (Adv. Math. 332 (2018) 349-402 [bib: GJR [27], Brini [51]; not on the
  arXiv page]) s.0.2, eq. (5): the formal quintic invariants `N~_(g,d) = Int_[Mbar_g(P^4,d)]^vir e(R pi_*(S^-5))`, the Euler class being
  defined after localisation, then specialised at `lambda_i = zeta^i` (or `zeta^i lambda_0`); s.0.3, eq. (6): the twisted theories on
  `P^m`; s.0.6 Theorems 3-4 (formal quintic: finite generation and holomorphic anomaly equations; proved in 1803.01409). Read [A].
  LIT s.1.1 confirmed. New key LP18.
* **L15. Holomorphic anomaly equations.** Lho-Pandharipande, arXiv:1803.01409v2 (2020; no journal reference on the abstract page; a
  search result suggests Peking Math. J., DOI 10.1007/s42543-018-0008-0, to be confirmed): s.0.2 eq. (5) (the same definition; footnote:
  `S^-5 = (S^*)^(x)5`); Theorem 1 (finite generation in `C[L^(+-1)][A_2, A_4, A_6, C_0^(+-1), C_1^(-1), K_2]`); Theorem 2 (the
  holomorphic anomaly equations, matching the conjectural ones of the true quintic [BCOV-type, their ref. (2.52)]); s.0.4 "Constants":
  the equations determine `F~_g^B` up to finitely many constants of integration, which "can be effectively computed via the
  localization formula, but whether there exists a closed formula determining the constants is an interesting open question". The
  conifold point is not discussed (no occurrence of "conifold" or "gap" in the text). Read [A], [P]. New key LPfq (s.3; no year in the
  key, since the journal data are unconfirmed). LIT s.1.1 confirmed.
* **L16. Localisation vertex.** 1803.01409v2 s.3.2, eq. (26): `H_g^(p_i) = e(E_g^* (x) T_(p_i)P^4)/e(T_(p_i)P^4) . (5 lambda_i)/
  e(E_g^* (x) (5 lambda_i))`. Read [A]. LIT s.1.1 confirmed; LIT s.1.4(a) derives `h_(w_i)` from it (project, PROVED given quotes).
* **L17. Zagier-Zinger.** arXiv:0710.0889v1 (in: Modular forms and string duality, Fields Inst. Commun. 54, AMS, 2008, 163-177 [bib:
  GJR [35]]): eq. (1) `F(w,x)`; Lemma 1; the expansion (5) `F ~ e^(mu w) Sum_s Phi_s w^(-s)` coefficientwise in `x`; Theorem 3, eq. (6)
  (`Phi_0 = L`, `Phi_1 = ((n-2)(n+1)/(24n))(L - L^n)`); Theorem 4 (i) the ODEs (12) with `Phi_s(0) = delta_(0,s)`, (ii) `Phi_s in L Q[L]`,
  (iii) `Phi_s in Q(n)[L, L^-1, L^n]`; s.3, Conjecture (39). GJR s.1 credits ZZ with the specialisation `(lambda, zeta lambda, ...,
  zeta^4 lambda, 0)` that defines the formal quintic ("referred as formal quintic (see [27])"). Read [A], [P]. LIT_2 s.3.5 confirmed.
  The project cites Thm 4(ii) for `deg ph_n <= floor(4n/5)` (MECHANISM_STAGE4 s.3). New key ZZ08.
* **L18. GJR 2018 on the formal quintic.** arXiv:1812.11908v1 s.3.1: the "lambda-twisted theory" (`t = 0`, `lambda_i = zeta^i lambda`)
  is the formal quintic, eq. (9); s.4.1: the formal quintic CohFT, reconstructed by Givental-Teleman (their Theorem 4.1); Lemma 5.2 with
  eq. (17) (L27 below); Lemma 5.5 (`R_0 = Id`, degree subspaces), proved with Cor. 7.4; s.7 "A technical result of the formal quintic
  theory": Theorem 7.1 (regularity at `L_alpha = 0`, degree bound, polynomiality in `m`), Remark 7.2 (stated without proof; L39),
  Cor. 7.3 (ZZ's conjecture (39), stronger form), Cor. 7.4, Prop. 7.5 (formal stationary phase), Lemma 7.7 with (37) (L40), Lemmas
  7.8-7.9. No statement about the conifold point of the formal quintic. Read [A], [P].
* **L19. GJR 2017.** arXiv:1709.07392v1 (no journal reference): genus-two mirror theorem for the quintic (abstract); proof of Lemma 4.6:
  the auxiliary theory with `(C^*)^5` on the base, `lambda_i = xi^i lambda`, semisimple at `q = 0`; eq. (19) (normalisation of `S` by
  constant matrices, "by uniqueness of fundamental solutions of the quantum differential equation (see also [21, 13])"), eq. (20):
  `C(lambda, z) = diag(exp(Sum_(k>0, j != i) B_2k/(2k(2k-1)) z^(2k-1)/(lambda_i - lambda_j)^(2k-1)))`,
  `Gamma(t, z) = exp(Sum_(k>0) B_2k/(2k(2k-1)) z^(2k-1)/(t - 5H)^(2k-1))`, "these constant matrices together give exactly the constant
  term of the R-matrix". At `t = 0` this is (QRR). Read [A] (exact LaTeX). LIT s.1.2 confirmed. New key GJR17.
