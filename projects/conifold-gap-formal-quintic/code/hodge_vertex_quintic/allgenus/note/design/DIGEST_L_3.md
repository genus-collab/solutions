# DIGEST_L_3: s.1, what is known (part 2), and what is new

Continuation of `DIGEST_L_2.md`. Same format.

## 1.5 Givental-Teleman reconstruction and the q = 0 normalisation (QRR)

* **L20. Givental 2001a.** arXiv:math/0008067v4 (IMRN 2001, no. 23, 1265-1286 [bib: GJR [16]; the conifold note's Giv01]): s.2.2
  "Compensating constants", Theorem 1: for torus-equivariant GW theory with isolated fixed points the higher genus potential is given by
  formula (1) with the fundamental solution `S` replaced by `newS = S exp[-Sum z^(2k-1) N_(2k-1)(1/chi^i)/(2k-1) . B_2k/(2k)]`
  (`N` = Newton sums of the inverse tangent characters `chi^i` at the `i`-th fixed point); Remark after it: `newR(z)` tends to the
  diagonal matrix of these constants "in the limit of classical equivariant cohomology". Proof by localisation (s.2.3 Lemma, Mumford's
  GRR). Read [A]. LIT s.1.3 confirmed. Key Giv01.
* **L21. Givental 2001b.** arXiv:math/0108100v2 (Mosc. Math. J. 1 (2001), no. 4, 551-568, 645; root key MR1901075): 4.1 Theorem (the
  Hodge potential `E_X = exp[Sum B_2k/(2k)! s_k (z^(2k-1))^] D_X`), 4.2 Remark (a reformulation of Faber-Pandharipande and Mumford);
  9.1 Theorem (ancestor potential = semisimple formula provided `R|_(Q=0) = exp diag(b_1..b_N)`, `b_i = Sum N^(i)_(2k-1) B_2k/(2k)
  z^(2k-1)/(2k-1)`, `N^(i)` over COTANGENT characters); 9.2 Remarks (a): equivalent to Theorem 2 of [22] (= math/0008067), proved by
  fixed-point localisation; 9.3 Conjecture; s.10: L38 below. Untwisted. Read [A], [P]. LIT s.1.3 and STAGE5_SYNTHESIS_2 K9 confirmed.
* **L22. Coates-Givental.** arXiv:math/0110142v2 (Ann. of Math. (2) 165 (2007), no. 1, 15-53; root key MR2276766): Theorem 1, eq. (7)
  (quantum Riemann-Roch for `c(.) = exp(Sum s_k ch_k(.))`); Remarks (1)-(2) (for `E = C` it is Givental's 4.1 Theorem, i.e. Mumford and
  Faber-Pandharipande; proof in Appendix 1); s.4, eq. (9) and Corollary 1 (Euler-class twist acts by
  `prod_i exp{Sum_(m>0) B_2m/(2m(2m-1)) ((z/(lambda + rho_i))^(2m-1))^}` times `1/z` terms); s.7, eq. (11) (L37); Theorem 2 (the
  hypergeometric modification `I_E` lies on the twisted cone; genus 0). Read [A], [P]. LIT s.1.3 confirmed.
* **L23. Faber-Pandharipande.** arXiv:math/9810173v1 (Invent. Math. 139 (2000), no. 1, 173-199; root key MR1728879): s.0 Theorem 1
  (reconstruction of Hodge integrals), Theorem 4 (`Int_(Mbar_g) lambda_(g-1)^3 = |B_2g|/(2g) |B_(2g-2)|/(2g-2) 1/(2g-2)!`); s.1.1
  "Theorem (Mumford)" (unnumbered; `ch(E) = g + Sum_l B_2l/(2l)! (kappa_(2l-1) + boundary term)`; Bernoulli convention
  `t/(e^t - 1)`); Proposition 1 (virtual version); Proposition 2 (`D_(2l-1) Z^X_E = 0`); s.2.1 Mumford's formula `c(E) c(E^*) = 1` and
  its consequence (22). Read [A], [P]. LIT s.1.3 confirmed.
* **L24. Teleman.** arXiv:0712.0160v3 (2012-02-17; Invent. Math. 188 (2012), no. 3, 525-588; root key MR2917177). Numbers as in the v3
  PDF (they differ from ar5iv, which renumbers: s.7 item 1):
  * Theorem (2), in paragraph (1.15) "Givental's loop group conjecture" [ar5iv: 1.9]: the CohFTs with underlying semisimple Frobenius
    algebra `A` are the `Sp_L^+ x| H^++` orbit of the trivial theory; theories with flat vacuum form the orbit of `T_z Sp_L^+ T_z^-1`;
  * 8.4 Proposition ([D, G1]) [ar5iv: 8.3]: solutions `R == Id mod z` of the canonical-frame ODE with the symplectic condition are
    unique up to right multiplication by `exp(H_1 z + H_3 z^3 + ...)`, `H_(2i+1)` constant diagonal; in the homogeneous case a unique
    Euler-invariant solution; 8.5 Corollary [ar5iv: 8.4];
  * (8.6) "Ambiguity for inhomogeneous theories" [ar5iv: s.8.2]: the Hodge CohFT `Z[h]`; 8.7 Proposition [ar5iv: 8.5]: (i)
    `Z[h] = T_z^-1 exp h(z) T_z (I_A)`; (ii) all CohFTs with flat vacuum on a fixed semisimple pointed Frobenius manifold are
    classified by `E exp h(z)` with arbitrary `h` and the same `E`.
  Read [A] (s.1.9, s.8.1-s.8.2) and [P] (v2 and v3; v2 numbers these (8.7) and 8.8 Proposition, so the version must be named). The project's "Theorem (2) (s.1.9), Prop. 8.5(ii) (s.8.2)" (MECHANISM_STAGE5 s.2;
  REFEREES_SYNTHESIS_2) is ar5iv numbering; the conifold note's LABELS ("(1.15); Prop 8.7(ii)") is v3 numbering. LIT_3 s.5 ("Teleman
  not opened") is superseded: now read verbatim.
* **L25. PPZ.** Pandharipande-Pixton-Zvonkine, arXiv:1303.1043v3 (2015-03-18; J. Amer. Math. Soc. 28 (2015), no. 1, 279-309 [bib:
  GJR [31], CGYZ [PPZ15]]); the id is correct. s.1 Theorem 5 (Teleman, citing [20, Theorem 1]); s.2: Definition 2.2 (R-matrix action:
  vertices `Omega`, legs `R^-1(psi)`, edges `(eta^-1 - R^-1(psi') eta^-1 R^-1(psi'')^t)/(psi' + psi'')`; "the same convention is used
  by Givental and Teleman"); Proposition 2.3 (`R Omega` is a CohFT); Proposition 2.4 (left group action); Definition 2.5 (translation
  `T Omega`); Propositions 2.7-2.8; Proposition 2.9 (`T_a R Omega = R T_b Omega` if `T_a = R T_b`, `R` symplectic); Definition 2.13
  (unit-preserving action `R.Omega = R T Omega`, `T(z) = z 1 - z R^-1(z) 1`); Proposition 2.14. Read [A], [P]. Key PPZ15.
* **L26. Localisation.** Graber-Pandharipande, arXiv:alg-geom/9708001v2 (Invent. Math. 135 (1999), no. 2, 487-518; root key MR1666787):
  virtual localisation, used by LP (eq. (26)) and in the project's calibration. Read [abs] only.
* **L27. (QRR) for the formal quintic is in the literature as a cited input.** GJR arXiv:1812.11908v1, proof of Lemma 5.2: the
  Birkhoff factorisation `S_(i alpha) C_alpha^-1 = e^(u_alpha/z) R_(i alpha)` "by the results of [16, 12]" (= Givental math/0008067,
  Coates-Givental), with, exactly (authors' LaTeX, [A]),
  `C_alpha(z) = e^(Sum_(k>0, beta != alpha) B_2k/(2k(2k-1)) z^(2k-1)/(zeta^alpha lambda - zeta^beta lambda)^(2k-1))
              . e^(Sum_(k>0) B_2k/(2k(2k-1)) z^(2k-1)/(-5 h_alpha)^(2k-1))`   (17),
  `h_alpha = zeta^alpha lambda`. So the twist weight `-5 h_alpha` is confirmed in the authors' source (LIT s.1.2 had reconstructed the
  sign from structure; STAGE5_SYNTHESIS_2 O15 asked for a reading of the typeset sign: the LaTeX source settles it). `R(q = 0) =
  C_alpha^-1 = h_(w_alpha)` is LIT's reading (LIT s.1.2). Neither GJR paper proves (17); GJR 2017 (19)-(20) is the second statement
  (L19). The step "the R-matrix of the full CohFT restricts at `q = 0` to the vertex R-matrices" uses the uniqueness of L24 (8.4, 8.7).

## 1.6 Topological recursion and the model curve (used ONLY on `x^4 + v/x`)

* **L28. Eynard-Orantin.** arXiv:math-ph/0702045v4 (Commun. Number Theory Phys. 1 (2007), no. 2, 347-452; root key Eynard:2007kz):
  s.4.1 Definition 4.1 (vertex (4-1), propagator (4-2), arrow propagator `dE_q(p) = (1/2) Int_q^(q-bar) B(xi, p)` (4-3), with the
  footnote that this is opposite to earlier notation; root (4-4)); Definition 4.2 (recursion (4-9)); Definition 4.3 (`F^(g)`, (4-10));
  s.4.4: Theorem 4.1 (`W_3^(0)`), 4.2 (poles only at branch points), 4.3, **4.4 (the linear loop equation (4-23)-(4-24))**, 4.5
  (quadratic loop equation), 4.6 (symmetry), Corollary 4.1 ((4-28)-(4-31)), 4.7 (dilaton, (4-32)); s.5.1 **Theorem 5.1** (variations:
  `D_Omega W_k^(g) = Int_(dOmega) W_(k+1)^(g) Lambda` (5-37), `D_Omega F^(g) = -Int W_1^(g) Lambda` (5-38), for
  `Omega = Int_(dOmega) B Lambda` with the path avoiding the branch points; hard edges excluded, footnote 8). Read [A], [P]. The
  project's "linear loop equation of its s.4" (Q) is Theorem 4.4 (STAGE2_SYNTHESIS M6: "statement number not re-checked"; now pinned).
* **L29. DNOPS.** Dunin-Barkowski-Norbury-Orantin-Popolitov-Shadrin, arXiv:1509.06954v2 (J. Inst. Math. Jussieu 18 (2019), no. 3,
  449-497, DOI 10.1017/S147474801700007X [abs]; root key Dunin-Barkowski:2015caa, whose root data are outdated: s.3): s.8, superpotential
  (8-1) `f(p,a) = p^n + a_2 p^(n-2) + ... + a_n + a_(n+1)/(p - a_1)`, `n >= 1`, semisimple Frobenius manifold ([6, Section 5] there);
  **Theorem 8.1**: the CohFT of this Frobenius manifold coincides with the one reconstructed from the spectral curve
  `(CP^1, f(p,a), p, dp_1 dp_2/(p_1 - p_2)^2)`; proof parallel to Theorem 7.1, via Theorem 5.1; `d = (n-2)/n`; Remark 8.2 (hypermaps).
  The model is `n = 4` on the line `a = (0,0,0,0,v)` (STAGE2_SYNTHESIS_2 s.6.2 item 6), `d = 1/2`. Read [A], [P].
* **L30. ABDKS.** Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin, arXiv:2212.00320v3 (J. Eur. Math. Soc. (2025), DOI
  10.4171/JEMS/1615 [abs]): **Theorem 7.1**: for `x(z)` an arbitrary rational function and `y = z`, the closed formula (7.9) for all
  `omega^(g)_(m,0)` (a corollary of Theorem 1.14, the universal x-y swap formula); Remark 7.3 (limits at discriminants). Used by the
  project for L-regularity and for the second route to the model's gap (STAGE2_SYNTHESIS M7, M8). Read [A], [P].
* **L31. Hock.** arXiv:2502.18115v2 (Commun. Math. Phys. 406, 192 (2025), DOI 10.1007/s00220-025-05373-8 [abs]): Theorem 1.1
  (`F_g - F_g^vee` for compact spectral curves, `g >= 2`, logarithmic TR); **Corollary 1.2**: explicit `F_g`, `g >= 2`, when `y` is
  unramified (`y = z`: formula (1.3); `y = log z`: (1.4)). Optional route to the model's gap (M8). The project's remark that the printed
  example of its s.3.2 has a sign slip (TR_MODEL s.4.1, via STAGE2_SYNTHESIS M8) was not re-checked by me. Read [H], [P].
* **L32. Carlet-van de Leur-Posthuma-Shadrin.** arXiv:2211.12259v1 (Commun. Number Theory Phys. 17 (2023), no. 3, 643-708 [abs]):
  the Hurwitz space `H_(0;1,N-1)` (two poles, of orders 1 and `N-1`; the model is `N = 5`); Theorem 2 (restriction of the descendant
  potential to the times `t^1_k` enumerates rooted hypermaps; proof via TR on `x = z^(N-1) + 1/z`, `y = -z`); Theorems 14 and 17
  (ancestor and descendant Hirota equations); Propositions 22, 24, 25 (Lax form: the rationally `(N-1)`-constrained KP hierarchy and
  its extension); some properties of the extra flows are conjectural (their Remark 28, cited in their introduction). Project use: MECHANISM
  s.1 (ML on the `H^2` tower, not in the chain), NATURAL_FORM s.2 (not used). Read [A], [P].
* **L33-L34. Background.** DOSS arXiv:1211.4021 (CMP 328 (2014) 669-700; root DuninBarkowski:2012bw); DNOPS arXiv:1605.07644 (Proc.
  Sympos. Pure Math. 100 (2018) 297-331; root dunin2016primary); Eynard arXiv:1110.2949 (root Eynard:2011ga). [abs] only; optional.

## 1.7 Mellin-Barnes integrals, stationary phase, Gamma factors

* **L35. Iritani.** arXiv:0903.1463v3 (Adv. Math. 222 (2009), no. 3, 1016-1079; root key MR2553377): Theorem 1.1 (= 4.11, the
  Gamma-integral structure matches the integral local system; weak Fano toric orbifolds; under Conjecture 4.3 and Assumption 2.7(c));
  Theorem 1.2 (= 4.14, central charge of `O_X` = oscillatory integral over the real Lefschetz thimble); Remark 4.15; Theorem 4.17,
  Remark 4.18; **Lemma 4.20** (exact `q`-expansion, to any order, of the thimble integral at `z = -1` (Iritani's `e^(-W/z)`) near a
  fixed point, in a chart, by termwise Gamma integrals). Read [A]. LIT_2 s.3.3 confirmed, except "at `z = 1`" should read `z = -1`.
* **L36. Galkin-Golyshev-Iritani.** arXiv:1404.6407v4 (Duke Math. J. 165 (2016), no. 11, 2005-2077): Theorem 1.3.1 (= 5.0.1, Gamma
  conjectures I and II for `P^(N-1)`); (5.2.3) `Psi(t) = (2 pi i)^-1 Int Gamma(s)^N t^(-Ns) ds`; Proposition 5.2.4:
  `Psi(t) ~ C t^(-(N-1)/2) e^(-Nt)(1 + O(1/t))`, `C = N^(-1/2)(2 pi)^((N-1)/2)`, and the connection formula `Psi = Int Gamma^ . Pi`;
  Corollary 5.2.6. Read [A]. LIT_2 s.3.4 corrected: the function is `Psi` (not `phi`) and the exponent is `-(N-1)/2`.
* **L37. Coates-Givental s.7.** eq. (11) `b_rho(z)`; "(11) coincides with the stationary phase asymptotics of the integral
  `(2 pi z (lambda + rho))^(-1/2) Int_0^oo e^((-x + (lambda+rho) ln x)/z) dx`" (asserted, with a reference to a table of functions);
  Theorem 2; s.8 eq. (13). Read [A]. LIT_2 s.3.2 confirmed.
* **L38. Givental s.10.** arXiv:math/0108100v2 s.10: Theorem 10.1 (mirror oscillating integral of `CP^n`); at `Q = 0` the R-matrix from
  the chart `x_0 = Q/(x_1...x_n)`; the integral factorises into `prod Int e^(x/z) x^((lambda_i - lambda_0)/z) dx` ~ `Gamma(1 + .)`;
  10.4 Lemma (Stirling, "well-known"); "For generic `lambda`, the solution depends continuously on `Q` up to `Q = 0`" (asserted);
  extension to toric Fano manifolds asserted; 10.5 Theorem, 10.6 Corollary. Untwisted; outline. Read [A]. LIT_2 s.3.1 confirmed.
* **L39. GJR Remark 7.2 (stated without proof).** Exact display (authors' LaTeX): with `W = Sum_i (x_i + lambda_i ln x_i) -
  (q^-1 prod_i x_i)^(1/m)` and the thimble integrals `I_alpha` of `e^(W/z) dx_0..dx_(m-1)/(x_0...x_(m-1))`,
  `(-2 pi z)^(-m/2) e^(-h_alpha log q/z) C_alpha(z) I_alpha(q,z) = I~_alpha(q,z) := z Sum_d q^d prod_(k=0)^(md-1)(m h_alpha + kz) /
  prod_(k=1)^d ((h_alpha + kz)^m - lambda^m)`. Critical points `(x_i)_alpha = L_alpha - lambda_i`, `y_alpha = m L_alpha`. GJR use it for
  Cor. 7.3-7.4 and hence Lemma 5.5. The three items LIT_3 s.5 could not read (the display of `W`, the sign of the exponent of `2 pi z`,
  the factor `m` in (37)) are now read exactly. Read [A].
* **L40. GJR Lemma 7.7, eq. (37) (proved).** At `L_alpha = 0`: `(1 + Sum r_(alpha,k;l) z^l)|_(L_alpha = 0) =
  exp[m Sum_(j>=1) (-1)^(mj+1) B_(mj+1)(k/m)/(mj(mj+1)) z^(mj)/lambda^(mj)]`, from the general display with `lambda_i = zeta^i lambda`.
  For `m = 5`, `k = 1`, `lambda = 1`, since `(-1)^(5j+1) = (-1)^(j+1)`, this is the project's `Ginf(z) = exp(Sum_l (-1)^(l+1)
  B_(5l+1)(1/5) z^(5l)/(l(5l+1)))` (LIT_2 s.3.6 confirmed, factor `m` present). Proof: the integrand at `L_alpha = 0` is a product of
  one-variable integrals, compared with `ln Gamma(z + a)` expanded with Bernoulli polynomials [Hermite 1895; Nemes 2013, (1.8)] =
  DLMF 5.11.8; the continuation to `L_alpha = 0` rests on Prop. 7.5 (rational coefficients regular there). Read [A].
* **L41. DLMF** (Version 1.2.8, 2026-09-15): 5.11.1 (Stirling's series, `|ph z| <= pi - delta`); 5.11.8 (`Ln Gamma(z + h)` with
  `(-1)^k B_k(h)/(k(k-1) z^(k-1))`; `h` anywhere in C since Version 1.0.11; source Rowe 1931); 24.4.1 (`B_n(x+1) - B_n(x) =
  n x^(n-1)`, used in the formal proof of 5.11.8, REFEREES_SYNTHESIS A7); 2.3.8 (Watson's lemma); 2.3.13-2.3.18 (Laplace's method);
  s.2.3(iv) (stationary phase). Read in full.
* **L42. Formal stationary phase.** Statement numbers of Etingof arXiv:2409.03117v2 (Thms 2.16-2.17, 3.1, 3.5), Kirwin
  arXiv:0810.1700v2 (Thm 1.1), Polyak arXiv:math/0406251 (Prop. 2.17 etc.), Bar-Natan-Garoufalidis-Rozansky-Thurston
  arXiv:math/9801049v4 (Defs 2.7-2.8, Prop. 2.13) are LIT_2 s.2's (summariser reads), not re-checked; MECHANISM_STAGE5 s.3.2: none states
  the formal lemmas for complex data with numbers, and the project's proofs are self-contained. Context citations only.

## 1.8 What is new in the note, against L1-L42

* **L43.** Conjecture N(g) for all `g >= 2` (NATURAL_FORM s.1 item 7; MECHANISM_STAGE5 s.2): no statement of this graph-local kind
  for the formal quintic in any paper read. New (conditional on (Q), (QRR); (U) is proved in the note).
* **L44.** The all-genus conifold gap of the formal quintic (the leg-free graph sum, i.e. the sum of the graphs of L43): open per CGYZ
  Remark 1.9 (v2); one of the two "key statements" of the GJR approach there. New as a theorem, conditional; its derivation from N(g)
  must be written in the note (s.7 item 4).
* **L45.** The factorisation `R = R' M_0(tau(u) zhat)` through the Landau-Ginzburg model `x^4 + v/x` with primitive form `dx`
  (NATURAL_FORM s.0): not found in the literature read. A factorisation of different nature (master space, `R = R^reg R^prin`) is CGYZ
  Prop. 3.1 (L8); cite it as related.
* **L46.** Theorem Gamma (no pole at `L = 1`, value `h_(w_0)`): no proof found; GJR Remark 7.2 states an equivalent identity for their
  LG integral without proof; Givental s.10 is the untwisted outline (L38). New as a written proof; wording modest (STAGE5_SYNTHESIS_2 s.8).
* **L47.** `Ginf`: GJR Lemma 7.7 (37) proves the formula for GJR's LG series; the note's two-variable chart at `L = 0` identifies it
  with the note's series. Not new as a formula; cite GJR (STAGE5_SYNTHESIS_2 K7).
* **L48.** (QRR) is quoted, not proved (L27, L24); (CYC), Theorem U, Theorem P(b), the glue and the formal lemmas are the project's own
  (no statement found in the papers read; no wider search was made for them).
