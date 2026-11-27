# Stage 5, task `lit`, part 2: formal stationary phase (s.2) and prior art for the chart mechanism (s.3)

Continuation of `LIT.md` (s.0 bottom line, s.1 QRR). Part 3 (`LIT_3.md`): consequences for Stage 5, what was not read, record.
Same labels and notation. (L1)-(L3) are the lemmas named in `../../MECHANISM_STAGE4.md` s.5 item 1 and
`../../stage4/STAGE4_SYNTHESIS.md` s.1.4.

## 2. Formal stationary phase: what can be cited, with statement numbers

"Read" means opened in full text. Where only an HTML rendering was available (marked [H]), it was read through the fetch tool's
summariser, and statement numbers were confirmed by a second, targeted query. "PDF" means read in my own text extraction of the
arXiv PDF.

### 2.1 Laplace / stationary phase for actual integrals (the analytic base)

* **DLMF s.2.1(iii), eqs. (2.1.13)-(2.1.16)** (online, read). Poincare asymptotic expansions; such expansions are unique (the
  sentence after (2.1.16)). Sources there: Olver 1997b pp. 16-22.
* **DLMF s.2.3(iii), eqs. (2.3.13)-(2.3.18)** (read). This is Laplace's method for `Int_a^b e^(-x p(t)) q(t) dt`, with the minimum
  of p at the endpoint a.
  * Hypotheses: p', q continuous near a; the minimum is approached only at a; the expansions (2.3.14) hold with `mu, lambda > 0`,
    `p_0 > 0`; the integral converges absolutely for large x.
  * Conclusion (2.3.15), with coefficients from (2.3.16), explicitly (2.3.17), and the residue formula (2.3.18)
    `b_s = (1/mu) res_(t=a)[q(t)/(p(t)-p(a))^((lambda+s)/mu)]`.
  * An interior minimum is two endpoint integrals with `mu = 2`. This splitting is standard; I did not see it stated in
    s.2.3(iii) itself.
  * s.2.3(iv), eqs. (2.3.19)-(2.3.23): stationary phase, one variable. Sources: Olver 1997b pp. 80-88; Olver 1974.
* **Etingof arXiv:2409.03117v2** ([H]; arXiv's own HTML of v2, 2026-08-28).
  * Theorem 2.16 (several variables, steepest descent). Hypotheses: D compact with smooth boundary; f, g continuous on D and smooth
    inside; f has a unique global minimum at an interior point c with `f''(c) > 0`. Conclusion:
    `Int_D g e^(-f/hbar) dx = hbar^(d/2) e^(-f(c)/hbar) I(hbar)`, where I extends to a smooth function on `[0, oo)` and
    `I(0) = (2 pi)^(d/2) g(c)/sqrt(det f''(c))`.
  * Theorem 2.17: stationary phase, same type, for a non-degenerate critical point, with g vanishing to all orders on the boundary.
  * Theorem 2.18 and Corollary 2.19: the Morse lemma.
  * s.2.9: the proofs are said to parallel the one-variable case. Remark 2.22: several critical points add.
* **Kirwin arXiv:0810.1700** ([H]).
  * Theorem 1.1 (d variables). Hypotheses: `f in C^(N+2)`, `g in C^N`, a unique non-degenerate minimum at 0, convergence for some
    `k_0`. Conclusion: the full expansion of `Int e^(-k f) g`, whose coefficients are explicit finite sums of Taylor derivatives at 0
    in normalised coordinates, times `(2 pi)^(d/2)/sqrt(det Hf(0))`.
  * Corollary 3.6 is a Stirling-type series from the case d = 1. Theorem 2.1 is the degenerate case.
  * Only real f, and no parameters.

### 2.2 Formal Gaussian integrals, Wick, Feynman

* **Etingof v2** ([H]).
  * Theorem 3.1: Wick.
  * Theorem 3.5: Feynman's theorem. For smooth S with a non-degenerate minimum at 0, the normalised integral of
    `l_1...l_N e^(-S/hbar)` equals `(2 pi)^(d/2)/sqrt(det B)` times a sum over graphs, with vertices `B_i = d^i S(0)` and
    propagator `B^-1`. This is stated as an identity of asymptotic expansions in hbar.
  * Remark 3.6: the sum converges hbar-adically. Theorem 3.8: the weighted version. Theorem 3.11: connected graphs.
  * So the coefficients of the asymptotic expansion are, explicitly, universal polynomials in finitely many Taylor coefficients and
    in `B^-1`, times `det^(-1/2)`.
* **Polyak arXiv:math/0406251** ([H]).
  * Prop. 2.2: Gaussian integral for positive-definite A. Remark 2.3: complex A with positive-definite imaginary part.
  * Prop. 2.4: a linear source term. Thms 2.7-2.8: Wick. Prop. 2.17 (eq. (19)): Feynman expansion. Prop. 2.20: connected graphs.
  * The paper says explicitly that it sets aside questions of measures and convergence.
* **Bar-Natan-Garoufalidis-Rozansky-Thurston, Aarhus II, arXiv:math/9801049** ([H]). This is formal Gaussian integration of
  diagram-valued integrands `P exp(Q/2)`.
  * Def. 2.7-2.8: integrable iff the strut matrix is invertible; then `Int^FG = <exp(-Q^-1/2), P>`.
  * Ex. 2.11: compatibility with LINEAR changes of variables (an exercise).
  * Prop. 2.13: Fubini. If the X-block is invertible, the integral over `Z = X U Y` is the iterated integral.
  * Prop. 2.15: integration by parts.
  * Weight systems specialise this to ordinary formal Gaussian integrals with polynomial perturbations.

### 2.3 The lemmas needed by the chart argument, against these sources

| lemma | what the chart argument uses | read statements that cover it | gap |
|---|---|---|---|
| (L1a) changes of variables | formal SP invariant under an analytic (or formal) change of coordinates at a non-degenerate critical point, complex data | real phases: Etingof Thm 2.16 + uniqueness DLMF s.2.1(iii) (actual integrals are invariant, hence so are expansions); linear: Aarhus II Ex. 2.11 | no numbered FORMAL nonlinear statement for complex data; route below |
| (L1b) Fubini | iterated SP = total SP at a non-degenerate point (Schur complement) | Aarhus II Prop. 2.13 (formal Gaussian integrals, perturbation expanded at a fixed point) | the chart argument integrates the inner variables at a MOVING critical point; that needs (L1a) translations as well |
| (L1c) parameters | coefficients analytic in parameters; value at a parameter point is the SP of the specialised integrand | coefficients are universal polynomials in Taylor data and `B^-1` (Etingof Thm 3.5; Kirwin Thm 1.1) plus the implicit function theorem | elementary, but no numbered statement read; Taylor-in-q of the PHASE (critic step (e)) is a separate formal identity in `Q[[q]]((z))`, not covered |
| (L2) delta lemma | integrating `(sigma, x_0)` with phase linear in sigma restricts to `G = 0` exactly, all orders | none numbered. Polyak s.3.2 (unnumbered): `delta(F) = (2 pi)^(-l) Int e^(i <xi, F>) d xi`, xi called Lagrange multipliers, physics-level. Givental math/0108100 s.10: constrained critical points with a Lagrange multiplier, and he integrates directly over `{x_0...x_n = Q}` | write out the Wick proof (critic s.1.4(b)) using Aarhus II Def. 2.8 (invertible, possibly indefinite, quadratic part) |
| (L3) Stirling as SP of Euler's integral, also with `x^m` | `Int_0^oo e^((-x + a log x)/z) x^m dx/x` has normalised formal SP = Stirling factor of `Gamma(a/z + m)` | DLMF (5.2.1), (2.3.13)-(2.3.18), s.2.1(iii), (5.11.1), (5.11.3), (5.11.8) (Bernoulli polynomials, shift h in C); Etingof s.2.6; Kirwin Cor. 3.6 | none for real positive `a/z`; complex a by polynomial identity in `1/a`, coefficient by coefficient |

**Route for (L1a)-(L1c) (SKETCH, mine; not in any source I read).**
* Each coefficient of a normalised formal SP series is a fixed polynomial in finitely many Taylor coefficients of phase and
  amplitude at the critical point, in the entries of `(-Hess)^-1`, and in `det(-Hess)^(-1/2)` (Etingof Thm 3.5).
* For real-analytic phases with a non-degenerate MINIMUM, cut-off amplitudes and real data, these are the coefficients of the
  asymptotic expansion of an actual integral (Etingof Thm 2.16, Kirwin Thm 1.1). Actual integrals satisfy change of variables and
  Fubini exactly, and expansions are unique (DLMF s.2.1(iii)). So the polynomial identities hold on the real-definite locus.
* That locus is Zariski dense in the space of jets, so the identities hold for all complex jets.
* Caveats:
  * Fubini for actual integrals with a moving inner critical point needs Laplace's method uniformly in the outer variables. That is
    standard (Hormander s.7.7, NOT read) but not covered by the statements above.
  * `det^(-1/2)` has a sign or branch ambiguity for indefinite or complex Hessians (e.g. `det = -1/x_0^2` in the `(sigma, x_0)`
    block). It cancels in series normalised by `c_0 = 1` provided the same branch is used on both sides.

**(L3) status.** For `a/z > 0` real, Laplace's method at the interior minimum of `-x + a log x` (DLMF (2.3.15)-(2.3.18) after
splitting) gives an asymptotic expansion. It must equal Stirling's (5.11.1)/(5.11.3), shifted by (5.11.8) for `x^m`, by uniqueness.
The formal SP coefficients are those Laplace coefficients (Etingof Thm 2.16/3.5). Each is a polynomial in `1/a`, so the identity
extends to complex a. PROVED modulo writing it out; every input is read.
* Prior assertions of the same: Coates-Givental math/0110142 s.7 after eq. (11), the series `b_rho` is asserted to equal the
  stationary-phase asymptotics of `(2 pi z (lambda + rho))^(-1/2) Int_0^oo e^((-x + (lambda+rho) ln x)/z) dx`, asserted with a citation to
  Jahnke-Emde-Losch. Givental math/0108100 Lemma 10.4, Stirling for `ln Gamma(1+s)`, called well known, same tables. GJR 1812.11908
  Lemma 7.7 proof, via Hermite and Nemes (= DLMF (5.11.8)).

## 3. Prior art for the chart mechanism

### 3.1 Givental arXiv:math/0108100v2, s.9-10 (PDF; the passage read line by line)

* **What is stated, and with what proof.**
  * Thm 9.1 (stated; s.1.3): its proof is by citation, via Remark 9.2(a), to math/0008067 and to localisation.
  * Conjecture 9.3: the non-equivariant limit of the Bernoulli-normalised R.
  * Thm 10.1 (proved, two lines): the oscillating integral over the cycle in `{x_0...x_n = Q}`, with the form
    `dx_0 ^ ... ^ dx_n / d(x_0...x_n)`, satisfies the equivariant quantum differential equation of CP^n.
  * Remark 10.2 and Example 10.3: the asymptotic solution `R(z) exp(U/z)` comes from stationary phase at the critical points,
    illustrated in one variable.
* **The chart mechanism (outline).**
  * The critical points are constrained extrema with a Lagrange multiplier P. The R-matrix at Q = 0 is read off the integrals
    `L_j(z d/dt) I`, which carry an extra factor proportional to `Q/x_j`.
  * The integral is rewritten in the chart `(x_1..x_n)`, `x_0 = Q/(x_1...x_n)`. That chart still contains the critical point as
    `Q -> 0`.
  * The phase becomes `Q/(x_1...x_n) + Sum x_i + Sum (lambda_i - lambda_0) ln x_i`. At Q = 0 the integral vanishes unless j = 0. For
    j = 0 it factorises into `prod_i Int e^(x/z) x^((lambda_i - lambda_0)/z) dx` (measure dx, not dx/x).
  * Each factor's SP equals that of `Gamma(1 + (lambda_i - lambda_0)/z)`. Lemma 10.4 then gives
    `ln R^00|_(Q=0) = Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_(alpha != 0) (lambda_alpha - lambda_0)^(1-2k)`, which is `b_0` of Thm 9.1.
* **What is asserted without proof.** Continuity: for generic lambda the asymptotic solution depends continuously on Q up to Q = 0;
  it also has the non-equivariant limit for Q != 0. The extension to toric Fano manifolds, asserted to carry over without essential change.
  Theorem 10.5 and Corollary 10.6 rest on this; the abstract presents the Virasoro result as an outline with details to follow.
* **Hypotheses:** CP^n with `T^(n+1)` and generic lambda; toric Fano asserted. **Twisted / formal quintic: NO.**
  * Same mechanism as the critic's (b)-(d): the chart `x_0 = -q y^5/prod x_j`, factorisation into Euler integrals at q = 0, Stirling.
  * Differences: no O(5) factor, i.e. no y-variable and no `Gamma(1 + 5/z)`; the constraint is built into the cycle instead of being
    imposed by the sigma-integral; the continuity in Q is asserted for asymptotic solutions instead of being reduced to
    coefficients polynomial in L.

### 3.2 Coates-Givental arXiv:math/0110142v2, s.7-8 (PDF)

* **s.7.** The twisted genus-0 cone for the Euler class is the untwisted cone multiplied by `prod_rho b_rho(z)`, eq. (11). The
  identification of `b_rho` with the SP of the Euler integral is asserted (s.2.3 above).
* Theorem 2: the hypergeometric modification `I_E`, eq. (12), lies on the twisted cone. Hypotheses: E is a sum of line bundles
  with integer `rho_i(d)`; the fibrewise S^1-parameter lambda; genus 0.
* **s.8, proof of Thm 2.**
  * eq. (13): an Euler-integral transform of `J_X` equals `I_E sqrt(e(E)) prod_i (2 pi z (lambda + rho_i))^(-1/2) Int_0^oo
    e^((-x_i + (lambda+rho_i) ln x_i)/z) dx_i`. This is termwise integration in Q (as in `Gamma(x+1) = x Gamma(x)`).
  * Then they conclude that the asymptotic expansion of (13) equals `I_E sqrt(e(E)) prod b_(rho_i)`, and a Lemma (unnumbered)
    places it on the untwisted cone.
* **Covers:** the twist half of the mechanism (the Euler integral of a twist = the Stirling / QRR factor; termwise Gamma integrals
  = the I-function, as in the critic's step (e)). It does so at genus 0, in the fibre-equivariant setting, with the asymptotic
  identification asserted.
  * NOT the formal quintic's colour-0 stationary phase at a moving critical point. There is no chart at q = 0.

### 3.3 Iritani arXiv:0903.1463v3 (PDF)

* Theorem 4.14 (= 1.2). Hypotheses: `rho-hat in cl(C~_X)` and Conjecture 4.3, the mirror theorem, proved later. The quantum
  cohomology central charges of `O_X` and `O_pt` are oscillatory integrals over the real Lefschetz thimble and a compact cycle,
  eqs. (70)-(71). Remark 4.15: this solves a connection problem for the QDE in z.
* s.4.4.
  * Thm 4.17: the equivariant oscillatory integral (72) over the real thimble equals `Int H(q, e^(pi i) z) Td~`, for real q and
    z > 0. Remark 4.18: otherwise the right side is the asymptotic expansion as `q -> +0`.
  * Lemma 4.20: near a fixed point sigma, for equivariant parameters with large real part, the thimble integral at `z = 1` has an
    expansion in q to any order M. It is obtained in the chart `{w_j : j not in I_sigma}` (the other `w_i` are monomials in q and
    `w_j`) by Taylor-expanding the q-terms of the exponent and integrating termwise into Gamma functions.
  * Lemma 4.21 compares with the H-function and fixes the connection coefficients.
* **Covers:** the analytic, exact (z = 1) counterpart of the chart elimination and of the critic's step (e), for compact toric
  orbifolds (weak Fano in the main theorems). It contains no `z -> 0` stationary-phase statement and no twist. **Formal quintic: NO.**

### 3.4 Galkin-Golyshev-Iritani arXiv:1404.6407v4 (PDF)

* Thm 1.3.1 (= Thm 5.0.1): Gamma conjectures I and II for `P^(N-1)`.
* s.5.2, eq. (5.2.3): the Mellin solution `phi(t) = (2 pi i)^-1 Int Gamma(s)^N t^(-Ns) ds`.
* Prop. 5.2.4: `phi(t) ~ C t^((N-1)/2) e^(-Nt)` as `t -> oo` in `|arg t| < pi/2`. The proof follows Dubrovin and Meijer: stationary
  phase with Stirling, at LEADING order. The connection formula `phi = Int Gamma-hat cup J` comes from residues.
* Cor. 5.2.6: Gamma conjecture I.
* Non-equivariant Fano; large t (not q = 0); leading order only. **Formal quintic: NO.**

### 3.5 Zagier-Zinger arXiv:0710.0889v1 (PDF)

* eq. (1): `F(w, x) = Sum_d x^d prod_(r<=nd)(nw + r) / prod_(r<=d)((w + r)^n - w^n)`.
  * For n = 5, `lambda_j = zeta^j lambda`, `w = lambda_i/z`, `x = q`, this is EXACTLY the restriction to `p_i` of the formal
    quintic's I-function (checked: `prod_j (lambda_i - lambda_j + kz) = (lambda_i + kz)^5 - lambda^5`).
* Lemma 1: if `M^k F = F`, every x-coefficient of `log F` is `O(w)` as `w -> oo`.
* eq. (5): `F ~ e^(mu w) Sum_s Phi_s w^(-s)` as `w -> oo`, with `mu, Phi_s in Q[[x]]`. It is obtained from Lemma 1, i.e.
  coefficientwise in x.
* Thm 3 (eq. (6)): `mu = Int_0^x (L-1) du/u`, `Phi_0 = L`, `Phi_1 = ((n-2)(n+1)/(24n))(L - L^n)`.
* Thm 4: (i) a first-order ODE system (12) with `Phi_s(0) = delta_(0,s)`; (ii) `Phi_s in L Q[L]`; (iii) `Phi_s in Q(n)[L, L^-1, L^n]`.
* No Stirling factor and no analytic statement at fixed x (s.3 of ZZ is conjectural). So `Phi^ZZ` is normalised by the
  I-function (`Phi_s(0) = delta`), not by QRR. This is the normalisation behind (Gamma3) and LP's R_k (`LIT.md` s.1.1).

### 3.6 Guo-Janda-Ruan arXiv:1812.11908v1, s.7, "A technical result of the formal quintic theory" (PDF): the closest prior art

* **Setting.**
  * The theory is the `(C^*)^m`-equivariant, O(m)-twisted theory of `P^(m-1)`. For m = 5, `lambda_i = zeta^i lambda`, t = 0 (s.3.1)
    it is the formal quintic.
  * `I_alpha(q,z)` is the integral over the m Lefschetz thimbles `gamma_alpha` in `(C^*)^m` of
    `e^(W/z) dx_0 ^...^ dx_(m-1)/(x_0...x_(m-1))`. GJR follow Givental ([15], [12] there): this is Givental's equivariant mirror.
  * The display of W is garbled in my extraction. s.7.1 gives its critical points: `(x_i)_alpha = L_alpha - lambda_i`,
    `y_alpha = m L_alpha`, with `prod_i (L_alpha - lambda_i) = (m L_alpha)^m q` and `L_alpha = lambda_alpha + O(q)`. At
    `lambda_i = zeta^i lambda` this is `L_alpha = zeta^alpha lambda L`.
  * These are the critic's critical data (`x_j = zeta^j - L`, `y = 5L`) up to `x -> -x`.
* **Thm 7.1** (at `lambda_i = zeta^i lambda`). The normalised `z -> 0` expansion (34) of `I_alpha` has three properties:
  (1) regularity at `L_alpha = 0`; (2) a degree bound; (3) polynomiality in m. Its proof combines Prop. 7.5, Lemma 7.8 and Lemma 7.9.
* **Prop. 7.5** (proved; formal stationary phase).
  * The expansion is Wick's at the critical point. Each coefficient is a polynomial in `(Q^-1)_ij`, `(x_i)_alpha^-1` and the higher
    derivatives of W.
  * All of these are rational in `L_alpha` and regular at `L_alpha = 0`. Cor. 7.6 gives the same for `(z q d/dq)^k I_alpha`.
* **Lemma 7.7** (proved, at the level of formal stationary phase). This gives the values at `L_alpha = 0`, i.e. `q^-1 = 0`, the
  orbifold point.
  * At `L_alpha = 0` the integrand of `L_alpha^-k (z q d/dq)^k I_alpha` becomes a product of one-variable Euler integrals. The
    value is therefore the Stirling series with Bernoulli polynomials (Hermite; Nemes; = DLMF (5.11.8)):
    `exp(Sum_i Sum_(j>=1) (-1)^(j+1) B_(j+1)(k/m) z^j/(j(j+1) lambda_i^j))`, as I read the display.
  * At `lambda_i = zeta^i lambda` this is eq. (37): `exp(m Sum_(j>=1) (-1)^(mj+1) B_(mj+1)(k/m) z^(mj)/(mj(mj+1) lambda^(mj)))`.
  * **For m = 5, k = 1, lambda = 1 this is `exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, exactly the project's `Ginf`.**
    Stage 4 lists `Ginf` as CONJECTURAL, VERIFIED for `m <= 100` at two primes and `m <= 25` exactly.
  * The factor m in front of the sum in (37) is a separate glyph in my extraction; with it, (37) is exactly the specialisation of
    the general display. This should be confirmed against the typeset PDF.
  * The continuation to `L_alpha = 0` is of the same kind as the critic's step (c): the coefficients are rational in `L_alpha` and
    regular at 0, so they are evaluated there and compared with the SP of the limiting integrand. The parameter lemma (L1c) is used
    without being stated.
* **Lemma 7.9** (s.7.2). The degree bound comes from the Picard-Fuchs recursion for the `r_k`, with initial data (37) at `L = 0`
  (`X = 1`). This is the analogue of the project's (HPF) recursion with the free coefficients fixed at the orbifold point.
* **Remark 7.2** (NOT proved). It states
  `(2 pi z)^(+-m/2) e^(h_alpha log q/z) C_alpha(z) I_alpha(q,z) = Itilde_alpha(q,z)`, where Itilde is ZZ's series (their (23)) and
  `C_alpha` are the QRR constants of their (17). The sign of the exponent is not legible in my extraction.
  * No argument is given, and no reference beyond the opening attribution of the integrals to Givental [15, 12].
* **Cor. 7.3** ("a direct consequence"; ZZ's conjecture [35, eq. (39)], in a stronger form) and **Cor. 7.4** (m = 5, the structure
  of `[z^k] R_(j0)`).
  * The identification of the SP coefficients with R-matrix entries, `[z^l] R_(4 alpha) = r_(alpha,0,l)` and
    `[z^l] R_(0 alpha) = r_(alpha,1,l)`, is made through Remark 7.2.
  * Lemma 5.5 (`R_0 = Id`; `R_k` in the degree-k subspace) uses Cor. 7.4.
* **Assessment (mine).**
  * Read coefficientwise in q, Remark 7.2 says that `C_alpha` times the normalised SP series of `I_alpha` is the normalised ZZ
    series, i.e. `Phi^SP_LG = h Phi^ZZ`. This is the project's (Gamma3), with GJR's LG integral in place of the full-Stirling
    Mellin-Barnes integral.
  * The two SP series should agree by (L1) + (L2): both are stationary phases of the same mirror integral, and the critic's
    7-variable form contains both. This identification is written neither in GJR nor in the project.
  * Read analytically (fixed q, `z -> 0`), Remark 7.2 is the uniformity statement (Gamma6)/(A2) of Stage 4, because ZZ's (5) is
    coefficientwise in q.
  * Hence GJR s.7 PROVES the orbifold-point normalisation (q = oo; Lemma 7.7): the analogue at `L = 0` of the critic's step (d) at
    `L = 1`. It STATES WITHOUT PROOF the large-radius normalisation (q = 0; Remark 7.2), which is what the chart argument would prove.
  * Their Cor. 7.3-7.4 and Lemma 5.5 depend on Remark 7.2. v1 is the only version; the abstract page lists no journal.
  * **Covers the O(5)-twisted formal quintic: YES** (m = 5, `lambda_i = zeta^i lambda`, t = 0).

### 3.7 Checked and set aside

* Coates-Iritani arXiv:1203.4193 ([H]). Convergence of Givental's formula, in the conformal, non-equivariant case, where R is fixed
  by Euler homogeneity. There is no q = 0 Stirling normalisation.
* Lho-Pandharipande 1702.06096 and 1803.01409: no stationary phase; the ZZ normalisation (`LIT.md` s.1.1).
* Iritani Thm 4.17 and the Gamma-integral structure: exact identities at fixed z, not the `z -> 0` expansion at q = 0.
