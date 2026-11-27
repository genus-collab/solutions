# Stage 5 / lemmas: the formal stationary-phase lemmas (2026-11-14)

Folder `allgenus/stage5/lemmas/` (task key `lemmas`). New files only; no project module is imported. AI-written, unrefereed; no human
has read it. Labels: PROVED (argument written out here, or a quoted statement with its number, read by me), VERIFIED (exact
computation, range and field stated), CONJECTURAL, SKETCH, REFUTED. Finite checks are evidence, never proof.
Parts: this file (s.0-s.2), `FORMAL_SP_LEMMAS_2.md` (s.3-s.5: FL1, translation, parameters L1c, FL4, changes of variables L1a),
`FORMAL_SP_LEMMAS_3.md` (s.6-s.8: Fubini L1b, delta lemma L2, the interface for the chart argument), `FORMAL_SP_LEMMAS_4.md`
(s.9-s.10: Stirling L3, the residue FL2-FL3), `FORMAL_SP_LEMMAS_5.md` (s.11-s.12: checks, caveats, record).
Users: the chart argument (STAGE4_SYNTHESIS s.1.4, written out by `stage5/chart/`) and Theorem P(b) (POLY s.6 with the poly referee's
REFEREE s.3 and item C6). Stage-4 notation (z, zeta, u, L, b_k, h, Gamma_St, c(z; L)) as in the Stage-5 brief.

## 0. Bottom line

* **Every lemma requested is PROVED here, algebraically, over an arbitrary commutative Q-algebra of coefficients.** The formal SP is
  defined by Wick contraction of Taylor data (s.1). The proofs use only: the Wick recursion (= Gaussian integration by parts),
  the Leibniz rule, the Piola identity, formal implicit functions, and ONE exact Gaussian integral (s.2, (G5), via the
  polynomial-identity principle). No asymptotic analysis, no uniformity estimate and no Laplace-method theorem is used, except in the
  second (classical) proof of L3, which quotes DLMF.
* **Status table** (all PROVED; checks are exact, own code, each with failing controls, s.11):

| label | statement (short) | where | exact check (own code) |
|---|---|---|---|
| (U) | each z^k-coefficient of SP^norm is a universal polynomial in finitely many Taylor coefficients and in the entries of -Hess^-1 | s.1.5 | (used by all) |
| (G1)-(G5) | Gaussian toolkit: generating function, linear substitution, integration by parts, block factorisation, reweighting | s.2 | via L1b, T-q |
| FL1 | formal SP kills exact forms (also convergent infinite sums of them) | s.3.1 | `chk_change`: n = 1, 2 to z^8, n = 3 to z^6 |
| (T) | translation invariance (by `p = O(q)` or `p = O(z)`); moving the base point (T') | s.3.2 | `chk_param` (T(b) with `p = -z`; T-q, T-z) |
| L1c | parameters: (i) analytic / formal dependence, (ii) Taylor expansion and evaluation, (iii) perturbation in q (T-q), (iv) z-dependent part of the phase moved into the amplitude (T-z) | s.4 | `chk_param`: weight <= 8, n = 1, 2 |
| FL4 | the full formal SP commutes with d/dq (and D_q = q d/dq) | s.4.5 | `chk_param` |
| L1a | invariance under analytic, formal and parameter-dependent changes of variables | s.5 | `chk_change` (random maps, n = 1, 2 to z^8, n = 3 to z^6) |
| FL3 (SP) | invariance under the CFU map termwise in u (L1a over Q[[u]]) | s.5.4 | `chk_param` (moving critical points) |
| L1b | Fubini: total SP = iterated SP, with the inner critical point, Gaussian factor and series entering the outer integral as stated | s.6 | `chk_fubini`: (1,1), (1,2), (2,1) to z^8, (2,2) to z^6, 3 random phases each |
| L2 | delta lemma: restriction to G = 0 times the Gaussian factor, exactly | s.7 | `chk_delta`: 2, 3 variables to z^8, 4 to z^6 |
| L3 | Stirling: formal SP of Euler's integral = z^(a/z+m) Gamma_St(a/z+m), all m in Z, both measures; Gamma_St(y+1) = y Gamma_St(y) | s.9 | `chk_stirling`: z^40 and z^80 (a = 1), z^12 (a = 3/7, -2), m = -2..5 |
| FL2 | the formal residue at X = 0 (termwise in u) kills exact forms | s.10 | `chk_res`: u^6, w^-14 |
| FL3 (residue) | the residue is invariant under X -> X(1 + O(u)) termwise in u | s.10 | `chk_res` |

* **Who uses what.** R5 (Theorem P(b), POLY s.6 with REFEREE s.3): the poly referee's four lemmas (C6) are FL1, FL2, FL3 (s.5.4 for
  the SP, s.10.3 for the residue) and FL4; Theorem P(a) also uses `Gamma_St(y+1) = y Gamma_St(y)` (= L3(iv)) and the shift
  `sigma -> sigma + z` (Lemma T(b), or FL1 with a convergent sum). The chart argument uses L1a, L1b, L1c, L2, L3 and (G4) (s.8).
  Step 6 of Theorem P(b) also relates the residue in the original and the hatted variables; that relabelling is a SKETCH (s.10.4).
* **Interface** (s.8): for each step (a)-(e) of the chart argument, the lemma used, its hypotheses, and where they can fail (in
  particular: the 7-variable critical point needs `L != 0, 1, zeta^j`; at `L = 1` only the chart integral, not the 7-variable one,
  is defined; step (e) needs L3 with NEGATIVE m).
* **One bug of my own, found and fixed** (s.11.3): the first FL4 check failed; the cause was a truncation error in my test harness
  (an amplitude of negative weight with pre-truncated phase data), not the lemma.
* **Not done here**: the chart argument itself (other agent); Theorem P(b) itself; Theorem U; (QRR).

## 1. Setting and definitions

### 1.1 Coefficients and variables
`A` is any commutative Q-algebra (in the applications: Q, Q(zeta), Q(zeta)(L), C, Q[[u]], Q(zeta)[[sigma - L]], rings of convergent
power series in a parameter). `x = (x_1, .., x_n)`; `A[[t]]` formal power series in `t = x - c` for a point `c in A^n`. An analytic
germ at `c` enters only through its Taylor series. `z` is a formal variable; `A[[t]]((z))` = Laurent series in `z` (finitely many
negative powers) with coefficients in `A[[t]]`. Nothing is ever evaluated at `z != 0`.

### 1.2 Weight and the completed ring
`wt(t^alpha z^j) := |alpha|/2 + j`. `R_n(A)` := the formal sums `Sum c_(alpha,j) t^alpha z^j` (`j in Z`) such that for every `W` only
finitely many terms with `c != 0` have weight `<= W`.
**Lemma 1.1 (PROVED).** (a) `R_n(A)` is a ring containing `A[[t]]((z))`. (b) `d/dt_i` and multiplication by `t_i` act termwise and
continuously. (c) If every term of `P in R_n(A)` has weight `>= 1/2`, then `exp(P) = Sum_m P^m/m!` converges in `R_n(A)`.
*Proof.* (a) The finiteness condition implies a minimal weight. In a product the coefficient of `t^gamma z^l` is a sum over pairs of
terms whose weights add up to `wt(t^gamma z^l)`; each weight is bounded below, so only finitely many pairs occur, and only
finitely many products have weight `<= W`. (b) clear. (c) `P^m` has weight `>= m/2`, so below any `W` only `m <= 2W` contribute. QED
**Parameters.** For `A = A_0[[q]]` we use `R_n(A_0)[[q]]` (power series in `q` with coefficients in `R_n(A_0)`); everything is done
coefficientwise in `q`. Alternatively we give `q` weight 1 (ring `R_(n;q)(A_0)`, same finiteness condition); this is what makes
the substitution `q -> z` legitimate in L1c(iv) (s.4).

### 1.3 Wick expectation
For a symmetric `C in A^(n x n)` put `<t^alpha>_C := d_xi^alpha exp(xi^T C xi/2)|_(xi=0)`, extended A-linearly to polynomials. It is
zero for `|alpha|` odd and equals the sum, over perfect matchings of the multiset `{1^(alpha_1), .., n^(alpha_n)}`, of the products of
the `C_ij` of the matched pairs (expand the exponential). It satisfies the **Wick recursion**
(W) `<t_i t^beta>_C = Sum_j C_ij beta_j <t^(beta - e_j)>_C`
(apply `d_xi_i` first: `d_xi_i e^(xi^T C xi/2) = (C xi)_i e^(..)`, then `d_xi^beta` at 0 by Leibniz; `(C xi)_i` is linear).
The **z-scaled functional** `E_zC : R_n(A) -> A((z))`, `E_zC[t^alpha z^j] := <t^alpha>_C z^(j + |alpha|/2)`, extended linearly, is
well defined (it preserves weight; finitely many terms per output power) and A-linear. It is the formal Gaussian integral with
covariance `zC`, normalised to `E_zC[1] = 1`.

### 1.4 The formal stationary phase
**Data.** A phase `F in A[[x - c]]` with `grad F(c) = 0` and `H := Hess F(c) in GL_n(A)` (a *nondegenerate critical point over A*);
an amplitude `a in A[[x - c]]((z))` (more generally `a in R_n(A)` in `t = x - c`). Write
`F(c + t) = F(c) + t^T H t/2 + F_3(t)`, `F_3 in (t)^3`, and `C := -H^-1`. Every term of `F_3/z` has weight `>= 1/2`.
**Definition 1.2 (normalised SP).** `SP^norm_c[F; a] := E_zC[ exp(F_3(t)/z) a(c + t) ] in A((z))`.
**Definition 1.3 (full SP, bookkeeping).** `SP_c[e^(F/z) a dx] := [e^(F(c)/z)] (2 pi z)^(n/2) [Delta^(-1/2)] SP^norm_c[F; a]`,
`Delta := det(-H)`. The bracketed factors are formal symbols (`[e^(al/z)][e^(be/z)] = [e^((al+be)/z)]`; `[Delta^(-1/2)]` a chosen
square root in an extension of `A`). Every identity between full SPs below is an identity between normalised series together with
an explicit identity between these symbols, under a stated rule for the square-root branches. The project normalises series by
their constant term (`c_0 = 1`), so only the z-dependence of the symbols matters there, and that is always explicit.
*Meaning (not used in proofs).* For real-analytic real data with `H` negative definite, the full SP is the Laplace asymptotic
expansion, as `z -> 0+`, of the integral of `e^(F/z) a dx` over a small neighbourhood of `c`, with the positive root (s.9.3 in one
variable). The phase convention of the Stage-4 notes `e^(Phi/z)` (gamma) and `e^(-W/w)` (poly) is `F = Phi`, resp. `F = -W`, `z = w`.

### 1.5 Universality
**Lemma U (PROVED).** Write `a = Sum_(j >= j0) a_j z^j`. Then
`[z^k] SP^norm_c[F; a] = Sum (1/m!) [t^alpha](F_3^m a_j) <t^alpha>_C`, the sum over `m >= 0`, `j0 <= j <= k` and `|alpha| = 2(k - j + m)`.
Since `F_3^m in (t)^(3m)`, only `m <= 2(k - j)` and `|alpha| <= 6(k - j)` occur. So `[z^k]` is a polynomial with rational coefficients
in finitely many Taylor coefficients of `F` at `c` (orders `3..6(k - j0)`), of the `a_j` (orders `<= 6(k - j0)`) and in the entries of
`C = -H^-1`. *Proof:* expand Definition 1.2. QED
**Consequences.** (U1) Every ring homomorphism `rho: A -> B` commutes with `SP^norm` (apply `rho` to all data; `rho(H)` stays
invertible). (U2) If the Taylor data at `c` and `C` depend holomorphically (resp. formally) on parameters, so does every coefficient.

## 2. The Gaussian toolkit (all PROVED)

**(G1) Generating function.** `<e^(xi . t)>_C = e^(xi^T C xi/2)` coefficientwise in `xi` (the definition in s.1.3).
**(G2) Linear substitution.** For `J in A^(n x m)` and a polynomial `P`: `<P(J y)>_(C') = <P(t)>_(J C' J^T)`. *Proof:*
`<e^(xi . J y)>_(C') = <e^((J^T xi) . y)>_(C') = e^(xi^T J C' J^T xi/2)` by (G1); compare coefficients. For invertible `J` the substitution
`t = J y` preserves weight, so (G2) extends to `R_n(A)`, also for `E_zC`.
**(G3) Gaussian integration by parts.** If `C` is invertible, `H := -C^-1`, then for every `g in R_n(A)`:
`E_zC[ d_i g ] = -z^-1 E_zC[ (H t)_i g ]`, i.e. `E_zC[ d_i g + z^-1 (H t)_i g ] = 0`.
*Proof.* Multiplying by `C`, the claim is `E_zC[t_k g] = z Sum_j C_kj E_zC[d_j g]`. For `g = t^beta z^l` both sides equal
`z^(l + (|beta|+1)/2) Sum_j C_kj beta_j <t^(beta - e_j)>_C` by (W). Both sides are linear and involve finitely many terms of `g` per
output power, so the identity extends to `R_n(A)`. QED
**(G4) Block factorisation.** If `C = diag(C_x, C_y)` then `E_zC[P] = E_(zC_y)[ E_(zC_x)[P] ]` for `P in R_(n+m)(A)`, the inner functional
acting on the x-monomials with y and z as coefficients. In particular, for phases and amplitudes in separate variables,
`SP^norm_((c1,c2))[F_1(x) + F_2(y); a_1(x) a_2(y)] = SP^norm_c1[F_1; a_1] SP^norm_c2[F_2; a_2]`. *Proof:* on monomials no matching
pairs an x with a y; continuity. QED
**(G5) Reweighting.** Let `C_0 in Sym_n(A)` and let `k` be a symmetric matrix of indeterminates. In `A[[k]]`:
`<e^(t^T k t/2) t^alpha>_(C_0) = det(I - C_0 k)^(-1/2) <t^alpha>_(C_1)`,  `C_1 := Sum_(r>=0) (C_0 k)^r C_0 = (I - C_0 k)^-1 C_0`,
with `det(I - C_0 k)^(-1/2) := exp(Sum_(r>=1) tr((C_0 k)^r)/(2r))`.
*Proof.* For every monomial in the entries of `k`, both coefficients are polynomials over Q in the entries of `C_0` (left:
`<(t^T k t/2)^r t^alpha>_(C_0)/r!`; right: product of the two series). A polynomial in the `n(n+1)/2` entries of a symmetric matrix that
vanishes on the nonempty open set of real positive-definite matrices vanishes identically; an identity of polynomials over Q then
holds after substitution in any Q-algebra. So let `C_0` be real positive definite and `k` real symmetric with `||k|| ||C_0|| < 1`.
With the Gaussian measure `N_C` (covariance `C`) and Fubini/dominated convergence (`e^(|t^T k t|/2)` is `N_(C_0)`-integrable):
`int e^(t^T k t/2) t^alpha dN_(C_0) = det(C_0)^(-1/2) det(C_0^-1 - k)^(-1/2) int t^alpha dN_(C_1) = det(I - C_0 k)^(-1/2) int t^alpha dN_(C_1)`
(complete the square: `C_1^-1 = C_0^-1 - k`), with the positive root. The moments of `N_C` are `<t^alpha>_C` ((G1) is its
moment generating function). Both sides are real-analytic in `k` near 0 and their Taylor series are the two formal sides (on the
left by termwise integration, justified by the same domination). QED This is the only real integral used in s.1-s.8; it is exact.
**(G5') With z.** For `K_1 in Sym_n(I)`, `I` an ideal of topologically nilpotent elements (entries in `(y)A[[y]]`, in `qA_0[[q]]`, or of
weight `>= 1/2`), `H_0 := -C_0^-1` (`C_0` invertible), `K := H_0 + K_1`, and `P in R_n(A)` (resp. `R_n(A_0)[[q]]`):
`E_(zC_0)[ e^(t^T K_1 t/(2z)) P ] = (det K/det H_0)^(-1/2) E_(-zK^-1)[ P ]`.
*Proof.* Rescaling, `E_(zC_0)[e^(t^T K_1 t/(2z)) t^alpha] = z^(|alpha|/2) <e^(t^T K_1 t/2) t^alpha>_(C_0)`; substitute `k = K_1` in (G5) (legitimate:
`K_1` is topologically nilpotent); `I - C_0 K_1 = H_0^-1 K` and `C_1 = -K^-1`. QED
