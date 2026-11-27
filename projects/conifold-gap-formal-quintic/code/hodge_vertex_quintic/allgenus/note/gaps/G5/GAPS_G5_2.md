# Gap closing, cluster G5, part 2: Corollary B (`cor:mai-gap`) and the constant of [CGYZ26, Conj. 1.5]

Continuation of `GAPS_G5.md` (labels, notation, s.0 summary). Same rules.

## 2. Item 2: `cor:mai-gap`

### 2.1 What was missing, and what was found
`OUTLINE_8.md` gives the derivation in one line with `\draftnote{not stated as a corollary in the sources; each step is}`; the review
(M2) asks to cite the group property of the unit-preserving action; DL s.7 item 4 asks for (a) the reconstruction (Teleman), (b) the
conifold-frame map, (c) the comparison theorem; DA A15 marks the corollary [M]. Searching HVQ:
* **FOUND (a written proof of the composition step):** `natural_form/genus4_referee/review/REVIEW.md` s.2.1-s.2.3 (Steps 0-5) proves,
  in the Feynman (decorated-graph) form used by the outer sum, for `R = R' M` with `R'`, `M` symplectic and `= 1` at `z = 0`: the
  point-frame total of `R` equals the sum of the outer terms with vertices `M.omega^` and edges, leaves of `R'` (edge splitting
  `E_R = E_M + M^-1 E_R' M^-T`, `T_R = T_M + M^-1 T_R'`, orbit-stabiliser on marked graphs, the exact sequence
  `1 -> prod_v Aut(gamma_v) -> Aut(Gt, S) -> Stab(d) -> 1`, the leaf bookkeeping with the vanishing of the unstable genus-0 pieces by
  the dimension count). Labelled ARGUED there; I checked it line by line and find it complete for `n = 0` (its Remark (ii): no property of
  `tau` is used). It is the same statement as [PPZ15, Prop. 2.14] specialised to `n = 0`, so the note may cite PPZ and need not
  transcribe it. S15 s.5 item 3 says that this regrouping is "needed only to identify the TOTAL with `F_g`".
* Nothing in HVQ states the corollary for all `g` (NF s.6 verifies the total at genus 5; FQ2 s.4 at genus 2).

### 2.2 Definitions (to be printed in s02/s10)
* **The ring.** `NN` is the presented differential ring of `def:set-confexp`: the localised polynomial ring
  `Q[L^(+-1), Lam, c_0, .., c_3, i_0, i_1, i_2][c_0^-1, W^-1]`, `W := c_0^2 + c_0 i_1 - c_1 i_0` (`= C_0^2 C_1`), with the derivation `D`
  fixed on the generators (`D c_j = c_(j+1)`, `D c_3` by PF4, `D i_j = i_(j+1)`, `D i_2 := i_3` solved from the concomitant relation, which is
  linear in `i_3` with unit coefficient `u c_0` up to sign, `D Lam = 1`, `D L = (L^6 - L)/5`). Its evaluation `ev: NN -> ` (functions of `q`),
  `c_j -> D^j C_0`, `i_j -> D^j I~_1`, `Lam -> log q`, is a differential homomorphism. (`Lam` is NOTATION's "formal `log q`", R8 (h).)
* **The conifold expansion** in the frame `(alpha, beta)`: `Cf_(alpha,beta) : NN -> Q((s))`, the ring homomorphism with
  `c_j -> D^j varpi_0`, `i_j -> D^j varpi_1 - j D^(j-1) varpi_0 - Lam' D^j varpi_0`, `Lam -> Lam' := log(1 - u)` (`D Lam' = 1`; any primitive of
  `1` gives the same value on `Lam`-free elements), `L -> s^-1`. `c_0` and `W` go to units (`varpi_0(0) = 1`, `D t_c(0) = -1`). It commutes
  with `D` by `lem:gau-confexp` (PF4 for `varpi_0`; the concomitant `CC(varpi_0, varpi_1) = 0`, `lem:gau-iso`). `C_0 -> varpi_0`, `I_1 -> varpi_1`,
  `T -> t_c`, `C_1 -> D t_c`. No relation beyond PF4 and the concomitant is imposed, so `Cf` needs no admissibility hypothesis.
* **The potential in `NN`.** By `prop:set-cyclicLR`, `thm:gau-R0` (iii) and `cor:gau-rowsNN` (cluster G1, R1), every row coefficient
  `Q^(i)_(a,n)` lies in `NN (x) Q(zeta)` (`Lam`-free). Put `Psi := Q(z = 0)`, `R~ := Psi^-1 Q = Psi^T eta^-1 Q` (Teleman form; `Psi^T eta^-1 Psi = 1`,
  s.3.2), `Delta_i^(1/2) = zeta^(4i) C_0/L`. Define
  `F_g^NN := Sum_(Gamma in G_(g,0)) |Aut Gamma|^-1 Cont_Gamma[E_(R~), q_(R~); omega]` in `NN (x) Q(zeta)`,
  the genus-`g` leg-free graph sum of `def:set-graphsum` ([PPZ15] Defs. 2.2, 2.13), with the edge written without division:
  `E_X(psi, psi') := (N_X(psi, psi') - N_X(-psi', psi'))/(psi + psi')`, `N_X := 1 - X^-1(psi) X^-1(psi')^T` (a power series for every `X`; equal to
  PPZ's edge when `X` is symplectic, [PPZ15] Rem. 2.1). At genus `g` it is a polynomial with rational coefficients in finitely many row
  coefficients, in `Delta_i^(+-1/2)` and in `zeta`. By the colour covariance (`def:set-rows`) the Galois group of `Q(zeta)` permutes the
  colours, so `F_g^NN` lies in `NN`. By `thm:set-GT` ((T), with (QRR) `thm:set-qrr`): `ev(F_g^NN) = F_g`.
* **Definition (conifold expansion of `F_g`).** `F_g^cf := Cf_(alpha,beta)(F_g^NN)` in `Q((s))`. "The gap of the total" means
  `F_g^cf - c_g t_c^(2-2g)` lies in `Q[[u]]`; as `t_c = u + O(u^2)` lies in `Q[[u]]`, `Q[[u]] = Q[[t_c]]`, this is "no fractional power of `u`,
  no `t_c^k` with `3 - 2g <= k <= -1`, leading term `c_g t_c^(2-2g)`".

### 2.3 Statement (Corollary B; `cor:mai-gap`)
Let `g >= 2` and `(alpha, beta)` in `Q^2`, `varpi_1 = u + (7/10) u^2 + ..`, `t_c = varpi_1/varpi_0`. Then `F_g^cf = c_g t_c^(2-2g) + (element of Q[[u]])`,
`c_g = (-1)^(g-1) B_2g 5^(1-g)/(2g(2g-2))`. Inputs: Theorem A (`thm:mai-main`); (T) `thm:set-GT`; [PPZ15] Def. 2.13 and Prop. 2.14;
`lem:flh-frame`; `def:set-confexp`, `lem:gau-confexp`, `cor:gau-rowsNN`; `thm:gau-R0` (iv), `thm:gam-D`, `prop:gau-symp`, `def:set-M0`.

### 2.4 Proof
*Step 1 (expand the graph sum).* `Cf := Cf_(alpha,beta)` (extended `Q(zeta)`-linearly) is a ring homomorphism and `F_g^NN` is a polynomial
expression in finitely many elements of `NN (x) Q(zeta)` (the division-free edges, the leaves `q_(R~) = R~^-1(z)(-z 1_omega)`,
`1_omega = Sum_i Delta_i^(-1/2) e~_i`, the ancestor correlators of `omega`, which are rational multiples of powers `Delta_i^((2h-2+n)/2)`,
`2h - 2 + n > 0`). Hence `F_g^cf = Sum_Gamma |Aut Gamma|^-1 Cont_Gamma[E_(Cf R~), q_(Cf R~); Cf omega]`. By `thm:gau-R0` (iv) and `D = 1`
(`thm:gam-D`), `Cf(Q) = P_can Q^m(v(u))`, a product of symplectic matrices (`prop:gau-symp`; the model rows, `def:set-M0`); `Cf(Psi)` is
orthonormal; so `Cf(R~)` is symplectic and the division-free edge of `Cf(R~)` is PPZ's edge.
*Step 2 (regroup).* Let `M := M_0(tau(u) zhat)`, `zhat = z/s^4`: as a series in `z` over `Q(zeta)((s))`, `M(0) = 1` and `M` is symplectic.
Then `X := Cf(R~) M^-1` is symplectic with `X(0) = 1`, and `Cf(Q) M^-1 = Cf(Psi) X =: R'` is the outer matrix of `def:set-outer` (as
corrected by R6). By [PPZ15, Prop. 2.14] (the unit-preserving action of Def. 2.13 is a left group action),
`Cf(R~).(Cf omega) = X.(M.(Cf omega))`. Take the genus-`g`, `n = 0` parts and integrate over `Mbar_g`: the left side is the sum of Step 1;
the right side is `Sum_Gamma |Aut Gamma|^-1 Cont_Gamma[E_X, q_X; M.(Cf omega)]` (Def. 2.2 with the translation of Def. 2.13 written as
dilaton-shifted ancestor potentials, `def:flh-outer`). Since `Cf(Psi)` is an isometry onto `(H, eta)` mapping `1_omega` to `1`, `lem:flh-frame`
(for `n = 0`: `E_(Psi X) = E_X`, `q_(Psi X) = q_X`) turns each term into `Cont_Gamma[E_(R'), q_(R'); M.omega] =: T_Gamma`, with `omega` the TFT
expanded in the frame as in `def:set-outer`. So `F_g^cf = Sum_Gamma |Aut Gamma|^-1 T_Gamma`.
(If the note's `T_Gamma` already carries `|Aut Gamma|^-1`, drop the factor.) The written proof found in s.2.1 is this step.
*Coefficient ring.* PPZ work over `Q`; the proofs of their Props. 2.4, 2.8, 2.9, 2.12, 2.14 (and of REVIEW s.2.3) are identities between
finite sums, polynomial in the coefficients of the R-matrices and of the CohFT, valid over any commutative `Q`-algebra; here
`Q(zeta)((s))`.
*Step 3 (Theorem A).* For every `Gamma` with an edge, `T_Gamma` lies in `Q[[u]]`; the vertex graph gives `c_g t_c^(2-2g) + Q[[u]]`.
*Step 4.* Summing the finitely many graphs: `F_g^cf` lies in `c_g t_c^(2-2g) + Q[[u]]`. QED.

### 2.5 Frame scope
* **As proved:** the frames of `def:set-frames`, `(alpha, beta)` in `Q^2` (D14), with `varpi_0(0) = 1` and `varpi_1 = u + O(u^2)`, exactly the
  scope of Theorem A; nothing else is used.
* **Remark (complex frames; proof).** Each `s`-coefficient of `Cf_(alpha,beta)(F_g^NN) - c_g t_c^(2-2g)` is a polynomial in `(alpha, beta)` with
  rational coefficients (the images of `c_j`, `i_j` are linear in `(alpha, beta)`; `c_0`, `W` go to series with constant terms `1`, `-1` and
  coefficients in `Q[alpha, beta]`; `t_c/u` is a unit of `Q[alpha, beta][[u]]`). The coefficients of the negative and the non-integral powers
  of `u` vanish on `Q^2`, hence identically: Corollary B holds for `(alpha, beta)` in `C^2`. (The same applies to each `T_Gamma`.)
* **Remark (scalings; proof).** For `lambda, mu` in `C^x` take the pair `(omega_0, omega_1) := (mu varpi_0, lambda varpi_1)`, `tau_c := omega_1/omega_0
  = (lambda/mu) t_c`, and `Cf'` the conifold expansion with `varpi_0, varpi_1` replaced by `omega_0, omega_1`. Then
  `Cf'(F_g^NN) - c_g lambda^(2g-2) tau_c^(2-2g)` lies in `C[[u]]`; the constant does not depend on `mu` and is `c_g` iff `lambda^(2g-2) = 1`.
  *Proof.* The endomorphism `sigma` of `NN` with `c_j -> mu c_j`, `i_j -> lambda i_j + (lambda - mu)(Lam c_j + j c_(j-1))` (i.e. `C_0 -> mu C_0`,
  `I_1 -> lambda I_1`), `Lam -> Lam`, `L -> L`, commutes with `D` (PF4 is linear, the concomitant bilinear) and `Cf' = Cf o sigma`.
  `sigma(C_1) = (lambda/mu) C_1`, so by the row recursion of `prop:set-cyclicLR` the rows scale by `diag(mu^-1, lambda^-1, lambda, mu, 1)`, an
  isometry of `eta^-1` (`nu`-pairs `(0,3)`, `(1,2)`, `(4,4)`); `Psi` scales the same way, so `sigma(R~) = R~`, while `Delta_i^(1/2) -> mu Delta_i^(1/2)`,
  i.e. `omega -> mu.omega` (TFT rescaling). By `lem:flh-comm`, `R~.(mu.omega) = mu.(R~.omega)`, whose `(g,0)` part is `mu^(2g-2) F_g`: so
  `sigma(F_g^NN) = mu^(2g-2) F_g^NN` and `Cf'(F_g^NN) = mu^(2g-2) F_g^cf = c_g lambda^(2g-2) tau_c^(2-2g) + mu^(2g-2) Q[[u]]`. QED.
  In particular the constant `c_g` refers to the normalisation `varpi_1 = u + O(u^2)`; `lambda = sqrt 5` gives `(-1)^(g-1) B_2g/(2g(2g-2))`, the
  form of [HKQ09, (1.2)] as DIGEST_L_4 s.4 transcribes it (not re-read by me; FQ2 s.4 observed the factor `5` at genus 2).

### 2.6 The constant of [CGYZ26, Conj. 1.5]: how `c_g` relates to it (CITED, read)
*Read:* arXiv:2609.13921v2 (current version, 2026-10-31), s.1 (Defs. 1.1, 1.2, Ex. 1.3, Rem. 1.4, Conj. 1.5, Def. 1.8, Thm 2, Rem. 1.9),
and the restatements in s.4.2 (Thm 5, Ex. 4.7) and s.5.1, s.5.3 ((5.2), (5.5)), in the arXiv HTML rendering through a text-extraction tool
(quotations limited to short fragments; the formulas of Ex. 1.3 and Conj. 1.5 were confirmed in two independent readings, those of
Defs. 1.1-1.2 and Rem. 1.4 in one; `\draftnote` for the writer: check them once in the PDF). Their content, in our words:
* `K = Q(q)`, `D = q d/dq`; `I_0 = 1 + O(q)`, `I_1 = I_0 log q + O(q)`, `tau = I_1/I_0`; `R` = the differential `K`-algebra generated by
  `I_0^(+-1)`, `(D tau)^(+-1)` (Def. 1.1); `F_g(Q(q))` the genus-`g` potential at `Q = e^tau`, an element of `R` for `g >= 2` (Thm 2).
* An admissible pair `(omega_0, omega_1)` of holomorphic periods at the conifold: `omega_0(0) != 0`, `omega_1` vanishing simply, the
  logarithmic period of the form `omega_1 log Delta + O(1)`, and the substitution `I_0 -> omega_0`, `D tau -> D tau_c`, `tau_c := omega_1/omega_0`,
  extended to derivatives and inverses, must preserve the relations of `R`; it defines the differential homomorphism
  `iota : R -> R^coni` (inside `C((Delta))`), and `f_c := iota(f)` is the conifold expansion (Def. 1.2), `F_g^c(tau_c) := iota(F_g)`.
* Ex. 1.3 (quintic): `Delta = 1 - 5^5 q`, `omega_0 = 1 + O(Delta^3)`, `omega_1 = Delta + (7/10) Delta^2 + O(Delta^3)`; Rem. 1.4: with `omega_1`
  fixed, `omega_0 + c_1 omega_1 + c_2 omega_2` (`c_i` in `Q`, `omega_2 = Delta^2 + O(Delta^3)`) is again admissible (same Picard-Fuchs and
  Wronskian identities); Thm 2 is stated for "normalized `omega_1`" (Ex. 4.7: `Delta + O(Delta^2)`).
* Conj. 1.5 (proved for the quintic as Thm 2): for every admissible pair and `g >= 2`,
  `F_g^c(tau_c) - (-1)^(g-1) B_2g/(5^(g-1) 2g(2g-2) tau_c^(2g-2))` lies in `C[[tau_c]]`.
*Dictionary.* The formal quintic has the quintic's PF4 (`C_0 = Sum (5d)!/(d!)^5 q^d`), so `R` is literally the same ring of periods, and:
`Delta = u`; `omega_0 = f_0` (the only power-series solution `1 + O(u^3)`); `omega_2 = f_2`; `omega_1 = varpi_1` (the only solution
`u + (7/10) u^2 + O(u^3)`); hence Ex. 1.3 is the frame `(alpha, beta) = (0, 0)`, `tau_c = t_c`, and Rem. 1.4's family is `(alpha, beta) = (c_1, c_2)`
in `Q^2`, exactly `def:set-frames`. `iota` and `Cf_(alpha,beta)` are differential homomorphisms for the same `D`, both expand `K` at `Delta = 0`
(`q = (1 - L^-5)/5^5 -> (1 - u)/5^5`), and `iota(I_0) = omega_0`, `iota(D tau) = D tau_c` are `Cf(C_0) = varpi_0`, `Cf(C_1) = D t_c`. So on the
common subring, the differential ring generated by `C_0^(+-1)`, `C_1^(+-1)` over `Q[q, (1 - 5^5 q)^-1] = Q[L^5, L^-5]` (it lies in `R` and in `NN`),
the two expansions coincide wherever `iota` is defined: two differential homomorphisms that agree on the generators agree there.
*The constant.* `c_g = (-1)^(g-1) B_2g 5^(1-g)/(2g(2g-2))` is the constant of [CGYZ26, Conj. 1.5] exactly: same sign, same power of 5, same
normalisation of the vanishing period (`omega_1 = Delta + O(Delta^2)` against `varpi_1 = u + O(u^2)`, `Delta = u`). By s.2.5, for every
admissible pair of CGYZ's family with normalised `omega_1` (any `omega_0 = mu(f_0 + c_1 varpi_1 + c_2 f_2)`) the formal quintic has the gap
with the same constant. So Corollary B is the formal-quintic analogue of [CGYZ26, Conj. 1.5], with the same constant, in every
admissible pair with normalised vanishing period.
*What is not claimed, and the differences (for `rmk:mai-constant`).*
1. Two different theories: the formal quintic's `F_g` is not the quintic's (`F_2(0) = -19/9600` against `-5/144`, appendix C); the equality
   of the constants is a statement about two different potentials. Nothing about the quintic follows.
2. The formal quintic's rows involve `L` (eigenvalues `zeta^i L`), and `L` is not in `K = Q(q)`; `NN` contains `L^(+-1)` and `Cf(L) = s^-1 =
   u^(-1/5)`, outside `C((Delta))`. The total `F_g^cf` is nevertheless a Laurent series in `u` (Corollary B). Whether `F_g^NN` lies in CGYZ's
   `R` (the `L`-free part) is neither needed nor proved here.
3. `iota` is defined on functions under the admissibility assumption; `Cf` is defined on the presented ring `NN` (no assumption,
   `lem:gau-confexp`), and the conifold expansion of `F_g` is that of its graph-sum representative `F_g^NN`. They agree on the common
   subring above wherever `iota` is defined. The injectivity of `ev` on `NN` (algebraic independence of `D^j C_0` (`j <= 3`), `D^j I~_1` (`j <= 2`), `L`, `log q` over `Q`
   modulo the concomitant relation) would make `F_g^cf` independent of the representative; it is not needed and not proved here.
4. The degree-0 term of `F_g` is a constant and does not affect the statement.
5. [CGYZ26, Rem. 1.9] names the conifold gap of the formal quintic as one of the two key statements of the approach of Guo-Janda-Ruan,
   without stating it; do not claim that Corollary B is literally that statement.
6. Compare with Conj. 1.5 only: Def. 1.8 and Thm 5 of CGYZ26 print `(-1)^g` with a constant `kappa` (as DIGEST_L_4 s.4 notes).

### 2.7 For the writer of s10 (and s01)
* `cor:mai-gap`: print 2.2 (the definition of `F_g^NN` and of the conifold expansion; it may go to `def:set-confexp`) and 2.3-2.4; replace
  `\draftnote{not stated as a corollary in the sources; each step is}` by the proof. Cite [PPZ15] Def. 2.13 and Prop. 2.14 (with Prop. 2.4
  for the plain action), as the review's M2 asks; list in `sec:int-inputs`: Theorem A, (T), [PPZ15] Prop. 2.14, `lem:gau-confexp`.
* The remark on complex frames and scalings (2.5) is optional; if `rmk:mai-constant` mentions "the same normalisation of the vanishing
  period", the scaling remark is its proof (`lambda`), and it also covers CGYZ's `omega_0(0) != 0` (`mu`).
* `rmk:mai-constant`: replace `\draftnote{compare def:set-confexp with CGYZ26 Defs. 1.1-1.2, Rem. 1.4}` by the dictionary and items 1-6 of
  2.6; keep one `\draftnote{CGYZ26 cited from arXiv v2; re-check the latest version before submission}` (as in s01 item 2).
* Conditionality: Corollary B is conditional exactly as Theorem A, plus (T) and [PPZ15]; it inherits R1 (`prop:set-cyclicLR`,
  `cor:gau-rowsNN`) through the definition of `F_g^NN` and of `T_Gamma`.
