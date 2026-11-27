# REF_1: referee report on `lem:set-Runique` (flatness and the value at q = 0 determine R)

Date: 2026-11-18. Scope: `sections/s02_setup_4.tex` (l.64-75 the lemma and proof; l.7-21 `def:set-rows`; l.29-48
`thm:set-qrr`), `sections/s01_intro_3.tex` items (QRR), (T) (l.10, l.13, l.21), `sections/s02_setup_5.tex` l.40-47
(`thm:set-GT`). Read-only; no source file was changed.

## Verdict: SOUND AFTER FIXES (wording only)

The statement is true and the proof is complete in substance: every step was re-derived by hand, and the recursion
and the symplecticity argument were also run exactly in a generic truncated model. No error and no mathematical gap
was found. Five wording findings follow (W1-W5); W1-W3 make explicit what the proof uses silently, W4-W5 concern the
sentences that use the lemma. Nothing elsewhere claims more than the lemma gives.

## What was checked, step by step

1. **Form of the flatness equation.** With `w = R e~_i = sum_j Y_ji e~_j`, `partial_T = C_1^{-1} D_q` acting on
   components in the constant flat basis, and `C_1 H* e~_j = zeta^j L e~_j`: the `e~_j`-component of
   `C_1 (z partial_T w + mu_i w - H*w) = 0` is `z (D Y)_ji + z (Theta Y)_ji = (zeta^j - zeta^i) L Y_ji`, which is
   `z D Y + z Theta Y = [DU, Y]`, `DU = diag(zeta^j L)`. Coefficient of `z^0`: `0 = [DU, Y_0]`, true for `Y_0 = 1`;
   coefficient of `z^(k+1)`: the displayed recursion. Correct, index placement included.
2. **Theta.** `Theta_ji = eta(e~_j, D e~_i)`; differentiating `eta(e~_i, e~_j) = delta_ij` (the Gram matrix of `eta` in
   the flat basis is constant, `def:set-fq`) gives `Theta_ji + Theta_ij = 0`. The flat components of `e~_i` lie in
   `Q(zeta)[[q]]` (see 7), so those of `D e~_i` lie in `q Q(zeta)[[q]]`, hence so does `Theta_ji`. Correct.
3. **Recursion.** Off-diagonal entries of `Y_(k+1)` from all of `Y_k`; then the diagonal entries of `Y_(k+1)` from its
   off-diagonal entries. The induction is well founded. Off-diagonal entries lie in `q Q(zeta)[[q]]` since both
   `D Y_k` and `Theta Y_k` do and `((zeta^j - zeta^i) L)^(-1)` is in `Q(zeta)[[q]]`. Correct.
4. **Existence.** It IS proved, though tersely. The diagonal equation at level `k+1` is
   `D (Y_(k+1))_ii = - sum_(l != i) Theta_il (Y_(k+1))_li`; the right side lies in `q Q(zeta)[[q]]` (factor `Theta`), and
   `D_q : Q(zeta)[[q]] -> q Q(zeta)[[q]]` is surjective with kernel `Q(zeta)` (characteristic 0). So it is solvable at
   every order, for every choice of constant term; the level-0 diagonal equation is `Theta_ii = 0`; the `z^0` equation
   holds. Hence solutions <-> families of diagonal constants, i.e. `Y -> Y diag(c_i(z))`, `c_i` in `1 + z Q(zeta)[[z]]`,
   which is the expected freedom (a constant diagonal matrix commutes with `DU`). No obstruction at any order.
5. **The condition at q = 0.** The proof does not assume that `R(q=0)` is diagonal: it proves it (off-diagonal entries
   lie in `q Q(zeta)[[q]]`). So the off-diagonal half of the quoted value in `thm:set-qrr` is automatic, and the
   condition only fixes the constants `(Y_k)_ii(0) = [z^k] h_{w_i}`; `h_{w_i}` is in `1 + z Q(zeta)[[z]]`, consistent
   with `Y_0 = 1`. `R(0,z) e~_i(0) = sum_j Y_ji(0) e~_j(0)`, so the condition of the lemma is exactly
   `Y(q=0) = diag(h_{w_i})`.
6. **Symplecticity.** Transposing the equation at `-z` with `Theta^T = -Theta`:
   `z D Y(-z)^T = z Y(-z)^T Theta - Y(-z)^T DU + DU Y(-z)^T`; adding `Y(-z)^T (z D Y)` the `Theta` terms and the
   `Y(-z)^T DU Y` terms cancel, giving `z D N = [DU, N]`. Off-diagonal: if `N_ji != 0` with lowest coefficient at
   `z^m`, the coefficient of `z^m` gives `(zeta^j - zeta^i) L [z^m]N_ji = 0`; so `N_ji = 0`. Diagonal: `D N_ii = 0`, so
   `N_ii` in `Q(zeta)[[z]]` equals its value at `q = 0`, which is `Y_ii(0,-z) Y_ii(0,z)` because `Y(q=0)` is diagonal,
   `= h(-z) h(z) = 1` as the exponent in `eq:set-hw` has only odd powers of `z`. `Y(-z)^T Y(z) = 1` is (a): with `P` the
   matrix of the `e~_i` in the flat basis, `P^T eta P = 1` and `R = P Y P^(-1)`. Correct. (Only `h(-z)h(z) = 1` is used
   here, not the form of `h`.)
7. **Unstated assumptions.** (i) `e~_i` in `H (x) Q(zeta)[[q]]`: by `cor:set-cyclicPF5` (a), `e_i` has coefficients
   `(1/5) prod kappa^LR_b mu_i^(-a)`, with `kappa^LR` in `Q[[q]]` and `mu_i^(-1) = zeta^(-i) C_1 / L`; and
   `Delta_i^(1/2) = zeta^(4i) C_0 / L` (`def:set-frob`) needs no root outside `Q(zeta)[[q]]`:
   `(zeta^(4i))^2 = zeta^(8i) = zeta^(3i)`. Re-derived `Delta_i = C_0^2 zeta^(3i) / L^2` from
   `eta(e_i,e_i) = eta(1,e_i) = kappa_0 kappa_1 kappa_2 mu_i^(-3)`. (ii) The `e~_i` are a basis over `Q(zeta)[[q]]`:
   `P^T eta P = 1` makes `P` invertible. (iii) `L = (1 - 5^5 q)^(-1/5)` is in `1 + q Q[[q]]`, a unit. (iv)
   `zeta^j - zeta^i != 0` in the field `Q(zeta)`. All hold; (i)-(iii) are not said in the proof (W1, W2).
8. **Citations and circularity.** The proof uses only: `cor:set-cyclicPF5` (a) (assumes (M), (G0); its proof uses
   `prop:set-cyclicLR` (i), proved in `a1_cyclicLR.tex`, which contains no reference to `R`, the rows, or (QRR));
   `def:set-frob` (orthonormality); `def:set-genus0` (`D_q`, `C_1`, `L`); `eq:set-hw`. It does not use
   `prop:set-cyclicLR` (iii), `cor:set-cyclicPF5` (b), `eq:set-qrr`, `lem:set-hrot`, (T) or (QRR). No circularity. The
   forward reference to the lemma inside `def:set-rows` is harmless (the lemma is about an arbitrary series with (b)).
9. **Exact check.** Script `chk.py` in the session scratchpad folder `ref1` (Python, fractions): generic antisymmetric
   `Theta` in `q Q[[q]]`, three colours, generic unit `L`, generic odd-exponent `h_i`; recursion solved through `z^5`,
   `q^6`; all diagonal equations hold and `Y(-z)^T Y(z) = 1` exactly (0 violations). This checks the algebra of the
   proof, not the formal-quintic data.
10. **Use elsewhere.** The lemma is cited once (`def:set-rows`, l.15: "(b) and the value at q = 0 determine R, and (a)
    follows from them"), which is exactly its content. (QRR), (T), `thm:set-qrr`, `thm:set-GT` do not rely on it and do
    not overstate it. Symplecticity is now both quoted in (T) and proved from (QRR) (W4): redundant, not wrong.

## Findings

### W1 (wording) `s02_setup_4.tex` l.70, first two sentences of the proof: basis not stated

Old:
```
Write $R(z)\tilde e_i=\sum_jY_{ji}(z)\tilde e_j$, $Y=1+\sum_{k\geq1}Y_kz^k$ with entries in $\bbQ(\zeta)\fps q$, and $\Dq\tilde e_i=\sum_j\Theta_{ji}\tilde e_j$ ($\Theta$ a local letter). As $\eta(\tilde e_i,\tilde e_j)=\delta_{ij}$ is constant, $\Theta$ is antisymmetric, and its entries lie in $q\bbQ(\zeta)\fps q$, because $\Dq$ maps $\bbQ(\zeta)\fps q$ into $q\bbQ(\zeta)\fps q$.
```
New:
```
The $\tilde e_i$ have flat components in $\bbQ(\zeta)\fps q$ (\cref{cor:set-cyclicPF5}~(a), \cref{def:set-frob}: $\kappa^{\mathrm{LR}}_b$, $C_1/L$ and $C_0/L$ lie in $\bbQ\fps q$) and are orthonormal, so they form a basis of $\bbH\otimes\bbQ(\zeta)\fps q$ over $\bbQ(\zeta)\fps q$, and $x=\sum_j\eta(\tilde e_j,x)\tilde e_j$. Write $R(z)\tilde e_i=\sum_jY_{ji}(z)\tilde e_j$, $Y=1+\sum_{k\geq1}Y_kz^k$ with entries in $\bbQ(\zeta)\fps q$, and $\Dq\tilde e_i=\sum_j\Theta_{ji}\tilde e_j$, $\Theta_{ji}=\eta(\tilde e_j,\Dq\tilde e_i)$ ($\Theta$ a local letter; $\Dq$ acts on the flat components). As $\eta(\tilde e_i,\tilde e_j)=\delta_{ij}$ and the pairing of the flat basis are constant, $\Theta$ is antisymmetric, and its entries lie in $q\bbQ(\zeta)\fps q$, because $\Dq$ maps $\bbQ(\zeta)\fps q$ onto $q\bbQ(\zeta)\fps q$, with kernel $\bbQ(\zeta)$.
```

### W2 (wording) `s02_setup_4.tex` l.74: "is a unit" and the existence half of "exactly one"

Old:
```
For $j\neq i$ this determines the off-diagonal entries of $Y_{k+1}$ from $Y_k$, as $(\zeta^j-\zeta^i)L$ is a unit; they lie in $q\bbQ(\zeta)\fps q$, because the left side does. For $j=i$, with $k+1$ in place of $k$, it reads $\Dq(Y_{k+1})_{ii}=-\sum_{l\neq i}\Theta_{il}(Y_{k+1})_{li}\in q\bbQ(\zeta)\fps q$, which determines $(Y_{k+1})_{ii}$ up to its constant term; the equation with $j=i$ and $k=0$ holds, as $\Theta_{ii}=0$. So the solutions $Y$ correspond to the families of constant terms of the diagonal entries, $Y(q=0)$ is diagonal, and the condition at $q=0$ fixes it as $\diag(h_{\bw_i})$.
```
New:
```
(The coefficient of $z^0$ is $[\Dq U,Y_0]=0$.) For $j\neq i$ this determines the off-diagonal entries of $Y_{k+1}$ from $Y_k$, as $(\zeta^j-\zeta^i)L$ is a unit of $\bbQ(\zeta)\fps q$ ($L\in1+q\bbQ\fps q$); they lie in $q\bbQ(\zeta)\fps q$, because the left side does. For $j=i$, with $k+1$ in place of $k$, it reads $\Dq(Y_{k+1})_{ii}=-\sum_{l\neq i}\Theta_{il}(Y_{k+1})_{li}$. The right side lies in $q\bbQ(\zeta)\fps q$, so this equation has a solution in $\bbQ(\zeta)\fps q$ for every prescribed constant term, and exactly one; the equation with $j=i$ and $k=0$ holds, as $\Theta_{ii}=0$. So, by induction on $k$, for every family of constants $c_{k,i}\in\bbQ(\zeta)$ ($k\geq1$) there is exactly one solution $Y$ with $(Y_k)_{ii}(q=0)=c_{k,i}$, and every solution has $Y(q=0)$ diagonal. The condition at $q=0$ is $Y(q=0)=\diag(h_{\bw_i})$, that is $c_{k,i}=\coef{z^k}h_{\bw_i}$; this proves existence and uniqueness, $R(z)$ being the endomorphism with the matrix $Y(z)$ in the basis $\tilde e_i$.
```

### W3 (wording) `s02_setup_4.tex` l.74, the diagonal of N: say why the value at q = 0 is a single product

Old:
```
and $\Dq N_{ii}=0$, so $N_{ii}$ is its value $h_{\bw_i}(-z)h_{\bw_i}(z)=1$ at $q=0$, the exponent of $h_{\bw_i}$ being odd in $z$ \eqref{eq:set-hw}.
```
New:
```
and $\Dq N_{ii}=0$, so $N_{ii}\in\bbQ(\zeta)\fps z$ is its value at $q=0$, which is $h_{\bw_i}(-z)h_{\bw_i}(z)=1$ because $Y(q=0)$ is diagonal and the exponent of $h_{\bw_i}$ is odd in $z$ \eqref{eq:set-hw}.
```

### W4 (wording) `def:set-rows`, l.15 last sentence; `thm:set-GT`: status of (a) and of the definite article

"The R-matrix of (T)" presupposes that (T) names one series, while (T) only asserts that `Omega = R.TFT` "for an
R-matrix R". With the lemma the definite article is justified without any uniqueness statement for the classification
over `Q(zeta)[[q]]`: (T) and (QRR) together say that the unique series of the lemma satisfies `Omega = R.TFT`. Also
(a) is now both quoted and proved; say which one the note relies on.

Old (`s02_setup_4.tex` l.15, last sentence):
```
By \cref{lem:set-Runique}, (b) and the value at $q=0$ determine $R$, and (a) follows from them.
```
New:
```
By \cref{lem:set-Runique}, which uses only (M) and (G0), (b) and the value at $q=0$ determine $R$, and (a) follows from them. So (QRR) and (T) together say that the graph sum of the unique series of \cref{lem:set-Runique} is the cohomological field theory of the formal quintic; no uniqueness of the $R$-matrix in (T) is used, and (a), although part of the quoted statement (T), is not an independent input.
```
Old (`s02_setup_5.tex` l.42):
```
symplectic as in \cref{def:set-rows}~(a) \cite[Theorem~(2) in (1.15)]{MR2917177},
```
New:
```
symplectic as in \cref{def:set-rows}~(a) (which also follows from (QRR), \cref{lem:set-Runique}) \cite[Theorem~(2) in (1.15)]{MR2917177},
```

### W5 (wording) `s01_intro_3.tex` l.13 and l.21: "R(0)=1" clashes with "R(0,z)" and "R(q=0)"

In (QRR) and `thm:set-qrr` the argument 0 of `R(0,z)` and `R(q=0)` is `q = 0`, where `R` is `diag(h)`, not 1; in (T)
and in the conditionality paragraph `R(0)=1` means `z = 0`. `def:set-rows` and `def:set-frob` write `R|_{z=0}=1`.

Old (l.13): `symplectic and with $R(0)=1$ \cite[Theorem~(2) in (1.15)]{MR2917177}`

New: `symplectic and with $R|_{z=0}=1$ \cite[Theorem~(2) in (1.15)]{MR2917177}`

Old (l.21): `only through $R(0)=1$, its coefficient ring`

New: `only through $R|_{z=0}=1$, its coefficient ring`

## Not findings (checked, no change needed)

* The lemma sits after `lem:set-hrot` but does not use it; the placement only makes `eq:set-hw` a backward reference.
* `thm:set-qrr` states the value "in the basis e~_i(q=0)" including diagonality; by step 5 diagonality is automatic,
  so the quoted statement is slightly redundant, not inconsistent.
* The hypothesis "Assume (M) and (G0)" of the lemma is exactly what `cor:set-cyclicPF5` (a) needs; neither (QRR) nor
  (T) is assumed, as it should be.
* The intro sentences "(T) identifies that R-matrix" and "`thm:int-D` uses the R-matrix that (T) identifies only
  through ... and (QRR)" remain accurate; with the lemma they could be strengthened (the rows are determined by (QRR)
  alone), but they do not overclaim.
* Not checked: the cited literature for (QRR) and (T) themselves (outside the scope of this piece).
