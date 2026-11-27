# Referee report on cluster G2, part 2: item 5 (s.2.9-s.2.10), corrections (s.3)

Continuation of `REFEREE.md`. Part 3 (`REFEREE_3.md`): the referee's own checks (s.4), the dependency on cluster G1 (s.5),
refutation attempts (s.6), files (s.7).

### 2.9 Item 5: `cor:gau-colours` and `lem:gau-colmod` (GAPS_G2_5 s.9.1-s.9.3). Accepted with corrections

**The three-line computation (9.2).** Each equality was checked.
* The first is R0 for colour `i`, i.e. the definition of `D_i`.
* `QQ^(m,i) = Rot_i QQ^(m,0)` (C1).
* `P_can Rot_i(Y) = Rot_i(P_can Y)`, by (C3) and the multiplicativity of 6.4 (ii).
* R0 for colour 0.
* The scalar rule `Rot_i(Y d(z)) = Rot_i(Y) d(zeta^-i z)`.
* (C2).

The ring `R_zeta = Q(zeta)((s))[[z]]` is a power-series ring over a field, hence a domain, and every object lies in it:
* `P_can`, with `u = s^5`;
* the model rows. `varsigma = s tau^(-1/4)`, so `Fs_b(zeta^-i z/varsigma^4)` has a Laurent series in `s` at each order of `z`;
* the actual rows after the conifold expansion (`L -> 1/s`; Step 6 of 9.4).

The non-vanishing of the row-0 entry is right. Its `z^0` term is `zeta^i L/varpi_0`, the conifold image of
`eta(1, e~_i) = zeta^i L/C_0`; cluster G1's `GAPS_G1_4` s.3.1 has the same computation.

**Recommended simplification.** The same chain, read without its first equality, proves R0 for colour `i`:
* `P_can QQ^(m,i) = Rot_i(QQ^(0) D_0) = QQ^(i) D_0(zeta^-i z)`.
* So `thm:gau-R0` (iv) needs to be proved for colour 0 only. Then `D_i := D_0(zeta^-i z)` satisfies `QQ^(i) D_i = P_can QQ^(m,i)`.
* The zero-divisor argument is then only the uniqueness of `D_i`.
* This removes the need for formal WKB uniqueness in colours `i != 0` inside `thm:gau-R0`. It is still used in
  `lem:gau-colrows` Step 2.

**Exponents.** 9.1 writes R0 with UR's `e^((u_i - u_i^c)/z)`, while OUTLINE_2 s.0.4 item 12 fixes the wording "the same
`e^(u_i/z)` stripped on both sides; UR's `e^((u_i - u_i^c)/z)` is not used".
* On amplitudes the two describe the same identity `QQ^(i) D_i = P_can QQ^(m,i)`. The model column with exponential `e^(U_i/z)`
  times the constant `e^(u_i^c/z)` has exponential `e^(u_i/z)` (`lem:set-canon`).
* So this is wording only; the editor should use the design's sentence (PS-C8).

**Fields.** `[z^n] D_i = zeta^(-in) [z^n] D_0`, so `D_i` lies in `Q(zeta)[[z]]`. For `i != 0` it lies in `Q[[z]]` iff `D_0` lies in `Q[[z^5]]`.
* Before `thm:gam-D` nothing gives `D_0` in `Q[[z^5]]`. That statement is UR's purity "(H)", which UR calls equivalent to `D` being a
  series in `z^5`.
* `D_0 = Phi^can(z,1)/h(z)` (PS s.4.2 (e)), and `h = h_(w_0)` has the term `-(3/20) z`.
* So the correction to `thm:gau-R0` (iv), DB D14, PS s.4.1 and DD1 D1.8 (G2 s.11 item 3) is right and necessary.
* "Gradedness is essential" is correct: `D'_i = D_i c(z)` for `P_can c(z)`.

**`lem:gau-colmod`.**
* The scaling identity `calT^(i)_k [G(zeta^-i t)] = zeta^i [calT^(0)_k G](zeta^-i t)` was re-derived, using `t d/dt` commuting with
  dilations.
* The five relations of the T-system were checked one by one. The closing relation uses `zeta^(-5i) = 1`.
* The uniqueness pivot was re-derived. The part of `-prod_k calT^(i)_k` linear in `t`, acting on `t^n`, is
  `t^n + 4n zeta^(4i) t^(n+1)`, because the five indices `-1, -2, 2, 1, 0` sum to 0. So the pivot is `4n zeta^(4i)`, not 0 for
  `n >= 1`.
* (3) follows.
* If `def:set-M0` stays the definition of the colour-`i` model rows, then (C1) is a one-line identity. The content of the lemma is
  that these rows are WKB solutions of the model system, which `thm:gau-R0` (i) needs. Both are written. Accepted.

### 2.10 Item 5': `lem:gau-colrows` (GAPS_G2_5 s.9.4). Accepted with corrections

* **Step 1.**
  * `PF5 = z^5 (D^5 - q prod(5D + k)) - 1 = z^5 D PF4 - 1` (`lem:set-pf5`), so every `z`-dependent monomial has degree 5.
  * The conjugation `zD -> zD + Du` is right.
  * `(zeta^-i z) D + L = zeta^-i (zD + zeta^i L)` turns `PF5^(u_0)` into `PF5^(u_i)`.
  * Correct.
* **Step 2.** The claim is right but the parenthetical is not. "The substitution `z -> zeta^i z` maps colour-`i` amplitudes to
  colour-0 ones" is the INVERSE of Step 1's substitution. "Step 1 with `i` replaced by `-i`" is a different statement: it maps
  colour 0 to colour `-i`. Correction C3.
* **Step 2, alternative.** A shorter route exists. The transport operator of PF5 in colour `i` is `T(L g) = 5 zeta^(4i) D g`.
  * Hand derivation:
    * `[z^1] PF5^(u)[f] = zeta^(4i)(5 L^-1 D f + (L^-1 - L^4) f)`;
    * use `L^5 (1 - 5^5 q) = 1`, `D L = (L^6 - L)/5` and `5^4 q = (1 - L^-5)/5`;
    * insert `f = L g`.
  * VERIFIED in `ref_wkb_run_6_10.log`.
  * So `lem:gau-wkb` holds verbatim in every colour, over any field containing `zeta`: the kernel is `Q(zeta) L`, the free
    constants sit in `g_n(q = 0)`, and no logarithm occurs.
* **Step 3.**
  * `d = 1` follows from `h_(w_0)(zeta^-i z) = h_(w_i)(z)` (`lem:set-hrot`) and `h_(w_i)` being a unit.
  * The input `psi_i(q = 0, z) = h_(w_i)(z)` is proved, for every colour, in cluster G1 (`GAPS_G1_4` s.3.1:
    `Phi_i(z, 1) = h_(w_i)(z)`), given the quoted (QRR).
* **Step 4.**
  * The induction needs `kappa'_a` independent of `z` and of `i`.
  * G1's `prop:set-cyclicLR` gives exactly the forward relations `X_(a+1) = (C_1 kappa_a)^(-1) (zD + zeta^i L) X_a`, with
    `X_a := zeta^(-i) Q^(i)_a` and `kappa^LR = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)`, for every colour (`GAPS_G1_2` s.2.0).
  * Correct.
* **Step 5.** Correct, with the conventions `Delta_i^(1/2) = zeta^(4i) C_0/L` and `QQ^(i)_a = zeta^i X^(i)_a`. These agree with G1's
  `X_a := zeta^(-i) Q^(i)_a`. The phase is `zeta^(i(a+1-n))`, as in CS s.1.
* **Step 6.** Correct. The membership of the coefficients in `NN` is G1's `cor:gau-rowsNN` (`GAPS_G1_3` s.2.10 (c)).
* **Notation.** `psi_i` is the sources' letter. With G1's recommended definition `Phi_i := C_0 Q^(i)_0/(zeta^i L)`, `psi_i = L Phi_i`.
  * Steps 1-3 then prove the clause `Phi_i(zhat, u) = Phi_0(zeta^-i zhat, u)` of `def:set-rows` (OUTLINE_3 s02 item 5). G1 s.3.1
    lists that clause as a "consistency" remark without proof.
  * `lem:gau-colrows` is exactly that missing proof. `def:set-rows` should cite it for both covariance clauses.
* **Evidence.** `ref_wkb_run_6_10.log` (s.4.3) gives Steps 1-3 for `i = 1..4`, `z^0..z^6`, `q^0..q^8`, exact over `Q(zeta)`, with
  the constants computed from the weights `w_i`. Two failing controls:
  * the wrong normalisation `h_(w_0)`;
  * a perturbation `+ q z^2` of PF5, which breaks the invariance of Step 1.
* **Status** (s.5): PROVED given cluster G1's items 1, 2, 3(a) and `lem:gau-wkb` stated over `Q(zeta)`.

## 3. Corrections

**Required.**
* **C1 (status of item 5; G2 s.0 table row 5, s.0 text, s.9.4 Status, s.11 item 2, s.12).**
  * "(i) is the open item R1 of cluster G1" is out of date. Cluster G1 reports R1 as PROVED (`GAPS_G1.md` s.0, `GAPS_G1_2`-`_3`).
    It is unrefereed at the time of this report.
  * "Item 5 inherits one dependency" understates the inputs. List all of them with their status:
    * `prop:set-cyclicLR` (G1 item 2);
    * `D u_i = zeta^i L` (G1 s.2.10 (a));
    * (M) (G1 item 1, CITED: Givental alg-geom/9603021, Cor. 11.7);
    * the colour-`i` value at `q = 0` (G1 item 3(a), PROVED given (QRR));
    * `cor:gau-rowsNN` (G1 s.2.10 (c));
    * `lem:gau-wkb` over `Q(zeta)` (s03).
  * The R7 conditionality list should keep `cor:gau-colours` until G1's referee has accepted those items.
* **C2 (`lem:gau-wkb`).** The outline states it over `Q`, colour 0, with `d` in `1 + zQ[[z]]`. `lem:gau-colrows` needs it over `Q(zeta)`.
  State it for every colour over any field `k` containing `Q(zeta)` (proof unchanged: s.2.10, transport operator).
* **C3 (s.9.4 Step 2).** Replace "(Step 1 with `i` replaced by `-i`)" by "(the inverse of the substitution of Step 1)". Or use the
  shorter route of s.2.10.
* **C4 (`thm:gau-R0` (iv)).** Adopt G2 s.11 item 3 (`D_i` in `Q(zeta)[[z]]`). Also adopt the simplification of s.2.9: R0 for
  colours `i != 0` is a consequence of R0 for colour 0, (C1)-(C3).
* **C5 (notation clashes not registered in G2 s.11 item 10).** Checked against `NOTATION.md` and `NOTATION_2.md`:
  * `gamma`. NOTATION table 3 has `gamma` = the number `7/10` (`lem:set-frobbasis`). G2 uses `gamma` as the free constant of
    `lem:gau-scal` (b)-(d), `lem:gau-lead`, 5.1, 5.3, s.7 and s.8. Rename it, e.g. to sans `c_n` (it IS the coefficient `c_n` of
    the freedom) or to `mu`.
  * `beta := u Bprod` (s.2, s.3 (f), (g), `lem:gau-scal`). It clashes with the frame parameter `beta` of `(alpha, beta)` (NOTATION
    table 3). Rename it, e.g. `Bprod = b_U/u` with a local letter.
  * `k`. G2 uses it both for the field and for the `z`-order, in the same formula ("`P = Sum_k z^k P_k`, `P_k` in `gl_5(k[[u]])`").
    NOTATION reserves `k, m, n` for orders. Use `F` (as G2 s.3 already does) or `\Bbbk` for the field.
  * `S := diag(zeta^-a)` (s.6.4). It clashes with Givental's `S ~ R e^(u/z)` (`def:set-rows`; NOTATION_2 s.10 item 1). G2 itself
    recommends keeping at most one sentence on it; write `zeta^(-gr)` instead.
  * `V(A, A~)` (solution module). It clashes with `V = tau^(-5/4)` (s08, NOTATION table 3) and with `VV = \cV`. Choose another
    symbol, e.g. `\mathrm{Gau}(A, \tilde A)`.
  * Local letters to confine to one proof or rename:
    * `e_c` (idempotents `e_i`);
    * `h_j` (vertex genus `h`; `h_n = [z^n] h`);
    * `phi_j` (CYC's `phi_i`; legs `phi^b_k`);
    * `rho_a` and `rho^U_a`, two different uses inside G2 (`rho_3`);
    * `psi_i` (write `L Phi_i`; `psi` is the cotangent class).
* **C6 (summary of the G2 run; not the note).** The summary JSON attributes the reproduction of `pi(z)` and `p_(1,b)(0)` to
  `g2_gauge`. These values are printed by `g2_G.py` (`g2_G.log`); `g2_gauge.log` does not contain them.

**Recommended.**
* **R1.** Index shifts by `m = 1..5`, with `gl^(5) = gl^(0)` (UR's convention). Then `lem:gau-comm` (d), (e) and `lem:gau-lift`
  need no special reading for `m = 0`.
* **R2.** In 5.1 (b) say "for any `P^nat` in `V(A, A~)` with `P^nat_0 = frak d Theta_U`". This is what the proof uses. It makes
  consequence 2 (`P^nat = 1` for `(A~, A~)`) literal. As written, consequence 2 is also correct (s.2.4).
* **R3 (`rmk:gau-Uremarks` (a)).** The role of the simple pole can be stated sharply. For two systems WITHOUT a pole (all weights
  units, (EP)), the same proof gives five free constants per order: `kappa^U_j` loses `j/(5u)`, so every `j` resonates. The
  solution module then has rank 5. VERIFIED here on random systems (s.4.2), as the failing control of the uniqueness clause.
  Label it VERIFIED, with the one-line argument as SKETCH, unless the editor writes the five-constant version of `lem:gau-lead`.
* **R4 (evidence lines).**
  * The two `P_can` lines of `g2_colours.log` ("colour covariance of `P_can`" and its control) cannot fail. The phase computed in
    `g2_colours.py` is identically 0 and the control is a constant expression. Delete them or label them "not a test". The
    statement is proved by 6.4 (i).
  * `g2_G.py` prints `pi(z)` and `p_(1,b)(0)` but does not compare them with the sources. They agree with GI s.5 and UR s.5 (this
    referee's comparison):
    * frame (0,0): `(127/3150) z^3 + (1416472237/725359556250000) z^8` and `(8/45, 0, 8/45, -1/15, -1/15)`;
    * frame `(1/3,-2/7)`: `-(239/450) z^3` and `(-7/45, 0, -7/45, -1/15, -1/15)`. The latter matches UR's `Xhat(2,4) = 7/225`.
  * Write "printed; agrees with" in s.10.
  * `g2_lemma2.py`'s "trace form" uses the shift decomposition, so it is only partly independent of the closed formula. The fully
    independent linear-algebra check is `ref_proj_run.log`.
* **R5.** The exponent sentence of 9.1: use the design's wording (s.2.9).
