# Outline, part 10: s02 set-up, items 10-17 (continuation of `OUTLINE_3.md`)

Same conventions as `OUTLINE_3.md`. Subsections: `sec:set-frames` items 10-12; `sec:set-model` items 13-15; `sec:set-hat` item 16;
`sec:set-N` item 17.

10. `def:set-frames` (defn) with `lem:set-frobbasis` (lem; DB S2): PF4 maps `Q[[u]]` to itself, local exponents `0, 1, 1, 2`, solutions
    `f_k = u^k + O(u^3)`; the log solution `varpi_1 log u + rho_3`, `rho_3` in `u^3 Q[[u]]`; obstruction matrix `[[-2,0,0],[17/5,0,0],
    [-7/5,-7/5,2]]`, kernel `Q(0, 1, 7/10)`, so `varpi_1 = f_1 + (7/10) f_2` up to scale. Proof DB s.3 S2 steps 1-6 [W] (CF s.2.2-s.2.3; CA
    s.3); RS2 item 1 (the 3x3 matrix replaces "computed"). The conifold frames: `varpi_0 = f_0 + alpha varpi_1 + beta f_2`, `(alpha, beta)` in
    `Q^2` (DEC D14), `t_c = varpi_1/varpi_0`; then `varpi_0(0) = 1`, `t_c'(0) = 1`, `D t_c(0) = -1`, `[u^1] varpi_0 = alpha`, `[u^2] t_c +
    alpha = 7/10` (DB D4; DA D11); scale of `varpi_1`: `u + O(u^2)` (DA C17). In these frames the normalisation (N0) of `lem:gau-diag` holds
    (G2 s.2, last paragraph). Printed name: "the conifold frames of `def:set-frames`" (never "programme family", REV_2 M9).
11. `def:set-confexp` (defn; REV R1 item 3, R8 (h); G5_2 s.2.2 [G] with G5R F2.1, F2.5). The differential ring `\Nring := Q[L^(+-1),
    \logq, c_0, .., c_3, i_0, i_1, i_2][c_0^-1, \mathsf W^-1]`, `\mathsf W := c_0^2 + c_0 i_1 - c_1 i_0`, with the derivation `D` fixed on the
    generators: `D c_j = c_(j+1)`, `D c_3` by PF4, `D i_j = i_(j+1)`, `D i_2 := i_3` solved from the concomitant relation (linear in `i_3`
    with a unit coefficient), `D \logq = 1`, `D L = (L^6 - L)/5`. Evaluation `ev : \Nring -> ` functions of `q` (`c_j -> D^j C_0`, `i_j ->
    D^j I~_1`, `\logq -> log q`), a differential homomorphism (`ev(\mathsf W) = C_0^2 C_1`). The conifold expansion `\Cexp_(alpha,beta) :
    \Nring -> Q((s))`: `c_j -> D^j varpi_0`, `i_j -> D^j varpi_1 - j D^(j-1) varpi_0 + \Xi D^j varpi_0`, `\logq -> log(1 - u) = -\Xi`, `L ->
    s^-1`; `c_0` and `\mathsf W` go to units; `C_0 -> varpi_0`, `I_1 -> varpi_1`, `T -> t_c`, `C_1 -> D t_c`. No relation beyond PF4 and
    the concomitant is imposed. A change of the primitive of `1` composes `\Cexp` with the automorphism `i_j -> i_j - k c_j`, `\logq ->
    \logq + k` of `\Nring`, which fixes `C_0`, `I_1`, `\mathsf W`, `C_1`, `Delta_i^(1/2)`, the rows and `F_g^\Nring` (G5R s.2.3 F2.1; never
    "any primitive gives the same value on `\logq`-free elements"). That `\Cexp` is a ring map commuting with `D` is `lem:gau-confexp`
    (s03); that the rows lie in `\Nring (x) Q(zeta)` is `cor:gau-rowsNN` (s03), not FQ2 (FQ2 s.4 goes to appendix C). "Regular at the
    conifold point": image in `Q[[u]]` (DA D11, D18). Never the code name. The comparison with [CGYZ26] Defs. 1.1-1.2 is `rmk:mai-constant`.
12. `def:set-tau` (defn): `tau(u) = 1/2F1(1, 4/5; 9/5; u) = 1/Sum_n 4u^n/(5n+4)`, `v(u) = u tau^(-5/4)`, `varsigma = v^(1/5) = s tau^(-1/4)`
    (DA D14; DB D7). `lem:set-canon` (lem): `u_i - u_i^c = -(5/4) zeta^i u^(4/5)/tau(u)`, so `tau zhat = -(5/4) zeta^i z/(u_i - u_i^c)` (DB D6;
    UR s.1; RN s.3(a): termwise `Int_0^u t^(-1/5)(1-t)^(-1) dt` [W*]); orientation remark NFS s.2 item 7 (immaterial). Uses `D u_i = zeta^i L`
    (`cor:set-cyclicPF5` (a)); the additive constants are fixed here; the colour argument needs only `u_i - u_i^c = zeta^i (u_0 - u_0^c)`
    (G2_6 s.11 item 4).
13. `def:set-model` (defn): `W_v(p) = p^4/4 + v/p`; formal integrals `I_m(v,z) = Int p^m e^(-W_v(p)/z) dp`; rows `Y = (I_0, -I_(-1), -I_3,
    I_2, I_1)`, i.e. `Y_b = eps_b I_(nu_b)`; dictionary `H^b <-> eps_b p^(nu_b) dp`; homogeneity `p -> kp, v -> k^5 v, z -> k^4 z` (DB D8; the
    sources' `lambda = -(x^4/4 + v/x)` is `-W_v`; DA C9 only in appendix B, DEC D13).
    * `prop:set-modrows` (prop; DB S18 (a)-(b)): (a) `z dY_b/dv = bmod_b Y_(b+1)`, `bmod = (1, -1/v, 1, -1, -1)`, from (Mod1) `z dI_m/dv =
      -I_(m-1)` and (Mod2) `m z I_(m-1) - I_(m+3) + v I_(m-2) = 0`; (b) with `Y_b = e^(U/z) varsigma^(nu_b - 1) Fs_b(t)`, `t = z/varsigma^4`,
      `U = -(5/4) zeta^i varsigma^4`, the system is the T-system `Fs_1 = \mathcal T_(-1) Fs_0, .., Fs_0 = -\mathcal T_0 Fs_4`, `\mathcal T_k =
      -zeta^i + (t/5)(k - 4t d/dt)`, with a unique solution `Fs_0(0) = 1` (pivot `4n zeta^(4i)`, G2R_2 s.2.9); `Fs_b(0) = eps_b`; homogeneous
      WKB columns with the same exponential differ by a constant, homogeneous only if constant ("uniqueness half", no hypothesis on the
      quintic). Proof DB s.3 S18 steps 1-4 [W] (UR s.1-s.2; PS s.1 (ii)-(v)). `Fs_0 = 1 - (3/20)t + (9/800)t^2 + ..` only as an illustration.
14. `def:set-M0` (defn): `Rhat^mod_b^(i)(t) = zeta^(i(b+1)) Fs_b(zeta^-i t)`, `Psi_0 = Rhat^mod(0)` (entries `eps_b zeta^(i(b+1))`,
    `Psi_0^T eta^-1 Psi_0 = 1`), `M_0(t) := Psi_0^-1 Rhat^mod(t)`, the model R-matrix in Teleman's colour-colour form (DB D15, D18; DD1 D1.5;
    LG s.0 item 4).
    * `lem:set-M0sympl` (lem; DEC Q3; one label for G1_5 s.5.1 = G4_2 s.5): `M_0(-t)^T M_0(t) = 1`, hence `M_0(t) M_0(-t)^T = 1`, i.e. `M_0`
      lies in Givental's group. Proof G4_2 s.5 [G] (the pairing `Y(-z)^T eta^-1 Y(z)` of a solution and an opposite solution of the
      self-dual model system is constant, `bmod_(a-1) = bmod_(3-a)`; then homogeneity), with G4R_2 s.5 R-5 (`eta(H^a, H^b) = 5 delta_(a+b =
      3 mod 5)`; indices of `bmod`, `eps`, `nu` mod 5) and local letters for `P_ij`, `B`, `S_ij` (G1R_2 s.6). The earlier `\draftnote` goes.
15. `def:set-modcohft` (defn; GS s.2.10): `Omega^mod_v := (Psi^m_v)_*(M_0(z/varsigma^4).\TFT^m_v)`, `Psi^m_v := Q^m(v; z = 0)`; the model TFT
    `\TFT^m_v` is DEFINED by `Deltahat_i^(1/2) := zeta^(4i) varsigma` in the normalised canonical basis, unit `Sum_i Deltahat_i^(-1/2) e~_i`
    (the identification with a Hurwitz Frobenius manifold is not used, G4_2 s.7; its name only in appendix B); the action is [PPZ15] Def.
    2.13 (with Def. 2.5). Ancestor correlators `< .. >_h(v)`; dilaton shift, base point `-z 1`; `F^m_h = c_h v^(2-2h)` is a theorem (s06).
    * `lem:set-modhom` (lem): every ancestor correlator `<prod_i psi^(k_i) H^(a_i)>_h(v)` is `C v^e`, `C` in `Q`, `e = 2 - 2h + Sum wt(k_i, a_i)`,
      `wt(k,a) = (4k + nu_a - 4)/5`; regular iff `e >= 0` (DC1 D3, D4; DC2 D6; S15 (I1)). Proof G1_5 s.5.2 [G] (degree count; independence
      of the fifth root; the Z/5 selection), with the naturality sentence and the unit check of G1R_2 s.7; local letters for G1R's `Lam`,
      `Pi`. The earlier `\draftnote` goes.
    * `lem:set-dilaton` (lem): statement G1_4 s.4 (i) and (iii) with G1R_2 s.5 C4-1, C4-2: (i) for a CohFT `Omega` with flat unit ([PPZ15]
      Def. 0.3 (i)-(iii)), `<psi 1 . prod_j psi^(k_j) a_j>_h = (2h - 2 + n) <prod_j psi^(k_j) a_j>_h`; (ii) `Omega^mod_v` has the flat unit
      `H^0` (G1R_2 s.7 (U)); (iii) for a complete local `Q`-algebra with maximal ideal `m` (e.g. `Q[[u]]`; a local letter, not `\Lambda`) and
      `\varkappa` in `1 + m`, with `Fbar_h(q) := Sum_n <(q + z 1)^n>_h/n!` a formal function of `y = q + z 1` (`y` formal: a hypothesis of the
      statement), `d^n Fbar_h(q/\varkappa) = \varkappa^(2h-2+n) d^n Fbar_h(q)` for `2h - 2 + n > 0`, and `Fbar_h(\varkappa q) = \varkappa^(2-2h)
      Fbar_h(q)` for `n = 0`, `h != 1`. Proof G1_4 s.4 [G] with [PPZ15] Def. 0.3 (iii), Prop. 2.12, Def. 2.13 and [AC96] (1.7), `kappa_0 = 2g -
      2 + n` (p.5). One sentence: as identities of Taylor series at `-z 1` they need `\varkappa(0) = 1` (G3R_2 s.13 K4). The genus-1
      logarithmic term `-(5/24) log \varkappa` is not printed (not used; it needs `Int_(Mbar_(1,1)) psi_1 = 1/24`, not pinned). Wherever the
      lemma is applied (`lem:tra-F`, `lem:tra-conv` (f), `lem:flh-dil`), state `d(0) = 1` (`lem:gau-diag`; GS s.2.3; DEC Q5).
16. `def:set-hat` (defn): `Hhat^a = s^(1-nu_a) H^a`, `z = s^4 zhat`, `Qhat_a = s^(1-nu_a) Q_a`; `xhat = varpi_0` (`Deltahat_i^(1/2) = zeta^(4i)
    xhat`), unit `1hat = s 1`, rescaling `Omega_(g,n) -> s^-(2g-2+n) Omega_(g,n)`; the hatted TFT `\hat\TFT`; `eta` preserved (`nu_a + nu_b = 2`
    on eta-pairs). DA D12; OU2 s.1; DC2 D14. The `u`-series property of `Qhat` is NOT assumed (`cor:mai-natfact`).
    * `def:set-outer` (defn; REV R6): `M := M_0(tau(u) zhat)`; `R' := Q M^-1`, `Q` the matrix of the lowered rows of `def:set-rows`
      expanded in the frame (`def:set-confexp`; the rows lie in `\Nring (x) Q(zeta)`, `cor:gau-rowsNN`), as formal series in `zhat` with
      coefficients Laurent series in `s`; the outer graph sum `R'.(M.\TFT)` (vertices `(M.\TFT)_(g,n)`, edges and leaves of `R'` by [PPZ15];
      the constant part `R'_0` is not unital, `def:flh-outer`; leaf `T' = z y + q_(R')`, FR C-list); the outer term `T_Gamma :=
      Cont_Gamma[E_(R'), q_(R'); M.\TFT]` of a leg-free stable graph (the full vertex tensors contain the colour sum, DA D13, C19); the vertex
      graph = the one-vertex graph with its `R'`-leaf decorations (DA D13); the hatted term `G_Gamma` := the same contraction with the
      hatted data of `def:set-hat` (rows `Qhat` at `zhat`, `\hat\TFT`, `1hat`); `lem:flh-hat` (s07) proves `G_Gamma = u^(2g-2) T_Gamma`; `T_Gamma`
      is regular iff it lies in `Q[[u]]` (DA D18, A8; DC2 D16). `M_0` is the model's (`def:set-M0`; OUTLINE_2 s.0.4 item 18).
    * `rmk:set-reversion` (rmk; DA A8): the `t_c`-polar part from `[u^k] G_Gamma` by reversion of `t_c = u(1 + ..)`, triangular with unit
      diagonal [W, NF s.6].
17. `def:set-N` (defn) "Property N(g)": the statement of `thm:mai-main` for one `g` (per outer graph; every conifold frame). The values of
    `c_g` (`1/1200, 1/25200, 1/180000, 1/660000`) and the equality of its four forms: `rmk:set-cg` (DA A16; S01 C2; digest arithmetic
    with `B_4 .. B_10`); no computation is cited.

Corrections to apply in s02: DA C4/C5 (one letter for the Teleman form `R~`; the I-function-normalised matrix is not used); DA s.7
item 7 (`Psi(u)` only in `cor:mai-natfact`; `R' = Q M^-1` suffices here); DA s.7 item 8 (`xhat = varpi_0 = 1 + alpha u + ..`, not
`1 + O(u^3)`); G1R s.2.3 and G1R_2 s.5-s.7 (local letters). Length of s02: 10 pp.
