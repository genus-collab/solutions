# Referee report `cyc_formal`, part 2: C4, C5, the Theorem, the Corollary, the frames (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict, s.1-s.5). Same labels and abbreviations. Part 3 (`REFEREE_3.md`): own code, errata,
record.

## 6. Item (6): C4, formal Picard-Lefschetz: PROVED as written

**`tau` is well defined.** `R = K[l]` is a polynomial ring, so `l -> l + 1` defines a K-algebra endomorphism `tau`. It is an
automorphism, with inverse `l -> l - 1`.

**`tau` commutes with D.** Both `D tau` and `tau D` are additive and satisfy `X(ab) = X(a) tau(b) + tau(a) X(b)`, i.e. they are
tau-twisted derivations. They agree on K, where `tau` is the identity, and on `l`: `D(l + 1) = D l = 1 - 1/u = tau(D l)`. Hence they
agree on R. (CYC's one-line reason "D l in K" is this argument.)

**`tau` preserves S and B.** `tau f_i = f_i` and `tau y_3 = g(l + 1) + h = y_3 + g`, so `tau(S) = S`; this needs no statement that S is
the whole solution space. `P(f, f')` has coefficients in K, so `tau P(f, f') = P(tau f, tau f')`. `B` takes values in Q, which `tau`
fixes. Hence `B(f, y_3) = B(f, y_3 + g)`, and `B(f, g) = 0` for every power-series solution `f`.

**VERIFIED (`cyc_ref_sol` S6).**
* `tau D = D tau` on six random elements of R (l-degree 3, exact);
* `tau y_3 = y_3 + g` exactly;
* `B(tau y_i, tau y_j) = B(y_i, y_j)`, computed directly on the transformed elements (not through the matrix);
* `B(f, g) = 0` for five random power-series solutions `f`;
* `B(y_3, g) = 1`, as CYC s.5 says.

**Controls.** The map `l -> l + u` does not commute with D, and it does not map `y_3` to a solution.

**Remark.** C4 is why the frame condition `B(varpi_0, g) = 0` is vacuous. It explains `fq_conifold`'s recorded `pairing_vacuous`.

## 7. Item (7): C5, special geometry: PROVED as written

**The adapted symplectic basis exists.** Put `Lam := span(y_0, y_1)`; it is isotropic by C4 and alternation. The map
`S -> Lam*`, `x -> B(., x)|Lam`, has kernel `Lam^perp = Lam` (B is non-degenerate and `dim S = 4`), so it is onto. This gives:
* `y_3'` with `B(y_0, y_3') = 1` and `B(y_1, y_3') = 0`;
* `y_2` with `B(y_0, y_2) = 0` and `B(y_1, y_2) = 1`;
* finally, replace `y_2` by `y_2 - B(y_2, y_3') y_0`, so that `B(y_2, y_3') = 0`.

In this basis Bmat is the signed permutation `J`, and `J^-1 = -J`. CYC says "linear algebra gives"; this is that linear algebra.
VERIFIED: the basis was constructed over Q in five frames, and its Bmat is `J` (`cyc_ref_frame` F1).

**The Leibniz identities.** I checked them by hand. On `y_j = y_0 f_j`, `T_k = y_0 Sum_m C(k,m) D^(k-m) y_0 s_m` with `s_0 = 0`. With
`e = Dt`:
* `D = e th`;
* `D^2 = (De) th + e^2 th^2`;
* `D^3 = (D^2 e) th + 3 e De th^2 + e^3 th^3`.

These hold for arbitrary series. VERIFIED on random data, where `tau_1` and `tau_2` are not zero, so the identities are tested
non-trivially (`cyc_ref_leib` L1-L2). The controls (a wrong binomial coefficient; a missing factor 3) fail.

**The conclusions.** `T_1 = T_2 = 0` and `T_3 = 1/u` (C3), and `y_0` and `e` are units of `Q[[u]]`. Hence `tau_1 = tau_2 = 0` and
`tau_3 = 1/(u y_0^2 e^3)`. Two further facts:
* The adapted-basis expansion `tau_k = -th^k f_3 + f_3 th^k f_0 - t th^k f_2 + f_2 th^k f_1` holds, and `th^k 1 = th^k t = 0` for
  `k >= 2`. VERIFIED on random data (L3).
* `tau_2 = 0` gives `th^2 f_3 = -tY`. Then `tau_3 = -th(th^2 f_3) - t th Y = (th t) Y = Y`.

VERIFIED in five frames below `u^76`: `tau_0 = tau_1 = tau_2 = 0`, `u y_0^2 Dt^3 tau_3 = 1`, `Y = th^2 f_2 = tau_3` and
`th^2 f_3 = -tY` (F1).

**Where isotropy enters (clarification, E7).** The quantities `T_k`, `s_r`, `tau_k` are basis-independent: `Sum Binv_ij y_i (x) y_j` is
the inverse bivector of B. So `tau_1 = tau_2 = 0` and `tau_3 = 1/(u y_0^2 Dt^3)` hold for every unit `y_0` and every `t` with `Dt` a
unit, isotropic or not.
* VERIFIED for the non-isotropic `t = f_1/y_0` (F1, "observation").
* Isotropy is used only in the next step: expanding `tau_k` in a symplectic basis whose first two vectors are `y_0` and `y_1 = t y_0`.
* For `t = f_1/y_0` no such basis exists, and the conclusion `th^2(y/y_0) in Y span(1, t)` fails (control, F1).

CYC's proof uses isotropy at exactly this point, so it is correct. One sentence would save the reader the search.

## 8. Item (8): the Theorem and the Corollary

### 8.1 Theorem CYC: PROVED as written
* **Operator identities.** `D (1/Dt) D = Dt th^2` and `M = (Dt/y_0) th^2 o (u y_0^2 Dt^3) o th^2 o (1/y_0)`, checked by hand.
  VERIFIED as operator identities for random `y_0`, `t` (`cyc_ref_leib` L4). The control, with `u y_0^2 Dt^2` in the th-form, fails.
* **M kills S.** For `y = Sum c_i y_i` in the adapted basis, C5 gives `u y_0^2 Dt^3 th^2(y/y_0) = c_2 - c_3 t`, and `th^2` kills it.
  VERIFIED on `y_2` and `y_3'` (value 1 and `-t`) in five frames (F2).
* **Leading coefficient.** It is the product of the factors: `(1/y_0)(1/Dt)(u y_0^2 Dt^2)(1/Dt)(1/y_0) = u`. VERIFIED (L4, F2).
* **Wronskian argument.** `E = M - PF4` has order at most 3. Its coefficient row `e` satisfies `e Wr = 0`, and `det Wr = -2/u^2` is
  not zero, so `e = 0`. Correct.

**VERIFIED, and stronger.**
* The full operator expansion of `M` (all coefficients of `D^0..D^4`) equals PF4 below `u^77` in five frames: the project's `(0,0)`,
  `(1/3,-2/7)`, `(7/11,-3/5)`, the rescaled frame `(-5/2, 4)` with `mu = -3`, `nu = 7/2`, and the critic's frame (F2).
* The same holds identically in `(alpha, beta)` below `u^27`: the generic frame `y_0 = f_0 + alpha g + beta f_2`, with coefficients in
  `Q[alpha, beta]` (`cyc_ref_sym` Y2). Every `u`-coefficient of `M - PF4` is the zero polynomial, while the series of `t` carries
  polynomials of degree up to 28.

**Controls (all fail as they must).**
* `t = f_1/y_0` and `t = (f_1 + f_2/2)/y_0` (non-isotropic pairs);
* the middle factor `u y_0^2 Dt^3`;
* `y_1 = f_1` in the generic frame.

**Scale invariance.** Under `y_0 -> mu y_0` and `g -> nu g`, `t` and `Dt` scale by `nu/mu`. The five factors of M then scale by
`(1/mu)(mu/nu)(nu^2)(mu/nu)(1/mu) = 1`, so "every normalisation of `g`" (CYC s.0) holds. VERIFIED with `mu = -3`, `nu = 7/2`.

### 8.2 The Corollary: PROVED as written
**PF5 = z^5 D PF4 - 1.** Since `z` is a constant, `(zD)^5 = z^5 D^5` and `5zD + kz = z(5D + k)`. The identity is therefore
equivalent to `D^5 - q o prod_(k=1..5)(5D+k) = D o PF4`. That is exact in `Q[u][D]` (`cyc_ref_alg` A2), and it follows from
`q(5D + 5) = 5Dq`.

Controls: with `q` on the right, or with `k = 0..4`, the identity fails.

**The coefficients `c_a`.**
* `fq_rows` (docstring of Rows) has `kappa = (1, L^5/(C0^2 C1^3), 1, C0/C1, C0/C1)`.
* UR s.1 then gives `c = C_1 kappa = (C1, L^5/(C0^2 C1^2), C1, C0, C0)`, with `prod c = L^5 = 1/u`.
* With `C0 -> y_0` and `C1 -> Dt` this is CYC's `c`.
* The chain `(1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D o (1/y_0) = (1/y_0) D M` holds, checked by hand.

VERIFIED:
* `O o (1/y_0) = (1/y_0) o D o PF4` below `u^76` in five frames (`cyc_ref_frame` F3), and identically in `(alpha, beta)` below `u^26`
  (Y3);
* UR/GI's `b_a = -c_a/(1-u)` and `prod b = -1/(u(1-u)^5)` (F3).

Controls, which fail as they must: using `kappa_1` in place of `c_1` (i.e. `Dt^3`), and swapping `c_2` with `c_3`.

**Both directions are used, and both are proved.**
* "Cyclic system => PF5" is used in two places. In P(b) step 5 (PB s.2.8): `Phi^can`, from the gauge, solves HPF. In the glue (PB s.3):
  `A^can_n` solves the WKB recursion of PF5.
* "PF5 => closing relation" is also used in two places. For the quintic rows (UR s.1; `fq_rows`: `X_0 = psi_i/C0`, with
  `e^(u_i/z) z psi_i` the WKB solution of PF5, and the four definitions). For `Q^SP`, whose rows are built "by the cyclic system" from
  `Phi^SP` (PB s.2.1, s.2.8).
* The Corollary's "iff" gives both directions. The closing relation needs only that `Psi` solves PF5; whether the `psi_n` are Laurent
  polynomials in L (no logs) plays no role here.

**The ring (E8).** The identities hold in `Q[[u]][D]`, so they act on every differential module over `Q[[u]]` in which `z` is a
constant. That includes:
* the s-series of UR s.1 (`D s = -(1 - s^5)/(5 s^4)`, the unique extension of D);
* after conjugation by `e^(u_i/z)` (the automorphism `D -> D + zeta^i L/z`), the rows `X_a` of `fq_rows`.

CYC does not say this. One line should be added.

## 9. Item (9): the frames: agree

**The project's conifold frames are all conifold frames in CYC's sense.**
* `fq_conifold` (lines 85-87) sets `varpi0 = f_0 + alpha g + beta f_2` and `tc = g/varpi0`. Here `g = f_1 + gamma f_2` with the fixed
  scale `g = u + (7/10)u^2 + (41/75)u^3 + ...`.
* `pb_core.set_frame` uses the same formula.
* GI s.1 and UR s.1 use "frame `(alpha, beta)`: `varpi_0 = varpi_0^cf`, `t_c = u + O(u^2)`".
* POLY uses `C0 = varpi_0` in frame `(0,0)`.
* Every one of these has `varpi_0(0) = 1` and `Dt(0) = -1`. CYC also allows every non-zero scale of `y_0` and of `g` (s.8.1).

VERIFIED:
* `cyc_ref_xproj` X3: the project's own code returns the same `varpi_0` and `t_c` in `(0,0)` and `(1/3,-2/7)`;
* `cyc_ref_frame` F3, frame `(0,0)`: `varpi_0 = 1 + O(u^3)` and `t_c = u + (7/10)u^2 + (41/75)u^3`, as in GI;
* PB s.2.8's hypotheses `t_c'(0) = varpi_0(0) = 1` hold in every such frame.

**The critic's frame labels (E3).** `cr5_cyc_run` and CYC s.8 parametrise frames as `y_0 = f_0 + a f_1 + b f_2`, not
`f_0 + alpha g + beta f_2`. Since `f_1 = g - (7/10) f_2`:
* the critic's "(1/3, -2/7)" is the project's `(alpha, beta) = (1/3, -109/210)`;
* the critic's "(7/11, -3/5)" is the project's `(7/11, -23/22)`.

`pb_cyc` uses the project's parametrisation. So the "3 frames" of pb and of the critic (STAGE5_SYNTHESIS s.1 table, s.4.1) are
different sets, apart from `(0,0)`. This is harmless, because the Theorem covers every frame and my run covers both labellings, but
the labels should be corrected.

**"Every frame" (E5).** CYC s.0's last line reads "a theorem, in every frame". The Theorem as proved covers every conifold frame. It
does not cover the large-radius frame `(I_0, I_1)` or the abstract ring N of `fq_ring`. An analogous argument at `q = 0` (the
isotropy of `(I_0, I_1)` at the point of maximal unipotent monodromy) is classical. It is not given in CYC and not checked here
(SKETCH at most), and nothing in the chain needs it.

What the chain needs is the closing relation for the conifold-frame (s-series) rows of UR s.1. Two facts identify the N-rows of
`fq_rows` with these series under `Cf`; both are properties of the project's construction, not part of CYC, and both are consistent:
* the rows are polynomials in `L`, `c_j` and `D^j I_1`, hence independent of `Lam`;
* `B(varpi_0, g) = 0` holds by C4.

The docstring's "elements that are Lam-free in N give Lam-independent images" should read "Lam-independent when written in `D^j I_1`":
`i_0` is Lam-free in N, but its image depends on the chosen value of `Lam`.
