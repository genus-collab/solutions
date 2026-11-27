# GAPS_G4, part 3: item 2, `lem:mod-xfun` (DIGEST item G3, M6) -- PROVED

Conventions: GAPS_G4.md s.1. Sources searched: SYN2 M6 ("standard; linear loop equation, statement number not re-checked"), XY
conventions ("by the linear loop equation"), LL s.1, LR F4 (a consistency check, not a control), DC2 s.7 G3 (the [DW] route via the trace
and `x = oo`, unchecked). No source writes the passage. EO v4 was read in the full text (PDF, s.4.4 and Appendix A): it contains the
statement and, inside the proof of its Cor. 4.1, exactly the deduction needed (s.3 below). The proof of record below is self-contained.

## 1. Statement (proposed wording)

**Lemma (`lem:mod-xfun`).** Let `2h - 2 + n > 0`, `n >= 1`, let `p_2, .., p_n` be points off the branch points, and let `f` be holomorphic
near the five critical values `x(varsigma_b)` (e.g. `f` in `Q(v)[X]`, `v != 0` fixed). Then for EACH `b`:
  `Res_(p = varsigma_b) f(x(p)) omega_(h,n)(p, p_2, .., p_n) = 0`.
Consequently every pairing of `omega_(h,n)` in which one slot is paired with `f(x)` by the residue pairing `Sum_b Res_(varsigma_b)` (the other
slots paired in any admissible way: residues at the branch points, circles, rays) vanishes, and in that sense the leg class `LL = Q[v][p] +
Q(v)[x]` acts modulo `Q(v)[x]`.
**Precision (finding F2, for `def:mod-pair`).** The slot carrying `f(x)` must be paired by residues at the branch points. If `deg f >= 1`,
`f(x)` has a pole at `p = 0` (`x = v/p + ..`), inside the circle `|p| = R` of `def:mod-pair`, and the circle integral equals
`Res_(p=0) f(x) omega`, in general non-zero: e.g. `Res_(p=0) x^5 omega_(1,1) = -15/4` at `v = 1` (`g4_w11.log`). The circle form is valid only
for legs holomorphic at `p = 0`, i.e. in `Q[v][p]`. The chain uses the residue pairing (LD, `prop:mod-LD`, pairs by `Sum_c Res_c`), so
nothing changes in the mathematics; only the wording of `def:mod-pair` and `thm:mod-Lreg` must say "residue pairing" for `LL`.

## 2. Proof (self-contained; uses only the recursion `conv:mod-kernel` and the curve)

*Lemma T (trace; this is EO v4 Thm 4.4, (4-23), for our curve).* For `2h - 2 + n > 0`, `n >= 1`, and `p_2, .., p_n` fixed off the branch
points, `Tr omega (X) := Sum_(i=1..5) omega_(h,n)(p^i(X), p_2, .., p_n) = 0` identically, where `p^1(X), .., p^5(X)` are the five points of
`x^-1(X)` (the roots of `T^5 - 4 X T + 4 v`), and `Tr` is the push-forward of 1-forms along `x : P^1 -> P^1`.
*Proof.* (i) For `X` not a critical value and `xi` not in `x^-1(X)`: `prod_i (xi - p^i(X)) = xi^5 - 4 X xi + 4 v = 4 xi (x(xi) - X)`, hence
`Sum_i log(xi - p^i(X)) = log(4 xi) + log(x(xi) - X)` locally, and since `B(p, xi) = d_p d_xi log(p - xi)`,
  `Sum_i B(p^i(X), xi) = d_X d_xi log(x(xi) - X) = dX dx(xi)/(X - x(xi))^2`.
(ii) Choose small discs `U_b` around the branch points with `sigma_b(U_b) = U_b`. For `p_1` outside `U := union U_b`, the recursion reads
`omega_(h,n)(p_1, J) = Sum_b (2 pi i)^-1 Oint_(dU_b) K(p_1, q) R_b(q, J)`, `K(p_1, q) = (1/2) [Int_(sigma_b q)^q B(p_1, .)]/((y(q) - y(sigma_b q))
dx(q))`, where `R_b(q, J)` (the bracket `omega_(h-1,n+1)(q, sigma_b q, J) + Sum' omega omega` of DNOPS (2-28)-(2-29), the `omega_(0,2)(q, p_l)`
factors included) does not depend on `p_1`, and the path from `sigma_b q` to `q` stays in `U_b`. (iii) For `X` outside the compact set
`x(closure U)`, all `p^i(X)` lie outside `U`, and by (i)
  `Sum_i [Int_(sigma_b q)^q B(p^i(X), .)] = dX [1/(X - x(q)) - 1/(X - x(sigma_b q))] = 0`  (`x(sigma_b q) = x(q)`),
so summing (ii) over `i` (a finite sum, exchanged with the contour integrals) gives `Tr omega (X) = 0` for these `X`; `Tr omega` is a
meromorphic 1-form in `X`, hence `Tr omega = 0`. QED (EO's own proof, App. A (A-1), (A-5), is this computation.)

*Proof of the Lemma.* Fix `b`, put `X_b := x(varsigma_b) = (5/4) zeta^b varsigma^4`. The fibre `x^-1(X_b)` is `varsigma_b` (local degree 2)
and three further points `q_1, q_2, q_3` (local degree 1). None of the `q_k` is a branch point, because the critical values `X_(b')` are
pairwise distinct. So `omega_(h,n)(., J)` is holomorphic at the `q_k` (poles only at branch points: EO Thm 4.2, or directly from (ii) by
shrinking the discs), and so is `f(x) omega_(h,n)(., J)`. *Trace formula for residues:* for a 1-form `eta` meromorphic near the fibre,
`Res_(X = X_b) Tr(eta) = Sum_(q in x^-1(X_b)) Res_q eta` (the preimage of a small circle around `X_b` is a union of small loops, one around
each point of the fibre, mapped onto the circle with the local degree; `Int_circle Tr(eta) = Sum_q Int_loop eta`). With `eta = f(x)
omega_(h,n)(., J)`: `Tr(eta) = f(X) Tr omega = 0` (Lemma T), so `0 = Res_(varsigma_b) f(x) omega + Sum_k Res_(q_k) f(x) omega =
Res_(varsigma_b) f(x) omega`. QED.
*Consequences.* `omega_(h,n)` is symmetric (EO Thm 4.6; the engine differentials are `(-1)^n` times EO's), so the slot carrying `f(x)` may
be taken to be the first. In a pairing with several slots, the other slots are residues at the branch points, circles or rays; all contours
can be taken disjoint from small discs around the branch points used in slot 1 and `omega_(h,n)` has no diagonal poles (EO Thm 4.2), so the
slot-1 residue can be taken first, for fixed points of the other contours, and it vanishes. For legs `phi, phi'` in `LL` with
`phi - phi'` in `Q(v)[x]`, all residue pairings agree. QED.
*Remarks.* (1) Without distinct critical values one still gets `Sum_b Res_(varsigma_b) f(x) omega = 0` (sum the trace formula over the
critical values). (2) The [DW] route of DC2 s.7 G3 (residue theorem on `P^1`, the poles of `f(x)` at `p = 0, oo` = the fibre over `X = oo`,
and Lemma T at `X = oo`) is also correct, but gives only the sum over `b` and needs the trace formula at the ramified point `p = oo`;
the per-`b` form above is the one EO state (Cor. 4.1 (4-26), (4-27) for `f = 1, X`).

## 3. The pinned literature statement (alternative citation; read in the full text)

EO = arXiv:math-ph/0702045v4 (PDF, "7 May 2007"), s.4.4 ("Properties of correlation functions", printed pp. 23-24), proofs in App. A:
* **Theorem 4.2**: for `(k, g) != (1, 0)`, `W^(g)_(k+1)(p, p_1, .., p_k)` "has poles (in any of its variables ...) only at the branch
  points".
* **Theorem 4.4**: "for every k and g", `Sum_i W^(g)_(k+1)(p^i, p_1, .., p_k)/dx(p^i) = delta_(k,1) delta_(g,0) dx(p_1)/(x(p) - x(p_1))^2`
  (4-23), and the same in a later slot for `k >= 1` (4-24); `p^i` are all the points with `x(p^i) = x(p)`. (With EO's `W^(0)_1 = 0`, (4-7).)
* **Corollary 4.1**: for every branch point `a_i`, `Res_(p -> a_i) W^(g)_(k+1)(p, ..) = 0` (4-26) and `Res_(p -> a_i) x(p) W^(g)_(k+1)(p, ..) = 0`
  (4-27); the proof (App. A, (A-25)) states the general case: for any rational `R(x)` with no pole at `x(a_i)`, the residue at `a_i` of
  `R(x(p)) W` equals half the residue of `R(x(p))` times the symmetrised `W(p) + W(pbar)`, which vanishes "due to theorem 4.4".
* *Hypotheses and translation.* EO's setting (compact curve, `x`, `y` meromorphic, simple branch points, `B` the Bergman kernel) holds;
  EO's kernel (4-3), (4-9) is minus the engine kernel, so `W^EO_(h,n) = (-1)^n omega_(h,n)` (`conv:mod-kernel`); vanishing statements are
  insensitive to the sign. EO's one-line deduction in (A-25) uses tacitly that the other points over `x(a_i)` are regular, which holds here
  (distinct critical values). The DC2 labels "(4-23)-(4-24) = linear loop equation" and "Cor. 4.1 (4-26)-(4-29)" are confirmed on the PDF
  (the ar5iv HTML rendering numbers some of these equations differently; cite the PDF numbers).

## 4. Evidence (VERIFIED; `g4_w11.py` -> `g4_w11.log`, exact over `Q(zeta5)`)

* With the closed form `omega_(1,1) = -(5/8) p^4 (p^10 + 18 v p^5 + 6 v^2)/(p^5 - v)^4 dp` (itself checked against the recursion, G4_5 s.4):
  `Res_(varsigma_b) x^k omega_(1,1) = 0` for all `b = 0..4`, `k = 0..6` (35 residues, exact); control `Res_(varsigma_b) p^k omega_(1,1)`, `k = 1..9`:
  40 of 45 non-zero; `Res_(p=0) x^5 omega_(1,1) = -15/4` (the circle pairing of `x^5` is not zero).
* This is a consistency check at `(1,1)` in the global coordinate; the lemma rests on s.2 (LR F4 applies: no failing control exists for a
  statement that is identically true in an `x`-adapted local coordinate).
