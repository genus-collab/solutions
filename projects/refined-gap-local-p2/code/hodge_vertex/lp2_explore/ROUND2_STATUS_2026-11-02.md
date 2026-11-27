<!-- workflow wf_1ab43e48-d73 (3 provers, 6 referees, synthesis); full reports in its journal.jsonl -->
# Round 2 status, 2026-11-02: H1 and the T_F continuity lemma

**Summary.** H1 is settled, in the form of a definition followed by theorems. It is not settled as a claim about the code. The T_F continuity lemma is proved, modulo two things: standard facts about formal Laplace expansions (F1–F3), and BS genus-0 inputs the user has allowed (the Givental mirror theorem / eq:asymp, (Normalisation), eq:canR). Six referee reports were filed, two for each of the three proofs. All six say "sound with minor fixes", and none raises a major objection. TR is not used anywhere.

## 1. H1

- **Definition.** The rows used by the theorem are R^{T_F}_{0i} := lim_{α→α*} R^α_{0i}, where R^α_{0i} are BS's canonical rows (eq:canR) on generic or rational tori.
  - Frame B multiplies by h_w(z)^{-1}.
  - Rows 1 and 2 are given by the buildRmatrix formulas at c1 = C1cf.
- **Theorems (all orders).**
  - The limit exists and equals 𝒜/N, the normalised saddle expansion of BS's oscillatory integral.
  - It solves eq:pf at every colour, and it is unique up to a factor c(z).
  - It has the closed form Δ_i^{-1/2} h_w(z) 𝔖(v_i,z), with 𝔖 = Σ r_l (vz)^l / Π_{j≤l+1}(1+jz/E). The z-only factor is exactly 1.
- **Status of the code.** The stored file cPP_z9 is an implementation, checked exactly:
  - against the closed form at z^0..z^9, at 7 weight pairs;
  - against eq:pf, exactly through z^10 at 6 points;
  - at z^0 and z^1 symbolically, by hand. Three independent computations agree on z^1.

  The claim that "genRmat computes R at all orders" is **not** proved, and has been downgraded to evidence. Referees flagged four gaps: the PF equation is solved only at sample points; the rational reconstruction in u is heuristic; the ansatz degree bounds are unproved; and levels 0–6 come from a checkpoint whose run was not logged. No proof step uses the code.

## 2. Continuity lemma: proved

The proof runs in four steps.
1. In log coordinates the critical points are q_j = α_j + Y, with L̂(Y) = 0, det Hess = W/3 and W = −Y L̂′(Y).
2. The saddle coefficients are analytic in (y,α) on the semisimple locus. That locus contains the diagonal for every y ∉ {0, −1/27}.
3. By WKB uniqueness, the ratio 𝒜/R^α does not depend on y. Its value at y→0 comes from a direct decoupling into three Gamma integrals: N = √(ε1ε2)Θ(−z/ε1)Θ(−z/ε2), independent of α and of the colour.
4. Hence eq:normNLR is derived, not cited, and BS's bare "by continuity" assertion (line ~4345) is now proved.

## 3. Referee objections and how they were handled

- **Universal constant.** The Gaussian square roots fix N only up to a universal constant c with c² = ±1. Read literally, the formula 𝒜(α*)/N could break symplecticity by a sign. *Fixed:* the rows are defined as the limit, and the z^0 term Δ^{-1/2} fixes c.
- **Conifold point.** At y = −1/27, for colours 1 and 2, only row 0 and the rows in the conifold frame are analytic. The large-radius rows 1 and 2 are not, because A has a logarithm there. *Restricted accordingly.*
- **Circular citation.** One proof cited BS 4337–4343 for the y→0 limit, but BS deduce that limit from eq:normNLR. *Replaced* by the direct decoupling argument.
- **F1–F3 not written out.** These are coordinate invariance, the vanishing of total derivatives, and differentiation under Asym. They are only cited or sketched. F1 carries real weight: the closed form is computed in BS's s-chart.
- **Smaller fixes.**
  - Fixing c(z) at v = 0 needs the extra hypothesis that the coefficients extend analytically to v = 0.
  - L1 is sharpened.
  - "Cf is the identity on G_{T_F}" is a convention, i.e. part of the definition of Cf.
  - BS's row-2 insertion identity in the proof of cfasym depends on the frame. It is no longer needed, because F2 now comes from the closed form.
- **BS errata.**
  - Line ~3383 has the wrong sign: W_i should be +(E/9)(3Y+E)².
  - In q^{1/3}-coordinates the determinant in eq:crituiWi is 243·W.
  - In "q_0^cr = c(1,ω,ω²)" the vector is the q^{1/3}-vector, not q.
  - Line drift: eq:Wadm is at 3225 and eq:Rcf at 5334.

## 4. Round-1 statements that become proofs

- **Closure lemma:** unconditional.
- **F2 (all orders):** follows from the closed form. It no longer depends on the content of eq:Rcf.
- **L1, sharpened:**
  - [t_c^k] of the frame-B row lies in span{λ_1..λ_k};
  - [t_c^k] of ρP̂0 lies in span{λ_1..λ_{k+1}}.
- **Stored form of the rows:** proved.
- **Analyticity (R) at colours 1 and 2:** proved.
- **Z∞:** lim_{z→∞} z v𝔖 = −t_c, i.e. round-1's lim zρP̂0 = E.
- **Symplecticity:** Π_ij = δ_ij for every choice of (c1, XX), so K = 1. The proof is BS (Symplecticity) on generic tori, passed to the limit, plus the drop-out identity. It no longer needs the round-1 Lagrange identity.
- **Regular-colour induction, step 1:** proved.
- **Residue lemma:** its inputs for row 0 (L1, F2, K = 1) are proved, given the round-1 σ-form.
- **Lemma H2:** its stage-1 inputs are proved. The hand identities (I1)/(I2) still need refereeing.
- **Frame-B skeleton:** N1, N2 and N5 are resolved.

## 5. What remains open

- F1–F3 in written form.
- Lemma H2's identities (I1)/(I2), and with them rows 1 and 2 of the residue lemma.
- Mixed edges E_0j: not touched in round 2.
- Givental–Teleman on T_F, needed for the gap corollary. It is proved conditionally, using BS lem:lho/thm:haegen and BS line 4056 (regularity of GW^{T_max} on the diagonal, which BS argue only briefly).

## 6. Revised next steps

1. Write out F1–F3 as formal Wick identities in n dimensions. Alternatively, prove PF for 𝔖 directly from a recurrence for r_l. That would make the regularity theorem independent of Laplace asymptotics, because the closed form could then serve as the definition.
2. Referee (I1)/(I2), which closes Lemma H2 (Z∞ is now available as an input).
3. Feed the sharpened L1 and the closed form into the residue lemma.
4. Mixed edges E_0j, using the closed forms at all three colours. This is the main new task.
5. Frame-B skeleton step 4: check BS line 4056.
6. Gap corollary, combining the regularity theorem with Proposition 11.
7. Report the BS errata.

