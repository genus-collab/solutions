<!-- workflow wf_74ed0422-50f (3 provers, 6 referees, synthesis); full reports in its journal.jsonl -->
**Round 3 status note (2026-11-02): H2, the residue lemma for all rows, mixed edges, and PF / F1-F5**

**What is now proved to all orders.** Every statement below is all-orders algebra. Exact computer checks were run as controls only and are never used as proof steps. Six referee reports came back, two per item, all "sound with minor fixes" and none with a major objection.

1. **Lemma H2 (singular colour, frame B, conifold polarisation).**
   - Put σ = ρ·P̂(-ψ) with ρ = -Ev/t_c. Then σ^(-1) = (0, E, -E²) and σ^(k) = 0 for k ≤ -2.
   - For k ≥ 0, σ_0^(k) and σ_1^(k) lie in Λ_{k+1}, and σ_2^(k) lies in Λ_k.
   - At no order is there a ψ-polynomial part, and no component has a constant at k ≥ 0.
   - Round 1's hand identities (I1) and (I2) are now proved. They follow from closed forms for rows 1 and 2 (Lemma 1) and a conjugation by ρ (Lemma 2).
   - The proof uses only round-2 results: Prop 8, and from Cor 9 the facts L1, F2, Z∞ and analyticity (R), plus K = 1. The constants of σ_2 are removed by the ψ → ∞ limit of K = 1.
   - The closure lemma, the rank-3 flow and round 1's version of H2 are not used.

2. **Residue lemma for all three rows (Theorem E00).**
   - Lemma R (algebraic, any rank) combined with H2 gives an explicit leaf formula for Ê_00.
   - This yields F3 (regularity), the L3 filtration (i-1)+(l-1) ≤ k, symmetry, leg bounds, and the rank-one value Ê_00|₀ = -(e1e2E⁴/27) λ_1 λ_1'.
   - It needs only Steps 1-2 of H2 and K = 1.

3. **Mixed and regular edges (Theorem M).**
   - E_0j and E_j0 are regular. E_j0 is the transpose of E_0j; this uses the off-diagonal part of Prop 10.
   - Explicit leaf formula: rank one at t_c = 0, leg bounds L2 in Λ_{k+1} with no constant, for both orientations.
   - E_jj' (j, j' = 1, 2), the dilatons Dil_j and Δ_j are regular, and Δ_j ≠ 0. This is F1 at the regular colours.
   - The E0j report treated H2 as an assumed hypothesis. Item 1 now discharges it, so Theorem M rests only on the round-2 base.

4. **Picard-Fuchs (PF) for the closed form, and F1-F5.**
   - Theorem A: conjugating L_pf by e^{u_i/z} gives (z³/27)(τ³ - 27σ³) in the basis V^n/(a)_n. So the closed form solves eq:pf at every colour with no Laplace input, and it is the unique PF solution in that ansatz, up to a factor depending on z only.
   - F1-F5 are proved in n dimensions by formal Wick algebra: coordinate invariance, vanishing of total derivatives, differentiation under Asym with a moving critical point, Fubini, and formal Stirling.
   - Corollary B: BS's Wick numbers satisfy the recurrence for all l.
   - Round-2 Lemma 3, Prop 5, Theorem 6, Prop 8, Cor 9 and Prop 10 now hold modulo only the BS genus-0 inputs C1-C3, which you have allowed.

**Referee objections, and where the reports disagree.**
- **E0j outcome label.** Both E0j referees objected that "proved modulo cited" was wrong while H2 was an unproved assumption. Item 1 removes the objection; I relabel the result "proved modulo the round-2 base".
- **Carrying K = 1 and orthogonality to the conifold point.** Several referees noted that K = 1 and (O) were proved only for y ≠ -1/27, and asked for an explicit argument at the conifold. It is one line: each z-coefficient of the pairing is a meromorphic germ at the conifold, and [z^n] commutes with [t_c^k]. I added this as Lemma CL.
- **Gauge drop-out.** This needs the same c1 at every colour. That holds, because c1 depends only on y. The "c1 = 1 gauge" used in the checks is a hybrid (c1 = 1 but XX = θ log C1cf).
- **Normalisation ρ(0) = 1.** This needs t_c = Δ_c + O(Δ_c²), equivalently C1cf(0) = -1. It is cited now; only the exact rank-one constants depend on it.
- **Theorem A rigour.** The two PF referees framed the infinite sums differently: one accepted V-adic continuity, the other preferred a z-adic argument. Both arguments are valid. I adopt the z-adic one, since each z-order is then a finite computation.
- **Corollary A5 normalisation.** The statement had a wrong factor: the solution is z·b_1(z)·vS, not b_1(z) times the closed form. Fixed.
- **Corollary B dependencies.** One referee flagged that Corollary B also uses the PF-under-the-integral identity of round-2 Lemma 3, step 2. Its formal content is F3 applied three times plus F2, but this is not yet written out. I also now state the l-independent sign cancellation in the s3 integral, and the convention for branches at α = -E.
- **Wording fixes.**
  - H2(d) now holds only within the family c1 = C1cf·f(Δ_c).
  - The off-shell control values are now indexed from k = 0.
  - (M-f) is demoted to a remark.
  - The χ_j evidence line, which had an inconsistent third entry, is dropped.
  - "Formal Stirling" replaces "proves BS exmpl:stirling".
  - The eq:pf line number is now about 2627-2632.

**What Theorem LC_P2 and Theorem R_P2 now depend on.**
- **Edge inputs T-F3, T-L2, T-L3: proved** (items 2-3), given:
  - (a) BS C1-C3 (allowed);
  - (b) the written-out step 2 of round-2 Lemma 3, which is routine from F2 and F3;
  - (c) the reading of BS 5257-5267 in round-2 Prop 8 step 2, which a referee re-checked.
- **Remaining inputs that are not edge statements:**
  - (d) L1 and F3 for Dil_0, already covered by round-2 Z∞ (to be restated);
  - (e) F1 at colour 0 (homogeneity Δ_0 = 3e1e2E³/(t_cρ)²): proved;
  - (f) Hodge side: the dilaton equation, and the ports of Theorem C and of the string identity with c → E;
  - (g) Lemma A, M_irr(h) for h ≥ 5. This is proved modulo literature for the conifold; we still need to check that it transfers to the weights (E, -e1, -e2);
  - (h) writing out the port of Theorem R itself (rank 3, c → E, x c³ → -3e1e2E³).
- **For the gap identification only:** (i) the frame-B formula Prop 13.2_P2 on T_F.
- **For identifying the closed form with the Gromov-Witten rows** (the z-only factor c(z) = 1): round-2 Prop 8 plus the transport from y → 0.

**Revised next steps.**
1. Write one consolidated document: rows on T_F defined by the closed form, with PF and F1-F5, then H2, E00 and Theorem M. Apply every referee fix, and write out (b), (d) and (h).
2. Port Theorem C and the string identity (f), and check that Lemma A transfers (g).
3. Prove the frame-B formula on T_F and derive the gap corollary.
4. Optional:
   - a direct proof that the Wick numbers satisfy the recurrence, without Laplace;
   - closed forms for the higher anti-diagonals of Ê_00;
   - the t_c = 0 values of the regular rows, which would explain the rank-2 relation.

