# REF_3: the membership argument (A3.4 (a)) and the "formal change of frame" edits (A3.1-A3.3)

Date 2026-11-18. Read-only on the paper. Verdict: **SOUND AFTER FIXES**.
The algebra of the membership argument is correct in every step; no error was found. The fixes are wording:
one sentence claims slightly more than is proved (F1), and a few sentences are inconsistent or ambiguous (F2-F9).

Read: `main.tex` (abstract); `s01_intro.tex`, `s01_intro_2.tex` l.24-75, `s01_intro_3.tex` l.19-36; `s02_setup_3.tex`
l.26-29; `s02_setup_6.tex`; `s03_gauge_8.tex`; `s10_main.tex`, `s10_main_2.tex`; `b_signs_2.tex` l.92-110; `b_signs.tex`
l.6; `a1_cyclicLR.tex` l.5; `AMEND_1_chain_2.md` (A3), `APPLIED_1.md`, `audit_2026-11-18/verify/E7.md`.
Script: `<user home>\AppData\Local\Temp\claude\d7c1b608d03d194a\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\ref3\check_membership.py`
(exact rational checks of the identities of step (1), (3) and of all zeta-weights of step (2); all pass).

---------------------------------------------------------------------------------------------------------------------

## 1. The membership argument, re-derived (`b_signs_2.tex` l.108)

Step (1). N = N_0[log q] with N_0 = Q[L^{+-1}, c_0^{+-1}, c_1, c_2, c_3, i_0, i_1, i_2][W^{-1}], because W is free of log q.
The substitution (i_0, i_1, i_2) -> (x, pi_1, pi_2), x = i_0/c_0, pi_j = i_j - x c_j, is triangular and invertible
(c_0 is a unit), and W = c_0^2 + c_0(pi_1 + x c_1) - c_1 c_0 x = c_0^2 + c_0 pi_1. So W lies in the x-free polynomial
ring, the localisation commutes with adjoining x, and N_0 = B[x], a polynomial ring over the domain B. CORRECT.
The shift a_k (i_j -> i_j - k c_j, log q -> log q + k) preserves N_0, sends x -> x - k, fixes pi_j and W. A polynomial
in x over a Q-algebra that is invariant under x -> x - k for one k != 0 is constant. CORRECT.
Inputs used: F_g^N in N (not only in N (x) Q(zeta)): `def:mai-Fcf`, "Second" (the sigma_k argument), which uses
`thm:con-glue`. F_g^N free of log q: it is a polynomial in the explicit representatives of `cor:gau-rowsNN` (i), which are
free of log q. F_g^N fixed by a_k: `def:set-confexp` and the end of the justification of `def:mai-Fcf`; this needs that
a_k commutes with D, i.e. a_k(i_3) = i_3 - k c_3, which holds (the relation Conc(c_0, I_1) = 0 determines the third jet
J_3 = log q c_3 + 3 c_2 + i_3 of I_1 from c_j and the jets J_0, J_1, J_2, all fixed by a_k). CORRECT.

Step (2). rho: L -> zeta L, fixing zeta and the other generators, is an automorphism of N (x) Q(zeta) (c_0 and W are fixed).
It commutes with D: D L = L(L^5 - 1)/5 has weight 1; D c_3 has coefficients in Q[L^5]; i_3 is defined by a relation
with coefficients in Q[u] = Q[L^{-5}] and leading coefficient -L^{-5} c_0. CORRECT.
Weights: by Theorem C (`thm:gam-D`, D_0 = 1 in every frame) and `thm:gau-R0` (iii) (injectivity of L -> (1-5^5q)^{-1/5} on
Laurent polynomials), A^act_n = A^can_n = L^{4n+1} ph_n(L^{-5}) as Laurent polynomials, of weight 4n+1 = 1-n mod 5. The
recursion Q_{a+1,n} = coef (D Q_{a,n-1} + L Q_{a,n}) has rho-fixed coefficients (c_0^2/W, W^2 L^{-5}/c_0^2, c_0^2/W,
c_0^{-1}), and both terms have weight a+2-n if Q_{a,m} has weight a+1-m. So rho Q^{(0)}_{a,n} = zeta^{a+1-n} Q^{(0)}_{a,n}, and
with colour i := zeta^{i(a+1-n)} colour 0 and rho(zeta) = zeta: rho Q^{(i)} = Q^{(i+1)}. Delta_i^{1/2} = zeta^{4i} c_0/L ->
zeta^{4i-1} c_0/L = Delta_{i+1}^{1/2}. Psi_Q: m_a has L-weights (1, 2, -2, -1, 0) = a+1 mod 5, so (Psi_Q)_{aj} ->
(Psi_Q)_{a,j+1}. Hence rho(R~_{ij}) = R~_{i+1,j+1} (eta^{-1} is rational), rho(1_TFT) is 1_TFT relabelled, and the TFT
correlators are rational multiples of powers of Delta_i. CORRECT.
"rho permutes the colours and fixes F_g^N": F_g^N is a sum over all colourings of a contraction built only from
R~_{ij}, Delta_i^{+-1/2} and rational numbers; rho acts as the relabelling i -> i+1 of the summation index. This is the
argument for sigma_k in `def:mai-Fcf`, and it is simpler here because rho fixes zeta. CORRECT.
"therefore contains only powers of L^5": B = sum_{k=0..4} L^k B_5, a direct sum (W is free of L), rho = zeta^k on the
k-th summand; the invariants of B (x) Q(zeta) are B_5 (x) Q(zeta), and its intersection with B is B_5. CORRECT; the
direct-sum decomposition is not said (F6).
Theorem C is really needed: without D_0 = 1 one has rho(Q^{(0)}) = Q^{(1)} D_0(zeta^{-1}z)/D_0(z), not a colour permutation.

Step (3). c_j = D^j c_0; W = c_0^2 (W/c_0^2) and W^{-1} lie in S; pi_1 = (W - c_0^2)/c_0; D W = 2 c_0 c_1 + c_0 i_2 - c_2 i_0
= 2 c_0 c_1 + c_0 pi_2 (Leibniz only; no relation of N is used), so pi_2 = (D W - 2 c_0 c_1)/c_0. S is D-stable over
Q[L^{+-5}] as D L^5 = L^5 (L^5 - 1). CORRECT.

Remark (optional shortening, not a defect). Steps (1) and (3) can be avoided: the representatives of `cor:gau-rowsNN` (i)
and Delta_i^{+-1/2} visibly lie in S[L] (x) Q(zeta) (row 0 is A_n(L)/c_0, the recursion uses D, L, (W/c_0^2)^{-1},
(W/c_0^2)^2 c_0^2 L^{-5}, c_0^{-1}), S[L] = sum_k L^k S is direct, and the rho-invariants in N are S. The printed route
through B[x] is equally valid.

---------------------------------------------------------------------------------------------------------------------

## 2. The recollection of the definitions of [CGYZ26]

No extract of [CGYZ26] is stored under `audit_2026-11-18/work/` (the 26 hits in `work/U08/main.txt` are the text of the
note), and E7 says "CGYZ26 itself not consulted". So before this report the recollection at l.100 had not been compared
with the source by the audit. `refs.bib` pins version 2; the two remarks carry a draft note (HTML rendering).
For this report the HTML rendering of arXiv:2609.13921v2 was queried through an automatic page reader (first 100000
characters; answers are paraphrases, so this is a consistency check, not a line-by-line reading). Result: l.100 agrees
with it on every point:
* Def. 1.1: K = Q(q); the differential K-algebra generated by I_0^{+-1} and (D tau)^{+-1}, tau = I_1/I_0, D = q d/dq.
* Def. 1.2: Delta in K a local parameter, expansion in C((Delta)); omega_0(0) != 0, omega_1 vanishing simply; the pair
  is admissible if I_0 -> omega_0, D tau -> D tau_c (extended to derivatives and inverses) preserves the relations; the
  expansion is then a differential homomorphism iota; "conifold regular" = image in C[[Delta]].
* Ex. 1.3: Delta = 1 - 5^5 q, omega_0 = 1 + O(Delta^3), omega_1 = Delta + (7/10) Delta^2 + ... (= (f_0, varpi_1)).
* Rem. 1.4: omega_0' = omega_0 + c_1 omega_1 + c_2 omega_2, c_1, c_2 in Q, omega_2 = Delta^2 + O(Delta^3) (= `def:set-frames`);
  admissibility is asserted with a short argument.
* Conj. 1.5: F_g^c - (-1)^{g-1} B_{2g} / (5^{g-1} 2g(2g-2) tau_c^{2g-2}) in C[[tau_c]], "for any admissible pair".
* Def. 1.8: (-1)^g with a constant kappa, omega = Delta + O(Delta^2). Rem. 1.9, Prop. 3.1, Thm. 2: as described in the note.
Two points where the note says more than this reading supports:
(a) The membership paragraph asserts "F_g^cf is its conifold expansion in the sense of Def. 1.2" without "as recalled"
    and without the admissibility hypothesis (F1).
(b) Def. 1.2 asks only that omega_1 vanish simply, and Conj. 1.5 is stated for any admissible pair; the normalisation
    omega = Delta + O(Delta^2) appears as a hypothesis in Def. 1.8 (2) and in Prop.-Def. 4.5 (1), and in the pair of
    Ex. 1.3. By the note's own rescaling paragraph the constant depends on it. The pin at l.100 should say that
    Conj. 1.5 is read with this normalisation (F3). To be confirmed on the PDF.

---------------------------------------------------------------------------------------------------------------------

## 3. Does the conclusion follow from membership alone?

Write iota for the expansion of [CGYZ26, Def. 1.2] on their algebra R (functions). What is wanted: iota(F_g) = F_g^cf.
* ev(S) lies in R: ev(L^5) = (1-5^5 q)^{-1} in Q(q), ev(c_0) = C_0, ev(W/c_0^2) = C_1, and ev is a differential homomorphism.
  So membership gives F_g = ev(F_g^N) in R. No further input.
* iota o ev and C_{alpha,beta} are differential homomorphisms on S which agree on L^{+-5} (-> Delta^{-+1} = u^{-+1}), c_0^{+-1}
  (-> varpi_0^{+-1}), (W/c_0^2)^{+-1} (-> (D t_c)^{+-1}); hence they agree on S, and iota(F_g) = iota(ev(F_g^N)) = C(F_g^N).
  Injectivity of ev is NOT needed, and no "evaluation maps agree" statement beyond this is needed.
* What IS needed and is not proved in the note: that iota exists, i.e. that the pair of the frame is admissible. In the
  note's terms admissibility is the statement that every relation among the evaluated generators is respected by
  C_{alpha,beta}; the note cannot prove it without the injectivity of ev (true by E7, with an undeclared external input).
  It is asserted in [CGYZ26, Ex. 1.3, Rem. 1.4]. l.100 says "wherever theirs is defined" and l.110 says "under
  admissibility", but the membership paragraph states the conclusion without it. (F1)
* The membership itself is conditional on the inputs of Theorem C. The paragraph says "This holds" flatly. (F1)
* The equality of the constants is a comparison of two printed formulas; it holds as recalled.
* l.106 draws the conclusion ("Consequently ... has the conifold gap with the constant of Conj. 1.5") one paragraph
  before the membership on which its reading in the sense of [CGYZ26] depends. (F2)
With F1, F2 the wording claims exactly what is proved.

---------------------------------------------------------------------------------------------------------------------

## 4. Findings and replacement texts

**F1 (gap in wording; the only substantive one).** `sections/b_signs_2.tex` l.108.
Old:
```latex
then $F_g=\ev(F_g^{\Nring})$ lies in the algebra of \cite[Def.~1.1]{CGYZ26}, and $F_g^{\mathrm{cf}}$ is its conifold expansion in the sense of \cite[Def.~1.2]{CGYZ26}. This holds (the letters
```
New:
```latex
then $F_g=\ev(F_g^{\Nring})$ lies in the algebra of \cite[Def.~1.1]{CGYZ26} as recalled above, since $\ev$ maps $L^5$, $c_0$ and $\mathsf W/c_0^2$ to $(1-5^5q)^{-1}$, $C_0$ and $C_1$; and for every frame whose pair $(\varpi_0,\varpi_1)$ is admissible, which is asserted in \cite[Ex.~1.3, Rem.~1.4]{CGYZ26} and is not proved here, $F_g^{\mathrm{cf}}$ is the conifold expansion of $F_g$ in the sense of \cite[Def.~1.2]{CGYZ26}: that expansion composed with $\ev$, and $\Cexp_{\alpha,\beta}$, are homomorphisms of differential rings on $\mathscr S$ which agree on its generators. The injectivity of $\ev$ is not used. The membership holds, under the inputs of \cref{thm:gam-D} (the letters
```

**F2 (wording).** `sections/b_signs_2.tex` l.106, last sentence.
Old:
```latex
the formal quintic has the conifold gap with the constant of \cite[Conj.~1.5]{CGYZ26}.
```
New:
```latex
the formal quintic has the conifold gap in the sense of Corollary~B, with the constant of \cite[Conj.~1.5]{CGYZ26}; by the membership statement below this is also the gap for the conifold expansion of \cite[Def.~1.2]{CGYZ26}, when the pair is admissible.
```

**F3 (wording; pin of the recollection).** `sections/b_signs_2.tex` l.100, two replacements.
Old: `In our notation, \cite[Def.~1.1]{CGYZ26} works over`
New: `In our notation, and as read in version~2 of the preprint, \cite[Def.~1.1]{CGYZ26} works over`
Old:
```latex
and $\varpi_1=u+O(u^2)$ is the normalisation of \cite[\S4.2, Prop.-Def.~4.5\,(1)]{CGYZ26}.
```
New:
```latex
and $\varpi_1=u+O(u^2)$ is the normalisation of the pair of \cite[Ex.~1.3]{CGYZ26} and of \cite[Def.~1.8\,(2); \S4.2, Prop.-Def.~4.5\,(1)]{CGYZ26}; \cite[Conj.~1.5]{CGYZ26} is read with this normalisation, on which its constant depends (see the rescaling below).
```
(The pins Def. 1.8 (2) and Ex. 1.3 come from the automatic reading of the HTML; confirm on the PDF, as the existing
draft notes already require.)

**F4 (wording; reads as a statement that [CGYZ26, Thm. 2] is the gap of the formal quintic).** `sections/s01_intro_3.tex` l.30.
Old: `the formal quintic is a different theory, and its conifold gap is \cite[Theorem~2]{CGYZ26}.`
New: `the formal quintic is a different theory. The conifold gap of the quintic threefold is \cite[Theorem~2]{CGYZ26}.`

**F5 (wording; the only remaining sentence that reads as a statement about functions near u = 0).** `sections/s01_intro.tex` l.16.
Old: `Near the conifold point the lowered rows`
New: `In the conifold expansion defined below, where $L=s^{-1}$, the lowered rows`

**F6 (wording; a step not said).** `sections/b_signs_2.tex` l.108, step (2).
Old: `and fixes $F_g^{\Nring}$, which therefore contains only powers of $L^5$.`
New:
```latex
and fixes $F_g^{\Nring}$. As $\mathsf W$ is free of $L$, $\mathsf B=\bigoplus_{k=0}^4L^k\mathsf B_5$, where $\mathsf B_5$ (local) is the ring displayed in (3), and $\varrho$ acts on the $k$-th summand by $\zeta^k$; so $F_g^{\Nring}\in\mathsf B_5$.
```
and in (3), old: `(3) $\bbQ[L^{\pm5}][c_0^{\pm1},c_1,c_2,c_3,\pi_1,\pi_2][\mathsf W^{-1}]\subset\mathscr S$, as`
new: `(3) $\mathsf B_5:=\bbQ[L^{\pm5}][c_0^{\pm1},c_1,c_2,c_3,\pi_1,\pi_2][\mathsf W^{-1}]\subset\mathscr S$, as`
(add `$\mathsf B_5$` to the list of local letters).

**F7 (wording; cross-reference).** `sections/s02_setup_6.tex` l.20, last sentence of `def:set-confexp`.
Old: `The comparison with the conifold expansion of \cite{CGYZ26} is \cref{rmk:mai-constant}.`
New: `The comparison with the conifold expansion of \cite{CGYZ26} is in \cref{sec:sgn-lit} (see also \cref{rmk:mai-constant}).`

**F8 (wording; consistency of the phrase).** `sections/s01_intro_3.tex` l.32. The abstract and `def:set-points` say "formal
change of frame"; here "isotropic frame" is used, a term the note defines only for planes, in the appendix (`lem:gau-iso`).
Old: `is a formal change of isotropic frame on the ring $\Nring$`
New: `is a formal change of frame on the ring $\Nring$`

**F9 (wording).** `sections/s10_main_2.tex` l.40 (`rmk:mai-constant`) and `sections/b_signs.tex` l.6.
Old (s10_main_2): `for the dictionary and the rescaling of the two periods of a frame, \cref{sec:sgn-lit}.`
New:
```latex
for the dictionary, the rescaling of the two periods of a frame, and the statement that $F_g$ lies in the algebra of \cite[Def.~1.1]{CGYZ26}, so that for an admissible pair $F_g^{\mathrm{cf}}$ is its conifold expansion in the sense of \cite[Def.~1.2]{CGYZ26}, see \cref{sec:sgn-lit}.
```
Old (b_signs): `The comparison with \cite{CGYZ26} at the end is proved there and enters no proof of the body.`
New: `The comparison with \cite{CGYZ26} at the end is proved there, for the definitions of \cite{CGYZ26} as recalled there, and enters no proof of the body.`

---------------------------------------------------------------------------------------------------------------------

## 5. Consistency of the "formal change of frame" wording (check 4)

| place | text now | judgement |
|---|---|---|
| abstract (`main.tex` l.19) | "expand ... in a conifold frame: ... substituted ... This is a formal change of frame; no analytic continuation is involved"; "the expanded $R$-matrix"; "the conifold gap for this expansion" | consistent. The phrase "the holomorphic limit in the conifold frame" of A3.1 is NOT in the file (APPLIED_1, doubt 7, still mentions it); its absence is the better state, as the phrase is undefined in the note |
| intro l.16 | "Near the conifold point the lowered rows ... grow with $L$" | F5 |
| intro l.18, l.28 | "expands at $u=0$ in the frame"; "gives its conifold gap ... (Corollary B)" | acceptable with item 3 of "What is not claimed" |
| intro l.26 | B-model prediction, literature | not a claim of the note |
| `cor:int-gap` and l.36 | "the image under $\Cexp$ of the genus-$g$ graph sum written in $\Nring$" | consistent |
| `rmk:int-constant` | "nothing about the quintic threefold follows" | consistent |
| "What is not claimed" 1 | ambiguous "its" | F4 |
| "What is not claimed" 3 | new item | consistent; F8 for the word "isotropic" |
| `def:set-points` | "the note uses their conifold expansions ... formal change of frame, and no analytic continuation is involved" | consistent |
| `def:set-confexp` | defines the map; last sentence points to the wrong place | F7 |
| `lem:set-canon` | "constants $u_i^{\mathrm c}$, their values at $u=0$" | formal (values of the series); the constants are said to play no role. Acceptable |
| Theorem C / abstract "two normalisations" | identity of Laurent polynomials in $L$ (`thm:gau-R0` (iii)) | algebraic in $L$; no continuation. Consistent |
| `cor:mai-gap`, `def:mai-Fcf` | "the conifold expansion of $F_g$" defined as $\Cexp(F_g^{\Nring})$ | consistent |
| `rmk:mai-constant` | "nothing about the quintic follows" | consistent; F9 |
| `sec:sgn-lit` l.100-110 | comparison, rescaling, membership, differences | F1, F2, F3, F6 |
| `a1_cyclicLR.tex` l.5 | "no ... analytic continuation or conifold frame is used" | consistent |

No sentence asserts a statement about the quintic threefold, and after F4, F5 none suggests analytic continuation. The
cautious form "no analytic continuation is involved" is used everywhere; the stronger "is not an analytic continuation"
(which would need the irreducibility of the monodromy) appears nowhere.

## 6. Not checked

The proofs of Theorem C, `thm:con-glue`, `cor:gau-rowsNN`, `lem:mai-Psi` (taken as stated); the last 48000 characters
of the HTML of [CGYZ26] and its PDF; the admissibility claims of [CGYZ26].
