# Working title and abstract (design; rules 5 and 8; drafts for the authors, DEC D12)

## Working title
**Conifold regularity of the formal quintic, graph by graph**

Alternatives for the authors: "The formal quintic at the conifold point: a graph-local conifold gap"; "Graph-local conifold regularity
for the formal quintic". `\draftnote{working title}` in `main.tex`. Short title for the running head: "The formal quintic at the
conifold point".

## Abstract (plain text; the LaTeX form is in `NOTE/main.tex`)

The formal quintic of Lho and Pandharipande is the O(5)-twisted equivariant Gromov-Witten theory of P^4 specialised at the fifth roots
of unity. Its genus-g potential is a Givental-Teleman sum over stable graphs of genus g, with the R-matrix normalised at q = 0 by
quantum Riemann-Roch. Near the conifold point u = 1 - 5^5 q = 0 we factorise the R-matrix as R = R' M_0(tau(u) zhat), where M_0 is the
R-matrix of the Landau-Ginzburg model p^4/4 + v/p with primitive form dp and tau(u) zhat is the natural variable defined by the
canonical coordinate. We prove that, for every g >= 2 and in every conifold frame of a two-parameter rational family, the term of
every leg-free stable graph of genus g with at least one edge is regular at the conifold point, and the single-vertex graph contributes
c_g t_c^(2-2g) plus a regular series, where t_c is the flat coordinate of the vanishing period and c_g = (-1)^(g-1) B_2g /
(2g(2g-2) 5^(g-1)). Summing over graphs gives the conifold gap for the formal quintic in every genus. The proof reduces the statement,
graph by graph, to a lemma about the model, proved by topological recursion on the model curve, and to the equality of two
normalisations of the R-matrix: the canonical one at the conifold point and the quantum Riemann-Roch one at q = 0. We prove this
equality through an identity of formal stationary phase for a Mellin-Barnes integral, which gives the value at q = 0 of its saddle
series. The results are conditional on the quantum Riemann-Roch normalisation of the R-matrix at q = 0 as stated in the literature, and
they use cited theorems: the genus-0 mirror theorem and the genus-0 axioms of the twisted theory, Teleman's classification, the
R-matrix action of Pandharipande, Pixton and Zvonkine, and results on topological recursion, the last only for the model curve. We do
not treat the quintic threefold itself.

## Checks against the rules (author-facing)
* Rule 5: no stage, workflow, agent, referee, status tag, file name or date; "Conjecture N" is not named (the note calls the property
  N(g) only in `def:set-N`; the abstract states the theorem).
* Rule 8: conditionality exactly as OUTLINE s.0.2: (QRR) as the only quoted input about the formal quintic's R-matrix; (M) and (G0) as
  "the genus-0 mirror theorem and the genus-0 axioms of the twisted theory"; (T) as "Teleman's classification"; (Q) as "the R-matrix
  action" ([PPZ15]) and "results on topological recursion" ([Eynard:2007kz], [DuninBarkowski:2012bw] as recalled in
  [Dunin-Barkowski:2015caa], [ABDKS25]; [Hock25] is no longer an input, DEC Q4); the classical formula [AC96, (1.7)] is not named in the
  abstract (it is named in `sec:int-inputs`); TR only on the model curve; nothing about the quintic threefold.
* REV_2 M5: (a) "leg-free stable graph"; (b) the model in the note's normalisation `p^4/4 + v/p` (DEC D13); (c) "in every conifold frame
  of a two-parameter rational family" (DEC D14: `(alpha, beta)` in `Q^2`); (d) the genus-0 mirror theorem stays in the list (R1, DEC D2).
* No overclaim: "the conifold gap for the formal quintic" is Corollary B (`cor:mai-gap`), now proved with its definition of the conifold
  expansion of `F_g` (G5_2 s.2.2-s.2.4); no novelty claim for the stationary-phase lemmas, `Ginf` or the chart mechanism; Theorem Gamma
  is described, not claimed as first.
