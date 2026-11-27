# Stage 5, task `lit`, part 3: consequences for Stage 5 (s.4), what was not read (s.5), record (s.6)

Continuation of `LIT.md` (s.0-s.1) and `LIT_2.md` (s.2-s.3). Same labels and notation.

## 4. Consequences for Stage 5 (my assessment; nothing here is a proof)

1. **(QRR) can now be cited for the formal quintic itself.**
   * The citation is GJR arXiv:1812.11908v1 eq. (17), in the proof of Lemma 5.2, with `R = e^(-u/z) S C_alpha^-1`, so
     `R(q=0) = diag(h_(w_i))` to all orders. A second, independent statement is GJR arXiv:1709.07392v1 eqs. (19)-(20), taken at t = 0.
   * Both rest on citations of Givental arXiv:math/0008067 Thm 1 (with the Remark after it) and of Coates-Givental Thm 1 / Cor. 1.
   * The step "the R-matrix of the full CohFT restricts at q = 0 to the vertex R-matrices" uses uniqueness of R for a semisimple
     CohFT, i.e. Teleman arXiv:0712.0160. I did not open that paper; it should be added to (Q) with a statement number.
   * The sign of the twist weight is settled (`LIT.md` s.1.4-1.5). The human check of the `z^(5k)` part (Stage 4, part 2, s.6.2
     item 3) can be done by reading GJR (17), math/0008067 s.2.2-2.3 and FP Prop. 2.
2. **The closed form of the free coefficients.**
   * Stage 4 lists `Sum_m f_m z^(5m) = Ginf(z)` as CONJECTURAL. GJR Lemma 7.7, eq. (37) (m = 5, k = 1, lambda = 1) PROVES the same
     formula for the orbifold-point values of THEIR normalised SP series of the LG integral, at formal-SP level.
   * To transfer it to the project's `f^SP` (and, by Theorem P(b), to `f^can`), one must identify GJR's LG SP series with the
     project's full-Stirling Mellin-Barnes SP series `c(z; L)`. Both are normalised stationary phases of the same mirror integral,
     so this is a formal-Fubini statement, (L1) + (L2): the same lemmas the chart argument needs.
   * Stage 5 should record it as "PROVED modulo (L1)-(L2), via arXiv:1812.11908 Lemma 7.7", once those lemmas are written. It
     should cite GJR rather than claim novelty.
   * The `Ginf` identity alone does not give (N). (N) is an identity at `L = 1`.
3. **(Gamma) itself.** I found no proof in the literature.
   * GJR Remark 7.2 states the equivalent identity for the formal quintic's LG integral, coefficientwise in q (`LIT_2.md` s.3.6),
     without proof, and builds on it.
   * Givental math/0108100 s.10 is the untwisted outline of the same mechanism, with the continuity in Q asserted.
   * So a written, refereed chart argument would be new as a proof, as far as I can tell. It would also supply the missing step
     in GJR s.7. Stage 5 should say so modestly: GJR may regard Remark 7.2 as standard via Givental's equivariant mirror theorem
     (their [15]), which I did not open.
4. **A second chart is available.** GJR's Lemma 7.7 is the same mechanism as the critic's step (d), but at the orbifold point `L = 0`
   (`q = oo`), where the integrand also factorises into Euler integrals.
   * Stage 5 could run both charts on the critic's 7-variable integrand: at `L = 1` it gives (Gamma); at `L = 0` it gives `Ginf`.
   * They would share one proof of (L1)-(L2). The `L = 0` chart is a free consistency test of the bookkeeping (z-powers, Gaussian
     prefactors, `x^m` insertions).
5. **Formal lemmas: write them out.**
   * No source I read gives numbered FORMAL versions of nonlinear changes of variables, parameters, or the delta lemma for complex
     data.
   * Suggested route (`LIT_2.md` s.2.3): reduce polynomial identities in the jets to the real-definite Laplace case. The inputs are
     Etingof Thms 2.16 and 3.5, Kirwin Thm 1.1 and DLMF s.2.1(iii), then Zariski density.
   * Fubini: Aarhus II Prop. 2.13, together with translation invariance for moving critical points. (L2): the Wick argument with an
     invertible quadratic part (Aarhus II Def. 2.8).
   * (L3) is covered by DLMF (5.2.1), (2.3.13)-(2.3.18), (5.11.1), (5.11.8) and uniqueness (s.2.1(iii)).
6. **Pitfalls seen in the sources.**
   * Normalisations.
     * LP eq. (18) and ZZ (5) are normalised by `Phi_s(0) = delta_(0,s)` (the I-function), not by QRR.
     * Pandharipande arXiv:1712.02528 uses `R^-1` relative to Givental-Teleman (footnote 11).
     * GJR use `x_i = L - lambda_i`; the critic uses `x_j = zeta^j - L`.
     * Givental math/0108100 Thm 9.1 uses COTANGENT characters with a plus sign; math/0008067 Thm 1 uses TANGENT characters with a
       minus sign. These agree.
   * Gamma-function normalisation. The Gamma integrals in Givental s.10 have measure `dx`, those in Iritani Lemma 4.20 `dw/w`, and
     those in GJR s.7 `dx/x` with `x^(k/m)` inserted. The critic's control "measure dy/y fails at z^2" shows that this matters.

## 5. Not read, failures, limitations (reported plainly)

* **Books and other sources I could not open** (no legitimate online text). None of these can serve as PROVED input under the rules.
  * Hormander, The Analysis of Linear PDO I, s.7.7: statement numbers such as "Thm 7.7.5" are commonly cited; I have NOT verified
    any of them.
  * Olver, Asymptotics and Special Functions (1974/1997). I only saw DLMF's source notes: pp. 16-22 for s.2.1(iii), pp. 80-88 for
    Laplace, pp. 87-88 and 293-295 for Stirling.
  * de Bruijn; Wong (1989); Bleistein-Handelsman.
  * Mumford (1983), Towards an enumerative geometry. Used only as quoted by FP, Givental and CG.
  * Jahnke-Emde-Losch (the Stirling source cited by CG and Givental); Hermite; Nemes (cited by GJR).
* **arXiv papers not opened.**
  * Teleman arXiv:0712.0160 (uniqueness of R).
  * Givental "Equivariant GW invariants" (IMRN 1996) and "Elliptic GW invariants" (math/9803053), cited by GJR and Givental.
  * Coates-Corti-Iritani-Tseng arXiv:1606.07254 and 1310.4163.
  * Lho arXiv:1807.05503 (equivariant HAE). It is possibly the study of the dependence on specialisation announced in LP
    1803.01409 footnote 4; not checked.
* **Read only through an HTML rendering, via the fetch tool's summariser** (not verbatim; numbers confirmed by a second query):
  arXiv:1702.06096; math/9810173; 1712.02528; math/0406251; 2409.03117v2; 0810.1700; math/9801049; 1203.4193; q-alg/9706004
  (Aarhus I: no numbered Fubini statement there).
* **Read in my own text extraction of the arXiv PDF.** The extraction loses minus signs and Greek letters, so signs were taken from
  the ar5iv LaTeX where it was available, or from the structure of the formula. Versions served on 2026-11-14:
  * 1803.01409v2; 0710.0889v1; 0903.1463v3; 1404.6407v4;
  * math/0108100v2; math/0110142v2; math/0008067v4;
  * 1709.07392v1; 1812.11908v1.
* **Partial reads and truncations.**
  * ar5iv renderings of 0903.1463 (from s.4.3) and 1812.11908 (from s.6) were truncated, so the PDFs were used.
  * In 1812.11908 s.7, the display of W and the sign of the exponent of `(2 pi z)` in Remark 7.2 were not legible, and the factor m
    in (37) is inferred from glyph order. All three should be read in the typeset PDF by a human.
* **Statement numbers that are mine to confirm, not the papers'.**
  * DLMF s.2.3(iii) treats the endpoint minimum. The interior case by splitting is my addition.
  * LP's `R_1`, `R_2` in 1803.01409 eq. (18): only their vanishing at q = 0 is used; their signs were not read reliably.
* **Earlier claims checked.**
  * Stage-4 K1 is confirmed: 1702.06096 eq. (5) and 1803.01409 eq. (5) both define the formal quintic.
  * Stage-4 O11 is confirmed in its numbers: 1803.01409 s.2.4.4 eq. (18) and Prop. 5 exist, but they carry the ZZ normalisation.
  * "arXiv:math/0108100 Thm 9.1 (the Bernoulli form)": confirmed, with the cotangent convention. In that paper Thm 9.1 is proved by
    citation (Remark 9.2(a)).

## 6. Record

* **Files** (all new, in `allgenus/stage5/lit/`; nothing else in the project was created, modified or deleted):

  | file | content |
  |---|---|
  | `LIT.md`, `LIT_2.md`, `LIT_3.md` | this note, each under 250 lines |
  | `lit_vertex_check.py` | exact check of s.1.5: LP (26) vs the project's vertex convention, `N_(2,0)`, `N_(3,0)`, control |
  | `logs/lit_vertex_check.log` | its output; ends ALL CHECKS PASSED; no paths printed |
  | `pdftext.py` | my stdlib helper, kept for reproducibility: text extraction from PDF content streams. No paper text is stored in the project |

* **Runs.** `py -I -B lit_vertex_check.py`, working directory the session scratch folder, under 1 s. It was the only computation,
  and at most one of my processes ran at a time.
  * My first version compared the two vertex classes as polynomials in FREE lambda classes and reported FAIL. They agree only
    modulo Mumford's relation, so the comparison now reduces both classes to a normal form first.
  * The FREE comparison is kept as line (A0), "differ (expected)". The integrals (B) were unaffected.
* **Web access.**
  * arxiv.org (pdf, abs, html), ar5iv.labs.arxiv.org and dlmf.nist.gov.
  * Four web searches with generic queries (paper titles and topics), containing no personal data.
  * The PDFs were fetched through the fetch tool, which stored them in its own cache. I copied each into a separate empty folder
    in the session scratch area and read it there with `pdftext.py` (`py -I -B`, file paths passed as arguments). Extracted texts
    stayed in the scratch area.
* **Compliance.**
  * No existing file was modified.
  * There is no `__pycache__`.
  * Python 3.14 stdlib only; no wolframscript and no Mathematica.
  * No process of anyone else was touched.
