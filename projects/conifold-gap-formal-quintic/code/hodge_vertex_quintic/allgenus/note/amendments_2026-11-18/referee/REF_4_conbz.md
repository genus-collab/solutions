# REF 4: the repaired proof of `eq:con-bz` (Theorem P(a) (i), Step 3)

Scope: `sections/s08_connection_4.tex` l.61-103 (proof of `thm:con-Pa` (i)), with `def:gam-c` and its Taylor
data (`s08_connection_3.tex`), `def:fsp-stirling`, `lem:fsp-stirfe` (`a_formalsp_5.tex`), `lem:fsp-transl`,
`lem:fsp-generic` (`a_formalsp_3.tex`). Read-only; no source file was changed.

## Verdict: SOUND

The new argument is correct and complete. Every step was re-derived by hand and the three identities involved
were checked by exact rational series arithmetic. There is no error and no gap. Six wording findings follow;
none of them affects the mathematics.

## 1. Step-by-step re-derivation

Notation: k := `\Bbbk` = A_L (x)_Q Q(zeta), sigma = L + t, rho = rho(L+t) = rho_0 + alpha t, w = z/rho.

**Rings.**
* k = Q(zeta)[L, L^{-1}, (L^5-1)^{-1}]. A_L -> k is injective (Q(zeta) is free over Q), hence so is
  A_L[[t,z]] -> k[[t,z]]. 5L is a unit; zeta^j - L is a unit because prod_j (zeta^j - L) = 1 - L^5.
  The embedding claim of l.83 is correct.
* rho_0 in k^x, so rho is a unit of k[[t]]; l_rho = log(1 + alpha t/rho_0) in t k[[t]] (logarithm of a
  series with constant term 1); U_rho = rho l_rho - alpha t = alpha^2 t^2/(2 rho_0) + ... in t^2 k[[t]];
  w in z k[[t,z]]; m_rho = sum_k b_k w^{2k-1} in z k[[t,z]] (converges z-adically).
* E in (t,z) A_L[[t,z]] (constant term S_0'(L) = 0); log a_0 = S_1(L+t) - S_1(L) + sum_k z^{2k-1} S_2k(L+t)
  in (t,z) A_L[[t,z]]. All exponentials are of series without constant term and all logarithms of series
  with constant term 1. T_z preserves the ideal (t,z), is (t,z)-adically continuous, hence commutes with
  exp and log; T_z f - f is divisible by z.
* (sigma+z)^5 - 1 has constant term L^5 - 1, a unit, so the quotient in the last display exists in A_L[[t,z]].

**The three "Taylor part" claims (l.88).** rho log rho - rho = U_rho + (rho_0 log rho_0 - rho_0) +
alpha t log rho_0, and log rho = log rho_0 + l_rho. Correct. mu(rho/z) = sum_k b_k z^{2k-1} rho^{1-2k} = m_rho.
Correct.

**`eq:con-factor`.**
* T_z rho = rho + alpha z = rho (1 + alpha w). Correct.
* T_z l_rho = log((rho + alpha z)/rho_0) = l_rho + log(1 + alpha w). Correct.
* T_z U_rho - U_rho = alpha z l_rho + rho (1 + alpha w) log(1 + alpha w) - alpha z; divided by z = rho w this
  is alpha l_rho + (w^{-1} + alpha) log(1 + alpha w) - alpha. Correct.
* T_z w = w/(1 + alpha w), so T_z m_rho = sum_k b_k w^{2k-1} (1 + alpha w)^{1-2k}. Correct.
* Sum: alpha l_rho + [(1/w + alpha - 1/2) log(1 + alpha w) - alpha + sum_k b_k w^{2k-1}(1+alpha w)^{1-2k}]
  - sum_k b_k w^{2k-1}. The bracket is the exponent of St_alpha(w) exactly as printed in `def:fsp-stirling`,
  and sum_k b_k w^{2k-1} is the exponent of St_0 by the same definition (m = 0). So the left side is
  alpha l_rho + log St_alpha(w) - log St_0(w). Correct.
* `lem:fsp-stirfe` as printed ends with "Hence St_m = Poch_m St_0" for every integer m, in Q[[w]]. This is
  the form used (m = alpha, an integer). The printed statement and the form used agree; nothing else of the
  lemma is used. The substitution w -> z/rho is a continuous homomorphism Q[[w]] -> k[[t,z]] and commutes
  with log. Correct.

**Six-factor decomposition (l.95), checked against the Taylor data printed in `def:gam-c`.**
* [t^{k+1}] U_rho = (-1)^{k+1} alpha^{k+1} rho_0^{-k}/(k(k+1)). For rho_N = 5 sigma: 5 (-1)^{k+1} L^{-k}/(k(k+1)).
  For rho_j = zeta^j - sigma (alpha = -1): (zeta^j - L)^{-k}/(k(k+1)), sum over j = p_k/(k(k+1)). This is the
  printed expansion of S_0(L+t) - S_0(L). Correct.
* [t^k] l_rho = (-1)^{k+1} alpha^k rho_0^{-k}/k. So (1/2) l_N - (1/2) sum_j l_j has coefficient
  ((-1)^{k+1} L^{-k} + p_k)/(2k): the printed S_1(L+t) - S_1(L). Correct.
* m_N + sum_j m_j = sum_k z^{2k-1} b_k [(5 sigma)^{1-2k} + sum_j (zeta^j - sigma)^{1-2k}]. Correct.
* Hence E = sum over the six rho of (T_z U_rho - U_rho)/z and
  log a_0 = (1/2) l_N - (1/2) sum_j l_j + m_N + sum_j m_j. The sum of `eq:con-factor` over the six factors
  carries -(1/2) l_N where log a_0 has +(1/2) l_N; the difference is (T_z - 1) l_N = log(1 + 5 w_N). The
  extra term of l.95 is exactly right.

**Multipliers.**
* Numerator: e^{5 l_N} Poch_5(w_N)(1 + 5 w_N) = (sigma/L)^5 prod_{k=1}^5 (1 + k z/(5 sigma))
  = prod_k (5 sigma + k z)/(5L)^5. (Poch_5(w) = prod_{i=0}^4 (1 + i w).) Correct.
* Each rho_j: e^{-l_j} Poch_{-1}(w_j) = (zeta^j - L)/(zeta^j - sigma) * (zeta^j - sigma)/(zeta^j - sigma - z).
  Product over j: prod_j (zeta^j - L)/(1 - (sigma+z)^5). Correct.
* prod_j (zeta^j - L)/(5L)^5 = (1 - L^5)/(5^5 L^5) = -q, an identity in A_L (no logarithm). Total:
  e^E T_z(a_0)/a_0 = q prod_k (5 sigma + k z)/((sigma+z)^5 - 1). Multiplying by
  T_z(sigma^5 - 1) a_0 gives b_z = e^E T_z((sigma^5-1) a_0) = q prod_k (5 sigma + k z) a_0. This is `eq:con-bz`.

**No circularity, no branches.** The proof uses: the Taylor data of `def:gam-c` (rational, in A_L),
`def:fsp-stirling` (definitions of Poch_M, St_m), `lem:fsp-stirfe` (Appendix A; its proof uses Bernoulli
numbers only), `lem:fsp-transl` (b) (Appendix A). `lem:con-form` (d), `lem:fsp-refl`, a branch of a
logarithm, the number lambda and anything proved after `thm:con-Pa` are not used. The only mention of
lambda (l.95, parenthesis) is a remark; see W5.

**Step 1 rewriting (chi A_St/A_St(L,0) -> chi a_0).** Correct: A_St(sigma,0)^2 = 5 sigma/(1 - sigma^5)
= -r(sigma), so A_St(sigma,0)/A_St(L,0) = (r(sigma)/r(L))^{1/2} with constant term 1, and
chi a_0 = (r(sigma)/r(L))^{1/2} b with b = chi exp(sum_k z^{2k-1} S_2k(sigma)) in A_sigma[[z]]
(sigma^{-1} = (sigma^5-1)/(5 delta), (1-sigma^5)^{-1} = -sigma/(5 delta)); S_0'' = -1/delta. So
`eq:fsp-generic` applies and gives `eq:con-shiftD` with c<chi> = SP°_L[S_0; chi a_0], which is the
definition of c<chi> in l.61.

**Uses of `eq:con-bz`.** One use only, l.103 of the same file (no other `\eqref`/`\cref` in `sections/`,
`main.tex`). It is used as b_z = q prod_k(5 sigma + kz) a_0 in A_L[[t,z]], with `lem:fsp-transl` (b)
(p = z, amplitude a = (sigma^5-1) a_0 in A_L[[t]][[z]], phase with S_0''(L) = -1/delta(L) a unit):
c<q prod_k(5 sigma+kz)> = SP°[S_0; b_z] = SP°[S_0; a] = c<sigma^5 - 1>. This is the form now proved, and it is
what Step 2 needs (q in A_L, A_L[z]-linearity of chi -> c<chi>).

## 2. Exact series check

Scripts (Python standard library, `fractions`, run with `py -I -B`):
`<user home>\AppData\Local\Temp\claude\d7c1b608d03d194a\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\ref4\conbz.py`
and `<user home>\AppData\Local\Temp\claude\d7c1b608d03d194a\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\ref4\control.py`.
Truncated power series in (t,z) over Q, exact.

| Identity | What was checked | Result |
|---|---|---|
| St_m = Poch_m St_0 (`lem:fsp-stirfe`) | m = -6..7, all coefficients through w^30 | holds |
| `eq:con-factor` | (alpha, rho_0) = (5,10), (-1,3), (-1,-7/2), (2,1/3), (-3,5), (1,2), (0,4): all monomials t^i z^j with i+j <= 8; (5,15/2), (-1,-2/3): i+j <= 13 | holds |
| `eq:con-bz` | built from the Taylor data as printed in `def:gam-c` (p_k from the generating function), at L = 2, 3, -2, 1/2, -3/2, 5/3, 7, -1/3, 4/7, -5, 9/4, 11: all monomials i+j <= 8; at L = 2, 3/5, -4: i+j <= 13 | holds |
| Controls | b_2 set to 0; b_1 perturbed | both identities fail, as they must |

Limits of the check: `eq:con-bz` was tested at rational values of L, not symbolically in L; `eq:con-factor`
with rational rho_0 (the identity is universal in rho_0, so this loses nothing but is not a symbolic check).
The hand derivation of section 1 is the proof; the computation is a confirmation.

## 3. Findings (all severity: wording)

**W1. T_z is defined on A_L[[t,z]] but applied in k[[t,z]].** l.74 defines T_z as an endomorphism of
`\AL\fps{t,z}`; l.85-93 apply it to l_rho, U_rho, m_rho, w, which lie in k[[t,z]] only. Also the two
facts used silently (T_z commutes with exp and log) deserve half a line.
* Location: `s08_connection_4.tex` l.83.
* Old: `in which $5L$ and the $\zeta^j-L$ are units; we prove \eqref{eq:con-bz} there.`
* New: `in which $5L$ and the $\zeta^j-L$ are units; we prove \eqref{eq:con-bz} there. We write $\mathsf T_z$ also for the endomorphism $f(t,z)\mapsto f(t+z,z)$ of $\Bbbk\fps{t,z}$; it is continuous and preserves the ideal $(t,z)$, so it commutes with $\exp$ and $\log$.`

**W2. The tensor product, the injectivity, and the letter `\Bbbk`.** The tensor product is unadorned and the
reason for injectivity is not given. `\Bbbk` denotes a field of characteristic 0 elsewhere
(`s03_gauge_2.tex` l.9, `a_formalsp_5.tex` l.88); here it is a ring that is not a field, and it is missing
from the list of local letters.
* Location: `s08_connection_4.tex` l.83.
* Old: `The ring $\AL\fps{t,z}$ embeds in $\Bbbk\fps{t,z}$, $\Bbbk:=\AL\otimes\bbQ(\zeta)$, in which`
* New: `The ring $\AL\fps{t,z}$ embeds in $\Bbbk\fps{t,z}$, where $\Bbbk:=\AL\otimes_{\bbQ}\bbQ(\zeta)=\bbQ(\zeta)[L,L^{-1},(L^5-1)^{-1}]$ (in this proof only; $\bbQ(\zeta)$ is free over $\bbQ$, so $\AL\to\Bbbk$ is injective), in which`
* Location: l.61. Old: `$\mathsf m_\rho$, $w$ are local to this proof.` New:
  `$\mathsf m_\rho$, $w$, $\Bbbk$ are local to this proof.`

**W3. `\rho` stands for the polynomial and for its expansion.** The display uses `\rho^{1-2k}` and `z/\rho`,
which need rho(L+t) to be a unit of k[[t]]; this is true but not said.
* Location: l.83.
* Old: `so that $\rho(L+t)=\rho_0+\alpha t$, and put`
* New: `so that $\rho(L+t)=\rho_0+\alpha t$ is a unit of $\Bbbk\fps t$, again written $\rho$, and put`

**W4. `w^{-1}` in k[[t,z]].** w is not a unit of k[[t,z]], so `(w^{-1}+\alpha)\log(1+\alpha w)` is the
power series, as in `def:fsp-stirling`.
* Location: l.93.
* Old: `$(\mathsf T_z\mathsf U_\rho-\mathsf U_\rho)/z=\alpha\,\mathsf l_\rho+(w^{-1}+\alpha)\log(1+\alpha w)-\alpha$, while`
* New: `$(\mathsf T_z\mathsf U_\rho-\mathsf U_\rho)/z=\alpha\,\mathsf l_\rho+(w^{-1}+\alpha)\log(1+\alpha w)-\alpha$ (here $w^{-1}\log(1+\alpha w)$ is the power series $\sum_{n\geq1}(-1)^{n+1}\alpha^nw^{n-1}/n$, as in \cref{def:fsp-stirling}), while`

**W5. The parenthesis about lambda.** In this proof S_0 is the formal phase in t^2 A_L[[t]] given by the
Taylor data (l.61), and each U_rho lies in t^2 k[[t]]: there are no linear terms to cancel, and lambda
does not occur. The parenthesis is true for the holomorphic function of `def:gam-c` but is not the
argument here, and it is the only place where the proof seems to touch the branch data. The comparison
is coefficient by coefficient with the printed Taylor data.
* Location: l.95.
* Old: `By \cref{def:gam-c}, $\mathsf S_0(L+t)-\mathsf S_0(L)=\mathsf U_{\rho_{\mathrm N}}+\sum_j\mathsf U_{\rho_j}$ (the terms linear in $t$ cancel by the definition of $\lambda$, that is, $\mathsf S_0'(L)=0$), $\mathsf S_1(L+t)-\mathsf S_1(L)=`
* New: `Since $[t^{k+1}]\,\mathsf U_\rho=(-1)^{k+1}\alpha^{k+1}\rho_0^{-k}/(k(k+1))$ and $[t^k]\,\mathsf l_\rho=(-1)^{k+1}\alpha^k\rho_0^{-k}/k$ for $k\geq1$, the Taylor data of \cref{def:gam-c} give, coefficient by coefficient, $\mathsf S_0(L+t)-\mathsf S_0(L)=\mathsf U_{\rho_{\mathrm N}}+\sum_j\mathsf U_{\rho_j}$, $\mathsf S_1(L+t)-\mathsf S_1(L)=`

**W6. `a_0` is declared local to the proof of (i) but used in the proof of (ii).** (Already noted as doubt 9
of `APPLIED_1.md`.)
* Location: `s08_connection_5.tex` l.4.
* Old: `and the amplitude is $a_0$, both with Taylor data in $\AL$, mapped to $A$ by $\jmath$. The letters $\phi$, $J$, $a$ are local to this proof.`
* New: `and the amplitude is $a_0:=A_{\mathrm{St}}/A_{\mathrm{St}}(L,0)$, both with Taylor data in $\AL$, mapped to $A$ by $\jmath$. The letters $\phi$, $J$, $a$, $a_0$ are local to this proof.`

## 4. Points checked and found in order (no change needed)

* `\alpha\in\bbZ` is what `def:fsp-stirling` requires for Poch_alpha and St_alpha; alpha = 5 and -1 are used.
* The sign pattern -(1/2) l_rho in `eq:con-factor` against +(1/2) l_N in log a_0 is handled by the extra term.
* `lem:fsp-transl` (b) is applied within its hypotheses (p = z in z A[[z]], a in A[[t]]((z))), and its
  translate b_p is literally e^E T_z(a).
* The comparison with the audit's repair (`audit_2026-11-18/verify/E1.md`, seven logarithm symbols): the
  note's version takes the difference (T_z - 1) before any symbol appears, so only series without constant
  term occur. It proves the same identity with less machinery, and it treats the first Mellin-Barnes
  summand (factors Gamma((zeta^j - sigma)/z), alpha = -1) directly, so the passage through
  `lem:con-form` (d) disappears.
* `lem:con-form` (d) is no longer cited in the proof of `thm:con-Pa`; it is still cited in
  `a_formalsp.tex` l.6 and `b_signs_2.tex` l.69 as a statement, which is consistent.
