# Mixed-scale vertex functionals (2026-10-30)

Setting: CY triple Hodge class with weights (c, -x, -1), c = 1 + x (mode B, eps1 = x, eps2 = 1), normalised as in hv_lib4.m hodgeNormU4.
Leaves lam_j(psi) = 1/(c - j psi) (j = 0: the constant 1/c; j = 1: phi of Theorem C), dilaton leaves T_c = -psi^2/(c - psi).
Xi_h(j_1..j_r) := Sum_m 1/m! < lam_{j_1} ... lam_{j_r} T_c^m >_h = d_{lam_j1}...d_{lam_jr} F_h(q_0), q_0 = c - c^2 lam_1 (dilaton-shifted input).

## Single probe of arbitrary scale cp (hvc_mixed_scale.wls)
Y[h,n](cp) = < phi_cp phi_c^n e^T >_h, phi_cp = 1/(cp - psi). Generic cp: NONZERO. Exact zeros: cp = c/k, k = 1..2h-2.
  genus 1: no zeros; genus 2: factor (cp - c)(2cp - c); genus 3: (cp - c)(2cp - c)(3cp - c)(4cp - c); remaining factor of degree h in cp.
  Y[h,n] = Y[h,0] prod_{i=0}^{n-1} (c/cp + i)/c^2 (string equation; this is Y[h,0]/(c cp) only for n = 1). The product never vanishes at cp = c/k, k >= 1, so the zeros are the same for every n. [Corrected 2026-10-31; the earlier version said Y[h,0]/(c cp)^n.]

## Conjecture M (hvc_multiprobe.wls, output hvc_multiprobe_out.txt)
For h >= 2 and integers j_i >= 0:  Xi_h(j_1..j_r) = 0  <=>  1 <= J := Sum j_i <= 2h - 3 + r   (equivalently J >= 1 and D := Sum (j_i - 1) <= 2h - 3).
Verified on all lattice points with J <= 2h - 1 + r: genus 2 r = 1, 2, 3; genus 3 r = 1, 2, 3 (zero exactly in the predicted range, nonzero outside).
For full x the polynomial Xi_h(j) does NOT factor through J (genus 2, r = 2 is irreducible-looking): the vanishing is a lattice phenomenon.

## Proved pieces
- string + dilaton: Xi(1, j) = (J/c^2) Xi(j);  Xi(0, j) = ((D - 2h + 2)/c^2) Xi(j)  (D, J of the remaining leaves).
- lambda_h limit (x -> 0, all genera): lambda_g formula <lambda_g prod tau_a>_g = b_g <tau_0^{2g} prod tau_a>_0 reduces Xi to genus-0 KW:
  Xi_lead = b_h N! [s^N] e^{J s} (s/(e^s - 1))^{N+1} = b_h (J-1)(J-2)...(J-N),  N = 2h - 3 + r  (Norlund B_N^{(N+1)}(J));
  proof of the identity: w = e^s - 1 turns it into Res_w (1+w)^{J-1}/w^{N+1} = C(J-1, N). Checked: ratio = 7/5760 (g = 2), 31/967680 (g = 3) = b_g.

## Why it matters
Fact G: the order-t^k corrections of the singular column have poles only at z = -c/j, j <= k+1, i.e. leaves lam_j with D = j - 1 <= k.
Conjecture V at genus h needs vanishing through order t^(2h-3), i.e. total D <= 2h - 3: exactly the range of Conjecture M.
If (i) Conjecture M holds and (ii) every order-t^k correction of the normalised dilaton input rho(t) q_t, of the legs and of the loops is a
combination of leaf products of total degree D <= k with J >= 1, then Conjecture V, hence (Reduction Lemma) the all-genus theorem, follows.
