Parity of the closed-form relations (main session, 2026-10-31). No TR anywhere.

Inputs
 - Route D, Theorem D5: the known conifold polar part (hypothesis (H), i.e. Brini-Schuler Prop. 3.1) gives one relation per level D at genus G,
       T_D(G) = Sum_{D(R)=D} (1/|Aut R|) prod_i ((-1)^{j_i} Pi_{j_i}^2 / j_i) Xi_G(R) = 0      (eps_j = Pi_j^2),
   given Lemma A in lower genera and for levels <= D-2 at genus G.
 - Route A, Corollary E (uses the KP hierarchy of the Marino-Vafa generating function): given Lemma A in lower genera and for levels < D at genus G,
       Xi_G(T) = (J-1)(J-2)...(J-r+1) Pi_{D+1} / prod_{t in T} Pi_t * Xi_G(D+1),     Pi_j = prod_{a=1}^{j-1} (j tau + a)/(j-1)!.

Result
 Substituting Corollary E into T_D gives T_D(G) = Pi_{D+1}(tau) kappa_D Xi_G(D+1), with
       kappa_D = Sum_{T at level D} (1/|Aut T|) (J-1)...(J-r+1) prod_i (-1)^{t_i} Pi_{t_i}/t_i  =  [y^D] ( tau y + log((e^y - 1)/y) ).
 Proof: kappa_D = (1/D) [y^D] (1 - alpha(y))^(-D) with alpha(y) = Sum_{t>=2} (-1)^t (Pi_t/t) y^(t-1); by Lagrange inversion this is
 [y^D] log(w(y)/y) for the inverse of y = w - Sum_{t>=2} (-1)^t (Pi_t/t) w^t. Since Sum_{t>=1} (Pi_t/t) z^t = log B_{tau+1}(z)
 (generalised binomial series, B = 1 + z B^(tau+1)), y = -log B_{tau+1}(-w), so w(y) = e^(tau y) (e^y - 1).
 Hence kappa_1 = tau + 1/2, kappa_D = B_D/(D D!) for even D (1/24, -1/2880, ...), and kappa_D = 0 for every odd D >= 3.

Consequence
 Given Lemma A in lower genera (and KP), the closed form forces every even level at genus G. At odd levels D >= 3 it is
 automatically satisfied and gives nothing. So Lemma A(G) reduces to the one-leaf values Xi_G(4), Xi_G(6), ..., Xi_G(2G-2)
 (the odd levels), given Lemma A in genera < G. Level 2 example: T_2 and Corollary E together give Xi_G(3) = Xi_G(2,2) = 0 for all G.

Check: kappa_parity.pl, exact rationals (Math::BigRat), D <= 10, tau = 0, 2/7, 3/2, -5/3: every value matches; output kappa_parity_out.txt.
