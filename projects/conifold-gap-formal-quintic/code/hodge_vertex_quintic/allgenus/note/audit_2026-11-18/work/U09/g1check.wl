(* Independent sanity check of thm:mod-tr at (h,n)=(1,1): direct TR on x=p^4/4+v/p, y=p, v=1,
   versus (-5)^0 * sum_a <tau_d H^a>_1 dxi_{a,d}, model correlators from the Givental graph sum
   with R = M0(z) (varsigma=1), PPZ unit-preserving action (translation T(z)=z1-zR^{-1}(z)1). *)
x[p_] := p^4/4 + v/p;
ord = 7;
sig = -e + Sum[a[k] e^k, {k, 2, ord}];
ser = Series[(x[c + sig] - x[c + e]) /. v -> c^5, {e, 0, ord + 1}];
sol = {};
Do[
  eq = SeriesCoefficient[ser, k] /. sol;
  s1 = Solve[eq == 0, a[k - 1]][[1]];
  sol = Join[sol, s1 /. sol], {k, 3, ord + 1}];
sigS = Simplify[sig /. sol];
xp[p_] := D[x[pp], pp] /. pp -> p;
integrand = (1/2) (1/(p0 - c - e) - 1/(p0 - c - sigS)) D[sigS, e]/((e - sigS)^3 xp[c + e]) /. v -> c^5;
res = SeriesCoefficient[Series[integrand, {e, 0, 0}], -1];
res = Simplify[res];
zeta = Exp[2 Pi I/5];
br[b_] := zeta^(-b);
omegaTR[pp_] := Sum[res /. {c -> br[b], p0 -> pp}, {b, 0, 4}];

(* model side *)
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[al[j] t^j, {j, 1, 4}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3]; F0c = -Tk[0, F4];
eqs = Table[SeriesCoefficient[F0c - F0, {t, 0, j}] == 0, {j, 1, 4}];
alsol = Solve[eqs, Table[al[j], {j, 1, 3}]][[1]];
Fs = Normal[Series[{F0, F1, F2, F3, F4} /. alsol, {t, 0, 2}]];
Print["Fs = ", Fs];
f1 = Table[SeriesCoefficient[Fs[[b + 1]], {t, 0, 1}], {b, 0, 4}];
eps = {1, -1, -1, 1, 1};
Psi0inv = Table[(1/5) eps[[b + 1]] zeta^(-i (b + 1)), {i, 0, 4}, {b, 0, 4}];
Rhat1 = Table[zeta^(i b) f1[[b + 1]], {b, 0, 4}, {i, 0, 4}];
M1 = Psi0inv . Rhat1;
Print["M1 symmetric? ", Chop[N[Max[Abs[M1 - Transpose[M1]]]]]];
sv = Table[zeta^(-i), {i, 0, 4}];  (* hat Delta_i^{1/2} = zeta^{4i} varsigma *)
tau0 = Table[(1/24) Sum[-M1[[j + 1, i + 1]] sv[[j + 1]], {j, 0, 4}] + (1/2) M1[[i + 1, i + 1]] sv[[i + 1]] + (1/24) sv[[i + 1]]^2 Sum[M1[[i + 1, k + 1]]/sv[[k + 1]], {k, 0, 4}], {i, 0, 4}];
tau0m = Table[(1/24) Sum[-M1[[j + 1, i + 1]] sv[[j + 1]], {j, 0, 4}] + (1/2) M1[[i + 1, i + 1]] sv[[i + 1]] - (1/24) sv[[i + 1]]^2 Sum[M1[[i + 1, k + 1]]/sv[[k + 1]], {k, 0, 4}], {i, 0, 4}];
tau1 = Table[sv[[i + 1]]/24, {i, 0, 4}];
xit[i_, p_] := -1/(sv[[i + 1]] (p - sv[[i + 1]]));
Dx[f_, p_] := p^2/(p^5 - 1) D[f, p];
omegaModel[pp_, t0_] := Module[{p}, (Sum[t0[[i + 1]] D[xit[i, p], p] + tau1[[i + 1]] D[Dx[xit[i, p], p], p], {i, 0, 4}]) /. p -> pp];
Do[
  pp = pt;
  Print["p0=", pt, "  TR=", N[omegaTR[pp], 20], "  model(PPZ sign)=", N[omegaModel[pp, tau0], 20], "  model(opposite dilaton sign)=", N[omegaModel[pp, tau0m], 20]],
  {pt, {3/10 + 7/10 I, 2 - I/3, -1.3 + 0.4 I}}];
