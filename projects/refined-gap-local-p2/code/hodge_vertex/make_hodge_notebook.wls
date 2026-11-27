(* make_hodge_notebook.wls -- assemble the runnable notebook localp2CY5_hodge_genus3.nb (project folder + copy here) from the
   validated code files (inlined as Input cells) and the notebook driver cells nb_cell_*.m. ASCII output. *)
SetDirectory[DirectoryName[$InputFileName]];
rd[f_] := Import[f, "Text"];
strip[s_, lines_List] := StringReplace[s, Table[l ~~ ("\r\n" | "\n") -> "", {l, lines}]];
codeCell[s_] := Cell[StringTrim[s], "Input"];
(* split a driver file into cells: a new cell starts at every comment line *)
splitCells[s_] := Module[{ls = StringSplit[s, "\n"], groups = {}, cur = {}},
   Do[If[StringStartsQ[l, "(*"] && cur =!= {}, AppendTo[groups, cur]; cur = {}]; AppendTo[cur, l], {l, ls}];
   If[cur =!= {}, AppendTo[groups, cur]];
   codeCell[StringRiffle[#, "\n"]] & /@ Select[groups, StringTrim[StringRiffle[#, "\n"]] =!= "" &]];
txt[s_] := Cell[s, "Text"];
sec[s_] := Cell[s, "Section"];
sub[s_] := Cell[s, "Subsection"];

cells = Flatten[{
   Cell["Local P2 x C2: Hodge-vertex graph sum on the fibrewise torus, genus 2, 3 (and 4 from the cluster)", "Title"],
   txt["Successor of localp2CY5_genus3.nb. The T_F-equivariant Gromov-Witten theory of Z = K_P2 x C^2 is reconstructed from three copies of the Hodge CohFT (the equivariant GW theory of C^3, vertex classes prod_a Lambda_g^vee(w_a)/w_a) with the e_T weights of the normal bundle O(-3) + O + O of the fixed P^2, w = {eps1+eps2, -eps1, -eps2}, instead of three copies of the GW theory of the point: Omega_Z = R'.(Hodge_w)^{+3} with R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] the scalar Bernoulli series (h_w . point = Hodge_w by Mumford's formula / Givental's quantum Riemann-Roch). R_Z is the R-matrix of localp2CY5_genus3.nb (identity column cPP_z6.m, rows by the T_R ansatz), so the totals F_{2,0}, F_{3,0} are unchanged (group-action property) while their distribution over the stable graphs changes. Evaluate top to bottom. Data files (in the notebook directory or its subfolder hodge_vertex/): cPP_z6.m, StableGraph-2.bin, StableGraph-3.bin, ref/F20_F30_y.m, hv_orbits_g3.mx (decorated graphs of genus 3, recomputed if absent), hv3_samples_B.mx (the 780 genus-3 sample records of the stored run; delete it to recompute, several hours on one kernel), hodge_vertex_genus2_results.m (stored genus-2 run for the regression)."],
   txt["Conventions. Vertices carry idempotent labels i in {0,1,2}; on T_F the canonical coordinates are LLL[i][y] = -alpha/(1 + 3 omega^i t), t = y^(1/3), alpha = (eps1+eps2)/3. Graph factors as in localp2CY5_genus3.nb (Dilaton[i,n], Leaf, EdgeTerm, sqrtDelta[i] = Sqrt[eps1 eps2/e3](3 e3 + e2 LLL[i])) but built from the rows of R' = R_Z h^(-1) (coefficient lists of R_Z times the coefficient list of h^(-1)), and the vertex factor sqrtDelta^(2g_v-2+|H(v)|) <prod tau> is replaced by sqrtDelta^(2g_v-2+|H(v)|) Int_{M_{g_v,n_v}} prod_a Lambda^vee_{g_v}(w_a)/w_a^{g_v} prod psi^a, so every psi assignment of degree at most the dimension contributes (the lambda classes make up the difference). A decorated graph is an orbit of labelled psi/dilaton decorations under Aut(Gamma); each stable graph contributes 1/|Aut(Gamma)| times the sum over labelled decorations and 3^V colourings, with 1/nd_v! for the dilaton leaves."],
   txt["Hodge integrals (genus <= 3). The normalised class prod_a Lambda_g^vee(w_a)/w_a^g is expanded in lambda monomials with the Mumford relations lambda_1^2 = 2 lambda_2, lambda_2^2 = 2 lambda_1 lambda_3, lambda_3^2 = 0; the integrals against psi classes are Int lambda_1^k psi^a by k-fold Mumford's formula lambda_1 = (kappa_1 - Sum psi_i + delta)/12 (kappa_1 = pi_* psi^2, <kappa_1 X> = <tau_2 X>, lambda_1 restricts to the components of the boundary), the lambda_g formula, the lambda_g lambda_{g-1} closed form (2g-3+n)! |B_2g|/(2^(2g-1) (2g)! prod (2a_i-1)!!) for a_i >= 1 with zeros removed by the string equation, and lambda_1 lambda_2 lambda_3 on M_3 = |B_4 B_6|/(2 4 6 4!). Checks: known values, the recursion against the closed forms (lambda_1^5 = 16 lambda_2 lambda_3, lambda_1^6 = 16 lambda_1 lambda_2 lambda_3 with psi insertions, lambda_1^6 on M_3 = 1/90720), string and dilaton equations, and the R-action of the scalar h on the point theory in the PPZ conventions reproducing the normalised Hodge class on M_2 (symbolically) and M_3 (42 graphs, exact rational weights)."],
   txt["Numerics. As in localp2CY5_genus3.nb: rational sample points (t0, x0 = eps1/eps2) with 350-digit arithmetic, C1i = 1/c1 and C1p = c1' symbolic so that the XX-dependence is separated; exact rational results by Rationalize, reconstructed as polynomials in XX with coefficients rational in y (denominators y^a (1+27y)^b) and homogeneous rational functions of (eps1, eps2) with denominators (eps1 eps2)^a (eps1+eps2)^b, every fit verified on spare points. The sum over the decorated graphs is compared with F_{2,0}, F_{3,0} at every point and symbolically; the total conifold polar part is compared with the gap coefficient."],
   sec["Core routines"],
   sub["Stable graphs, decorations, Witten-Kontsevich numbers (graphs.m)"], codeCell[rd["graphs.m"]],
   sub["Decorated graphs as Aut-orbits (orbits.m)"], codeCell[rd["orbits.m"]],
   sub["R-matrix rows, colour sums, fits, conifold frame (lp2_lib.m)"], codeCell[strip[rd["lp2_lib.m"], {"Get[\"graphs.m\"]; Get[\"orbits.m\"];"}]],
   sub["Hodge integrals, Hodge vertices, R' = R_Z h^(-1) (hv_lib3.m)"], codeCell[strip[rd["hv_lib3.m"], {"Get[\"lp2_lib.m\"];"}]],
   sub["Genus-general Hodge evaluator used by the cluster package for genus 4 (hv_lib4.m): every lambda class as a polynomial in lambda_1 and ch_3, integrals by Mumford's formula"], codeCell[strip[rd["hv_lib4.m"], {"Get[\"lp2_lib.m\"];"}]],
   sec["Setup"], splitCells[rd["nb_cell_setup.m"]],
   sec["Checks of the Hodge integrals and of the scalar R-action"], splitCells[rd["nb_cell_checks.m"]],
   sec["Genus 2 (regression against the stored run)"], splitCells[rd["nb_cell_genus2.m"]],
   sec["Genus 3"], splitCells[rd["nb_cell_genus3.m"]],
   sec["Display"], splitCells[rd["nb_cell_display.m"]],
   sec["Genus 4 (cluster results)"], splitCells[rd["nb_cell_genus4.m"]]
   }];
nb = Notebook[cells, WindowSize -> {1400, 900}, StyleDefinitions -> "Default.nb"];
header = "(* Content-type: application/vnd.wolfram.mathematica *)\n\n(*** Wolfram Notebook File ***)\n(* http://www.wolfram.com/nb *)\n\n(* Beginning of Notebook Content *)\n";
out = header <> ToString[nb, InputForm] <> "\n(* End of Notebook Content *)\n";
Export["localp2CY5_hodge_genus3.nb", out, "Text"];
CopyFile["localp2CY5_hodge_genus3.nb", FileNameJoin[{ParentDirectory[Directory[]], "localp2CY5_hodge_genus3.nb"}], OverwriteTarget -> True];
(* validation: read back, count cells, parse every Input cell *)
nb2 = Get["localp2CY5_hodge_genus3.nb"];
inputs = Cases[nb2, Cell[s_String, "Input"] :> s, Infinity];
bad = Select[inputs, ! SyntaxQ[#] &];
Print["notebook written: ", Length[Cases[nb2, Cell[__], Infinity]], " cells, ", Length[inputs], " Input cells, syntax errors in ", Length[bad], If[bad =!= {}, ": " <> StringTake[First[bad], Min[200, StringLength[First[bad]]]], ""]];
