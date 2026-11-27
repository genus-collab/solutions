(* make_hvc_notebook.wls -- assemble the runnable notebook resolvedconifold_hodge.nb (project folder + copy here) from the validated
   code files (inlined as Input cells) and the notebook driver cells nb_hvc_*.m. ASCII output. *)
SetDirectory[DirectoryName[$InputFileName]];
rd[f_] := Import[f, "Text"];
stripGet[s_] := StringReplace[s, StartOfLine ~~ "Get[\"" ~~ Except["\n"] ... ~~ "\n" -> ""];   (* remove the Get[...] lines (the code is inlined) *)
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
   Cell["Resolved conifold x C2: Hodge-vertex graph sum on the fibrewise torus, genus 2 and 3", "Title"],
   txt["Companion of localp2CY5_hodge_genus3.nb. The equivariant Gromov-Witten theory of X x C^2, X = O(-1) + O(-1) -> P^1, on the fibrewise torus with opposite weights on each O + O(-1) pair (alpha0 = 0, alpha1 = eps1, alpha2 = eps2 in the conventions of resolvedconifold_genus3.nb: tangent weights {0, eps1, eps2, -eps1, -eps2} along the fixed P^1) is reconstructed from two copies of the SAME triple Hodge class as for local P2, the equivariant GW theory of C^3 with weights w = {eps1+eps2, -eps1, -eps2} (mode B; mode A is the opposite sign), instead of two copies of the GW theory of the point: Omega_Z = R'.(Hodge_w)^{+2} with R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] (h_w . point = Hodge_w by Mumford/QRR; a scalar, so the totals F_{g,0} are unchanged by the group-action property). R_Z is the conifold R-matrix of resolvedconifold_genus3.nb built from the Wick expansion of the superpotential; at alpha0 = 0, alpha_i = eps_i it needs no QRR normalisation constants (they are h_{(eps1, eps2, -eps1, -eps2)} = 1, the Bernoulli series of the actual normal bundle O(-1)+O(-1)+O+O of P^1, whose rank-4 Hodge class is trivial by c(E) c(E^vee) = 1). Everything is exact: Q symbolic, eps1 = x symbolic, eps2 = 1; polar parts are taken at Q = 1 in t = Log Q (the gap singularity P_g Li_{3-2g}(Q)/(eps1 eps2), where one canonical coordinate runs into the logarithmic singularity of the superpotential). Evaluate top to bottom. Data files (in the notebook directory or its subfolder hodge_vertex_conifold/): StableGraph-2.bin, StableGraph-3.bin, ref/resolvedconifold_results.m (F_{1,1}, F_{2,0}, F_{3,0} on the general torus), ref/conifold_graphterms_genus2_3.m (point-vertex terms of 2026-09), the stored hvc_results_g3_*.m (genus 3 is read from them unless recompute3 = True), and for the comparison ../hodge_vertex/hodge_vertex_genus{2,3}_results.m."],
   txt["Conventions. Vertices carry idempotent labels i in {1, 2} (the two roots U_1, U_2 of U^2 - S U + P = 0, S = (eps1 - eps2)/eps1, P = -eps2/(Q eps1), elements of the quadratic algebra C[a, b] = a + b U); graph factors as in resolvedconifold_genus3.nb (Delta_i^(g_v - 1) per vertex, dilaton leaves (-1)^n R'_{0,i}[z^(n-1)], leaves (-1)^n R'_{1,i}[z^n], edge terms from the rows of R' and eta^-1), with the Witten-Kontsevich vertex <prod tau> replaced by Int_{M_{g_v,n_v}} prod_a Lambda^vee_{g_v}(w_a)/w_a^{g_v} prod psi^a, so every psi assignment of degree at most the dimension contributes. A decorated graph is an orbit of labelled psi/dilaton decorations under Aut(Gamma); each stable graph contributes 1/|Aut(Gamma)| times the sum over labelled decorations and 2^V colourings, with 1/nd_v! for the dilaton leaves. Every orbit value must have zero U-component (Galois symmetry), which is checked."],
   txt["Hodge integrals (genus <= 4, hv_lib4.m): every lambda class is a polynomial in lambda_1 and ch_3(E) (Newton's identities with ch_2 = ch_4 = 0), and Int lambda_1^k1 ch_3^k3 psi^a is computed by Mumford's formula for ch_1 and ch_3. Checks: known values up to genus 4, vanishing above the top degree, Mumford relations, closed forms, string and dilaton equations, and the R-action of the scalar h_w^-1 on the point theory in the PPZ conventions reproducing the normalised Hodge class on M_2 (symbolically) and M_3 (42 graphs, exact rational weights). The exponent s in R' = R_Z h_w^s is fixed by the requirement that the sum of the decorated graphs equals F_{g,0}: s = -1 for both signs of the weights, as for local P2."],
   sec["Core routines"],
   sub["Stable graphs, decorations, Witten-Kontsevich numbers (graphs.m)"], codeCell[rd["graphs.m"]],
   sub["Decorated graphs as Aut-orbits (orbits.m)"], codeCell[rd["orbits.m"]],
   sub["Conifold Frobenius structure, quadratic algebra, Wick R-matrix, graph factors and graph sum (rc_lib.m)"], codeCell[stripGet[rd["rc_lib.m"]]],
   sub["Hodge integrals, Hodge vertices, twisted rows R' = R_Z h^s, references (hvc_lib.m)"], codeCell[stripGet[rd["hvc_lib.m"]]],
   sec["Setup"], splitCells[rd["nb_hvc_setup.m"]],
   sec["Checks"], splitCells[rd["nb_hvc_checks.m"]],
   sec["Genus 2"], splitCells[rd["nb_hvc_genus2.m"]],
   sec["Genus 3"], splitCells[rd["nb_hvc_genus3.m"]],
   sec["Comparison with local P2"], txt["The decorated-graph-by-decorated-graph comparison with local P2 in the same frame (genus 2, 3, 4) belongs to the local P2 project since 2026-11-27 (its folder hodge_vertex/compare_conifold/, with frozen copies of the mode-B results hvc_results_g<g>_Bm.m of this project). This notebook no longer reads any file of that project."],
   sec["Display"], splitCells[rd["nb_hvc_display.m"]],
   sec["Genus 4 (cluster results)"], splitCells[rd["nb_hvc_genus4.m"]]
   }];
nb = Notebook[cells, WindowSize -> {1400, 900}, StyleDefinitions -> "Default.nb"];
header = "(* Content-type: application/vnd.wolfram.mathematica *)\n\n(*** Wolfram Notebook File ***)\n(* http://www.wolfram.com/nb *)\n\n(* Beginning of Notebook Content *)\n";
out = header <> ToString[nb, InputForm] <> "\n(* End of Notebook Content *)\n";
Export["resolvedconifold_hodge.nb", out, "Text"];
CopyFile["resolvedconifold_hodge.nb", FileNameJoin[{ParentDirectory[Directory[]], "resolvedconifold_hodge.nb"}], OverwriteTarget -> True];
(* validation: read back, count cells, parse every Input cell *)
nb2 = Get["resolvedconifold_hodge.nb"];
inputs = Cases[nb2, Cell[s_String, "Input"] :> s, Infinity];
bad = Select[inputs, ! SyntaxQ[#] &];
Print["notebook written: ", Length[Cases[nb2, Cell[__], Infinity]], " cells, ", Length[inputs], " Input cells, syntax errors in ", Length[bad], If[bad =!= {}, ": " <> StringTake[First[bad], Min[200, StringLength[First[bad]]]], ""]];
