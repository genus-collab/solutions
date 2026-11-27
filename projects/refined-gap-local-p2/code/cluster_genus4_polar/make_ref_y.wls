(* make_ref_y.wls -- (run once, locally, from this directory) write ref/F20_F30_y.m = {F20, F30} in y, XX, eps1, eps2 from the
   stored results of the genus-3 project (F20_result.mx, F30_result.mx in the project folder), and check them against ref/F20_F30.m
   (whose first element is F20 in y and second element F30/(eps1 eps2) in LL = (1+27y)^(-1/3)). *)
SetDirectory[DirectoryName[$InputFileName]];
proj = ParentDirectory[Directory[]];
Get[FileNameJoin[{proj, "F20_result.mx"}]]; F20y = Fg;
Get[FileNameJoin[{proj, "F30_result.mx"}]]; F30y = Fg;
Put[{F20y, F30y}, "ref/F20_F30_y.m"];
old = Get["ref/F20_F30.m"];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
Print["F20 agrees with ref/F20_F30.m[[1]]: ", Together[F20y - old[[1]]] === 0];
ll = 2/5; y0 = (1 - ll^3)/(27 ll^3); pt = {y -> y0, XX -> 3/7, e1 -> 5, e2 -> 11};
Print["F30 at a point versus the LL-form element 2 (times eps1 eps2 and times 1): ", {Together[(F30y /. pt) - (old[[2]] /. LL -> ll /. pt) (e1 e2 /. pt)], Together[(F30y /. pt) - (old[[2]] /. LL -> ll /. pt)]}];
Print["written ref/F20_F30_y.m; LeafCounts ", LeafCount /@ {F20y, F30y}];
