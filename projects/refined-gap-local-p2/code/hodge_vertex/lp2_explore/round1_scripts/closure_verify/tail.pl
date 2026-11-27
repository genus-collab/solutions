# ===== independent closure check (closure_verify) =====
my $theta = sub { my ($f)=@_; smul(sscal(K(QQ('1/3')), sadd($sS, sconst(K(-$ONE)))), sder($f)) };
my $mode = $ENV{MODE} // 'onshell';
my $c1; if ($mode eq 'onshell') { $c1 = sscal(K(-$ONE), $wc); } elsif ($mode eq 'lin') { $c1 = spoly(0,[K($ONE),K(QQ($ENV{Q}//'1/10'))]); } elsif ($mode eq 'one') { $c1 = sconst(K($ONE)); }
elsif ($mode eq 'shift') { # c1 = -wc * exp(q * Int 3 ds/(s-1)) approximated: c1 = -wc*(1-s)^(-3q)?? use theta log c1 = XX+q  -> c1 = -wc*(1-s)^{3q}; q integer only
  my $q = $ENV{Q} // 1; my $f = sconst(K($ONE)); $f = smul($f, spoly(0,[K($ONE),K(-$ONE)])) for 1..(3*$q); $c1 = smul(sscal(K(-$ONE),$wc), $f); }
my $ic1 = sinv($c1);
my $XXm = smul($theta->($c1), $ic1);   # XX := theta log c1, computed directly
print "mode $mode; XX = theta log c1: ", join(', ', map { kstr(sget($XXm,$_)) } 0..4), "\n";
if ($mode eq 'onshell') { my $d = ssub($XXm,$XX); print "  XX(direct) - XX(edges_sing): valuation ", sval($d), "\n"; }
my $iDel = sinv($Del); my $Dm1 = sadd($Del, sconst(K(-$ONE)));
my $d2 = smul(sscal(K($E/2), sadd(smul($ic1,$ic1), sconst(K(-$ONE)))), $iDel);
my $Rd = sadd(sadd(smul($Del,$theta->($XXm)), smul($Del, smul($XXm,$XXm))), sadd(smul($Dm1,$XXm), sscal(K(QQ('2/9')),$Dm1)));
print "  Rd = ", join(', ', map { kstr(sget($Rd,$_)) } 0..3), "  (valuation ", sval($Rd), ")\n";
my $c3 = smul($ic1, smul($ic1,$ic1));
my @A0 = ( [ sconst($K0), sconst(K($ONE)), sconst($K0) ],
           [ sconst($K0), smul($ic1,$d2), smul($c3,$iDel) ],
           [ sconst(K(-$E**3/27)), sscal(K(-$ONE), smul($c1, sadd(sadd(sscal(K($E),$d2), smul($Del,smul($d2,$d2))), sconst(K($E*$E/3))))),
             sscal(K(-$ONE), smul($ic1, sadd(sscal(K($E),$iDel), $d2))) ] );
# self-adjointness: G A0 symmetric, G = pairing matrix of edges_sing
{ my @G = ([$et00,0,0],[0,$et11,$et12],[0,$et12,0]); my $ok=1;
  for my $a (0..2) { for my $b (0..2) { my $x = sconst($K0); my $y = sconst($K0);
    for my $c (0..2) { $x = sadd($x, sscal(K(QQ($G[$a][$c])), $A0[$c][$b])) if $G[$a][$c] != 0; $y = sadd($y, sscal(K(QQ($G[$b][$c])), $A0[$c][$a])) if $G[$b][$c] != 0; }
    my $g = sstrip(ssub($x,$y)); $ok = 0 if $g->{v} < $g->{N}; } }
  print "  G A0 symmetric (to known order): ", ($ok?"YES":"NO"), "\n"; }
my $bad = 0; my $badpf = 0; my $badrow = 0; my $minN = 99;
for my $i (0..2) {
  my $oi = $om[$i]; my $vi; my $ivi;
  if ($i == 0) { $ivi = spoly(-1, [ K($ONE) ]); } else { $vi = spoly(0, [ ksub(K($ONE),$oi), $oi ]); $ivi = sinv($vi); }
  my $ti = spoly(0, [ ksc(QQ('-1/3'),$oi), ksc(QQ('1/3'),$oi) ]); my $tov = smul($ti,$ivi);
  my $Li = sscal(K(-$alpha), $ivi);
  my $Dop = sub { my ($F)=@_; sadd(sadd(smul(sscal(K(QQ('1/3')), sadd($sS, sconst(K(-$ONE)))), sder($F)), sscal(K(QQ('-1/3')),$F)), smul($tov,$F)) };
  my $Th = sub { my ($F)=@_; [ map { my $m=$_; my $x = smul($Li,$F->[$m]); $x = sadd($x, $Dop->($F->[$m-1])) if $m>=1; $x } 0..$NZ ] };
  my $zs = sub { my ($f,$F)=@_; [ map { $_ >= 1 ? smul($f,$F->[$_-1]) : sconst($K0) } 0..$NZ ] };   # z f F
  my $z2s = sub { my ($f,$F)=@_; [ map { $_ >= 2 ? smul($f,$F->[$_-2]) : sconst($K0) } 0..$NZ ] };
  my $mul = sub { my ($f,$F)=@_; [ map { smul($f,$_) } @$F ] };
  my $add = sub { my ($F,$G)=@_; [ map { sadd($F->[$_],$G->[$_]) } 0..$NZ ] };
  my $sub = sub { my ($F,$G)=@_; [ map { ssub($F->[$_],$G->[$_]) } 0..$NZ ] };
  my $P0 = $rows[$i][0];   # frame-B P0 from the edges_sing construction
  my $T1 = $Th->($P0); my $T2 = $Th->($T1); my $T3 = $Th->($T2);
  # PF: Del T3 + (E + (Del-1) z) T2 + (E^2/3 + (2/9)(Del-1) z^2) T1 + E^3/27 P0
  my $pfc = $ENV{PFCTRL} ? QQ('1/5') : QQ('2/9');
  my $pf = $add->($add->($mul->($Del,$T3), $add->($mul->(sconst(K($E)),$T2), $zs->($Dm1,$T2))),
                  $add->($add->($mul->(sconst(K($E*$E/3)),$T1), $z2s->(sscal(K($pfc),$Dm1),$T1)), $mul->(sconst(K($E**3/27)),$P0)));
  for my $m (0..$NZ) { my $g = sstrip($pf->[$m]); $minN = $g->{N} if $g->{N} < $minN; if ($g->{v} < $g->{N}) { $badpf++; print "  PF fail colour $i z^$m at s^$g->{v}\n" if $badpf < 4; } }
  # my rows
  my $v0 = $P0; my $v1 = $mul->($ic1,$T1);
  my $v2 = $mul->(smul($Del,$c1), $sub->($sub->($T2, $zs->($XXm,$T1)), $mul->($d2,$T1)));
  if ($mode eq 'onshell') { for my $m (0..$NZ) { for my $pr ([$v1,1],[$v2,2]) { my $g = sstrip(ssub($pr->[0][$m], $rowsTrue[$i][$pr->[1]][$m])); $badrow++ if $g->{v} < $g->{N}; } } }
  my @v = ($v0,$v1,$v2);
  for my $a (0..2) { my $lhs = $mul->($ic1, $Th->($v[$a])); my $rhs = $mul->($A0[$a][0],$v0); $rhs = $add->($rhs, $mul->($A0[$a][1],$v1)); $rhs = $add->($rhs, $mul->($A0[$a][2],$v2));
    $rhs = $sub->($rhs, $z2s->(smul($c1,$Rd), $v1)) if $a == 2 && !$ENV{NOANOM};
    for my $m (0..$NZ) { my $g = sstrip(ssub($lhs->[$m],$rhs->[$m])); $minN = $g->{N} if $g->{N} < $minN; if ($g->{v} < $g->{N}) { $bad++; print "  closure fail colour $i row $a z^$m at s^$g->{v}: ", kstr(sget($g,$g->{v})), "\n" if $bad < 6; } } }
}
print "PF residual (", ($ENV{PFCTRL}?"CONTROL 6y->(27/5)y":"true"), "): ", ($badpf?"FAILED in $badpf":"ZERO"), " all colours z^0..z^$NZ\n";
print "my rows v1,v2 vs edges_sing rowsTrue: ", ($mode eq 'onshell' ? ($badrow?"DIFFER($badrow)":"EQUAL") : "n/a"), "\n";
print "closure C1i Th v = (A0 - z^2 c1 Rd E21) v", ($ENV{NOANOM}?" [anomaly term OMITTED]":""), ": ", ($bad?"FAILED in $bad":"EXACT"), " (min known s-order ", $minN-1, ")\n";
