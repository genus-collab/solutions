# ===== referee E0j_ref2 independent checks =====
my $fails=0; my $KM=$NS-3;
my @ET = ([K($et00),$K0,$K0],[$K0,K($et11),K($et12)],[$K0,K($et12),$K0]);
sub sg { ($_[0]%2)?-$ONE:$ONE }
sub solve { my ($A,$r)=@_; my $n = scalar(@$r); my @M = map { [ @{$A->[$_]}, $r->[$_] ] } 0..$n-1;
  for my $c (0..$n-1) { my $p = $c; $p++ while $p < $n && kz($M[$p][$c]); return undef if $p >= $n; @M[$c,$p] = @M[$p,$c];
    my $iv = kinv($M[$c][$c]); $M[$c] = [ map { kmul($iv,$_) } @{$M[$c]} ];
    for my $q (0..$n-1) { next if $q == $c; my $f = $M[$q][$c]; next if kz($f); $M[$q] = [ map { ksub($M[$q][$_], kmul($f,$M[$c][$_])) } 0..$n ]; } }
  [ map { $M[$_][$n] } 0..$n-1 ] }
# coefficient at s^e of a product of two series (direct convolution, own code)
sub cprod { my ($f,$g,$e)=@_; my $x=$K0; for my $e1 ($f->{v}..$e-$g->{v}) { my $u=sget($f,$e1); next if kz($u); $x=kadd($x,kmul($u,sget($g,$e-$e1))); } $x }
# P_ij(z,w) coefficient [z^p w^q] at s^e = -(-1)^(p+q) <Phat_i[p],Phat_j[q]>
sub Pc { my ($Ri,$Rj,$p,$q,$e)=@_; my $v=$K0; for my $a (0..2){for my $b(0..2){ next if kz($ET[$a][$b]);
  $v=kadd($v,kmul($ET[$a][$b],cprod($Ri->[$a][$p],$Rj->[$b][$q],$e))); }} ksc(-sg($p+$q),$v) }
# (R1) independent difference quotient: (z+w)E = P(z,w)-P(-w,w)
sub Eind { my ($Ri,$Rj,$e,$N)=@_; my %P; for my $p(0..$N+1){for my $q(0..$N+1-$p){$P{$p}{$q}=Pc($Ri,$Rj,$p,$q,$e)}}
  my %Dd; for my $p(0..$N+1){for my $q(0..$N+1-$p){ $Dd{$p}{$q}=$P{$p}{$q}; } }
  for my $p(0..$N+1){for my $q(0..$N+1-$p){ $Dd{0}{$p+$q}=ksub($Dd{0}{$p+$q}, ksc(sg($p),$P{$p}{$q})); } }
  my %E; for my $d(0..$N){ $E{$d}{0}=$Dd{$d+1}{0}; for (my $m=$d;$m>=1;$m--){ $E{$m-1}{$d-$m+1}=ksub($Dd{$m}{$d-$m+1},$E{$m}{$d-$m}); }
     unless (kz(ksub($Dd{0}{$d+1},$E{0}{$d}))) { $fails++; print "  R1 inconsistency d=$d\n"; } }
  \%E }
for my $i(0..2){for my $j(0..2){ my $bad=0;my $c=0; for my $e(0..$KM){ my $E=Eind($rows[$i],$rows[$j],$e,$NE);
  for my $m(0..$NE){for my $n(0..$NE-$m){$c++; $bad++ unless kz(ksub($E->{$m}{$n},sget($Ed{$i}{$j}[$m][$n],$e)));}} }
  print "(R1) independent difference quotient vs code Ed_$i$j, s^0..s^$KM: ",($bad?"MISMATCH $bad/$c":"EQUAL $c"),"\n"; $fails+=$bad; }}
for my $j(1,2){ my $E=Eind($rows[0],$rows[$j],-1,$NE); my $nz=0; for my $m(0..$NE){for my $n(0..$NE-$m){$nz++ unless kz($E->{$m}{$n})}}
  my $E2=Eind($rows[0],$rows[$j],-2,$NE); my $nz2=0; for my $m(0..$NE){for my $n(0..$NE-$m){$nz2++ unless kz($E2->{$m}{$n})}}
  print "(R1b) [s^-1] E_0$j independent: nonzero entries $nz; [s^-2]: $nz2\n"; $fails+=$nz+$nz2; }
# (R2) Pihat_0j(u) = delta Delta_0
for my $j(0..2){ my $bad=0; for my $e(-2..$KM){ for my $n(0..$NZ){ my $acc=$K0; for my $q(0..$n){ my $v=$K0;
  for my $a(0..2){for my $b(0..2){ next if kz($ET[$a][$b]); $v=kadd($v,kmul($ET[$a][$b],cprod($rows[0][$a][$n-$q],$rows[$j][$b][$q],$e))); }}
  $acc=kadd($acc,ksc(sg($q),$v)); } $acc=ksub($acc,sget($Delta[0],$e)) if $j==0 && $n==0; $bad++ unless kz($acc); } }
  print "(R2) Pihat_0$j - delta Delta_0 = 0 (s^-2..s^$KM, u^0..u^$NZ): ",($bad?"FAIL $bad":"OK"),"\n"; $fails+=$bad; }
# leaf fit with free basis
sub fitfree { my ($G,$ls,$const)=@_; my @bf; push @bf, sub{ $_[0]==0?$ONE:$ZERO } if $const; for my $l(@$ls){ my $ll=$l; push @bf, sub{ QQ($ll)**$_[0]/($E**($_[0]+1)) } }
  my $n=scalar(@bf); return 'nofit' if $n>=scalar(@$G); return (grep {!kz($_)} @$G) ? 'FAIL' : 'OK' if $n==0;
  my @A; my @r; for my $m(0..$n-1){ push @A,[ map { K($_->($m)) } @bf ]; push @r,$G->[$m]; }
  my $x=solve(\@A,\@r); return 'sing' unless defined $x;
  for my $m($n..$#$G){ my $v=$K0; for my $c(0..$n-1){ $v=kadd($v,ksc($bf[$c]->($m),$x->[$c])) } return 'FAIL' unless kz(ksub($v,$G->[$m])); } 'OK' }
# (R3) legs
for my $j(1,2){ my (%res,%ctl,%ctl2); for my $k(0..$KM){ for my $b(0..$NE){ my @G=map{ sget($Ed{0}{$j}[$_][$b],$k) } 0..$NE-$b;
   next if scalar(@G) < $k+3; $res{fitfree(\@G,[1..$k+1],0)}++; $ctl{fitfree(\@G,[1..$k],0)}++; $ctl2{fitfree(\@G,[2..$k+1],0)}++; } }
  print "(R3) legs f_{$j,b} in Lam_{k+1}, no const: ",join(' ',map{"$_=$res{$_}"} sort keys %res),"; control Lam_k: ",join(' ',map{"$_=$ctl{$_}"} sort keys %ctl),"; control drop lam_1: ",join(' ',map{"$_=$ctl2{$_}"} sort keys %ctl2),"\n";
  $fails += ($res{FAIL}//0)+($res{sing}//0); }
for my $j(1,2){ my %res; for my $k(0..$KM){ for my $b(0..$NE){ my @G=map{ sget($Ed{$j}{0}[$b][$_],$k) } 0..$NE-$b; next if scalar(@G) < $k+3; $res{fitfree(\@G,[1..$k+1],0)}++; } }
  print "(R3b) E_${j}0 legs (singular end second slot) in Lam_{k+1}(w): ",join(' ',map{"$_=$res{$_}"} sort keys %res),"\n"; $fails += ($res{FAIL}//0); }
# (R4) rank one at s^0
for my $j(1,2){ my $bad=0; for my $m(0..$NE){for my $n(0..$NE-$m){ my $x=sget($Ed{0}{$j}[$m][$n],0); my $y=ksc($ONE/($E**$m),sget($Ed{0}{$j}[0][$n],0)); $bad++ unless kz(ksub($x,$y)); }}
  my $nzc=0; for my $n(0..$NE){ $nzc++ unless kz(sget($Ed{0}{$j}[0][$n],0)) }
  print "(R4) E_0$j|s^0 = lam_1(psi) chi(w): ",($bad?"FAIL $bad":"OK"),"; chi nonzero coefs $nzc\n"; $fails+=$bad; }
# (R5) closed form with wrong variants
my %S; for my $k(-1..$KM){ my $Kk=$k+1; for my $a(0..2){ my @G=map{ ksc(sg($_),sget($rows[0][$a][$_],$k)) } 0..$NZ;
  my @A; my @r; for my $m(0..$Kk){ my @row=($m==0?K($ONE):$K0); push @row,K(QQ($_)**$m/($E**($m+1))) for 1..$Kk; push @A,\@row; push @r,$G[$m]; } $S{$k}[$a]=solve(\@A,\@r); } }
for my $variant ('true','noL','wsign','dropk1top'){ for my $j(1,2){ my $bad=0;my $cnt=0;
  for my $k(0..$KM){ for my $m(0..$NE){for my $n(0..$NE-$m){ my $pred=$K0;
    for my $k1(0..$k){ next if $variant eq 'dropk1top' && $k1==$k; my $k2=$k-$k1; for my $l(1..$k1+1){ my $coef=$K0;
      for my $q(0..$n){ my $pv=$K0; for my $a(0..2){for my $b(0..2){ next if kz($ET[$a][$b]); $pv=kadd($pv,kmul(kmul($ET[$a][$b],$S{$k1}[$a][$l]),ksc(sg($q),sget($rows[$j][$b][$q],$k2)))) }}
        my $r=$n-$q; my $sgn = $variant eq 'wsign' ? $ONE : sg($r); $coef=kadd($coef,ksc($sgn*QQ($l)**$r/($E**($r+1)),$pv)); }
      my $lf = $variant eq 'noL' ? $ONE : QQ($l); $pred=kadd($pred,ksc(-$lf*QQ($l)**$m/($E**($m+1)),$coef)); } }
    $cnt++; $bad++ unless kz(ksub($pred,sget($Ed{0}{$j}[$m][$n],$k))); }} }
  print "(R5) closed form variant '$variant' E_0$j: ",($bad?"FAIL $bad/$cnt":"EQUAL $cnt"),"\n"; $fails+=$bad if $variant eq 'true'; } }
# (R6) control for (M-d): break (O)
{ my @Rp = ( [ @{$rows[1][0]} ], $rows[1][1], $rows[1][2] ); $Rp[0][1] = sadd($Rp[0][1], sconst(K($ONE)));
  my $E10=Eind(\@Rp,$rows[0],0,$NE); my $E01=Eind($rows[0],\@Rp,0,$NE); my $bad=0; for my $m(0..$NE){for my $n(0..$NE-$m){ $bad++ unless kz(ksub($E10->{$n}{$m},$E01->{$m}{$n})) }}
  print "(R6) control: (O) broken, E_10(n,m) vs E_01(m,n) at s^0 mismatches: $bad (expect >0)\n"; }
# (R7) Delta_j at s=0
for my $j(1,2){ my $o=$om[$j]; my $one_m=ksub(K($ONE),$o); my $pred=kmul(K($e1*$e2*$E**3/3), kmul(kmul($o,$o), kinv(kmul($one_m,$one_m))));
  print "(R7) Delta_$j(0) = ",kstr(sget($Delta[$j],0))," pred ",kstr($pred)," ",(kz(ksub($pred,sget($Delta[$j],0)))?"OK":"FAIL"),"\n"; }
# (R8) gauge independence of E_0j and E_jj'
{ my $bad=0; for my $pr([0,1],[0,2],[1,1],[1,2],[0,0]){ my ($i,$j)=@$pr; for my $e(-1..$KM){ my $A=Eind($rowsTrue[$i],$rowsTrue[$j],$e,$NE); my $B=Eind($rows[$i],$rows[$j],$e,$NE);
  for my $m(0..$NE){for my $n(0..$NE-$m){ $bad++ unless kz(ksub($A->{$m}{$n},$B->{$m}{$n})) }}}}
  print "(R8) edges identical in c1=1 and c1=C1cf gauges: ",($bad?"FAIL $bad":"OK"),"\n"; $fails+=$bad; }
# (R9) lowest s-order of rows at colour 0 and colours 1,2 (true gauge): (R) says >= 0 at colours 1,2
for my $i(0..2){ my $v=99; for my $a(0..2){ for my $m(0..$NZ){ my $g=sstrip($rowsTrue[$i][$a][$m]); $v=$g->{v} if $g->{v}<$g->{N} && $g->{v}<$v; } } print "(R9) colour $i true-gauge rows: lowest s-power $v\n"; }
print "REF TOTAL failures: $fails\n";
