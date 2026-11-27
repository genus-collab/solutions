# general_tail: basis {1, psi, lam_j}; exact difference-quotient identity valid with or without symplecticity / polynomial parts.
sub linfit2 { my ($f)=@_;
  for my $K (0..$NZ-3) { my @A = map { my $m = $_; [ ($m == 1 ? $ONE : $ZERO), (map { (QQ(-$_))**$m / $E**($m+1) } 1..$K), $f->[$m] ] } 1..$K+1;
    my $x = solve(\@A); my ($c1, @b) = @$x; my $c0 = $f->[0]; $c0 -= $b[$_-1]/$E for 1..$K;
    my $ok = 1; for my $m ($K+2..$NZ) { my $v = $ZERO; $v += $b[$_-1]*(QQ(-$_))**$m/$E**($m+1) for 1..$K; $ok = 0 if $v != $f->[$m]; }
    return { 0=>$c0, L=>-$c1, map { ($_ => $b[$_-1]) } 1..$K } if $ok; } return undef; }   # z-linear c1 z = -c1 psi
sub lamc { my ($i,$m)=@_; QQ($i)**$m/$E**($m+1) }   # [psi^m] lam_i(psi), i may be negative
sub qcoef { my ($al,$be,$m,$n)=@_;   # [psi^m psi'^n] of [f_al(psi) - f_al(-psi')] f_be(psi')/(psi+psi')
  return $ZERO if $al eq '0';
  my $fb = sub { my $n=shift; $be eq '0' ? ($n==0?$ONE:$ZERO) : $be eq 'L' ? ($n==1?$ONE:$ZERO) : lamc($be,$n) };
  if ($al eq 'L') { return $m == 0 ? $fb->($n) : $ZERO; }
  my $i = $al;
  if ($be eq '0') { return QQ($i)*lamc($i,$m)*lamc(-$i,$n); }
  if ($be eq 'L') { return lamc($i,$m)*(($n==0?$ONE:$ZERO) - $E*lamc(-$i,$n)); }
  my $l = $be; return QQ($i)*lamc($i,$m)*(QQ($i)*lamc(-$i,$n) + QQ($l)*lamc($l,$n))/($E*QQ($i+$l)); }
for my $set (['c1=1', $rows[0]], ['c1=C1cf', $rowsTrue[0]]) { my ($name,$R) = @$set; my (%cf, @va, @ea);
  for my $a (0..2) { my $v = 99; my $N = 99; for my $m (0..$NZ) { my $g = sstrip($R->[$a][$m]); $N = $g->{N} if $g->{N} < $N; $v = $g->{v} if $g->{v} < $v && $g->{v} < $g->{N}; }
    $va[$a] = $v; $ea[$a] = $N-1; for my $e ($v..$N-1) { my @f = map { realq(sget($R->[$a][$_],$e)) } 0..$NZ; my $h = linfit2(\@f); $cf{$a}{$e} = $h if $h; } }
  my $kmax = 99; for my $ab (keys %eta) { my ($a,$b) = split /,/, $ab; my $x = $va[$a] + $ea[$b]; $kmax = $x if $x < $kmax; }
  $kmax = $edN-1 if $kmax > $edN-1;
  for my $k (-2..$kmax) { my %C; my $miss = 0; my $lin = 0;
    for my $ab (keys %eta) { my ($a,$b) = split /,/, $ab;
      for my $e1 ($va[$a]..$ea[$a]) { my $e2 = $k-$e1; next if $e2 < $va[$b]; if (!$cf{$a}{$e1} || !$cf{$b}{$e2}) { $miss++; next; }
        for my $p (keys %{$cf{$a}{$e1}}) { for my $q (keys %{$cf{$b}{$e2}}) { my $v = $cf{$a}{$e1}{$p}*$cf{$b}{$e2}{$q}; next if $v == 0;
          $C{"$p|$q"} = ($C{"$p|$q"} // $ZERO) - $eta{$ab}*$v; $lin = 1 if ($p eq 'L' || $q eq 'L'); } } } }
    my $ok = 1; my $e00;
    for my $m (0..$NE) { for my $n (0..$NE-$m) { my $pr = $ZERO; for my $key (keys %C) { my ($p,$q) = split /\|/, $key; $pr += $C{$key}*qcoef($p,$q,$m,$n); }
      $e00 = $pr if $m == 0 && $n == 0; $ok = 0 if $pr != realq(sget($Ed00[$m][$n],$k)); } }
    print "  GENERAL $name [s^$k]: basis-pair identity vs code E_00 (all m+n <= $NE): ", ($ok ? "EQUAL" : "DIFFERENT"), ($miss ? " (missing fits: $miss)" : ""),
      "; psi-linear terms present: ", ($lin ? "YES" : "no"), "; E_00(0,0) = $e00\n"; } }
