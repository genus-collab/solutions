# diag: for off-shell rows, fit sigma_a[s^e](z) = c0 + c1 z + Sum_{j<=K} b_j/(E + j z) (c1 z = -c1 psi: psi-linear remainder)
sub linfit { my ($f)=@_;
  for my $K (0..$NZ-3) { my @A = map { my $m = $_; [ ($m == 1 ? $ONE : $ZERO), (map { (QQ(-$_))**$m / $E**($m+1) } 1..$K), $f->[$m] ] } 1..$K+1;
    my $x = solve(\@A); my ($c1, @b) = @$x; my $c0 = $f->[0]; $c0 -= $b[$_-1]/$E for 1..$K;
    my $ok = 1; for my $m ($K+2..$NZ) { my $v = $ZERO; $v += $b[$_-1]*(QQ(-$_))**$m/$E**($m+1) for 1..$K; $ok = 0 if $v != $f->[$m]; }
    return ($c0, $c1, $K, $NZ-$K-1) if $ok; } return (); }
for my $set (['c1=1', $rows[0]], ['c1=C1cf', $rowsTrue[0]]) { my ($name,$R) = @$set;
  for my $a (2) { for my $e (-1..4) { my @f = map { realq(sget($R->[$a][$_],$e)) } 0..$NZ; my @r = linfit(\@f);
    print "  $name sigma_$a [s^$e]: ", (@r ? "c0 = $r[0], z-linear coefficient = $r[1], leaves K = $r[2] ($r[3] extra checks)" : "no fit"), "\n"; } } }
