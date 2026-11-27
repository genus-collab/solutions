# t_mumford.pl h : report 0(f)(4) with the referee engine.  Xi^(h-1)_h = X_a - h Xi^(h)_h + X_b with
#   X_a = Int K lambda_{h-1} prod lam,  X_b = Int K lambda_h Sum_{n odd} n! ch_n prod lam   (dilaton leaves as points, so K implicit),
# from c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum_{n odd} n! ch_n  (Mumford: prod (1 - alpha_i^2) = 1).
# Also the (G) lambda_g formula Xi^(h)_h(R) = b_h (J-1)...(J-N) and Harer-Zagier Int Lambda^vee(1) K = chi(M_{h,r}) for the referee engine.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my $h = shift; my $maxd = 3*$h - 3;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
my @lam = RefHodge::lambdas($h, $maxd);
my $sumch = {}; for my $i (0..$#RefHodge::CHL) { my $l = $RefHodge::CHL[$i]; next if $l > $maxd; my @e = (0) x $RefHodge::NCH; $e[$i] = 1; $sumch->{join(',', @e)} = fact($l); }
my $Xb = RefHodge::cmul($lam[$h], $sumch, $maxd);
my $LV1 = {}; for my $k (0..$h) { $LV1 = RefHodge::cadd($LV1, $lam[$k], ($k % 2) ? Q(-1) : Q(1)); }
# component files for Xi^(h-1), Xi^(h)
my %CF; for my $r (1..2) { my @C; open my $fh, '<', "$FindBin::Bin/ref_engine_h${h}_r$r.txt" or next; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
    for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $m], Q($cf)]; } } close $fh; $CF{$r} = \@C; }
sub comp { my ($k, $p) = @_; my $C = $CF{scalar @$p} or die; my $s = Q(0); for my $t (@{$C->[$k] // []}) { my ($e, $cf) = @$t; my $x = $cf->copy; for my $i (0..$#$p) { $x *= Q($p->[$i]) ** $e->[$i] if $e->[$i]; } $s += $x; } return $s; }
my $K = $h + 1; my @sn = map { ((($_ % 2) ? Q(-1) : Q(1)) / (Q(4) ** $_) / fact(2*$_ + 1)) } 0..$K; my @inv = (Q(1));
for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s -= $sn[$m] * $inv[$n - $m]; } $inv[$n] = $s; } my $bh = $inv[$h];
my @pts; for my $j (2..2*$h) { push @pts, [$j]; } for my $a (2..2*$h) { for my $b ($a..2*$h) { next if $a + $b > 2*$h + 1; push @pts, [$a, $b]; } } push @pts, [0], [1], [0, 3], [1/2 => ()] if 0;
my ($bad, $badG) = (0, 0);
for my $p (@pts) { my $D = 0; $D += $_ - 1 for @$p; my $J = 0; $J += $_ for @$p; my $N = 2*$h - 3 + scalar(@$p);
  my $xa = RefHodge::xi_value(RefHodge::xi($h, $lam[$h-1], $p)); my $xb = RefHodge::xi_value(RefHodge::xi($h, $Xb, $p)); my $xg = RefHodge::xi_value(RefHodge::xi($h, $lam[$h], $p));
  my $ok = (comp($h - 1, $p) == $xa - $h * $xg + $xb); $bad++ unless $ok;
  my $eg = $bh; $eg *= ($J - $_) for 1..$N; my $okg = ($xg == $eg) && (comp($h, $p) == $eg); $badG++ unless $okg;
  printf "h=%d R={%s} D=%d: identity %s; X_a %s, X_b %s, X_a + X_b %s; (G) %s%s\n", $h, join(',', @$p), $D, ($ok ? "holds" : "FAILS"), ($xa == 0 ? "= 0" : "!= 0"), ($xb == 0 ? "= 0" : "!= 0"),
    ($xa + $xb == 0 ? "= 0" : "!= 0"), ($okg ? "ok" : "FAIL"), ($D <= 2*$h - 3 ? "   [Lemma-A range]" : ""); }
my $chi = RefHodge::bern(2*$h) / (2*$h * (2*$h - 2));
for my $r (0..3) { my $v = RefHodge::xi_value(RefHodge::xi($h, $LV1, [(0) x $r])); my $e = $chi; $e *= (2 - 2*$h - $_) for 0..$r-1; printf "(HZ) referee engine genus %d r = %d: Int Lambda^vee(1) K = chi(M_{h,r}): %s\n", $h, $r, ($v == $e ? "OK" : "FAIL $v vs $e"); }
printf "SUMMARY h=%d: Mumford identity failures %d, (G) failures %d\n", $h, $bad, $badG;
