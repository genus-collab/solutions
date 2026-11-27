package Pair;
# Pairing of Hodge classes (ch-polynomials, RefHodge format) with the K-twisted one-leaf moments refdata/moments_h<h>.txt:
#   P_C(j) = Int_{Mbar_{h,1}} C K lam_j = Sum_mu C_mu Sum_a S_mu(a) j^a      (lam_j = 1/(1 - j psi)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
our %MOM;   # MOM{h}{mu}[a]
sub load { my $h = shift; return $MOM{$h} if $MOM{$h}; my %M; open my $fh, '<', "$FindBin::Bin/refdata/moments_h${h}.txt" or die "no moments for h=$h";
  while (my $l = <$fh>) { chomp $l; next unless $l =~ /\|/; my ($mu, $rest) = split /\s*\|\s*/, $l, 2; $rest //= ''; my @v;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($a, $c) = split /:/, $t; $v[$a] = Q($c); } $M{$mu} = \@v; }
  close $fh; return $MOM{$h} = \%M; }
# j-polynomial (array of BigRat coefficients) of a ch-polynomial class
sub jpoly { my ($h, $C) = @_; my $M = load($h); my @p = map { Q(0) } 0..3*$h;
  for my $mu (keys %$C) { my $c = $C->{$mu}; next if $c == 0; my $v = $M->{$mu}; next unless $v; for my $a (0..$#$v) { next unless defined $v->[$a]; $p[$a] += $c * $v->[$a]; } }
  pop @p while @p > 1 && $p[-1] == 0; return \@p; }
sub peval { my ($p, $x) = @_; my $s = Q(0); my $xx = Q(1); for my $c (@$p) { $s += $c * $xx; $xx *= $x; } return $s; }
sub pstr { my $p = shift; my @t; for my $a (0..$#$p) { push @t, "$a:" . $p->[$a] if $p->[$a] != 0; } return @t ? join(' ; ', @t) : '0'; }
# divide polynomial p(j) by prod_{l in L} (j - l); returns (quotient, remainder-is-zero flag)
sub pdivroots { my ($p, @L) = @_; my @q = @$p; my $ok = 1;
  for my $l (@L) { my $lq = Q($l); my $n = $#q; if ($n < 0) { last } my @r = map { Q(0) } 0..$n-1; my $acc = Q(0);
    for (my $i = $n; $i >= 0; $i--) { $acc = $acc * $lq + $q[$i]; $r[$i-1] = $acc if $i >= 1; }
    $ok = 0 if $acc != 0; @q = @r; @q = (Q(0)) unless @q; }
  return (\@q, $ok); }
# classes
sub one { return RefHodge::cone(); }
sub maxd { my $h = shift; return $h == 1 ? 1 : 3*$h - 3; }
sub lambdas { my $h = shift; return RefHodge::lambdas($h, maxd($h)); }
sub mul { my ($h, $a, $b) = @_; return RefHodge::cmul($a, $b, maxd($h)); }
sub add { my ($a, $b, $c) = @_; return RefHodge::cadd($a, $b, $c // Q(1)); }
sub scal { my ($a, $c) = @_; my %r = map { $_ => $a->{$_} * $c } keys %$a; return \%r; }
sub chm { my ($l, $coef) = @_; my @e = (0) x $RefHodge::NCH; my ($ix) = grep { $RefHodge::CHL[$_] == $l } 0..$RefHodge::NCH-1; die "no ch_$l" unless defined $ix; $e[$ix] = 1; return { join(',', @e) => ($coef // Q(1)) }; }
# Lambda^vee(w)/w^h as ch-polynomial; Omega at numeric tau
sub lamvee { my ($h, $w) = @_; my @lam = lambdas($h); my $f = {}; for my $k (0..$h) { $f = add($f, $lam[$k], (($k % 2) ? Q(-1) : Q(1)) / (Q($w) ** $k)); } return $f; }
# components Omega^(k) = Lambda^vee(1) e_k(alpha + alpha^2), via exact interpolation in v = -1/sigma
sub omega_components { my $h = shift; my @taus = map { Q($_) } qw(1 2 3 1/2 1/3 -1/3 3/2 5 -7/2 7 1/5 -2/7 9 -9/4); @taus = @taus[0..$h];
  my @vs = map { Q(-1) / ($_ * (1 + $_)) } @taus; my @om = map { RefHodge::omega($h, [Q(1), $_, -1 - $_], maxd($h)) } @taus;
  my %mus; for my $o (@om) { $mus{$_} = 1 for keys %$o; } my @comp = map { {} } 0..$h;
  for my $mu (keys %mus) { my @y = map { $om[$_]{$mu} // Q(0) } 0..$h; my @C = interp(\@vs, \@y, $h); for my $k (0..$h) { $comp[$k]{$mu} = $C[$k] if $C[$k] != 0; } }
  return @comp; }
# e_m(beta), beta_i = alpha_i + alpha_i^2 (alpha = Chern roots of E), m = 0..h, from the lambda classes:
#   prod_i (1 + t alpha_i + t alpha_i^2) = c_{-r1}(E) c_{-r2}(E), r1 + r2 = -t, r1 r2 = t,
#   = Sum_{k,l} lambda_k lambda_l (-1)^{k+l} t^{min(k,l)} p_{|k-l|}(r1,r2)/2  (p_0 = 2, p_1 = -t, p_n = -t p_{n-1} - t p_{n-2}).
sub ek_beta { my $h = shift; my @lam = lambdas($h);
  my @pn = ([Q(2)], [Q(0), Q(-1)]);   # p_n as coefficient arrays in t
  for my $n (2..$h) { my @a = map { Q(0) } 0..$n; for my $i (0..$#{$pn[$n-1]}) { $a[$i+1] -= $pn[$n-1][$i]; } for my $i (0..$#{$pn[$n-2]}) { $a[$i+1] -= $pn[$n-2][$i]; } $pn[$n] = \@a; }
  my @E = map { {} } 0..$h;
  for my $k (0..$h) { for my $l (0..$h) { my $P = mul($h, $lam[$k], $lam[$l]); next unless %$P; my $mn = $k < $l ? $k : $l; my $d = abs($k - $l); my $sg = (($k + $l) % 2) ? Q(-1) : Q(1);
      for my $i (0..$#{$pn[$d]}) { my $c = $pn[$d][$i]; next if $c == 0; my $m = $mn + $i; next if $m > $h; $E[$m] = add($E[$m], $P, $sg * $c / 2); } } }
  return @E; }
sub bydeg { my $C = shift; my %B; for my $mu (keys %$C) { my $d = RefHodge::cdeg($mu); $B{$d}{$mu} = $C->{$mu}; } return \%B; }
sub interp { my ($v, $y, $d) = @_; my @C = map { Q(0) } 0..$d;
  for my $i (0..$d) { my @L = (Q(1)); my $den = Q(1);
    for my $m (0..$d) { next if $m == $i; my @N = map { Q(0) } 0..$#L+1; for my $k (0..$#L) { $N[$k+1] += $L[$k]; $N[$k] -= $L[$k] * $v->[$m]; } @L = @N; $den *= ($v->[$i] - $v->[$m]); }
    for my $k (0..$d) { $C[$k] += $y->[$i] * $L[$k] / $den; } }
  return @C; }
1;
